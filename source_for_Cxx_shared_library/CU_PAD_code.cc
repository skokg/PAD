#include "CU_kdtree_with_index.cc"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

// Planar equivalent of the spherical PAD engine. Rows are {x, y, amount}.
// Distances and cutoffs use the coordinate units; no spherical conversion or
// normalization is performed here. The wrapper supplies rng_state_global/ran2.
constexpr std::size_t PAD_dimensions = 2;
using PADPoint = kdtree::Point_str<PAD_dimensions>;
using PADKdTree = kdtree::KdTree<PAD_dimensions>;
using PADPoints = std::vector<std::vector<double>>;

struct PADResult
{
    double distance; // Euclidean, in the same units as x/y.
    double attributed_amount;
    std::size_t index1; // Original input order, including zero-valued points.
    std::size_t index2;
};
using PADResults = std::vector<PADResult>;

void validate_PAD_random_seed(const std::int64_t random_seed)
{
    if (random_seed < -1 || random_seed > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument("PAD random seed must be -1 (random) or an unsigned 32-bit value.");
}

void validate_PAD_cutoff(const double euclidian_attribution_distance_cutoff)
{
    if (!std::isfinite(euclidian_attribution_distance_cutoff) || euclidian_attribution_distance_cutoff < 0)
        throw std::domain_error("Euclidean PAD cutoff must be finite and non-negative.");
}

void validate_PAD_point(const std::vector<double> &point)
{
    if (point.size() != PAD_dimensions + 1)
        throw std::invalid_argument("Planar PAD points must contain x, y and an amount.");
    // The unchanged tree accumulates squared distances in float. Leave ample
    // headroom for differences, squares and their sum, not just coordinate casts.
    const double coordinate_limit = std::sqrt(std::min<double>(
        std::numeric_limits<float>::max(), std::numeric_limits<kdtree::PointType>::max())) /
        (4.0 * std::sqrt(static_cast<double>(PAD_dimensions)));
    for (std::size_t axis = 0; axis < PAD_dimensions; ++axis)
        if (!std::isfinite(point[axis]) || std::fabs(point[axis]) > coordinate_limit)
            throw std::domain_error("Planar coordinates must be finite and within the tree's safe squared-distance range; rescale or recenter them.");
    if (!std::isfinite(point.back()) || point.back() < 0)
        throw std::domain_error("PAD amounts must be finite and non-negative.");
}

std::vector<double> check_points(const PADPoints &points)
{
    if (points.empty())
        throw std::invalid_argument("PAD requires nonempty point arrays.");
    if (points.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) ||
        points.size() - 1 > std::numeric_limits<kdtree::IndexType>::max())
        throw std::length_error("Too many PAD grid points.");
    std::vector<double> values;
    values.reserve(points.size());
    bool positive = false;
    for (const auto &point : points)
    {
        validate_PAD_point(point);
        values.push_back(point.back());
        positive = positive || point.back() > 0;
    }
    if (!positive)
        throw std::domain_error("PAD requires at least one positive amount in each field.");
    return values;
}

PADPoint kdtree_Point_str_from_point(const std::vector<double> &point, const std::size_t index)
{
    PADPoint result;
    for (std::size_t axis = 0; axis < PAD_dimensions; ++axis)
        result.coords[axis] = static_cast<kdtree::PointType>(point[axis]);
    result.set_index(index);
    return result;
}

double squared_euclidian_distance(const PADPoint &a, const PADPoint &b)
{
    // Match findNearestNode_in_radius: products in PointType, sum in double.
    double sum = 0;
    for (std::size_t axis = 0; axis < PAD_dimensions; ++axis)
        sum += (a.coords[axis] - b.coords[axis]) * (a.coords[axis] - b.coords[axis]);
    return sum;
}

void construct_kdtree_with_shuffled_index_list(
    const PADPoints &points, const std::vector<double> &values,
    std::vector<std::size_t> &indices, PADKdTree &tree)
{
    if (points.size() != values.size() || !tree.is_empty())
        throw std::invalid_argument("Invalid PAD tree-construction state.");
    indices.clear();
    std::vector<PADPoint> active;
    for (std::size_t i = 0; i < points.size(); ++i)
        if (values[i] > 0)
        {
            active.push_back(kdtree_Point_str_from_point(points[i], i));
            indices.push_back(i);
        }
    std::shuffle(indices.begin(), indices.end(), rng_state_global);
    if (!tree.buildKdTree(active))
        throw std::runtime_error("PAD k-d tree construction failed.");
}

// False means the selected point has no acceptable neighbor. Its amount stays
// in the residual output, but it is removed from further attribution searches.
bool perform_one_PAD_iteration_with_attribution_distance_cutoff(
    const PADPoints &points1, std::vector<double> &values1, std::vector<double> &values2,
    std::vector<std::size_t> &indices1, std::vector<std::size_t> &indices2,
    PADKdTree &tree1, PADKdTree &tree2, long &idum, const bool first_field_turn,
    const double squared_euclidian_attribution_distance_cutoff, PADResult &result)
{
    const std::size_t index1 = indices1.back();
    const PADPoint point1 = kdtree_Point_str_from_point(points1[index1], index1);
    const auto *node = tree2.findNearestNode_in_radius(point1, squared_euclidian_attribution_distance_cutoff);
    const bool matched = node != nullptr;
    if (!matched)
    {
        tree1.deleteNode(point1);
        indices1.pop_back();
    }
    else
    {
        // Copy before deletion: deleting a node can replace its stored point.
        const PADPoint point2 = node->val;
        const std::size_t index2 = point2.index;
        const double amount = std::min(values1[index1], values2[index2]);
        values1[index1] -= amount;
        values2[index2] -= amount;
        if (values1[index1] == 0)
            tree1.deleteNode(point1);
        else
            std::swap(indices1.back(), indices1[static_cast<std::size_t>(
                std::floor(ran2(&idum) * static_cast<double>(indices1.size())))]);
        if (values2[index2] == 0)
            tree2.deleteNode(point2);
        result = PADResult{std::sqrt(squared_euclidian_distance(point1, point2)), amount,
            first_field_turn ? index1 : index2, first_field_turn ? index2 : index1};
    }
    while (!indices1.empty() && values1[indices1.back()] == 0) indices1.pop_back();
    while (!indices2.empty() && values2[indices2.back()] == 0) indices2.pop_back();
    return matched;
}

// Shared engine: called after validation/overlap preprocessing, with a positive
// residual in each field. Appends to any pre-existing zero-distance records.
// Shared RNG, like spherical PAD; the ctypes wrapper serializes calls.
void calculate_PAD_results_general(
    const PADPoints &points1, const PADPoints &points2, const double euclidian_attribution_distance_cutoff,
    std::vector<double> &values1, std::vector<double> &values2, PADResults &out,
    const std::int64_t random_seed = -1)
{
    validate_PAD_cutoff(euclidian_attribution_distance_cutoff);
    validate_PAD_random_seed(random_seed);
    std::uint32_t selected_seed;
    if (random_seed == -1)
    {
        std::random_device source;
        selected_seed = std::uniform_int_distribution<std::uint32_t>(
            0, std::numeric_limits<std::uint32_t>::max())(source);
    }
    else selected_seed = static_cast<std::uint32_t>(random_seed);
    rng_state_global.seed(selected_seed);
    std::cout << "PAD random seed: " << selected_seed << std::endl;
    const double squared_cutoff = euclidian_attribution_distance_cutoff > std::sqrt(std::numeric_limits<double>::max())
        ? std::numeric_limits<double>::infinity()
        : euclidian_attribution_distance_cutoff * euclidian_attribution_distance_cutoff;

    PADKdTree tree1, tree2;
    std::vector<std::size_t> indices1, indices2;
    construct_kdtree_with_shuffled_index_list(points1, values1, indices1, tree1);
    construct_kdtree_with_shuffled_index_list(points2, values2, indices2, tree2);
    const std::size_t available = out.max_size() - out.size();
    if (indices1.size() > available || indices2.size() > available - indices1.size())
        throw std::length_error("PAD result capacity exceeds vector::max_size().");
    out.reserve(out.size() + indices1.size() + indices2.size());
    long idum = 1;
    bool first_field_turn = true;
    while (!indices1.empty() && !indices2.empty())
    {
        PADResult result;
        const bool matched = first_field_turn
            ? perform_one_PAD_iteration_with_attribution_distance_cutoff(points1, values1, values2,
                indices1, indices2, tree1, tree2, idum, true, squared_cutoff, result)
            : perform_one_PAD_iteration_with_attribution_distance_cutoff(points2, values2, values1,
                indices2, indices1, tree2, tree1, idum, false, squared_cutoff, result);
        if (matched) out.push_back(result);
        first_field_turn = !first_field_turn;
    }
    // Trees own contiguous storage and release it on destruction, also on errors.
}

void validate_PAD_outputs(const PADPoints &points1, const PADPoints &points2,
                          const std::vector<double> &values1, const std::vector<double> &values2)
{
    if (&values1 == &values2)
        throw std::invalid_argument("PAD requires distinct residual vectors.");
    for (const auto &point : points1)
        if (&point == &values1 || &point == &values2)
            throw std::invalid_argument("PAD residual vectors must not alias input rows.");
    for (const auto &point : points2)
        if (&point == &values1 || &point == &values2)
            throw std::invalid_argument("PAD residual vectors must not alias input rows.");
}

PADResults calculate_PAD_results_assume_same_grid_and_remove_overlap(
    const PADPoints &points1, const PADPoints &points2, const double euclidian_attribution_distance_cutoff,
    std::vector<double> &values1, std::vector<double> &values2, const std::int64_t random_seed = -1)
{
    validate_PAD_outputs(points1, points2, values1, values2);
    validate_PAD_cutoff(euclidian_attribution_distance_cutoff);
    validate_PAD_random_seed(random_seed);
    if (points1.size() != points2.size())
        throw std::invalid_argument("Same-grid PAD requires matching point counts.");
    auto prepared1 = check_points(points1);
    auto prepared2 = check_points(points2);
    for (std::size_t i = 0; i < points1.size(); ++i)
        for (std::size_t axis = 0; axis < PAD_dimensions; ++axis)
            if (points1[i][axis] != points2[i][axis])
                throw std::invalid_argument("Same-grid PAD requires identical coordinates in identical order.");
    values1.swap(prepared1);
    values2.swap(prepared2);
    PADResults out;
    bool positive1 = false, positive2 = false;
    for (std::size_t i = 0; i < values1.size(); ++i)
    {
        const double overlap = std::min(values1[i], values2[i]);
        if (overlap > 0)
        {
            values1[i] -= overlap;
            values2[i] -= overlap;
            out.push_back(PADResult{0, overlap, i, i});
        }
        positive1 = positive1 || values1[i] > 0;
        positive2 = positive2 || values2[i] > 0;
    }
    if (positive1 && positive2)
        calculate_PAD_results_general(points1, points2, euclidian_attribution_distance_cutoff,
                                      values1, values2, out, random_seed);
    return out;
}

// Also usable on identical grids when overlap must NOT be removed beforehand.
PADResults calculate_PAD_results_assume_different_grid(
    const PADPoints &points1, const PADPoints &points2, const double euclidian_attribution_distance_cutoff,
    std::vector<double> &values1, std::vector<double> &values2, const std::int64_t random_seed = -1)
{
    validate_PAD_outputs(points1, points2, values1, values2);
    validate_PAD_cutoff(euclidian_attribution_distance_cutoff);
    validate_PAD_random_seed(random_seed);
    auto prepared1 = check_points(points1);
    auto prepared2 = check_points(points2);
    values1.swap(prepared1);
    values2.swap(prepared2);
    PADResults out;
    calculate_PAD_results_general(points1, points2, euclidian_attribution_distance_cutoff,
                                  values1, values2, out, random_seed);
    return out;
}

double calculate_PAD_from_PAD_results(const PADResults &results)
{
    double max_weight = 0, max_distance = 0;
    for (const auto &row : results)
    {
        if (!std::isfinite(row.distance) || row.distance < 0 ||
            !std::isfinite(row.attributed_amount) || row.attributed_amount < 0)
            throw std::domain_error("PAD results require finite, non-negative distances and weights.");
        max_weight = std::max(max_weight, row.attributed_amount);
        max_distance = std::max(max_distance, row.distance);
    }
    if (max_weight == 0)
        throw std::domain_error("PAD is undefined without positive attributed weight.");
    if (max_distance == 0) return 0;
    double weighted = 0, weights = 0;
    for (const auto &row : results)
    {
        const double weight = row.attributed_amount / max_weight;
        weighted += (row.distance / max_distance) * weight;
        weights += weight;
    }
    return std::min(weighted / weights, 1.0) * max_distance;
}

