// Linux build from this directory (planar ABI 2; old binaries must be rebuilt):
// g++ -std=c++11 -O2 -Wall -Wextra -pthread -shared -fPIC -o ../PAD_Cxx_shared_library.so CC_PAD_python_lib.cc
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <memory>
#include <mutex>
#include <random>
#include <stdexcept>
#include <vector>

// Only the utilities needed by PAD; failures throw rather than exiting Python.
namespace {
std::mt19937 rng_state_global;

double ran2(long *seed)
{
    if (*seed < 0)
    {
        if (*seed == std::numeric_limits<long>::min())
            throw std::invalid_argument("Random seed cannot be LONG_MIN.");
        rng_state_global.seed(-*seed);
        *seed = -*seed;
    }
    return std::uniform_real_distribution<double>(0.0, 1.0)(rng_state_global);
}
}

#include "CU_PAD_code.cc"

namespace {
thread_local char last_error[1024] = {};
// ctypes releases the GIL, so protect PAD's shared random generator.
std::mutex calculation_mutex;

void save_error(const char *message) noexcept
{
    std::snprintf(last_error, sizeof(last_error), "%s", message ? message : "Unknown PAD error");
}

PADPoints make_points(const double *x, const double *y, const double *values, std::size_t size)
{
    if (!x || !y || !values || size == 0)
        throw std::invalid_argument("PAD requires nonempty, non-null input buffers.");
    if (size > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) ||
        size - 1 > std::numeric_limits<kdtree::IndexType>::max())
        throw std::length_error("Too many PAD grid points.");
    PADPoints points;
    points.reserve(size);
    for (std::size_t i = 0; i < size; ++i)
        points.push_back({x[i], y[i], values[i]});
    // The core validates all rows, including zero-valued points, before building.
    return points;
}

double *pack_results(const PADResults &results, const std::vector<double> &remaining1,
                     const std::vector<double> &remaining2)
{
    const std::size_t limit = std::numeric_limits<std::size_t>::max() / sizeof(double);
    if (remaining1.size() > limit || remaining2.size() > limit - remaining1.size())
        throw std::length_error("PAD output is too large.");
    const std::size_t tail = remaining1.size() + remaining2.size();
    if (results.size() > (limit - tail) / 4)
        throw std::length_error("PAD output is too large.");
    std::unique_ptr<double[]> buffer(new double[results.size() * 4 + tail]);
    std::size_t position = 0;
    for (const PADResult &row : results)
    {
        buffer[position++] = row.distance;
        buffer[position++] = row.attributed_amount;
        buffer[position++] = static_cast<double>(row.index1);
        buffer[position++] = static_cast<double>(row.index2);
    }
    for (double value : remaining1) buffer[position++] = value;
    for (double value : remaining2) buffer[position++] = value;
    return buffer.release();
}

// The caller supplies accessible buffers of the stated lengths. No exception
// crosses the C boundary; nullptr indicates an error, NOT an empty result.
double *calculate(const double *x1, const double *y1, const double *values1, std::size_t size1,
                  const double *x2, const double *y2, const double *values2, std::size_t size2,
                  std::size_t *count, double euclidian_cutoff, std::int64_t random_seed,
                  bool same_grid) noexcept
{
    last_error[0] = '\0';
    if (count) *count = 0;
    try
    {
        if (!count) throw std::invalid_argument("Missing attribution-count output pointer.");
        std::lock_guard<std::mutex> lock(calculation_mutex);
        validate_PAD_cutoff(euclidian_cutoff);
        validate_PAD_random_seed(random_seed);
        auto points1 = make_points(x1, y1, values1, size1);
        PADPoints points2;
        if (same_grid)
        {
            if (!values2 || size2 != size1)
                throw std::invalid_argument("Same-grid fields require matching non-null input buffers.");
            points2.reserve(size2);
            for (std::size_t i = 0; i < size2; ++i)
                points2.push_back({points1[i][0], points1[i][1], values2[i]});
        }
        else points2 = make_points(x2, y2, values2, size2);
        std::vector<double> remaining1, remaining2;
        PADResults results = same_grid
            ? calculate_PAD_results_assume_same_grid_and_remove_overlap(
                points1, points2, euclidian_cutoff, remaining1, remaining2, random_seed)
            : calculate_PAD_results_assume_different_grid(
                points1, points2, euclidian_cutoff, remaining1, remaining2, random_seed);
        double *buffer = pack_results(results, remaining1, remaining2);
        *count = results.size();
        return buffer;
    }
    catch (const std::exception &error) { save_error(error.what()); }
    catch (...) { save_error("Unknown C++ exception during planar PAD calculation."); }
    return nullptr;
}
}

// Use a planar-specific ABI symbol so a spherical or legacy .so is rejected.
extern "C" int PAD_planar_wrapper_abi_version() noexcept { return 2; }
extern "C" const char *PAD_last_error() noexcept { return last_error; }
extern "C" void free_mem_double_array(double *buffer) noexcept { delete[] buffer; }

// ABI 2: four doubles per result, followed by size1 and size2 residual amounts.
// Cutoff is Euclidean in coordinate units; DBL_MAX is unrestricted. Seed -1
// chooses a random seed. Amounts are NOT normalized at this C/C++ layer.
extern "C" double *calculate_PAD_results_assume_same_grid_ctypes(
    const double *x, const double *y, const double *values1, const double *values2,
    std::size_t size, std::size_t *count, double euclidian_cutoff, std::int64_t random_seed) noexcept
{
    return calculate(x, y, values1, size, x, y, values2, size,
                     count, euclidian_cutoff, random_seed, true);
}

extern "C" double *calculate_PAD_results_assume_different_grid_ctypes(
    const double *x1, const double *y1, const double *values1, std::size_t size1,
    const double *x2, const double *y2, const double *values2, std::size_t size2,
    std::size_t *count, double euclidian_cutoff, std::int64_t random_seed) noexcept
{
    return calculate(x1, y1, values1, size1, x2, y2, values2, size2,
                     count, euclidian_cutoff, random_seed, false);
}
