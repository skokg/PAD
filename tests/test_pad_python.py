"""Core correctness tests for an installed pad_python; no optional dependencies.

Run explicitly after installation: python tests/test_pad_python.py
Do not add src to sys.path: wheel tests must exercise the installed module.
"""
from concurrent.futures import ThreadPoolExecutor
import ctypes as ct
import threading
import unittest
from unittest.mock import patch

import numpy as np
import pad_python as pad


class PADTests(unittest.TestCase):
    def assert_conservation(self, a, b, result, normalize=False):
        rows, remaining_a, remaining_b = result
        self.assertEqual(rows.ndim, 2)
        self.assertEqual(rows.shape[1], 6)
        self.assertEqual(remaining_a.shape, a.shape)
        self.assertEqual(remaining_b.shape, b.shape)
        self.assertTrue(np.all(np.isfinite(rows)))
        self.assertTrue(np.all(rows[:, :2] >= 0))
        self.assertTrue(np.all(remaining_a >= 0))
        self.assertTrue(np.all(remaining_b >= 0))
        used_a = np.zeros(a.shape)
        used_b = np.zeros(b.shape)
        for distance, amount, x1, y1, x2, y2 in rows:
            self.assertGreater(amount, 0)
            self.assertTrue(all(value == int(value) for value in (x1, y1, x2, y2)))
            self.assertTrue(0 <= x1 < a.shape[1] and 0 <= y1 < a.shape[0])
            self.assertTrue(0 <= x2 < b.shape[1] and 0 <= y2 < b.shape[0])
            self.assertAlmostEqual(distance, np.hypot(x2 - x1, y2 - y1), places=5)
            used_a[int(y1), int(x1)] += amount
            used_b[int(y2), int(x2)] += amount
        expected_a = a / a.sum() if normalize else a
        expected_b = b / b.sum() if normalize else b
        np.testing.assert_allclose(used_a + remaining_a, expected_a, atol=1e-13)
        np.testing.assert_allclose(used_b + remaining_b, expected_b, atol=1e-13)

    def test_native_abi(self):
        self.assertEqual(pad.libc.PAD_planar_wrapper_abi_version(), 2)

    def test_known_distance_and_endpoint_order(self):
        a, b = np.zeros((5, 7)), np.zeros((5, 7))
        a[1, 1], b[4, 5] = 2, 3
        result = pad.calculate_PAD_attributions(a, b, normalize=False, random_seed=42)
        np.testing.assert_allclose(result[0], [[5, 2, 1, 1, 5, 4]])
        self.assertEqual(result[2][4, 5], 1)
        self.assertEqual(pad.calculate_PAD_distance_from_attributions(result[0]), 5)
        self.assert_conservation(a, b, result)

    def test_identical_fields(self):
        a = np.array([[0., 2., 1.], [3., 0., 4.]])
        result = pad.calculate_PAD_attributions(a, a, random_seed=42)
        self.assertEqual(pad.calculate_PAD_distance_from_attributions(result[0]), 0)
        np.testing.assert_array_equal(result[1], np.zeros_like(a))
        np.testing.assert_array_equal(result[2], np.zeros_like(a))
        self.assert_conservation(a, a, result, normalize=True)

    def test_normalization_preserves_inputs(self):
        a, b = np.array([[2., 0.]]), np.array([[0., 8.]])
        before_a, before_b = a.copy(), b.copy()
        result = pad.calculate_PAD_attributions(a, b, random_seed=0)
        np.testing.assert_allclose(result[0], [[1, 1, 0, 0, 1, 0]])
        np.testing.assert_array_equal(a, before_a)
        np.testing.assert_array_equal(b, before_b)
        self.assert_conservation(a, b, result, normalize=True)

    def test_zero_cutoff_attributes_only_overlap(self):
        a, b = np.array([[3., 2., 0.]]), np.array([[1., 0., 4.]])
        for overlap in (True, False):
            with self.subTest(remove_overlap=overlap):
                result = pad.calculate_PAD_attributions(
                    a, b, normalize=False, remove_overlap=overlap,
                    distance_cutoff=0, random_seed=42,
                )
                np.testing.assert_allclose(result[0], [[0, 1, 0, 0, 0, 0]])
                self.assert_conservation(a, b, result)

    def test_cutoff_boundary_and_empty_result(self):
        a, b = np.zeros((4, 5)), np.zeros((4, 5))
        a[0, 0], b[3, 4] = 1, 1
        for cutoff, count in ((4.99, 0), (5., 1)):
            result = pad.calculate_PAD_attributions(a, b, distance_cutoff=cutoff, random_seed=1)
            self.assertEqual(result[0].shape, (count, 6))
            self.assert_conservation(a, b, result, normalize=True)
            if count == 0:
                with self.assertRaises(ValueError):
                    pad.calculate_PAD_distance_from_attributions(result[0])

    def test_nearest_neighbor_against_independent_distances(self):
        a, b = np.zeros((6, 8)), np.zeros((6, 8))
        a[1, 2] = 1
        b[1, 5], b[5, 7], b[0, 0] = 2, 2, 2
        # The engine starts with field one; its only point is fully consumed.
        result = pad.calculate_PAD_attributions(a, b, normalize=False, random_seed=7)
        candidates = np.argwhere(b > 0)
        expected = np.min(np.hypot(candidates[:, 0] - 1, candidates[:, 1] - 2))
        self.assertAlmostEqual(result[0][0, 0], expected, places=5)
        self.assertEqual(len(result[0]), 1)
        self.assert_conservation(a, b, result)

    def test_every_match_is_nearest_after_repeated_deletions(self):
        # No cutoff and no overlap preprocessing: every iteration emits a row,
        # so rows alternate the selected field, starting with field one.
        # Replay the returned amounts, checking every match against a brute-force
        # search of the still-positive opposite field. No native RNG is replicated
        # and any equally near target is accepted.
        generator = np.random.default_rng(834)
        for shape in ((1, 37), (29, 1), (12, 13)):
            a = generator.integers(0, 9, size=shape).astype(float)
            b = generator.integers(0, 9, size=shape).astype(float)
            for seed in (0, 1, 42, 0xFFFFFFFF):
                with self.subTest(shape=shape, seed=seed):
                    result = pad.calculate_PAD_attributions(
                        a, b, normalize=False, remove_overlap=False, random_seed=seed,
                    )
                    left, right = a.copy(), b.copy()
                    self.assertGreater(len(result[0]), 10)
                    for turn, row in enumerate(result[0]):
                        distance, amount, x1, y1, x2, y2 = row
                        i, j = (int(y1), int(x1)), (int(y2), int(x2))
                        self.assertGreater(left[i], 0)
                        self.assertGreater(right[j], 0)
                        if turn % 2 == 0:
                            candidates = np.argwhere(right > 0)
                            squared = (candidates[:, 0] - i[0]) ** 2 + (candidates[:, 1] - i[1]) ** 2
                        else:
                            candidates = np.argwhere(left > 0)
                            squared = (candidates[:, 0] - j[0]) ** 2 + (candidates[:, 1] - j[1]) ** 2
                        self.assertAlmostEqual(distance, np.sqrt(squared.min()), places=5)
                        self.assertEqual(amount, min(left[i], right[j]))
                        left[i] -= amount
                        right[j] -= amount
                    self.assertTrue(not np.any(left > 0) or not np.any(right > 0))
                    np.testing.assert_array_equal(result[1], left)
                    np.testing.assert_array_equal(result[2], right)
                    self.assert_conservation(a, b, result)

    def test_tied_targets_allow_either_nearest_point(self):
        a, b = np.zeros((3, 3)), np.zeros((3, 3))
        a[1, 1] = 1
        b[1, 0] = b[1, 2] = 1
        for seed in (0, 1, 42):
            result = pad.calculate_PAD_attributions(a, b, normalize=False, random_seed=seed)
            self.assertEqual(result[0].shape, (1, 6))
            self.assertEqual(result[0][0, 0], 1)
            self.assertIn(tuple(result[0][0, 4:]), ((0., 1.), (2., 1.)))
            self.assert_conservation(a, b, result)

    def test_diagonal_cutoff_at_adjacent_float_values(self):
        a, b = np.array([[1., 0.], [0., 0.]]), np.array([[0., 0.], [0., 1.]])
        boundary = np.sqrt(2.)
        for cutoff, count in ((np.nextafter(boundary, 0.), 0),
                              (boundary, 1), (np.nextafter(boundary, np.inf), 1)):
            with self.subTest(cutoff=cutoff):
                result = pad.calculate_PAD_attributions(a, b, distance_cutoff=cutoff, random_seed=42)
                self.assertEqual(len(result[0]), count)
                self.assert_conservation(a, b, result, normalize=True)

    def test_large_grid_nearest_neighbor_and_cutoff_precision(self):
        # Squared distances 100000000 and 100000001 round to the same float.
        # The old search retained the farther root and pruned the closer target,
        # then rejected the root at cutoff 10000, leaving a valid match unused.
        base_a, base_b = np.zeros((10001, 2)), np.zeros((10001, 2))
        base_a[10000, 1] = 1
        base_b[0, 0] = base_b[0, 1] = 1
        for transpose in (False, True):
            a, b = (base_a.T, base_b.T) if transpose else (base_a, base_b)
            expected_row = ([10000, 1, 10000, 1, 0, 1] if transpose
                            else [10000, 1, 1, 10000, 1, 0])
            for cutoff in (np.nextafter(10000., 0.), 10000.,
                           np.nextafter(10000., np.inf), None):
                for overlap in (False, True):
                    for seed in (0, 42):
                        with self.subTest(transpose=transpose, cutoff=cutoff,
                                          overlap=overlap, seed=seed):
                            result = pad.calculate_PAD_attributions(
                                a, b, normalize=False, remove_overlap=overlap,
                                distance_cutoff=cutoff, random_seed=seed,
                            )
                            if cutoff is not None and cutoff < 10000:
                                self.assertEqual(result[0].shape, (0, 6))
                                np.testing.assert_array_equal(result[1], a)
                                np.testing.assert_array_equal(result[2], b)
                            else:
                                np.testing.assert_array_equal(result[0], [expected_row])
                                self.assertEqual(
                                    pad.calculate_PAD_distance_from_attributions(result[0]), 10000,
                                )
                                expected_b = b.copy()
                                expected_b[(1, 0) if transpose else (0, 1)] = 0
                                np.testing.assert_array_equal(result[1], np.zeros_like(a))
                                np.testing.assert_array_equal(result[2], expected_b)
                            self.assert_conservation(a, b, result)

    def test_cutoff_leaves_no_admissible_residual_pair(self):
        generator = np.random.default_rng(907)
        a = generator.integers(0, 9, size=(7, 9)).astype(float)
        b = generator.integers(0, 9, size=(7, 9)).astype(float)
        for cutoff in (0., 1., 2.5):
            for seed in (0, 42, 99):
                for overlap in (True, False):
                    with self.subTest(cutoff=cutoff, seed=seed, overlap=overlap):
                        result = pad.calculate_PAD_attributions(
                            a, b, normalize=False, remove_overlap=overlap,
                            distance_cutoff=cutoff, random_seed=seed,
                        )
                        self.assert_conservation(a, b, result)
                        self.assertTrue(np.all(result[0][:, 0] <= cutoff))
                        left, right = np.argwhere(result[1] > 0), np.argwhere(result[2] > 0)
                        if left.size and right.size:
                            differences = left[:, None, :] - right[None, :, :]
                            squared = np.sum(differences ** 2, axis=2)
                            self.assertTrue(np.all(squared > cutoff ** 2))

    def test_normalization_with_extreme_amounts(self):
        # The first two cases would underflow/overflow in less stable schemes;
        # the last case's direct sum overflows even though each value is finite.
        cases = (
            (np.array([[1e-300, 3e-300, 0.]]), np.array([[0., 1e-300, 3e-300]]),
             np.array([[.25, .75, 0.]]), np.array([[0., .25, .75]])),
            (np.array([[1e307, 3e307, 0.]]), np.array([[0., 1e307, 3e307]]),
             np.array([[.25, .75, 0.]]), np.array([[0., .25, .75]])),
            (np.array([[1e308, 1e308, 0.]]), np.array([[0., 1e308, 1e308]]),
             np.array([[.5, .5, 0.]]), np.array([[0., .5, .5]])),
        )
        for a, b, expected_a, expected_b in cases:
            for overlap in (True, False):
                with self.subTest(a=a, overlap=overlap):
                    result = pad.calculate_PAD_attributions(a, b, remove_overlap=overlap, random_seed=3)
                    reference = pad.calculate_PAD_attributions(
                        expected_a, expected_b, normalize=False, remove_overlap=overlap, random_seed=3,
                    )
                    for actual, expected in zip(result, reference):
                        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-15)
                    self.assert_conservation(expected_a, expected_b, result)

    def test_unnormalized_extreme_and_unequal_amounts(self):
        for amount in (np.nextafter(0., 1.), 1e-300, 1e308):
            a, b = np.array([[amount, 0.]]), np.array([[0., amount]])
            result = pad.calculate_PAD_attributions(a, b, normalize=False, random_seed=3)
            np.testing.assert_array_equal(result[0], [[1., amount, 0., 0., 1., 0.]])
            np.testing.assert_array_equal(result[1], np.zeros_like(a))
            np.testing.assert_array_equal(result[2], np.zeros_like(b))
        a, b = np.array([[1e-300, 0.]]), np.array([[0., 1e300]])
        result = pad.calculate_PAD_attributions(a, b, normalize=False, random_seed=3)
        self.assertEqual(result[0][0, 1], 1e-300)
        self.assertEqual(result[2][0, 1], 1e300)  # Tiny subtraction rounds away.

    def test_inputs_preserved_for_all_options_and_layouts(self):
        base_a = np.array([[0., 1., 3.], [4., 0., 2.]])
        base_b = np.array([[2., 0., 1.], [1., 5., 0.]])
        layouts = {
            "c": lambda x: x.copy(),
            "fortran": np.asfortranarray,
            "transpose": lambda x: x.T,
            "negative_stride": lambda x: x[:, ::-1],
            "float32": lambda x: x.astype(np.float32),
            "integer": lambda x: x.astype(np.int32),
            "big_endian": lambda x: x.astype(">f8"),
        }
        for layout, convert in layouts.items():
            for readonly in (False, True):
                for overlap in (False, True):
                    for normalize in (False, True):
                        for cutoff in (None, 0., 1.):
                            with self.subTest(layout=layout, readonly=readonly, overlap=overlap,
                                              normalize=normalize, cutoff=cutoff):
                                a, b = convert(base_a.copy()), convert(base_b.copy())
                                before_a, before_b = a.copy(), b.copy()
                                a.setflags(write=not readonly)
                                b.setflags(write=not readonly)
                                result = pad.calculate_PAD_attributions(
                                    a, b, remove_overlap=overlap, normalize=normalize,
                                    distance_cutoff=cutoff, random_seed=4,
                                )
                                np.testing.assert_array_equal(a, before_a)
                                np.testing.assert_array_equal(b, before_b)
                                self.assertEqual(a.flags.writeable, not readonly)
                                self.assertEqual(b.flags.writeable, not readonly)
                                self.assert_conservation(a, b, result, normalize)
                                for output in result:
                                    self.assertFalse(np.shares_memory(output, a))
                                    self.assertFalse(np.shares_memory(output, b))
                                    output[...] = -999
                                np.testing.assert_array_equal(a, before_a)
                                np.testing.assert_array_equal(b, before_b)

    def test_results_survive_later_calls(self):
        a, b = np.array([[1., 0., 3.]]), np.array([[0., 4., 0.]])
        result = pad.calculate_PAD_attributions(a, b, random_seed=42)
        snapshots = [array.copy() for array in result]
        for seed in range(10):
            pad.calculate_PAD_attributions(b, a, random_seed=seed)
        for actual, expected in zip(result, snapshots):
            np.testing.assert_array_equal(actual, expected)

    def test_native_errors_are_reported_and_cleared(self):
        # Real buffers with truthful lengths: deliberately invalid VALUES only.
        valid, bad = np.array([1.]), np.array([np.nan])
        count = ct.c_size_t(999)
        pointer = pad.libc.calculate_PAD_results_assume_same_grid_ctypes(
            bad, valid, valid, valid, 1, ct.byref(count), 1., 42,
        )
        try:
            self.assertFalse(pointer)
            self.assertEqual(count.value, 0)
            self.assertIn(b"coordinates", pad.libc.PAD_last_error())
        finally:
            if pointer:
                pad.libc.free_mem_double_array(pointer)
        pointer = pad.libc.calculate_PAD_results_assume_same_grid_ctypes(
            valid, valid, valid, valid, 1, ct.byref(count), 1., 42,
        )
        try:
            self.assertTrue(pointer)
            self.assertEqual(count.value, 1)
            self.assertEqual(pad.libc.PAD_last_error(), b"")
        finally:
            if pointer:
                pad.libc.free_mem_double_array(pointer)

    def test_python_translates_native_failure(self):
        a = np.ones((1, 1))
        for message, expected in ((b"native test failure", "native test failure"),
                                  (b"", "Planar PAD calculation failed")):
            with patch.object(pad.libc, "calculate_PAD_results_assume_same_grid_ctypes",
                              return_value=ct.POINTER(ct.c_double)()), \
                    patch.object(pad.libc, "PAD_last_error", return_value=message):
                with self.assertRaisesRegex(RuntimeError, expected):
                    pad.calculate_PAD_attributions(a, a)

    def test_concurrent_seeded_calls_match_serial_results(self):
        a = np.array([[0., 1., 3.], [4., 0., 2.]])
        b = np.array([[2., 0., 1.], [1., 5., 0.]])
        seeds = (0, 1, 42, 0xFFFFFFFF) * 4
        serial = {seed: pad.calculate_PAD_attributions(a, b, random_seed=seed) for seed in set(seeds)}
        barrier = threading.Barrier(4)

        def calculate(seed):
            barrier.wait(timeout=30)
            return pad.calculate_PAD_attributions(a, b, random_seed=seed)

        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = [executor.submit(calculate, seed) for seed in seeds]
            for seed, future in zip(seeds, futures):
                actual = future.result(timeout=60)
                for output, expected in zip(actual, serial[seed]):
                    np.testing.assert_array_equal(output, expected)

    def test_native_errors_are_thread_local(self):
        barrier = threading.Barrier(2)

        def calculate(invalid):
            x = np.array([np.nan if invalid else 0.])
            values = np.array([1.])
            count = ct.c_size_t()
            pointer = pad.libc.calculate_PAD_results_assume_same_grid_ctypes(
                x, values, values, values, 1, ct.byref(count), 1., 42,
            )
            try:
                # Both calls complete before either reads its error buffer.
                barrier.wait(timeout=30)
                return bool(pointer), pad.libc.PAD_last_error()
            finally:
                if pointer:
                    pad.libc.free_mem_double_array(pointer)

        with ThreadPoolExecutor(max_workers=2) as executor:
            bad = executor.submit(calculate, True)
            good = executor.submit(calculate, False)
            bad_pointer, error = bad.result(timeout=60)
            good_pointer, no_error = good.result(timeout=60)
        self.assertFalse(bad_pointer)
        self.assertIn(b"coordinates", error)
        self.assertTrue(good_pointer)
        self.assertEqual(no_error, b"")

    def test_weighted_mean_extreme_distances_and_weights(self):
        for scale in (1e-300, 1., 1e308):
            rows = np.array([[scale, 1e308, 0, 0, 0, 0],
                             [scale / 2, 1e308, 0, 0, 0, 0]])
            answer = pad.calculate_PAD_distance_from_attributions(rows)
            self.assertTrue(np.isfinite(answer))
            np.testing.assert_allclose(answer, scale * .75, rtol=1e-14, atol=0)

    def test_zero_weight_rows_do_not_change_extreme_weighted_means(self):
        for distance in (0., np.nextafter(0., 1.), 1e-300, 1., 1e308):
            positive = np.array([[distance, 2., 0, 0, 0, 0]])
            baseline = pad.calculate_PAD_distance_from_attributions(positive)
            self.assertEqual(baseline, distance)
            ignored = np.array([[np.finfo(float).max, 0., 0, 0, 0, 0],
                                [0., 0., 0, 0, 0, 0]])
            for rows in (np.vstack((positive, ignored)), np.vstack((ignored, positive))):
                with self.subTest(distance=distance, rows=rows):
                    before = rows.copy()
                    self.assertEqual(pad.calculate_PAD_distance_from_attributions(rows), baseline)
                    np.testing.assert_array_equal(rows, before)

    def test_zero_weight_rows_still_require_valid_distances(self):
        # Filtering contributions must not silently accept invalid input rows.
        for distance in (-1., np.nan, np.inf, -np.inf):
            rows = np.array([[1., 1., 0, 0, 0, 0], [distance, 0., 0, 0, 0, 0]])
            with self.subTest(distance=distance):
                with self.assertRaises(ValueError):
                    pad.calculate_PAD_distance_from_attributions(rows)

    def test_mass_balance_across_options(self):
        a = np.array([[0., 1., 3.], [4., 0., 2.]])
        b = np.array([[2., 0., 1.], [1., 5., 0.]])
        for seed in (0, 42, 0xFFFFFFFF):
            for overlap in (True, False):
                for normalize in (True, False):
                    for cutoff in (None, 0., 1.):
                        with self.subTest(seed=seed, overlap=overlap,
                                          normalize=normalize, cutoff=cutoff):
                            result = pad.calculate_PAD_attributions(
                                a, b, random_seed=seed, remove_overlap=overlap,
                                normalize=normalize, distance_cutoff=cutoff,
                            )
                            self.assert_conservation(a, b, result, normalize)
                            if cutoff is not None:
                                self.assertTrue(np.all(result[0][:, 0] <= cutoff))

    def test_seed_repeats_within_environment(self):
        a = np.array([[0., 1., 3.], [4., 0., 2.]])
        b = np.array([[2., 0., 1.], [1., 5., 0.]])
        first = pad.calculate_PAD_attributions(a, b, random_seed=123)
        second = pad.calculate_PAD_attributions(a, b, random_seed=123)
        for left, right in zip(first, second):
            np.testing.assert_array_equal(left, right)

    def test_noncontiguous_integer_inputs(self):
        a = np.array([[0, 1], [3, 2]], dtype=np.int32).T
        b = np.array([[2, 0], [1, 4]], dtype=np.int32).T
        result = pad.calculate_PAD_attributions(a, b, normalize=False, random_seed=42)
        self.assert_conservation(a, b, result)

    def test_invalid_fields(self):
        good = np.ones((2, 2))
        invalid = [
            [[1., 1.], [1., 1.]], np.ones(4), np.empty((0, 2)),
            np.ones((2, 3)), np.zeros((2, 2)), -good,
            np.array([[np.nan, 1], [1, 1]]),
            np.array([[np.inf, 1], [1, 1]]), good.astype(complex),
            np.ma.array(good, mask=False), np.array([["a", "b"], ["c", "d"]]),
        ]
        for bad in invalid:
            for a, b in ((bad, good), (good, bad)):
                with self.subTest(value=repr(bad)):
                    with self.assertRaises((TypeError, ValueError)):
                        pad.calculate_PAD_attributions(a, b)

    def test_invalid_options(self):
        a = np.ones((1, 1))
        for options in (
            {"distance_cutoff": -1}, {"distance_cutoff": np.inf},
            {"distance_cutoff": np.nan}, {"random_seed": True},
            {"random_seed": 1.5}, {"random_seed": -2},
            {"random_seed": 0x100000000}, {"normalize": 1}, {"remove_overlap": 1},
        ):
            with self.subTest(options=options):
                with self.assertRaises((TypeError, ValueError)):
                    pad.calculate_PAD_attributions(a, a, **options)

    def test_weighted_mean_independent_answer(self):
        rows = np.array([[3., 2., 0, 0, 3, 0], [7., 6., 0, 0, 7, 0]])
        self.assertAlmostEqual(pad.calculate_PAD_distance_from_attributions(rows), 6)
        rows[:, 1] *= 1e307
        self.assertAlmostEqual(pad.calculate_PAD_distance_from_attributions(rows), 6)

    def test_invalid_weighted_means(self):
        for rows in (
            np.empty((0, 6)), np.ones((2, 4)), np.ones(6),
            [[-1, 1, 0, 0, 0, 0]], [[1, -1, 0, 0, 0, 0]],
            [[1, 0, 0, 0, 0, 0]], [[np.inf, 1, 0, 0, 0, 0]],
            [[1, np.nan, 0, 0, 0, 0]],
        ):
            with self.subTest(rows=repr(rows)):
                with self.assertRaises(ValueError):
                    pad.calculate_PAD_distance_from_attributions(rows)


if __name__ == "__main__":
    unittest.main(verbosity=2)
