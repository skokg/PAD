"""Planar PAD using the C++ wrapper (planar ABI 2; rebuild the shared library).

All results have six columns: Euclidean distance, amount, x1, y1, x2, y2.
The 2D-field API optionally normalizes each field (enabled by default) and always
return [attributions, remaining1, remaining2]. No cutoff is used by default.
"""
import ctypes as ct
import operator
from pathlib import Path

import numpy as np

_library_path = Path(__file__).resolve().parent / "PAD_Cxx_shared_library.so"
try:
    libc = ct.CDLL(str(_library_path))
except OSError as exc:
    raise ImportError(
        "Cannot load PAD_Cxx_shared_library.so. Rebuild it for this platform; "
        "see source_for_Cxx_shared_library/HOW_TO_COMPILE.txt."
    ) from exc
try:
    libc.PAD_planar_wrapper_abi_version.argtypes = []
    libc.PAD_planar_wrapper_abi_version.restype = ct.c_int
except AttributeError as exc:
    raise ImportError(
        "Rebuild PAD_Cxx_shared_library.so: planar wrapper ABI 2 is required."
    ) from exc
if libc.PAD_planar_wrapper_abi_version() != 2:
    raise ImportError("Incompatible PAD shared library: planar wrapper ABI 2 is required.")

ND_POINTER_1D = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags=("C_CONTIGUOUS", "ALIGNED")
)
libc.free_mem_double_array.argtypes = [ct.POINTER(ct.c_double)]
libc.free_mem_double_array.restype = None
libc.PAD_last_error.argtypes = []
libc.PAD_last_error.restype = ct.c_char_p
libc.calculate_PAD_results_assume_same_grid_ctypes.argtypes = [
    ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_1D,
    ct.c_size_t, ct.POINTER(ct.c_size_t), ct.c_double, ct.c_int64,
]
libc.calculate_PAD_results_assume_same_grid_ctypes.restype = ct.POINTER(ct.c_double)
libc.calculate_PAD_results_assume_different_grid_ctypes.argtypes = [
    ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_1D, ct.c_size_t,
    ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_1D, ct.c_size_t,
    ct.POINTER(ct.c_size_t), ct.c_double, ct.c_int64,
]
libc.calculate_PAD_results_assume_different_grid_ctypes.restype = ct.POINTER(ct.c_double)


def _check_array(array, name, ndim):
    if not isinstance(array, np.ndarray) or isinstance(array, np.ma.MaskedArray):
        raise TypeError(f"{name} must be an unmasked NumPy array.")
    if array.ndim != ndim or array.size == 0:
        raise ValueError(f"{name} must be a nonempty {ndim}-dimensional array.")
    if not np.issubdtype(array.dtype, np.number) or np.iscomplexobj(array):
        raise TypeError(f"{name} must contain real numeric values.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values.")


def _check_amounts(values, name):
    if np.any(values < 0) or not np.any(values > 0):
        raise ValueError(f"{name} must be nonnegative with at least one positive amount.")
    if values.size > np.iinfo(np.int32).max:
        raise ValueError("Too many grid points for the C++ tree.")


def _array(array, name):
    _check_array(array, name, 1)
    converted = np.require(array, dtype=np.float64, requirements=["C", "A"])
    if not np.all(np.isfinite(converted)):
        raise ValueError(f"{name} must be representable as finite float64 values.")
    return converted


def check_input_fields(fa, fb):
    """Validate equally shaped, nonnegative 2D fields; raise on invalid input."""
    for field, name in ((fa, "fa"), (fb, "fb")):
        _check_array(field, name, 2)
        _check_amounts(field, name)
    if fa.shape != fb.shape:
        raise ValueError("fa and fb must have identical shapes.")
    return True


def _seed(random_seed):
    if random_seed is None:
        return -1
    if isinstance(random_seed, (bool, np.bool_)):
        raise TypeError("random_seed must be an integer, not a boolean.")
    seed = operator.index(random_seed)
    if seed < -1 or seed > 0xFFFFFFFF:
        raise ValueError("random_seed must be None, -1, or an unsigned 32-bit integer.")
    return seed


def _cutoff(distance_cutoff):
    if distance_cutoff is None:
        return np.finfo(np.float64).max
    value = float(distance_cutoff)
    if not np.isfinite(value) or value < 0:
        raise ValueError("Euclidean distance cutoff must be finite and nonnegative, or None.")
    return value


def _normalized_values(values):
    """Normalize already validated float64 amounts without modifying inputs."""
    scaled = values / values.max()
    return np.ascontiguousarray(scaled / scaled.sum(), dtype=np.float64)


def calculate_PAD_attributions(
    fa, fb, remove_overlap=True, distance_cutoff=None, random_seed=None,
    normalize=True,
):
    """Calculate planar PAD from two equally shaped, equidistant 2D fields.

    Coordinates are generated from columns (x) and rows (y), shared by both
    fields. Separate coordinate arrays, 1D inputs and different grids are not
    supported by this Python API.

    remove_overlap=True attributes pointwise overlap at distance zero first.
    False skips preprocessing; zero-distance NN matches can still occur.
    normalize=True independently scales each field to sum to one before
    attribution; False keeps original amounts. Inputs are never modified.

    distance_cutoff=None means no cutoff. A finite nonnegative Euclidean cutoff
    is in grid-cell units, as are returned distances. Zero is valid.
    random_seed=None or -1 chooses and prints a seed when NN work is needed;
    uint32 seeds reproduce runs for identical inputs and library versions.

    Always returns [attributions, remaining1, remaining2]. Attributions is
    float64 with shape (N, 6), including (0, 6) when no matches are allowed:
    distance, amount, x1, y1, x2, y2. Endpoint coordinates are zero-based
    column (x) and row (y) positions, as in the original planar Python API.
    Residuals have the original 2D field shapes, in normalized units when
    normalize=True and original amount units otherwise.
    """
    check_input_fields(fa, fb)
    if not isinstance(remove_overlap, (bool, np.bool_)):
        raise TypeError("remove_overlap must be a boolean.")
    if not isinstance(normalize, (bool, np.bool_)):
        raise TypeError("normalize must be a boolean.")
    seed, cutoff = _seed(random_seed), _cutoff(distance_cutoff)
    values1 = _array(fa.ravel(order="C"), "fa")
    values2 = _array(fb.ravel(order="C"), "fb")
    # Check again after conversion, which may underflow extended-precision data.
    _check_amounts(values1, "fa")
    _check_amounts(values2, "fb")
    if normalize:
        values1 = _normalized_values(values1)
        values2 = _normalized_values(values2)
    x = np.tile(np.arange(fa.shape[1], dtype=np.float64), fa.shape[0])
    y = np.repeat(np.arange(fa.shape[0], dtype=np.float64), fa.shape[1])
    count = ct.c_size_t()
    if remove_overlap:
        result = libc.calculate_PAD_results_assume_same_grid_ctypes(
            x, y, values1, values2, values1.size, ct.byref(count), cutoff, seed,
        )
    else:
        # Use the general C++ engine on the SAME generated grid to skip overlap.
        result = libc.calculate_PAD_results_assume_different_grid_ctypes(
            x, y, values1, values1.size, x, y, values2, values2.size,
            ct.byref(count), cutoff, seed,
        )
    if not result:
        message = libc.PAD_last_error()
        raise RuntimeError(
            message.decode("utf-8", errors="replace") if message else "Planar PAD calculation failed."
        )
    try:
        n = count.value
        if n > values1.size + values2.size:
            raise RuntimeError("Invalid attribution count returned by PAD.")
        total = 4 * n + values1.size + values2.size
        # One owned copy instead of a list of Python float objects.
        packed = np.ctypeslib.as_array(result, shape=(total,)).copy()
        # The internal C++ buffer retains indices; expose original Python XY output.
        records = packed[:4 * n].reshape(n, 4)
        index1 = records[:, 2].astype(np.intp)
        index2 = records[:, 3].astype(np.intp)
        attributions = np.empty((n, 6), dtype=np.float64)
        attributions[:, :2] = records[:, :2]
        attributions[:, 2] = x[index1]
        attributions[:, 3] = y[index1]
        attributions[:, 4] = x[index2]
        attributions[:, 5] = y[index2]
        remaining1 = packed[4 * n:4 * n + values1.size]
        remaining2 = packed[4 * n + values1.size:]
        return [attributions, remaining1.reshape(fa.shape), remaining2.reshape(fb.shape)]
    finally:
        libc.free_mem_double_array(result)


def calculate_PAD_distance_from_attributions(PAD_attributions):
    """Return the amount-weighted mean distance; empty results are undefined."""
    rows = np.asarray(PAD_attributions, dtype=np.float64)
    if rows.ndim != 2 or rows.shape[1] != 6 or rows.shape[0] == 0:
        raise ValueError("PAD requires a nonempty six-column attribution array.")
    distances, weights = rows[:, 0], rows[:, 1]
    if (not np.all(np.isfinite(distances)) or not np.all(np.isfinite(weights))
            or np.any(distances < 0) or np.any(weights < 0) or not np.any(weights > 0)):
        raise ValueError("PAD requires finite nonnegative distances and positive total weight.")
    weights = weights / weights.max()
    scale = distances.max()
    if scale == 0:
        return 0.0
    fraction = np.sum((distances / scale) * weights) / np.sum(weights)
    return float(min(fraction, 1.0) * scale)
