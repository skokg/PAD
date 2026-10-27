# The PAD Python Software Package

#### Description:

Precipitation Attribution Distance (PAD) is a spatial verification measure for precipitation. It is based on a random nearest-neighbor attribution concept - it works by sequentially attributing randomly selected precipitation in one field to the closest precipitation in the other. Overall the PAD provides a good and meaningful estimate of precipitation displacement in two fields that tends to be in line with a subjective forecast evaluation. For a more detailed description of PAD, please refer to the paper listed in the References section.

The Python package can be used for efficient calculation of PAD value and related parameters. It uses a k-d tree implementation for fast computation, with O(n) memory complexity; runtime depends on the point distribution and tree shape. The input fields must be equally shaped 2D NumPy arrays on a shared rectangular, equidistant grid. For spherical geometry, see https://github.com/skokg/PAD_on_Sphere.

The underlying code is written in C++, and a precompiled shared library file is available for easy use with Python on Linux systems (the C++ source code is available in the `source_for_Cxx_shared_library` folder - it can be used to compile the shared library for other types of systems). Python ctypes library is used to access the functions in the shared library file. Place a compatible `PAD_Cxx_shared_library.so` in the same folder as `PAD_library.py`.

To build the shared library on Linux, run `make` in this folder (requires GNU Make and a C++11-compatible compiler, defaulting to `g++`). This creates `PAD_Cxx_shared_library.so` beside `PAD_library.py`. Rebuild any older shared library to match the current planar wrapper ABI 2. Run `make clean` to remove the shared library.

#### Usage:

To see how the package can be used please refer to the two examples: PAD_example_01.py and PAD_example_02.py. Run PAD_example_02.py from this folder so it can load the NetCDF fields from the `sample_data` subfolder.

`calculate_PAD_attributions(fa, fb)` returns `[attributions, remaining1, remaining2]`. Attributions have six columns: distance, amount, x1, y1, x2, y2, with zero-based column (x) and row (y) coordinates. Distances are Euclidean, in grid-cell units. Residual arrays retain the input shapes.

By default, `normalize=True` independently normalizes each field to sum to one, `remove_overlap=True` attributes pointwise overlap at zero distance first, and `distance_cutoff=None` applies no cutoff. A finite nonnegative cutoff may be supplied in grid-cell units. Attribution amounts and residuals use normalized units unless `normalize=False`. An optional `random_seed` gives reproducible results for identical inputs and library versions; otherwise a seed is randomly chosen.

Inputs must contain finite, nonnegative values with at least one positive value per field; masked arrays are not supported. Invalid inputs raise exceptions. If no matches satisfy the cutoff, attributions have shape `(0, 6)`; check for nonempty results before calculating the weighted PAD distance.

#### Author:

Gregor Skok, Faculty of Mathematics and Physics, University of Ljubljana, Slovenia

Email: Gregor.Skok@fmf.uni-lj.si

#### References:

Skok, G. (2023) Precipitation attribution distance. Atmospheric Research, 295. [ https://doi.org/10.1016/j.atmosres.2023.106998](https://doi.org/10.1016/j.atmosres.2023.106998)
