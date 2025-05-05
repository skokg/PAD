# The PAD Python Software Package

#### Description:

Precipitation Attribution Distance (PAD) is a spatial verification measure for precipitation. It is based on a random nearest-neighbor attribution concept - it works by sequentially attributing randomly selected precipitation in one field to the closest precipitation in the other. Overall the PAD provides a good and meaningful estimate of precipitation displacement in two fields that tends to be in line with a subjective forecast evaluation. For a more detailed description of PAD, please refer to the paper listed in the References section. 

The python package can be used for efficient calculation of PAD value and related parameters. It uses a k-d tree implementation for fast computation, and its time complexity for calculating the PAD value is O(n*log(n)) while the memory complexity is O(n). The input fields need to be available in a rectangular domain with an equidistant grid (a more general version of the package that can be used on a sphere is availabe on https://github.com/skokg/PAD_on_Sphere).  

The underlying code is written in C++, and a precompiled shared library file is available for easy use with python on Linux systems (the C++ source code is available in the `source_for_Cxx_shared_library` folder - it can be used to compile the shared library for other types of systems). Python ctypes library is used to access the functions in the shared library file. 

#### Usage:

To see how the package can be used please refer to the two examples: PY_PAD_example_01.py and PY_PAD_example_02.py.

#### Author:

Gregor Skok, Faculty of Mathematics and Physics, University of Ljubljana, Slovenia

Email: Gregor.Skok@fmf.uni-lj.si

#### References:

Skok, G. (2023) Precipitation attribution distance. Atmospheric Research, 295. [ https://doi.org/10.1016/j.atmosres.2023.106998]( https://doi.org/10.1016/j.atmosres.2023.106998)
