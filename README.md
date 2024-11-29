# The PAD Python Software Package

#### Description:

Precipitation Attribution Distance (PAD) is a spatial verification measure for precipitation. It is based on a random nearest-neighbor attribution concept - it works by sequentially attributing randomly selected precipitation in one field to the closest precipitation in the other. Overall the PAD provides a good and meaningful estimate of precipitation displacement in two fields that tends to be in line with a subjective forecast evaluation. For a more detailed description of PAD, please refer to the paper listed in the References section. 

The python package can be used for efficient calculation of PAD value and related parameters. It uses a k-d tree implementation for fast computation, and its time complexity for calculating the PAD value is O(n*log(n)) while the memory complexity is O(n). The input fields need to be available in a rectangular domain with an equidistant grid.  

The underlying code is written in C++, and a precompiled shared library file is available for easy use with python on Linux systems (the C++ source code is available in the `source` folder - it can be used to compile the shared library for other types of systems). Python ctypes library is used to access the functions in the shared library file. 

An example of how the package can be used is available (for more information, please refer to the Example section).

#### Usage:

There are three functions available: 

- function calculate_PAD(*fa, fb, max_number_of_nonzero_points*)

This function returns a single float value - the PAD value. If there is a problem with the calculation, the function displays an error warning and returns "None".

- function calculate_PAD_attribution_PDF(*fa, fb, max_number_of_nonzero_points*)

This function returns a two-dimensional numpy array representing the PAD attribution PDF (Probability Density Function). The array's dimensions are (number_of_bins, 2). Each row represents one bin, with the first value representing the bin center (expressed in grid points) and the second representing the corresponding PDF value. The bin width is one grid point; the first bin always starts at 0 grid points. If there is a problem with the calculation, the function displays an error warning and returns "None".

- function calculate_PAD_attributions(*fa, fb, max_number_of_nonzero_points*)

This function returns a two-dimensional numpy array representing the list of all PAD attributions. The array's dimensions are (number_of_all_attributions, 6). Each row represents one attribution, described with six values: attribution distance, attributed value, x coordinate in the fa field, y coordinate in the fa field, x coordinate in the fb field, and y coordinate in the fb field. If there is a problem with the calculation, the function displays an error warning and returns "None".

#### Function arguments:

- *fa* - a two-dimensional numpy array representing the first field. Only positive and zero values are allowed.

- *fb* - a two-dimensional numpy array representing the second field. Only positive and zero values are allowed. The array's dimensions need to be the same as the dimensions of *fa*. 

- *max_number_of_nonzero_points* - the maximum number of nonzero points to be used for the calculation (default value: *max_number_of_nonzero_points = 100000*). If *max_number_of_nonzero_points* is smaller than the number of all nonzero points in a field, the nonzero points will be randomly subsampled. If *max_number_of_nonzero_points < 0*, then all the nonzero points will always be used (there will not be any subsampling).

#### Example:

An example of how the package can be used is available in the file PY_PAD_example.py. 

#### Author:

Gregor Skok, Faculty of Mathematics and Physics, University of Ljubljana, Slovenia

Email: Gregor.Skok@fmf.uni-lj.si

#### References:

Skok, G. (2023) Precipitation attribution distance. Atmospheric Research, 295. [ https://doi.org/10.1016/j.atmosres.2023.106998]( https://doi.org/10.1016/j.atmosres.2023.106998)
