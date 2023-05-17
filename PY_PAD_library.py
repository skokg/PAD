
import numpy as np
import os
from ctypes import *

# search for the PAD C++ shared library file (PAD_Cxx_shared_library.so) in the same folder
libc = CDLL(os.path.abspath(os.path.expanduser(os.path.dirname(__file__)))+ os.path.sep + "PAD_Cxx_shared_library.so") 

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

def check_input_fields(fa, fb):
	
	# test if fields are numpy arrays
	if type(fa) is not np.ndarray or type(fa) is not np.ndarray:
		print("ERROR: the fa and fa are not numpy arrays of type numpy.ndarray - this is not permitted! Perhaps they are masked arrays, which is also not permitted. Returning \"None\" as result!")
		return(False)
	
	# compare dimensions of fields
	if fa.shape != fb.shape:
		print("ERROR: the fa and fb arrays do not have the same shape. Returning \"None\" as result!")
		return(False)
	
	# check dimensions of fields
	if fa.ndim != 2 or fb.ndim != 2:
		print("ERROR: the fa and fb arrays are not two-dimensional. Returning \"None\" as result!")
		return(False)
	
	# compare the array has some elements
	if fa.size == 0 or fb.size == 0:
		print("ERROR: the dimensions of fa or fb arrays are zero. Returning \"None\" as result!")
		return(False)
	
	# detect non-numeric values
	result=np.where(np.isfinite(fa) == False)
	if len(result[0]) > 0:
		print("ERROR: fa or fb arrays contain some non-numeric values. Returning \"None\" as result!")
		return(False)
	result=np.where(np.isfinite(fb) == False)
	if len(result[0]) > 0:
		print("ERROR: fa or fb arrays contain some non-numeric values. Returning \"None\" as result!")
		return(False)
	
	# detect masked array
	if isinstance(fa, np.ma.MaskedArray) or isinstance(fb, np.ma.MaskedArray) :
		print("ERROR: fa or fb arrays are masked arrays which is not allowed. Returning \"None\" as result!")
		return(False)
	
	# detect negative values
	if fa[fa<0].size > 0 or fb[fb<0].size > 0 :
		print("ERROR: fa or fb arrays contain some negative values which is not allowed . Returning \"None\" as result!")
		return(False)
	
	# check if any field is empty 
	if np.sum(fa) == 0 or np.sum(fb) == 0 :
		print("ERROR: fa or fb fields contain only zeroes - this is not allowed since PAD value cannot be defined in this case. Returning \"None\" as result!")
		return(False)
	
	return(True)

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

def check_if_max_number_of_nonzero_points_integer(max_number_of_nonzero_points):
	
	#check x is integer
	if isinstance(max_number_of_nonzero_points, int) != True:
		print("ERROR: Parameter max_number_of_nonzero_points is not an integer. Returning \"None\" as result!")
		return(False)
	
	return(True)



# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

ND_POINTER_2D = np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags="C")

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

libc.calculate_PAD_ctypes.argtypes = [ND_POINTER_2D, ND_POINTER_2D, c_size_t, c_size_t, c_size_t]
libc.calculate_PAD_ctypes.restype = c_double

def calculate_PAD(fa, fb, max_number_of_nonzero_points = 100000):
	
	# check input fields
	if check_input_fields(fa, fb) != True:
		return(None)
	
	#check max_number_of_nonzero_points is integer
	if check_if_max_number_of_nonzero_points_integer(max_number_of_nonzero_points) != True:
		return(None)
	
	# if needed make a copy and cast to float64 - just in case the input fields are integers
	if fa.dtype != np.float64:
		fa=fa.astype(float64, copy=True)
	if fb.dtype != np.float64:
		fb=fb.astype(float64, copy=True)
	
	PADtemp = libc.calculate_PAD_ctypes(fa, fb, fa.shape[1], fa.shape[0], max_number_of_nonzero_points)
	
	
	return(PADtemp)

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

libc.calculate_PAD_attributions_ctypes.argtypes = [ND_POINTER_2D, ND_POINTER_2D, c_size_t, c_size_t, c_size_t, POINTER(c_size_t)]
libc.calculate_PAD_attributions_ctypes.restype = POINTER(c_double)

libc.free_mem_double_array.argtypes = [POINTER(c_double)]
libc.free_mem_double_array.restype = None

def calculate_PAD_attributions(fa, fb, max_number_of_nonzero_points = 100000):
	
	# check input fields
	if check_input_fields(fa, fb) != True:
		return(None)
	
	#check max_number_of_nonzero_points is integer
	if check_if_max_number_of_nonzero_points_integer(max_number_of_nonzero_points) != True:
		return(None)
	
	# if needed make a copy and cast to float64 - just in case the input fields are integers
	if fa.dtype != np.float64:
		fa=fa.astype(float64, copy=True)
	if fb.dtype != np.float64:
		fb=fb.astype(float64, copy=True)
	
	c_number_of_attributions = c_size_t()
	
	results = libc.calculate_PAD_attributions_ctypes(fa,fb,fa.shape[1],fa.shape[0], max_number_of_nonzero_points, byref(c_number_of_attributions))
	
	number_of_attributions = c_number_of_attributions.value
	
	#print(number_of_attributions)
	
	attributions = np.asarray(results[0:number_of_attributions*6]).reshape(number_of_attributions,-1)
	
	libc.free_mem_double_array(results)
	
	return(attributions)

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------


def calculate_PAD_attribution_PDF(fa, fb, max_number_of_nonzero_points = 100000):
	
	PAD_attributions = calculate_PAD_attributions(fa, fb, max_number_of_nonzero_points=max_number_of_nonzero_points)
	
	if type(fa) is not np.ndarray:
		return(None)
	
	maxdistance=np.ceil(np.max(PAD_attributions[:,0])) + 1
	
	hist = np.histogram(PAD_attributions[:,0], bins=int(maxdistance), range=(0,maxdistance), density=True, weights=PAD_attributions[:,1])
	
	return( np.asarray([hist[1][:-1]+0.5, hist[0][:]]).transpose(1,0))

