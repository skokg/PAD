// Compile command for a linux system
// g++ -fopenmp -O2 -Wall -Wno-unused-result -Wno-unknown-pragmas -shared -o PAD_Cxx_shared_library.so -fPIC CC_PAD_python_lib.cc


#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <sys/resource.h>
#include <string.h>
#include <random>
#include <chrono>
#include <cstring>

using namespace std;

#define BAD_DATA_FLOAT -9999

// da prav dela error(..) - da prav displaya line number
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define FU __PRETTY_FUNCTION__
#define TT "\t"
#define tabt "\t"
#define ERRORIF(x) if (x) error(AT,FU, #x)

const rlim_t kStackSize = 1000 * 1024 * 1024;   // min stack size = 16 MB

long random_number_seed=-1;

//#include "CU_utils.cc"
#include "CU_utils_subset.cc"
#include "CU_PAD_code.cc"

vector <vector <double> > X_Y_val_arrays_to_XY_points(const double * const x, const double * const y, const double * const values, const size_t size)
	{
	vector <vector <double> > points;
	for (size_t il = 0; il < size; il++)
		{
		points.push_back({x[il],y[il],values[il]});
		}
	return(points);
	}

extern "C" void free_mem_double_array(double* a)
	{
	delete[] a;
	}

extern "C"  double * calculate_PAD_results_assume_different_grid_ctypes(const double * const x1, const double * const y1, const double * const values1, const size_t size1, const double * const x2, const double * const y2, const double * const values2, const size_t size2, size_t * const number_of_attributions)
	{
    vector <vector <double> > points1 =  X_Y_val_arrays_to_XY_points(x1,y1,values1, size1);
    vector <vector <double> > points2 =  X_Y_val_arrays_to_XY_points(x2,y2,values2, size2);

	vector <double> non_attributed_values1;
	vector <double> non_attributed_values2;

	vector <vector <double>> results = calculate_PAD_results_assume_different_grid(points1, points2, non_attributed_values1, non_attributed_values2);

	// SERIALIZE the output into a double vector
	vector <double> out;
	for (unsigned long il=0; il < results.size(); il++)
		{
		out.push_back(results[il][0]);
		out.push_back(results[il][1]);
		out.push_back(results[il][2]);
		out.push_back(results[il][3]);
		}
	out.insert(out.end(), non_attributed_values1.begin(), non_attributed_values1.end());
	out.insert(out.end(), non_attributed_values2.begin(), non_attributed_values2.end());

	double* out_arr = new double[out.size()];
	std::copy(out.begin(), out.end(), out_arr);

	*number_of_attributions = results.size();

	return(out_arr);
	}

