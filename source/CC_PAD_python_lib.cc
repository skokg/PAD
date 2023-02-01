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

#include "CU_utils_subset.cc"
#include "CU_CField_class_subset.cc"
#include "CU_PAD_code.cc"

CField array_to_Cfield(const double *indatav1, size_t dimx, size_t dimy)
	{
	CField f;
	f.allocatememory_and_set_dimx_dimy(dimx, dimy);
	f.allocatememory_lat_lon_array();
	size_t i;
	for (i = 0; i < dimx*dimy; ++i)
		f.data[i] = indatav1[i];

	return(f);
	}

double* CField_data_to_array(const CField &f)
	{
	double* out = new double[f.data.size()];
	std::copy(f.data.begin(), f.data.end(), out);
	return(out);
	}

extern "C" void free_mem_double_array(double* a)
	{
	delete[] a;
	}

extern "C"  double calculate_PAD_ctypes(const double *indatav1, const double *indatav2, size_t dimx, size_t dimy, size_t max_number_of_points)
	{
    CField f1 = array_to_Cfield(indatav1, dimx, dimy);
    CField f2 = array_to_Cfield(indatav2, dimx, dimy);
	return(calculate_PAD_with_subset_of_points(f1, f2, max_number_of_points));
	}


extern "C"  double * calculate_PAD_attributions_ctypes(const double *indatav1, const double *indatav2, size_t dimx, size_t dimy, size_t max_number_of_points, size_t *number_of_attributions)
	{
    CField f1 = array_to_Cfield(indatav1, dimx, dimy);
    CField f2 = array_to_Cfield(indatav2, dimx, dimy);

    //cout << f1.output() << endl;
    //cout << "---------------" << endl;
    //cout << f2.output() << endl;
    //cout << "---------------" << endl;

	vector <vector <double>> results = calculate_PAD_results_with_subset_of_points(f1,f2, max_number_of_points);

	//for (unsigned long il=0; il < results.size(); il++)
	//	cout << output_vector_as_string(results[il], " ") << endl;


	vector <double> out;
	for (unsigned long il=0; il < results.size(); il++)
		{
		out.push_back(results[il][0]);
		out.push_back(results[il][1]);
		int x1,y1;
		f1.DataLoc2XY(results[il][2], x1, y1);
		out.push_back((double)x1);
		out.push_back((double)y1);
		f2.DataLoc2XY(results[il][3], x1, y1);
		out.push_back((double)x1);
		out.push_back((double)y1);
		}

	//cout << output_vector_as_string(out, " ") << endl;

	double* out_arr = new double[out.size()];
	std::copy(out.begin(), out.end(), out_arr);

	*number_of_attributions = results.size();

	//cout << "aaa" << endl;

	return(out_arr);
	}
