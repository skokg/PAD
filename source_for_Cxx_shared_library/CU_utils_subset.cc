
// error function
void error(const char *location,const char *function, string txt)
	{
	cout << "ERROR: " << location << " " << function << ": " <<  txt << endl;
	exit(1);
	}

// this foreces the vector to free its memory - only using .clear() does not free the memory !
template<typename T>
void free_vector_memory(std::vector<T>& vec)
	{
	std::vector<T>().swap(vec);
	}

template<typename T>
string output_vector_as_string(const std::vector<T>& vec, string separator)
	{
	ostringstream s1;
	for (long il=0; il < (long)vec.size(); il++)
		{
		s1 << vec[il];
		if (il< (long)vec.size() - 1)
			s1 << separator;
		}
	return(s1.str());
	}



template<typename T>
T sum_vector(std::vector<T>& vec)
	{
	T sum=0;
	for (long il=0; il < (long)vec.size(); il++)
		sum+=vec[il];
	return(sum);
	}



double deg2rad(double kot)
	{
	return (kot*M_PI/180.0);
	}

void spherical_to_cartesian_coordinates(double lat, double lon, double r, double &x, double &y, double &z)
	{
	x=r*cos(lon) * cos(lat);
	y=r*sin(lon) * cos(lat);
	z=r*sin(lat);
	}


template<typename T>
double squared_euclidian_distance_multidimensional(const vector <T> &p1, const vector <T> &p2)
	{
	ERRORIF(p1.size() != p2.size());
	double sum=0;
	for (unsigned long il=0; il < p1.size(); il++)
		sum+=(p1[il]-p2[il])*(p1[il]-p2[il]);

	return(sum);
	}

// ----------------------------------
// random number generator (RNG) based on random.h library (requires C++11).
// --------------------------------------

std::mt19937 rng_state_global;  // global instance of RNG state

// wrapper function for random number generator based on random.h library - it is not threadsafe
double ran2(long *idum)
	{
	if (*idum < 0)
		{
		rng_state_global.seed(-*idum);
		*idum=-*idum;
		}
	return(std::uniform_real_distribution<double>(0.0, 1.0)(rng_state_global));
	}


void generate_random_binary_vector_with_fixed_number_of_1_values(long p,long N,vector <double> &v, long &seed)
	{
	v.clear();
	if (N < p)
		error(AT,FU, "N < p");

	long count=0;
	long pos;

	double to_set;
	double to_fill_on_start;
	long how_many_to_add;

	if ((double)p/(double)N < 0.5)
		{
		to_set=1;
		to_fill_on_start=0;
		how_many_to_add=p;
		}
	else
		{
		to_set=0;
		to_fill_on_start=1;
		how_many_to_add=N-p;
		}

	v.assign(N,to_fill_on_start);

	while (count < how_many_to_add)
			{
			pos = floor(ran2(&seed)*(double)N);
			if (pos > N - 1)
				error(AT,FU, "pos > N - 1");

			if (v[pos]==to_fill_on_start)
				{
				v[pos]=to_set;
				count++;
				}
			}

	// make a final check
	if (sum_vector(v) != p)
		error(AT,FU, "sum_vector(v) != p");
	}

