
// error function
void error(const char *location,const char *function, string txt)
	{
	cout << "ERROR: " << location << " " << function << ": " <<  txt << endl;
	exit(1);
	}


double squared_euclidian_distance(double x1, double y1, double x2, double y2)
	{
	return((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	}

template<typename T>
T sum_vector(std::vector<T>& vec)
	{
	T sum=0;
	for (long il=0; il < (long)vec.size(); il++)
		sum+=vec[il];
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


