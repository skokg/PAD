

// ******************************************************
// ******************************************************
// Class CField. Stores one 2-D field with dimension atributes and raw data.
//
class CField {
	public:
	int dimx, dimy;
	vector <double> data;
	vector <double> lat_data;
	vector <double> lon_data;
	// used optionally
	long unique_timestamp;

	// constructor
	CField() {freememory();}
	~CField() {freememory();}

	// allocates memory to data pointer using dimx, dimy
	void allocatememory()
		{
		data.assign((long)dimx*(long)dimy,0);
		}

	// allocates memory to data pointer using dimx, dimy
	void allocatememory_and_set_dimx_dimy(int dimx_, int dimy_)
		{
		dimx=dimx_;
		dimy=dimy_;
		allocatememory();
		}


	// allocates memory for lat and lon array using dimx, dimy
	void allocatememory_lat_lon_array()
		{
		lat_data.assign((long)dimy,0);
		lon_data.assign((long)dimx,0);
		}

	// frees memory
	void freememory()
		{
		data.clear();
		lat_data.clear();
		lon_data.clear();
		unique_timestamp=0;
		dimx=0;
		dimy=0;
		}

	string output()
		{
		if (!check_if_initilized())
			error(AT,FU, "!check_if_initilized()");

		ostringstream s1;
		s1.str("");

		long ix,iy;
		for (iy=0; iy < dimy; iy++)
			{
			s1 << iy << ": ";
			for (ix=0; ix < dimx; ix++)
				s1 << get(ix,iy) << " ";
			s1 << endl;
			}
		return(s1.str());
		}


	void output_Cfield_for_ctypes(vector <double> *v, int *dimx_, int* dimy_) const
		{
		*v=data;
		*dimx_=dimx;
		*dimy_=dimy;
		}


	bool check_if_initilized() const
		{
		if (lat_data.size()==0 || lon_data.size()==0 || data.size()==0) return (false);
		return(true);
		}

	bool check_if_data_is_initilized() const
		{
		if (data.size()==0 ) return (false);
		return(true);
		}

	bool compare_have_they_the_same_dimensions(const CField &f) const
		{
		if (dimx != f.dimx) return(false);
		if (dimy != f.dimy) return(false);
		return(true);
		}

	bool compare_have_they_the_same_dimensions_and_coordinates(const CField &f) const
		{
		if (lat_data != f.lat_data) return(false);
		if (lon_data != f.lon_data) return(false);
		if (dimx != f.dimx) return(false);
		if (dimy != f.dimy) return(false);
		return(true);
		}


	bool compare_are_they_the_same(const CField &f) const
		{
		if (unique_timestamp != f.unique_timestamp) return(false);
		if (lat_data != f.lat_data) return(false);
		if (lon_data != f.lon_data) return(false);
		if (data != f.data) return(false);
		return(true);
		}

	// overload operators == nad !=
	bool operator== (const CField& f) const
		{return(compare_are_they_the_same(f));}
	bool operator!= (const CField& f) const
		{return(!compare_are_they_the_same(f));}

	// area size of field
	long size() const
		{
		return((long)dimx*(long)dimy);
		}

	// gets value at x,y
	double get(int x, int y) const
		{
		return(data[(long)y*(long)dimx + (long)x]);
		}

	// calculates x,y location from the array offset
	void DataLoc2XY(long ix, int &x, int &y) const
		{
		x=ix%dimx;
		y=ix/dimx;
		}

	// calculates array offset from x,y
	long XY2DataLoc(int x, int y) const
		{
		return((long)y*(long)dimx + (long)x);
		}

	// sets all values of a field to value
	void Set_all_values(double value)
		{
		long ixl;
		for (ixl=0;ixl<size();ixl++)
			data[ixl]=value;
		}

	double get_sum() const
		{
		double sum=0;
		for (long il=0; il < size(); il++)
			if (data[il] != BAD_DATA_FLOAT)
				sum+=data[il];
		return(sum);
		}

 	bool check_if_field_has_no_missing_values() const
		{
		for (long il=0; il < size(); il++)
			if (data[il] == BAD_DATA_FLOAT)
				return(false);

		return(true);
		}

	void multiply_by_double(double x)
		{
		for (long il=0; il < size(); il++)
  			if (data[il] != BAD_DATA_FLOAT)
				data[il]*=x;
		}

};
