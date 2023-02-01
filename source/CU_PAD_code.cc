#include "CU_kdtree.cc"

void remove_point_from_the_kdtree(const CField &f_norm, long ind, kdtree::KdTree &tree)
	{
	// remove point from kdtree
	int ix, iy;
	f_norm.DataLoc2XY(ind, ix, iy);
	tree.deleteNode({(float)ix,(float)iy});
	}


vector <double> perform_one_PAD_iteration(CField &f1_norm, CField &f2_norm, const CField &f1_norm_orig, const CField &f2_norm_orig, vector <long> &list_nzp_f1, vector <long> &list_nzp_f2, kdtree::KdTree &tree1, kdtree::KdTree &tree2, long &idum, double maximal_reduction_factor, bool f1_is_fa)
	{
	// choose the last point from f1
	long ind1=list_nzp_f1.back();

	// find the closes non-zero point in f2
	// check same location first
	long ind2=-1;
	double min_squared_distance=1E100;
	if (f2_norm.data[ind1] > 0)
		{
		ind2=ind1;
		min_squared_distance=0;
		}
	// else search the border points of f2
	else
		{
		int x1,y1;
		f1_norm.DataLoc2XY(ind1, x1, y1);
		auto node = tree2.findNearestNode({(float)x1,(float)y1});
		ERRORIF(node == nullptr);
		ind2=f2_norm.XY2DataLoc(node->val[0], node->val[1]);
		min_squared_distance=squared_euclidian_distance(x1, y1, node->val[0], node->val[1]);
		}

	int xx1,yy1,xx2,yy2;
	f1_norm.DataLoc2XY(ind1, xx1, yy1);
	f2_norm.DataLoc2XY(ind2, xx2, yy2);
	//cout << ind1 << " " << ind2 << "; " << xx1 << " " << yy1 << " " << xx2 << " " << yy2 << endl;

	// value reduction - this depends on the original normalized values and the actual normalized values
	double value_reduction=min( min(f1_norm_orig.data[ind1],f2_norm_orig.data[ind2]) * maximal_reduction_factor, min(f1_norm.data[ind1],f2_norm.data[ind2]));
	f1_norm.data[ind1]-=value_reduction;
	f2_norm.data[ind2]-=value_reduction;
	//cout << "val.reduction " << value_reduction << endl;
	//cout << f1_norm.data[ind1] << endl;
	//cout << f2_norm.data[ind2] << endl;

	//cout << list_nzp_f1.size() << " " << list_nzp_f2.size() << endl;

	// if necessary remove the point from tree1
	if (f1_norm.data[ind1] == 0)
		remove_point_from_the_kdtree(f1_norm, ind1, tree1);
	// else swap with random point so this is not the last point anymore - so it is not automatically chosen the next time
	else
		swap(list_nzp_f1.back(),list_nzp_f1[ floor(ran2(&idum) * (double)list_nzp_f1.size()) ]);

	// if necessary remove the point from tree2
	if (f2_norm.data[ind2] == 0)
		remove_point_from_the_kdtree(f2_norm, ind2, tree2);

	//cout << list_nzp_f1.size() << " " << list_nzp_f2.size() << endl;

	vector <double> result;
	result.push_back(sqrt(min_squared_distance));
	result.push_back(value_reduction);
	if (f1_is_fa)
		{
		result.push_back(ind1);
		result.push_back(ind2);
		}
	else
		{
		result.push_back(ind2);
		result.push_back(ind1);
		}

	// remove all zero points at the end of vector
	while (list_nzp_f1.size() > 0 && f1_norm.data[list_nzp_f1.back()] == 0)
		list_nzp_f1.pop_back();

	// remove all zero points at the end of vector
	while (list_nzp_f2.size() > 0 && f2_norm.data[list_nzp_f2.back()] == 0)
		list_nzp_f2.pop_back();

	return(result);
	}

vector <vector <double>> calculate_PAD_results(const CField &fa, const CField &fb)
	{


	ERRORIF(!fa.compare_have_they_the_same_dimensions(fb));
	// check missing data
	ERRORIF(!fa.check_if_field_has_no_missing_values());
	ERRORIF(!fb.check_if_field_has_no_missing_values());

	// check for negative values
	for (long il=0; il < fa.size(); il++)
		ERRORIF( fa.data[il] < 0 || fb.data[il] < 0 );

	double fasum=fa.get_sum();
	double fbsum=fb.get_sum();

	ERRORIF(fasum == 0 || fbsum == 0);

	// normalize by sum
	CField fa_norm=fa;
	CField fb_norm=fb;
	fa_norm.multiply_by_double(1.0/(fasum));
	fb_norm.multiply_by_double(1.0/(fbsum));

	CField fa_norm_orig=fa_norm;
	CField fb_norm_orig=fb_norm;


	//WriteField("bbb1.nc",fa_norm);
	//WriteField("bbb2.nc",fb_norm);


	// generate lists of non-zero points
	vector <long> list_nzp_fa;
	vector <long> list_nzp_fb;
	for (long il=0; il < fa_norm.size(); il++)
		{
		if (fa_norm.data[il] > 0) list_nzp_fa.push_back(il);
		if (fb_norm.data[il] > 0) list_nzp_fb.push_back(il);
		}

	//cout << list_nzp_fa.size() << endl;
	//cout << list_nzp_fb.size() << endl;

	// shuffle points
	std::shuffle(list_nzp_fa.begin(), list_nzp_fa.end(), rng_state_global);
	std::shuffle(list_nzp_fb.begin(), list_nzp_fb.end(), rng_state_global);

	// generate points for kdtree
	vector<kdtree::Point> pointsa;
	vector<kdtree::Point> pointsb;
	for (unsigned long ip=0; ip < list_nzp_fa.size(); ip++)
		{
		int x1,y1;
		kdtree::Point p;
		fa_norm.DataLoc2XY(list_nzp_fa[ip], x1, y1);
		p={(float)x1,(float)y1};
		pointsa.push_back(p);
		}
	for (unsigned long ip=0; ip < list_nzp_fb.size(); ip++)
		{
		int x1,y1;
		kdtree::Point p;
		fb_norm.DataLoc2XY(list_nzp_fb[ip], x1, y1);
		p={(float)x1,(float)y1};
		pointsb.push_back(p);
		}
	// generate kdtrees
	kdtree::KdTree treea;
	treea.buildKdTree(pointsa);
	kdtree::KdTree treeb;
	treeb.buildKdTree(pointsb);

	//cout << list_bp_fa.size() << endl;
	//cout << list_bp_fb.size() << endl;



	//long idum=-1;
	long idum=1;
	double maximal_reduction_factor=1;

	/*for (long it = 0; it<10; it++)
		{
		perform_one_PAD_iteration(fa_norm, fb_norm, fa_norm_orig, fb_norm_orig, list_nzp_fa, list_nzp_fb, list_bp_fa, list_bp_fb, idum, maximal_reduction_factor);

		cout << list_nzp_fa.size() << endl;
		cout << list_nzp_fb.size() << endl;
		cout << list_bp_fa.size() << endl;
		cout << list_bp_fb.size() << endl;
		}
	*/


	vector <vector <double>> out;
	bool fa_turn=true;
	//long counter=0;
	while (list_nzp_fa.size() > 0 && list_nzp_fb.size() > 0)
		{
		vector <double> result;
		/*cout << "----" << endl;
		for (unsigned long il=0; il < list_nzp_fa.size(); il++)
			cout << list_nzp_fa[il] << "(" << fa_norm.data[list_nzp_fa[il]] << ") ";
		cout << endl;
		for (unsigned long il=0; il < list_nzp_fb.size(); il++)
			cout << list_nzp_fb[il] << "(" << fb_norm.data[list_nzp_fb[il]] << ") ";
		cout << endl;*/

		if (fa_turn)
			{
			result = perform_one_PAD_iteration(fa_norm, fb_norm, fa_norm_orig, fb_norm_orig, list_nzp_fa, list_nzp_fb, treea, treeb, idum, maximal_reduction_factor, fa_turn);
			fa_turn=false;
			}
		else
			{
			result = perform_one_PAD_iteration(fb_norm, fa_norm, fb_norm_orig, fa_norm_orig, list_nzp_fb, list_nzp_fa, treeb, treea, idum, maximal_reduction_factor, fa_turn);
			fa_turn=true;
			}
		/*for (unsigned long il=0; il < list_nzp_fa.size(); il++)
			cout << list_nzp_fa[il] << "(" << fa_norm.data[list_nzp_fa[il]] << ") ";
		cout << endl;
		for (unsigned long il=0; il < list_nzp_fb.size(); il++)
			cout << list_nzp_fb[il] << "(" << fb_norm.data[list_nzp_fb[il]] << ") ";
		cout << endl;
		*/
		//vector<std::shared_ptr<kdtree::KdTreeNode>> clustera = treea.findNearestNodeCluster({0,0}, 10000);
		//vector<std::shared_ptr<kdtree::KdTreeNode>> clusterb = treeb.findNearestNodeCluster({0,0}, 10000);
		//cout << counter << " | " << list_nzp_fa.size() << " " << list_nzp_fb.size() << " | " << clustera.size() << " " << clusterb.size() << endl;
		//cout << counter << " | " << list_nzp_fa.size() << " " << list_nzp_fb.size() << endl;

		out.push_back(result);
		//counter++;
		}

	//WriteField("ddd1.nc",fa_norm);
	//WriteField("ddd2.nc",fb_norm);

	//for (unsigned long ip=0; ip < list_bp_fa.size(); ip++)
	//	fa_norm.data[list_bp_fa[ip]]*=2;
	//for (unsigned long ip=0; ip < list_bp_fb.size(); ip++)
	//	fb_norm.data[list_bp_fb[ip]]*=2;
	//WriteField("ccc1.nc",fa_norm);
	//WriteField("ccc2.nc",fb_norm);


	//for (unsigned long il=0; il < out.size(); il++)
	//	cout << il << "\t" << output_vector_as_string(out[il],"\t") << endl;

	return(out);
	}

vector <vector <double>> calculate_PAD_results_with_subset_of_points(const CField &fa, const CField &fb, size_t max_number_of_nonzero_points)
	{
	// use all points if max_number_of_nonzero_points < 1
	if (max_number_of_nonzero_points < 1)
		return(calculate_PAD_results(fa, fb));

	ERRORIF(!fa.compare_have_they_the_same_dimensions(fb));
	// check missing data
	ERRORIF(!fa.check_if_field_has_no_missing_values());
	ERRORIF(!fb.check_if_field_has_no_missing_values());

	// check for negative values
	for (long il=0; il < fa.size(); il++)
		ERRORIF( fa.data[il] < 0 || fb.data[il] < 0 );

	double fasum=fa.get_sum();
	double fbsum=fb.get_sum();

	ERRORIF(fasum == 0 || fbsum == 0);

	// generate lists of non-zero points
	vector <long> list_nzp_fa;
	vector <long> list_nzp_fb;
	for (long il=0; il < fa.size(); il++)
		{
		if (fa.data[il] > 0) list_nzp_fa.push_back(il);
		if (fb.data[il] > 0) list_nzp_fb.push_back(il);
		}

	CField fa_subset=fa;
	CField fb_subset=fb;

	if (list_nzp_fa.size() > max_number_of_nonzero_points)
		{
		fa_subset.Set_all_values(0);
		vector <double> v;
		long idum=1;
		generate_random_binary_vector_with_fixed_number_of_1_values(max_number_of_nonzero_points,list_nzp_fa.size(),v, idum);
		for (unsigned long ip=0; ip < list_nzp_fa.size(); ip++)
			if (v[ip] == 1)
				fa_subset.data[list_nzp_fa[ip]] = fa.data[list_nzp_fa[ip]];
		cout << "Notice: Subsampling " << max_number_of_nonzero_points << " nonzero points from fa field." << endl;
		}

	if (list_nzp_fb.size() > max_number_of_nonzero_points)
		{
		fb_subset.Set_all_values(0);
		vector <double> v;
		long idum=1;
		generate_random_binary_vector_with_fixed_number_of_1_values(max_number_of_nonzero_points,list_nzp_fb.size(),v, idum);
		for (unsigned long ip=0; ip < list_nzp_fb.size(); ip++)
			if (v[ip] == 1)
				fb_subset.data[list_nzp_fb[ip]] = fb.data[list_nzp_fb[ip]];
		cout << "Notice: Subsampling " << max_number_of_nonzero_points << " nonzero points from fb field."  << endl;
		}

	//WriteField("eee1.nc",fa_subset);
	//WriteField("eee2.nc",fb_subset);

	return(calculate_PAD_results(fa_subset, fb_subset));
	}



double calculate_PAD_from_PAD_results(const vector <vector <double>> &results)
	{
	double sum_weights=0;
	double sum=0;
	for (unsigned long il=0; il < results.size(); il++)
		{
		sum+=results[il][0]*results[il][1];
		sum_weights+=results[il][1];
		}

	return(sum/sum_weights);
	}

double calculate_PAD_with_subset_of_points(const CField &fa, const CField &fb, size_t  max_number_of_nonzero_points)
	{
	vector <vector <double>> results = calculate_PAD_results_with_subset_of_points(fa,fb, max_number_of_nonzero_points);
	return(calculate_PAD_from_PAD_results(results));
	}


double calculate_PAD(const CField &fa, const CField &fb)
	{
	vector <vector <double>> results = calculate_PAD_results(fa,fb);
	return(calculate_PAD_from_PAD_results(results));
	}

