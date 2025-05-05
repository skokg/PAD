#include "CU_kdtree_with_index.cc"

kdtree::Point_str kdtree_Point_str_from_point(const vector <double> &point, const size_t ind)
	{
	kdtree::Point_str p = { {point.begin(), point.end() -1} , ind };
	return(p);
	}

vector <double> perform_one_PAD_iteration(const vector <vector <double>> &points1, const  vector <vector <double>> &points2,  vector <double> &values1, vector <double> &values2,  vector <size_t> &index_list1, vector <size_t> &index_list2, kdtree::KdTree &kdtree1, kdtree::KdTree &kdtree2, long &idum, const bool f1_is_fa)
	{
	// choose the last point from list1
	const auto ind1=index_list1.back();
	const kdtree::Point_str p1 = kdtree_Point_str_from_point(points1[ind1], ind1);

	// find the closes non-zero point in the other field
	//auto begin = std::chrono::high_resolution_clock::now();
	const auto node = kdtree2.findNearestNode(p1);
	//temp.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9);
	ERRORIF(node == nullptr);
	const kdtree::Point_str p2 = node->val;
	const auto ind2=p2.index;

	const double min_squared_distance=squared_euclidian_distance_multidimensional(p1.coords, p2.coords);

	//cout << ind1 << " " << values1[ind1] << endl;
	//cout << ind2 << " " << values2[ind2] << endl;

	vector <double> result;

	// value reduction
	const double value_reduction= min(values1[ind1],values2[ind2]);
	values1[ind1]-=value_reduction;
	values2[ind2]-=value_reduction;

	//cout << list_nzp_f1.size() << " " << list_nzp_f2.size() << endl;

	// if necessary remove the point from kdtree1
	if (values1[ind1] == 0)
		kdtree1.deleteNode(p1);
	// else swap with random point so this is not the last point anymore - so it is not automatically chosen the next time
	else
		swap(index_list1.back(),index_list1[ floor(ran2(&idum) * (double)index_list1.size()) ]);

	// if necessary remove the point from kdtree2
	if (values2[ind2] == 0)
		kdtree2.deleteNode(p2);

	//cout << list_nzp_f1.size() << " " << list_nzp_f2.size() << endl;



	//result.push_back(sqrt(min_squared_distance));
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

	//cout << index_list1.size() << " " << index_list2.size() << endl;

	// remove all zero points at the end of list
	while (index_list1.size() > 0 && values1[index_list1.back()] == 0)
		index_list1.pop_back();

	// remove all zero points at the end of vector
	while (index_list2.size() > 0 && values2[index_list2.back()] == 0)
		index_list2.pop_back();

		//cout << index_list1.size() << " " << index_list2.size() << endl;

	return(result);

	}



// check points for negative values
vector <double> check_points(const vector <vector <double> > &points)
	{
	// check and count non-zero points
	size_t count=0;
	vector <double> values;
	for (unsigned long il=0; il < points.size(); il++)
		{
		double val=points[il].back();
		ERRORIF( val < 0 );
		values.push_back(val);
		if (val > 0)
			{
			count++;
			}
		}

	// check if there are no non-zero points
	ERRORIF(count == 0);

	return(values);
	}


// from non-zero points construct kdtree and shuffled index list
void construct_kdtree_with_shuffled_index_list(const vector <vector <double>> &points,  const vector <double> &values, vector <size_t> &index_list, kdtree::KdTree &kdtree )
	{
	index_list.clear();
	ERRORIF(kdtree.is_empty() == false);

	vector<kdtree::Point_str> pointsX;
	for (size_t ip=0; ip < points.size(); ip++)
		if (values[ip] > 0)
			{
			pointsX.push_back(kdtree_Point_str_from_point(points[ip], ip));
			index_list.push_back(ip);
			}

	// shuffle list
	std::shuffle(index_list.begin(), index_list.end(), rng_state_global);

	// construct kdtree
	kdtree.buildKdTree(pointsX);
	}


vector <vector <double>> calculate_PAD_results_assume_different_grid(const vector <vector <double> > &points1, const vector <vector <double> > &points2, vector <double> &values1, vector <double> &values2)
	{
	auto begin = std::chrono::high_resolution_clock::now();
	// check points for negative values
	values1 = check_points(points1);
	values2 = check_points(points2);
	cout << "----- preprocessing: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9 << " s" << endl;

	//cout << sum_vector(values1) << endl;
	//cout << sum_vector(values2) << endl;

	//cout << "---" << endl;
	//cout << points1.size() << endl;
	//cout << points2.size() << endl;

	//cout << "---" << endl;

	//long temp_counter1=0;
	//long temp_counter2=0;
	//for (unsigned long il=0; il < values1.size(); il++)
	//	if (values1[il] > 0) temp_counter1++;
	//for (unsigned long il=0; il < values1.size(); il++)
	//	if (values2[il] > 0) temp_counter2++;
	//ERRORIF(temp_counter1 == 0 || temp_counter2 == 0);
	//cout << temp_counter1 << endl;
	//cout << temp_counter2 << endl;

	begin = std::chrono::high_resolution_clock::now();
	vector <size_t> index_list1;
	kdtree::KdTree kdtree1;
	construct_kdtree_with_shuffled_index_list(points1, values1, index_list1, kdtree1);
	vector <size_t> index_list2;
	kdtree::KdTree kdtree2;
	construct_kdtree_with_shuffled_index_list(points2, values2, index_list2, kdtree2);
	cout << "----- kdtree construction: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9 << " s" << endl;

	//kdtree1.printKdTree();
	//cout << "---" << endl;
	//cout << index_list1.size() << endl;
	//cout << index_list2.size() << endl;
	//cout << output_vector_as_string(index_list1, " ") << endl;
	//cout << values1[index_list1.back()] << endl;
	//cout << values2[index_list2.back()] << endl;

	vector <vector <double>> out;
	long idum=1;
	bool fa_turn=true;

	// iteratevely remove non-zero points
	begin = std::chrono::high_resolution_clock::now();
	while (index_list1.size() > 0 && index_list2.size() > 0)
		{
		vector <double> result;
		if (fa_turn)
			{
			result = perform_one_PAD_iteration(points1, points2, values1, values2, index_list1, index_list2, kdtree1, kdtree2, idum, fa_turn);
			fa_turn=false;
			}
		else
			{
			result = perform_one_PAD_iteration(points2, points1, values2, values1, index_list2, index_list1, kdtree2, kdtree1, idum, fa_turn);
			fa_turn=true;
			}
		//cout << output_vector_as_string(result, " ") << endl;
		//cout << index_list1.size() << " " << index_list2.size() << " " << output_vector_as_string(result, " ") << endl;

		if (result.size() > 0)
			out.push_back(result);
		}

	// to avoid memory leaks free kdtree memory since not all nodes were necesarrily removed from the tree
	kdtree1.free_memory();
	kdtree2.free_memory();

	cout << "----- attribution: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9 << " s" << endl;

	return(out);
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


