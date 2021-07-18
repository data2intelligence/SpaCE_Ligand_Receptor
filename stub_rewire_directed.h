/*
 * Author: Peng Jiang (pengj@alumni.princeton.edu)
 * Stub-rewiring network randomization. All edges may not be connected, but this procedure is always fast even for dense network, which will cause slow down for edge-swapping.
 * */

#ifndef STUB_REWIRE_DIRECTED_H_
#define STUB_REWIRE_DIRECTED_H_

#include <map>
#include <vector>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstring>
using namespace std;


class stub_rewire_directed
{
public:
	// input a list of directed edge with vertex indices
	stub_rewire_directed(const vector<pair<size_t, size_t> > &network);
	~stub_rewire_directed();

	// generate one randomized network
	double randomize();

	// output current network (in adjacency list) to vector of edge pairs
	size_t output(vector<pair<size_t, size_t> > &network) const;

private:
	vector<size_t> degree_vec_a, degree_vec_b, index_vec_a, index_vec_b;

	size_t Nvertex_a, Nvertex_b, Nedge,
		**adjacency_list_a, *adjacency_list_edge_index_a,
		*degree_cur_a, *degree_cur_b, *candidate_a, *candidate_b;

	bool *visited;

};

#endif /* STUB_REWIRE_DIRECTED_H_ */
