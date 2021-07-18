#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <climits>
#include <map>
#include <set>
#include <cmath>
using namespace std;

#include "stub_rewire_directed.h"


#define check(flag,msg) if(flag){cerr << msg << endl;exit(1);}
#define PROGRESS_STEP 100
#define EPS 1e-8


/*
 * Compute network sum
 */

void network_sum(
		const map<size_t, double*> &data_map,
		const vector<pair<size_t, size_t> > &network,
		const size_t N_edge, // the network may not have full size, after stub-rewiring

		double score[], const size_t N_sample)
{
	size_t i, j;

	const double *arr_a;
	const double *arr_b;

	memset(score, 0, N_sample * sizeof(double));

	vector<pair<size_t, size_t> >::const_iterator iter = network.begin();

	for(i=0; i<N_edge; i++, iter++)
	{
		arr_a = data_map.at(iter->first);
		arr_b = data_map.at(iter->second);

		for(j=0; j<N_sample; j++) score[j] += arr_a[j] * arr_b[j];
	}

	for(j=0; j<N_sample; j++) score[j] /= N_edge;
}


/*
 * Compute network score, using stub re-wiring (will always be fast even for dense network).
 */

void network_score(
		const map<size_t, double*> &data_map,
		vector<pair<size_t, size_t> > &network,

		const size_t N_rand,
		double score[], double pvalue[], const size_t N_sample,
		const size_t N_retry=10, const double ratio_thres=0.95)
{
	size_t i, j, N_edge, step = N_rand/PROGRESS_STEP;

	double ratio;
	double *score_rand = new double[N_sample], *average_rand = new double[N_sample];
	stub_rewire_directed *stub_rewirer = new stub_rewire_directed(network);

	memset(pvalue, 0, N_sample * sizeof(double));
	memset(average_rand, 0, N_sample * sizeof(double));

	// compute original network sum
	network_sum(data_map, network, network.size(), score, N_sample);

	for(i=0; i<N_rand; i++)
	{
		// print out a progress %
		if(i % step == 0){ cout << 100 * i / N_rand << '%' << endl;}

		// stub-rewiring randomization
		for(j=0; j<N_retry; j++)
		{
			ratio = stub_rewirer->randomize();
			if(ratio >= ratio_thres) break;
		}

		if(ratio < ratio_thres){
			cerr << "Warning, maximal reconnection ratio is " << ratio << endl;
		}

		// dump the network to edge vector
		N_edge = stub_rewirer->output(network);

		network_sum(data_map, network, N_edge, score_rand, N_sample);

		for(j=0; j<N_sample; j++)
		{
			// one-sided p-value
			if(score_rand[j] >= score[j]) pvalue[j] += 1;

			// average of network sum
			average_rand[j] += score_rand[j];
		}
	}

	for(j=0; j<N_sample; j++)
	{
		pvalue[j] /= N_rand;
		average_rand[j] /= N_rand;

		// smoothed ratio
		score[j] = (score[j] + EPS) / (average_rand[j] + EPS);
	}

	delete[] score_rand;
	delete[] average_rand;
	delete stub_rewirer;
}




int main(int argc, char *argv[])
{
	// gene name map
	vector<string>
		gene_names,	// gene name array, with index as the gene ID
		sample_names	// sample names in the ST matrix, i.e., the spot location
		;

	map<string, size_t> name_map;
	map<string, size_t>::iterator niter;

	// matrix
	size_t index;
	map<size_t, double*> data_map;
	map<size_t, double*>::iterator diter;

	// interaction, a -> b
	vector<pair<size_t, size_t> > network;

	// set to remove duplicated edges in the network
	set<uint> included;

	size_t i, j, N_sample, N_rand, key;

	double *data, *score, *pvalue;

	string line;
	ostringstream oss;
	istringstream iss;
	ifstream fin;

	check(argc != 5, "SpaCE_Ligand_Receptor matrix network output nrand")

	N_rand = atoi(argv[4]);

	// load matrix
	fin.open(argv[1]);
	check(fin.fail(), "Cannot open " << argv[1])

	getline(fin, line, '\n');
	iss.str(line);
	while(getline(iss, line, ' ')) sample_names.push_back(line);
	iss.clear();

	N_sample = sample_names.size();

	while(getline(fin, line, '\n'))
	{
		iss.str(line);
		getline(iss, line, ' ');

		niter = name_map.find(line);

		// ignore duplicated names
		if(niter != name_map.end()){
			cerr << "Jump duplicated gene name " << line << endl;
			continue;
		}

		// add a new gene name
		name_map[line] = index = gene_names.size();
		gene_names.push_back(line);

		data_map[index] = data = new double[N_sample];
		for(i=0;i<N_sample;i++)
		{
			iss >> data[i];
			check(iss.fail(), "Error: Fail reading number on line ID: " << line)
		}

		iss.clear();
	}

	fin.close();
	fin.clear();


	// load network
	fin.open(argv[2]);
	check(fin.fail(), "Cannot open " << argv[1])

	getline(fin, line, '\r');

	// read the header
	iss.str(line);
	while(getline(iss, line, '\t'));
	iss.clear();

	while(getline(fin, line, '\r'))
	{
		iss.str(line);

		// set i, j as empty first
		i = j = string::npos;

		getline(iss, line, '\t');
		getline(iss, line, '\t');

		niter = name_map.find(line);
		if(niter != name_map.end()) i = niter->second;

		getline(iss, line, '\t');
		getline(iss, line, '\t');

		niter = name_map.find(line);
		if(niter != name_map.end()) j = niter->second;

		iss.clear();

		// only keep an edge when both genes are identified in the expression matrix
		if(i== string::npos || j== string::npos) continue;

		key = i * gene_names.size() + j;

		// ignore duplicated edge
		if(included.find(key) != included.end())
		{
			cerr << "Jump duplicated edge\t" << gene_names[i] << '\t' << gene_names[j] << endl;
			continue;
		}

		// push in the edge
		network.push_back(pair<uint, uint>(i, j));

		included.insert(key);
	}

	fin.close();
	fin.clear();


	// Start the network score computation
	score = new double[N_sample];
	pvalue = new double[N_sample];

	network_score(data_map, network, N_rand, score, pvalue, N_sample);

	ofstream fout(argv[3]);
	fout << "Sample\tScore\tPvalue\n";

	for(i=0;i<N_sample;i++)
	{
		fout << sample_names[i] << '\t' << score[i] << '\t' << pvalue[i] << '\n';
	}

	fout.close();

	// free up allocated space
	for(diter=data_map.begin(); diter!=data_map.end(); diter++) delete[] diter->second;

	delete[] score;
	delete[] pvalue;

	return 0;
}
