#include "stub_rewire_directed.h"


stub_rewire_directed::stub_rewire_directed(const vector<pair<size_t, size_t> > &interactions)
: Nedge(interactions.size())
{
	size_t i, j;
	vector<pair<size_t, size_t> >::const_iterator iter;

	map<size_t, size_t> inxmap_a, inxmap_b;
	map<size_t, size_t>::iterator iter_a, iter_b;

	for(iter = interactions.begin(); iter != interactions.end(); iter++)
	{
		i = iter->first;
		j = iter->second;

		iter_a = inxmap_a.find(i);

		if(iter_a == inxmap_a.end())
		{
			iter_a = inxmap_a.insert(pair<size_t,size_t>(i, index_vec_a.size())).first;
			index_vec_a.push_back(i);
			degree_vec_a.push_back(0);
		}

		iter_b = inxmap_b.find(j);

		if(iter_b == inxmap_b.end())
		{
			iter_b = inxmap_b.insert(pair<size_t,size_t>(j, index_vec_b.size())).first;
			index_vec_b.push_back(j);
			degree_vec_b.push_back(0);
		}

		degree_vec_a[iter_a->second]++;
		degree_vec_b[iter_b->second]++;
	}

	Nvertex_a = degree_vec_a.size();
	Nvertex_b = degree_vec_b.size();

	candidate_a = new size_t[Nvertex_a];
	candidate_b = new size_t[Nvertex_b];

	degree_cur_a = new size_t[Nvertex_a];
	degree_cur_b = new size_t[Nvertex_b];

	adjacency_list_a = new size_t*[Nvertex_a];
	adjacency_list_edge_index_a = new size_t[Nedge];

	visited = new bool[Nvertex_b];

	for (i=j=0 ; i<Nvertex_a ; i++)
	{
		candidate_a[i] = i;
		adjacency_list_a[i] = adjacency_list_edge_index_a + j;
		j += degree_vec_a[i];
	}
}



stub_rewire_directed::~stub_rewire_directed()
{
	delete[] adjacency_list_a;
	delete[] adjacency_list_edge_index_a;
	delete[] degree_cur_a;
	delete[] degree_cur_b;
	delete[] candidate_a;
	delete[] candidate_b;
	delete[] visited;
}

double stub_rewire_directed::randomize()
{
	size_t i,j, inxa, inxb, offset, count, *neighbor_a, Ncandidate = Nvertex_b, connected = 0;

	// current degree of vertex b
	memset(degree_cur_a, 0, Nvertex_a*sizeof(size_t));
	memset(degree_cur_b, 0, Nvertex_b*sizeof(size_t));

	// visit flag of vertex b
	memset(visited, 0, Nvertex_b*sizeof(bool));

	// random order of vertex a in connecting
	random_shuffle(candidate_a, candidate_a + Nvertex_a);

	for(i=0;i<Nvertex_b;i++) candidate_b[i] = i;

	for(i=0;i<Nvertex_a;i++)
	{
		inxa = candidate_a[i];
		neighbor_a = adjacency_list_a[inxa];

		for(j=0; j<degree_vec_a[inxa]; j++)
		{
			// circular searching for unconnected vertex
			for(offset = rand() % Ncandidate, count=0; count < Ncandidate; offset = (offset+1) % Ncandidate, count ++)
			{
				inxb = candidate_b[offset];

				// not visited and not self interaction
				if(!visited[inxb] && index_vec_a[inxa] != index_vec_b[inxb]) break;
			}

			// cannot find any b to connect
			if (count == Ncandidate) break;

			degree_cur_a[inxa]++;
			degree_cur_b[inxb]++;

			neighbor_a[j] = inxb;
			visited[inxb] = true;
			connected++;

			// if vertex b is full, remove it from candidate list
			if(degree_cur_b[inxb] == degree_vec_b[inxb]) swap(candidate_b[offset], candidate_b[--Ncandidate]);
		}

		for(j=0; j<degree_cur_a[inxa]; j++) visited[neighbor_a[j]] = false;
	}

	return static_cast<double>(connected)/Nedge;
}



size_t stub_rewire_directed::output(vector<pair<size_t, size_t> > &network) const
{
	size_t i, j, inxa, inxb, N_edge=0;

	for(i=0;i<Nvertex_a;i++)
	{
		inxa = index_vec_a[i];

		for(j=0; j<degree_cur_a[i]; j++)
		{
			inxb = index_vec_b[adjacency_list_a[i][j]];

			// ignore all edges if the network is already full
			if(N_edge >= network.size())
			{
				cerr << "Warning: network size is already full" << endl;
				continue;
			}

			network[N_edge].first = inxa;
			network[N_edge].second = inxb;
			N_edge++;
		}
	}

	return N_edge;
}
