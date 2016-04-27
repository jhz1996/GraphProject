#include<iostream>
#include<string>
#include <vector>
#include <algorithm>
#include<list>
using namespace std;

struct node{
	int numMagi;
	vector<int> s;
	int subSeqSize;
	string charm;

};

class Graph{

public:
	int Dims;
	int incants;
	//node *planets;
	vector<node*> planets;
	int **g;
	//adjacency list


};

//Graph::Graph(int v){
//	this->Dims = v;
//	//this->planets = new list<int>[v]; //array of lists of type int inside them
//	g = new int *[v];
//}

//void Graph::addEdge(int from, int to, int weight){
//
//	//planets.at(from)->push_back(to);
//
//	g[from][to] = weight;
//
//}



int editDistance(const string& s1, const string& s2)
{
	const std::size_t len1 = s1.size(), len2 = s2.size();
	std::vector<std::vector<unsigned int>> d(len1 + 1, std::vector<unsigned int>(len2 + 1));

	if (s1 == s2){
		return 0;

	}

	d[0][0] = 0;
	for (unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
	for (unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

	for (unsigned int i = 1; i <= len1; ++i)
		for (unsigned int j = 1; j <= len2; ++j)
			// note that std::min({arg1, arg2, arg3}) works only in C++11,
			// for C++98 use std::min(std::min(arg1, arg2), arg3)
			d[i][j] = min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
	return d[len1][len2];
}



int CeilIndex(vector<int> &A, int l, int r, int key)
{
	while (r - l > 1)
	{
		int m = l + (r - l) / 2;
		if (A[m] >= key)
			r = m;
		else
			l = m;
	}
	return r;
}

vector<int> LongestIncreasingSubsequenceLength(int A[], int size)
{
	// Add boundary case, when array size is one

	int len; // always points empty slot
	vector<int> temp;
	//memset(tailTable, 0, sizeof(tailTable[0])*size);


	temp.push_back(A[0]);
	len = 1;
	for (int i = 1; i < size; i++)
	{
		if (A[i] < temp.front()){
			// new smallest value

			temp.pop_back();
			temp.push_back(A[i]);
		}
		else if (A[i] > temp.back()){
			//tailTable[len++] = A[i];
			temp.push_back(A[i]);

		}
		// A[i] wants to extend largest subsequence


		else
			// A[i] wants to be current end candidate of an existing
			// subsequence. It will replace ceil value in tailTable
			temp[CeilIndex(temp, -1, len - 1, A[i])] = A[i];
	}


	return temp;
}


int minDistance(int dist[], bool sptSet[], int realms)
{
	// Initialize min value
	int min = INT_MAX, min_index;

	for (int v = 0; v < realms; v++)
		if (sptSet[v] == false && dist[v] <= min)
			min = dist[v], min_index = v;

	return min_index;
}

void printSolution(int dist[], int n)
{
	printf("Vertex   Distance from Source\n");
	for (int i = 0; i < n; i++)
		printf("%d \t\t %d\n", i, dist[i]);
}

void dijkstra(int **graph, int src, int realms)
{
	int *dist = new int[realms];     // The output array.  dist[i] will hold the shortest
	// distance from src to i

	bool *sptSet = new bool[realms]; // sptSet[i] will true if vertex i is included in shortest
	// path tree or shortest distance from src to i is finalized

	// Initialize all distances as INFINITE and stpSet[] as false
	for (int i = 0; i < realms; i++)
		dist[i] = INT_MAX, sptSet[i] = false;

	// Distance of source vertex from itself is always 0
	dist[src] = 0;

	// Find shortest path for all vertices
	for (int count = 0; count < realms - 1; count++)
	{
		// Pick the minimum distance vertex from the set of vertices not
		// yet processed. u is always equal to src in first iteration.
		int u = minDistance(dist, sptSet, realms);

		// Mark the picked vertex as processed
		sptSet[u] = true;

		// Update dist value of the adjacent vertices of the picked vertex.
		for (int v = 0; v < realms; v++)

			// Update dist[v] only if is not in sptSet, there is an edge from 
			// u to v, and total weight of path from src to  v through u is 
			// smaller than current value of dist[v]
			if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
				&& dist[u] + graph[u][v] < dist[v])
				dist[v] = dist[u] + graph[u][v];
	}

	// print the constructed distance array
	printSolution(dist, realms);
}



int main() {
	int numRealms = 0;
	cin >> numRealms;

	Graph *graph = new Graph();

	for (int i = 0; i < numRealms; i++) {
		string charmOfRealm = "";
		cin >> charmOfRealm;
		int numMagi = 0;
		cin >> numMagi;
		int *magiArray = new int[numMagi];
		node*n = new node();

		for (int j = 0; j < numMagi; j++) {
			int powerMagi = 0;
			cin >> powerMagi;
			magiArray[j] = powerMagi;
		}

		n->numMagi = numMagi;
		vector<int> list = LongestIncreasingSubsequenceLength(magiArray, numMagi);
		n->s = list;
		n->subSeqSize = list.size();
		n->charm = charmOfRealm;
		graph->planets.push_back(n);


	}
	graph->g = new int *[numRealms];

	for (int i = 0; i < numRealms; ++i)
		graph->g[i] = new int[numRealms];

	for (int i = 0; i < numRealms; i++){
		//graph.g[i] = new int[numRealms]();
		for (int j = 0; j < numRealms; j++){
			//cout << graph->Dims;
			int s = editDistance(graph->planets.at(i)->charm, graph->planets[j]->charm);
			
			if ((graph->planets[i]->subSeqSize) >s && i != j){

				graph->g[i][j] = s;

			}
			else{
				graph->g[i][j] = 0;

			}


		}

	}

	for (int i = 0; i < numRealms; i++){

		for (int j = 0; j < numRealms; j++)

			cout << graph->g[i][j] << ", ";

		cout << endl;

	}


	dijkstra(graph->g, 0, numRealms);







}