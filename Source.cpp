#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <list>
#include <climits>
#include <queue>
#include <cstring>
using namespace std;

struct node{
	int numMagi;
	vector<int> s; //
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
	bool *path = new bool[Dims];
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
	std::vector<std::vector<unsigned int> > d(len1 + 1, std::vector<unsigned int>(len2 + 1));

	if (s1 == s2){
		return 0;

	}

	d[0][0] = 0;
	for (unsigned int i = 1; i <= len1; ++i){ 
		d[i][0] = i;
	}
	for (unsigned int i = 1; i <= len2; ++i){ 
		d[0][i] = i;
	}

	for (unsigned int i = 1; i <= len1; ++i){
		for (unsigned int j = 1; j <= len2; ++j){
			// note that std::min({arg1, arg2, arg3}) works only in C++11,
			// for C++98 use std::min(std::min(arg1, arg2), arg3)
			d[i][j] = min(d[i - 1][j] + 1, d[i][j - 1] + 1);
			d[i][j] = min( d[i][j] , d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) );
		}		
	}
			
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

	// int len; // always points empty slot
	// vector<int> temp;
	// //memset(tailTable, 0, sizeof(tailTable[0])*size);


	// temp.push_back(A[0]);
	// len = 1;
	// for (int i = 1; i < size; i++)
	// {
	// 	if (A[i] < temp.front()){
	// 		// new smallest value

	// 		temp.pop_back();
	// 		temp.push_back(A[i]);
	// 	}
	// 	else if (A[i] > temp.back()){
	// 		//tailTable[len++] = A[i];
	// 		temp.push_back(A[i]);

	// 	}
	// 	// A[i] wants to extend largest subsequence


	// 	else
	// 		// A[i] wants to be current end candidate of an existing
	// 		// subsequence. It will replace ceil value in tailTable
	// 		temp[CeilIndex(temp, -1, len - 1, A[i])] = A[i];
	// }


	// return temp;
	int totalSubSets = 1;
	std::vector<std::vector<int> > Subset;

	
	for(int i = 0; i < size; i++){
		cout << "i is " << i << endl;
		for(int k = 0; k < totalSubSets; k++){
			cout << "k is " << k << endl;
			cout << Subset.empty() << endl;
			if(Subset.empty()){
				cout << "Is Empty\n";
				Subset.push_back(A[i]);
			}
			else if(A[i] > Subset[k][i-1]){
				cout << "adding" << A[i] << "\n";
				Subset[k].push_back(A[i]); 
			}
			else if(A[i] < Subset[k][i-1]){
				cout << "less than" << A[i] << "\n";
				bool isInOtherSubSet = false;
				for(int h = 0; h < Subset.size(); h++){

					for(int s = 0; s < Subset[h].size(); s++){
						if(Subset[h][s] == A[i]){
							isInOtherSubSet = true;
						}
						
					}

				}
				if(!isInOtherSubSet){
					totalSubSets++;
					k = totalSubSets + 1;
				}
				else{
					//Do nothing, the number is already there
				}
			}
		}
	}
	int max = 0;
	int indexOfMax = 0;
	for(int i = 0; i < totalSubSets; i++){
		if (max < Subset[i].size()){
			max = Subset[i].size();
			indexOfMax = i;
		}
	}
return Subset[indexOfMax];


}


int minDistance(int dist[], int realms, bool visited[])
{
	// Initialize min value
	int min = INT_MAX;
	int min_index = -1; 

	for(int i = 0; i < realms; i++){
		if(dist[i] < min && dist[i] > 0 && visited[i] == false){ //might error if step is actually 0
			min = dist[i];
			min_index = i;
		}
	}
	// for (int v = 0; v < realms; v++)
	// 	if (sptSet[v] == false && dist[v] <= min){
	// 		min = dist[v]; 
	// 		min_index = v;
	// 	}



	return min_index;

	// using u, finds the smallest one that has not been visited
}

void printSolution(int dist[], int n)
{
	// cout << "Vertex   Distance from Source\n";
	for (int i = 0; i < n; i++){
		// printf("%d \t\t %d\n", i, dist[i]);
		// cout << i << " \t\t " << dist[i] << "\n";
	}
}

bool allvisited(bool visited[], int size){
	// cout << "ALL visited \n";
	for(int i = 0; i < size; i++){
		if(visited[i] == false){
			return false;
		}
	}
	return true;
}

std::vector<int> dijkstra(int **graph, int src, int realms, int end, int *&distoutput, int *&prevoutput)
{
	int *dist = new int[realms];     // The output array.  dist[i] will hold the shortest
	// distance from src to i
	int *prev = new int[realms];
	//int* path = new int[realms];
	bool *visited = new bool [realms];

	//bool *sptSet = new bool[realms]; // sptSet[i] will true if vertex i is included in shortest
	// path tree or shortest distance from src to i is finalized
	//std::vector<int> Q;
	// Initialize all distances as INFINITE and stpSet[] as false
	// cout << "Initializer \n";
	for (int i = 0; i < realms; i++){
		dist[i] = INT_MAX;
		prev[i] = -1; 
		visited[i] = false;
	}	//sptSet[i] = false;
		//Q.push_back(i);
	// Distance of source vertex from itself is always 0
	dist[src] = 0;
	prev[src] = -1;
	visited[src] = true;


	int u = src;
	//Q.erase(src);
	//order will not be there for Q
	bool first = true;
	while(!allvisited(visited, realms) ){
		if(first == true){
			// cout << "u is " << u << endl;
		}
		else{
			// cout << "minDistance to " << u <<" \n";

			cout << "dist: ";
			for(int q = 0; q < realms; q++){
				cout << dist[q] << " ";
			}
			cout << "prev: ";
			for(int q = 0; q < realms; q++){
				cout << prev[q] << " ";
			}
			cout << endl;
			u = minDistance(dist, realms, visited);
			// cout << "Traveled to " << u <<" \n";
			if(u == -1){
				//no more to be found or some are unreachable
				// cout << "u is -1 \n";

				break;
			}
			visited[u] = true;
		}
		// cout << "Checking for faster routes \n";
		for(int k = 0; k < realms; k++){
			//cout << "k is " << k << "\n";
			//cout << "graph[u][k] is " << graph[u][k] << endl;
			if(k != u && graph[u][k] != -1){ //Don't travel to your selected node
				int alt = dist[u] + graph[u][k];
				if(alt < dist[k] ){
					// cout << "made " << k << " have a distance of " << alt << endl;
					// cout << "it's  previous is now " << u << endl;
					dist[k] = alt;
					prev[k] = u;
				}



			}
		}
		if(first == true){
			first = false;
		}


	}

	// // Find shortest path for all vertices
	// for (int count = 0; count < realms - 1; count++)
	// {
	// 	// Pick the minimum distance vertex from the set of vertices not
	// 	// yet processed. u is always equal to src in first iteration.
	// 	int u = minDistance(dist, sptSet, realms);

	// 	// Mark the picked vertex as processed
	// 	sptSet[u] = true;

	// 	// Update dist value of the adjacent vertices of the picked vertex.
	// 	for (int v = 0; v < realms; v++)

	// 		// Update dist[v] only if is not in sptSet, there is an edge from 
	// 		// u to v, and total weight of path from src to  v through u is 
	// 		// smaller than current value of dist[v]
	// 		if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
	// 			&& dist[u] + graph[u][v] < dist[v]){
	// 			dist[v] = dist[u] + graph[u][v];
	// 			//path[v] = u;
	// 		}
	// }
	std::vector<int> rpath;
	int current = end;
	rpath.push_back(current);
	// cout << "pushing " << current << endl;
	while(current != -1){//tracing the path
		rpath.push_back(prev[current]);
		// cout << "pushing " << prev[current] << endl;
		current = prev[current];

	}
	//reverse the order because it's backwards
	int size = rpath.size();
	int routesize;
	std::vector<int> path;
	//int path[size -1] = new int[size - 1];
	//memset( path, 0, (size - 1)*sizeof(int) );
	int k = size -2;
	for(int h = 0; h < size -1; h++){
		path.push_back(rpath[k]);

		k--;
	} 
	cout << "Path: ";
	for(int h = 0; h < size - 1; h++){
		cout << path[h] << " ";
	}
	cout << endl;
	routesize = size -1;
	//*route = path;

	distoutput = dist;
	prevoutput = prev;


	// print the constructed distance array
	//printSolution(dist, realms);
	return path;
}



int main() {
	int numRealms = 0;
	cin >> numRealms;

	Graph *graph = new Graph();
	graph->Dims = numRealms;
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
		n->s = list; //list of the longest increasing subsequence
		n->subSeqSize = list.size();
		n->charm = charmOfRealm;
		graph->planets.push_back(n);//push the planet nodes into the greater graph


	}

	string start;
	int startNode;

	cin >> start;

	string end;
	int endNode;
	cin >> end;

	for (int i = 0; i < numRealms; i++){

		if (graph->planets[i]->charm == start){
			startNode = i;

		}

		else if (graph->planets[i]->charm == end){
			endNode = i;

		}

	}

	graph->g = new int *[numRealms];

	for (int i = 0; i < numRealms; ++i)
		graph->g[i] = new int[numRealms];

	for (int i = 0; i < numRealms; i++){
		//graph.g[i] = new int[numRealms]();
		// cout << "for i = " << i << endl;
		for (int j = 0; j < numRealms; j++){
			//cout << graph->Dims;
			// cout << "for j = " << j << endl;
			int s = editDistance(graph->planets.at(i)->charm, graph->planets[j]->charm);
			// cout << "the difference between the charms " << s << endl;
			// cout << "subSeqSize is for " << i << " is " << graph->planets[i]->subSeqSize << endl;
			if ((graph->planets[i]->subSeqSize) >=s && i != j){
				// cout << "takes " << s << " to get from " << i << " to " << j << endl;
				graph->g[i][j] = s;

			}
			else{
				graph->g[i][j] = -1; //-1 denotes it is visiting its own node or there is no path that exists
				// cout << "its own or doesn not exist\n";
			}


		}

	}

	for (int i = 0; i < numRealms; i++){

		for (int j = 0; j < numRealms; j++){

			cout << graph->g[i][j] << " ";

		}	
		cout << endl;

	}



	//graph->path=dijkstra(graph->g, startNode, numRealms, endNode);//needs to itereate through an array
	//int *route;
	//int *route = new std::vector<int>;
	int *s2fdist = new int [numRealms]; //start to finish with the optimal distances
	int *s2fprev = new int [numRealms]; //start to finish with the previous node

	std::vector<int> route = dijkstra(graph->g, startNode, numRealms, endNode, s2fdist, s2fprev);



	cout << "Printing the Subsets \n";
	for(int q = 0; q < numRealms; q++){
		cout << "For realm " << q << endl;
		for(int y = 0; y < graph->planets[q]->s.size(); y++){
			cout << graph->planets[q]->s[y] << " ";
		}
		cout << endl;
	}
	//int &rs = routesize;
	// cout << "Printing out route: ";
	if(route.empty()){
		cout << "IMPOSSIBLE \n"; //TODO, fix if the doesn't end up to start
	}
	else{
		for (int i = 0; i < route.size(); i++){

			 cout << route[i] << " ";

		}

		cout << endl;
		cout << "Incantations: " << s2fdist[endNode] << endl;
		int s2fgems = 0;
		for(int k = 0; k < route.size() - 1; k++){
			int depart = route[k];
			int destination = route[k+1];
			int incantations = graph->g[depart][destination];
			int localgems = 0;
			for(int p = 0; p < incantations; p++){
				cout << "adding to local: " << graph->planets[depart]->s[p] << endl;
				localgems = localgems + graph->planets[depart]->s[p];
			}
			cout << "adding to global: " << localgems << endl;
			s2fgems = s2fgems + localgems;
			cout << "totalgems: " << s2fgems << endl;
		}

		cout << "Gems: " << s2fgems << endl;
	}
	int *f2sdist = new int [numRealms]; //start to finish with the optimal distances
	int *f2sprev = new int [numRealms]; //start to finish with the previous node

	std::vector<int> returnroute = dijkstra(graph->g, endNode, numRealms, startNode, f2sdist, f2sprev);
	
	if(!route.empty() && returnroute.empty()){//TODO CHANGE
		cout << "IMPOSSIBLE \n";
	}
	else{
		for (int i = 0; i < returnroute.size(); i++){

			cout << returnroute[i] << " ";

		}

		cout << endl;
		cout << "Incantations: " << f2sdist[startNode] << endl;
		int f2sgems = 0;
		for(int k = 0; k < route.size() - 2; k++){
			int depart = returnroute[k];
			int destination = returnroute[k+1];
			int incantations = graph->g[depart][destination];
			int localgems = 0;
			for(int p = 0; p < incantations; p++){

				localgems = localgems + graph->planets[depart]->s[p];
			}
			f2sgems = f2sgems + localgems;
		}

		cout << "Gems: " << f2sgems << endl;
	}


	//system("pause");





}