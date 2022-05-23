#pragma once
#include <stdlib.h>
#include <fstream>
#include <list>
#include <vector>
#include <queue>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include "edge.h"
#include "node.h"


static int flag{ 0 };
static int MAXLISTSIZE;	//the length of inverted list
static float GAMMA;//the parameter gamma
static int BLOCKSIZE; // the size of blocks used in WH2H
static int MAXDEGREE;//the maximum degree of the graph

using namespace std;

/*
* The graph class is used to remove vertices by degree (H2H) or value (WH2H) while
* saving the weight of all the edges during processing. We create an inverted list
* to hasten the procedure of minimunDegree and record the iterator of all the vertices
* for the Convenience of updating.
*/
class graph
{
public:
	//construct function
	graph(const char* graph_file, int maxListSize, float gamma, int blocksize, int maxDegree, bool directed);

	//read the query frequency of high frequency vertices
	void readHFVertexData(const char* HFVertex_file);

	//Add a vertex with degree in indirected graph to the vertexList, and update vertexListIterator, degreePointer, and vertexIterator.
	void addVertexToVertexList(const int& nodeNum, const int degree);

	//Add a vertex with degree in directed graph to the vertexList, and update vertexListIterator, degreePointer, and vertexIterator.
	void addVertexToVertexList_directed(const int& nodeNum, const int degree);

	//Add a vertex with value in indirected graph to the vertexList, and update vertexListIterator, degreePointer, and vertexIterator.
	void addVertexToVertexList_HFVertices(const int& nodeNum, const int value);

	//Add a vertex with value in directed graph to the vertexList, and update vertexListIterator, degreePointer, and vertexIterator.
	void addVertexToVertexList_HFVertices_directed(const int& nodeNum, const int value);

	//Used for read data from indirected graph£¬add a vertex with degree to vertexList.
	void addVertexToVertexList_readData(const int& nodeNum, const int& degree);

	//Used for read data from directed graph£¬add a vertex with degree to vertexList.
	void addVertexToVertexList_readData_directed(const int& nodeNum, const int& degree);

	//Used for read graph data in WH2H£¬ add a vertex with degree to vertexList.
	void addVertexToVertexList_readData_HFVertices(const int& nodeNum, const int& value);

	//Used for read graph data in WH2H£¬ add a vertex with degree to vertexList.
	void addVertexToVertexList_readData_HFVertices_directed(const int& nodeNum, const int& value);

	//Remove a vertex with degree one in indirected graph.
	void removeDegreeOne(const int& number);

	//Remove a vertex with degree one in directed graph.
	void removeDegreeOne_directed(const int& number);

	//Remove a vertex with degree one in indirected graph.
	void removeDegreeOne_HFVertices(const int& number);

	//Remove a vertex with degree one in directed graph.
	void removeDegreeOne_HFVertices_directed(const int& number);

	//Remove a vertex by degree (>=2) in indirected graph.
	void removeByDegree(const int& number, const int& degree);

	//Remove a vertex by degree (>=2) in directed graph. Here degree means number of neighbors.
	void removeByDegree_directed(const int& number, const int& degree);

	//Remove a vertex by value (>=2) in indirected graph.
	void removeByValue(const int& number, const int& degree);

	//Remove a vertex by value (>=2) in directed graph.
	void removeByValue_directed(const int& number, const int& degree);

	//Calculate the value of a vertex in indirected graph.
	inline int calValue(const int& num);

	//Calculate the value of a vertex in directed graph.
	inline int calValue_directed(const int& num);

	//Calculate degree for directed graph.
	inline int calDegree(const int& num);

	//Remove vertices in indirected graph by H2H.
	void removeVertices_H2H();

	//Remove vertices in directed graph by H2H .
	void removeVertices_H2H_directed();

	//Remove vertices in indirected graph by WH2H.
	void removeVertices_WH2H();

	//Remove vertices in directed graph by WH2H.
	void removeVertices_WH2H_directed();

	//Update edges in directed graphs. 
	void updateEdges_out(const int& u, const int& outNum, const int& number, const int& calNum, const vector<edge>& edgeProcess);
	void updateEdges_in(const int& v, const int& inNum, const int& number, const int& calNum, const vector<edge>& edgeProcess);

	//test the correctness of queries
	int dijkstra(int x, int y);

	vector<node> nodes;//Used for saving nodes in tree decomposition.
	vector<int> level;//Record the removal sequence of vertices. 
	vector<int> allHFVertex;//Save high frequency vertices by query frequency.
	int numOfVertices; //number of vertices.
	int numOfEdges; //number of edges.
	vector<vector<int>> dis_dijkstra;

private:
	int sequence{ 0 };//Record the sequence of current vertex for processing.
	vector<pair<int, int>> HFVertex;//The first elements is high frequency vertex number, the second one is query frequency. 
	vector<list<edge>> es;//Save the indirected graph data, i.e., vertices and edges.
	vector<pair<list<edge>,list<edge>>> es_directed;//Save the directed graph data, i.e., vertices and edges. First is out-edges, and second is in-edges.
	vector<pair<list<edge>, list<edge>>> es_dijkstra;//Save the directed graph data, i.e., vertices and edges. First is out-edges, and second is in-edges.
	vector<bool> isRemoved;//Record whether a vertex is removed.
	list<list<pair<int, int>>*> vertexList;//inverted list, the first dimension is degree or value, the second dimension is vertex (ID, value).
	vector<list<list<pair<int, int>>*>::iterator*> vertexListIterator;//The iterator of inverted list saved by the first dimension (ID, value).
	vector<list<pair<int, int>>::iterator*> vertexIterator;//The iterator of vertices (ID, value).
	vector<bool> isAltered;//Record whether the degree or value of neighbors have been altered.
	queue<pair<int, int>> changedNeighbors; //Record vertices of which the degree or value has been altered. (vertex number, value or degree)
	vector<pair<int, int>> degree; //only used when processing a directed graph, degree.first is in-degree, degree.second is out-degree
	//vector<list<int>> neighbors; //only used when processing a directed graph, record neighbors of each vertex
	
};