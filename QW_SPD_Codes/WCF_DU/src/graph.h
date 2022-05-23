/*
 * @Descripttion: define the graph class and util functions of graph
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-09 15:37:10
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-10-16 15:04:38
 */

/*
This file referenced:
An Experimental Study on Hub Labeling based Shortest Path Algorithms [Experiments and Analyses]

Authors: Ye Li, Leong Hou U, Man Lung Yiu, Ngai Meng Kou
Contact: yb47438@umac.mo
Affiliation: University of Macau

The MIT License (MIT)

Copyright (c) 2016 University of Macau

All rights reserved. 
*/
#pragma once
#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <fstream>
#include <algorithm>
#include <cstdint> 
#include <iostream>
#include <assert.h>
#include <map>
#include <list>
#include "paras.h"

using namespace std;

//typedef unsigned int NodeID;
typedef int NodeID;
typedef long long EdgeID;
//used for h2h construction
struct edge {
	int u;
	int v;
	int weight;
};
//typedef unsigned int EdgeID;
typedef pair<NodeID, NodeID> Edge;
//typedef double EdgeWeight;
typedef  int EdgeWeight;
typedef pair<NodeID, EdgeWeight> NodeEdgeWeightPair;

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define INF_WEIGHT SP_Constants::INF_WEIGHT

// Unweighted Graph.
class Graph {
public:
	

	vector<vector<NodeID> > adj; // Adjacent Lists.
	vector<vector<NodeID> > r_adj; // Reverse adjacent lists for directed search.
	vector<EdgeID> vertices;
	vector<EdgeID> r_vertices;
	vector<NodeID> edges;
	vector<NodeID> r_edges;


	Graph() {
		numOfVertices = 0;
		numOfEdges = 0;
		adj.clear();
		r_adj.clear();
	}

	~Graph() {
		adj.clear();
		r_adj.clear();
	}

	// Load Graph file with format:
	// u	v
	// node id in the graph file should start with 0.
	// The graph file should imply whether it is directed. If it is directed, we will duplicaate (u ,v) with (v, u).
	bool load_graph(const char* graph_file) {
		
		ifstream ifs(graph_file);

		vector<Edge> es;


		for (NodeID u, v; ifs >> u >> v;) {
			numOfVertices = max(numOfVertices, max(u, v) + 1);
			if (u == v) continue; // Delete self loop.
			es.push_back(make_pair(u, v));
		}
		numOfEdges = es.size();
		if (DIRECTED_FLAG == true)
			numOfEdges += numOfEdges;

		if (ifs.bad()) return false;
		ifs.close();

		adj.resize(numOfVertices);
		if (DIRECTED_FLAG == true) 	r_adj.resize(numOfVertices);
		

		for (size_t i = 0; i < es.size(); ++i) {
			NodeID u = es[i].first, v = es[i].second;
			adj[u].push_back(v);
			if (DIRECTED_FLAG == false)
				adj[v].push_back(u);
			else {
				r_adj[v].push_back(u);
			}
		}

		for (NodeID v = 0; v < numOfVertices; ++v) {
			sort(adj[v].begin(), adj[v].end());
			if(DIRECTED_FLAG == true)
				sort(r_adj[v].begin(), r_adj[v].end());
		}

		long long sum_of_edges = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			sum_of_edges += adj[v].size();
		}

		vertices.resize(numOfVertices + 1);
		edges.reserve(sum_of_edges);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			vertices[v] = edges.size();

			for (NodeID i = 0; i < adj[v].size(); ++i)
				edges.push_back(adj[v][i]);
		}

		vertices[numOfVertices] = edges.size();

		if (DIRECTED_FLAG == true) {
			sum_of_edges = 0;
			for (NodeID v = 0; v < numOfVertices; ++v) {
				sum_of_edges += r_adj[v].size();
			}
			r_vertices.resize(numOfVertices + 1);
			r_edges.reserve(sum_of_edges);

			for (NodeID v = 0; v < numOfVertices; ++v) {
				r_vertices[v] = r_edges.size();

				for (NodeID i = 0; i < r_adj[v].size(); ++i)
					r_edges.push_back(r_adj[v][i]);
			}
			r_vertices[numOfVertices] = r_edges.size();
		}

		adj.clear();
		r_adj.clear();
		adj.shrink_to_fit();
		r_adj.shrink_to_fit();

		return true;
	}

};

// Weighted Graph.
class WGraph : public Graph {

public:
	vector< vector<EdgeWeight> > adj_weight; // Weights of adjacent lists respect to adj which is ordered by ids
	vector< vector<EdgeWeight> > r_adj_weight; // Weights of reverse adjacent lists;
	vector<NodeEdgeWeightPair> edges;
	vector<NodeEdgeWeightPair> r_edges;


	WGraph() {
		numOfVertices = 0;
		numOfEdges = 0;
		adj.clear();
		adj_weight.clear();
		r_adj.clear();
		r_adj_weight.clear();
	}

	~WGraph() {
		adj.clear();
		adj_weight.clear();
		r_adj.clear();
		r_adj_weight.clear();
	}

	bool load_graph(const vector<list<edge> >& es_list,const vector<int>& originalToNew,int numOfOriginalVertices,int graph_type=0){
		if(numOfVertices!=0){
			numOfVertices=0;
			numOfEdges=0;
		}
		vector<Edge> es;
		vector<EdgeWeight> es_weight;
		NodeID u, v;
		EdgeWeight w;
		for(int i=0;i<numOfOriginalVertices;++i){
			if (es_list[i].size() != 0) {
				for (list<edge>::const_iterator it = es_list[i].begin(); it != es_list[i].end(); it++) {
					if (originalToNew[it->u] != -1 && originalToNew[it->v] != -1){
						u=originalToNew[it->u];
						v=originalToNew[it->v];
						w=it->weight;
						numOfVertices = max(numOfVertices, max(u, v) + 1);
						if (u == v) continue; // Delete self loop.
						es.push_back(make_pair(u, v));
						es_weight.push_back(w);
						//cout<<u<<"-"<<it->u<<" "<<v<<"-"<<it->v<<" "<<w<<endl;//to be deleted
					}
				}
			}
		}
		cout<<"numOfVertices = "<<numOfVertices<<endl;
		
		adj.resize(numOfVertices);
		adj_weight.resize(numOfVertices);
		assert(adj.size() == numOfVertices && adj_weight.size() == numOfVertices && "No enough memory for adj lists!");
		if (DIRECTED_FLAG == true) {
			r_adj.resize(numOfVertices);
			r_adj_weight.resize(numOfVertices);
			assert(r_adj.size() == numOfVertices && r_adj_weight.size() == numOfVertices  && "No enough memory for adj lists!");
		}

		for (size_t i = 0; i < es.size(); ++i) {
			NodeID u = es[i].first, v = es[i].second;
			EdgeWeight w = es_weight[i];
			adj[u].push_back(v);
			adj_weight[u].push_back(w);
			if (DIRECTED_FLAG == false) {
				adj[v].push_back(u);
				adj_weight[v].push_back(w);
			}
			else {
				r_adj[v].push_back(u);
				r_adj_weight[v].push_back(w);
			}
		}

		for (NodeID v = 0; v < numOfVertices; ++v) {
			vector<NodeID>& adj_v = adj[v];
			vector<EdgeWeight>& adj_weight_v = adj_weight[v];
			vector<pair<NodeID, EdgeWeight> > adj_vw_v(adj_v.size());
			for (size_t i = 0; i < adj_v.size(); ++i) {
				NodeID u = adj_v[i];
				EdgeWeight w = adj_weight_v[i];
				adj_vw_v[i] = make_pair(u, w);
			}
			sort(adj_vw_v.begin(), adj_vw_v.end());
			for (size_t i = 0; i < adj_v.size(); ++i) {
				adj_v[i] = adj_vw_v[i].first;
				adj_weight_v[i] = adj_vw_v[i].second;
			}
			adj_vw_v.clear();

			if (DIRECTED_FLAG == true) {
				vector<NodeID>& r_adj_v = r_adj[v];
				vector<EdgeWeight>& r_adj_weight_v = r_adj_weight[v];
				vector<pair<NodeID, EdgeWeight> > r_adj_vw_v(r_adj_v.size());
				for (size_t i = 0; i < r_adj_v.size(); ++i) {
					NodeID u = r_adj_v[i];
					EdgeWeight w = r_adj_weight_v[i];
					r_adj_vw_v[i] = make_pair(u, w);
				}
				sort(r_adj_vw_v.begin(), r_adj_vw_v.end());
				for (size_t i = 0; i < r_adj_v.size(); ++i) {
					r_adj_v[i] = r_adj_vw_v[i].first;
					r_adj_weight_v[i] = r_adj_vw_v[i].second;
				}
				r_adj_vw_v.clear();
			}
		}

		long long sum_of_edges = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			sum_of_edges += adj[v].size();
		}

		vertices.resize(numOfVertices + 1);
		edges.reserve(sum_of_edges);
		EdgeID duplicated = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			vertices[v] = edges.size();

			NodeID last_v = numOfVertices;
			EdgeWeight last_weight = -1;
			for (NodeID i = 0; i < adj[v].size(); ++i) {
				if (adj[v][i] == last_v) {
					if (last_weight <= adj_weight[v][i]) {
						last_v = adj[v][i];
						last_weight = adj_weight[v][i];
						duplicated++;
						continue;
					}
				}
				edges.push_back(NodeEdgeWeightPair(adj[v][i], adj_weight[v][i]));
				last_v = adj[v][i];
				last_weight = adj_weight[v][i];
			}
		}

		vertices[numOfVertices] = edges.size();
		cout << sum_of_edges << " edges in forward adj." << endl;
		cout << edges.size() << " edges inserted." << endl;
		cout << duplicated << " edges deleted." << endl;

		if (DIRECTED_FLAG == true) {

			sum_of_edges = 0;
			for (NodeID v = 0; v < numOfVertices; ++v) {
				sum_of_edges += r_adj[v].size();
			}
			r_vertices.resize(numOfVertices + 1);
			r_edges.reserve(sum_of_edges);
			EdgeID duplicated = 0;
			for (NodeID v = 0; v < numOfVertices; ++v) {
				r_vertices[v] = r_edges.size();
				NodeID last_v = numOfVertices;
				EdgeWeight last_weight = -1;
				for (NodeID i = 0; i < r_adj[v].size(); ++i) {
					if (r_adj[v][i] == last_v) {
						if (last_weight <= r_adj_weight[v][i]) {
							last_v = r_adj[v][i];
							last_weight = r_adj_weight[v][i];
							++duplicated;
							continue;
						}							
					}
					r_edges.push_back(NodeEdgeWeightPair(r_adj[v][i], r_adj_weight[v][i]));
					last_v = r_adj[v][i];
					last_weight = r_adj_weight[v][i];
				}
			}
			r_vertices[numOfVertices] = r_edges.size();

			cout << sum_of_edges << " edges in backward adj." << endl;
			cout << r_edges.size() << " edges inserted." << endl;
			cout << duplicated << " edges deleted." << endl;
		}

		numOfEdges=edges.size();
		cout<<"numOfEdges = "<<numOfEdges<<endl;

		adj.clear();
		r_adj.clear();
		adj.shrink_to_fit();
		r_adj.shrink_to_fit();

		adj_weight.clear();
		r_adj_weight.clear();
		adj_weight.shrink_to_fit();
		r_adj_weight.shrink_to_fit();

		return true;
	}


	// Load Graph file with format:
	// u	v	w
	// node id in the graph file should start with 0.
	bool load_graph(const char* graph_file,int graph_type=0) {
		if(numOfVertices!=0){
			numOfVertices=0;
			numOfEdges=0;
		}
		ifstream ifs(graph_file);

		vector<Edge> es;
		vector<EdgeWeight> es_weight;

		NodeID u, v;
		EdgeWeight w;
		for (; ifs >> u >> v >> w;) {
			numOfVertices = max(numOfVertices, max(u, v) + 1);
			//cout<<"u="<<u<<" v="<<v<<" w="<<w<<endl;
			if (u == v) continue; // Delete self loop.
			es.push_back(make_pair(u, v));
			//es.push_back(make_pair(u-1, v-1));//modified by wanjingyi
			es_weight.push_back(w);
		}

		cout<<"numOfReadEdges="<<es.size()<<endl;
		cout<<"numOfVertices="<<numOfVertices<<endl;
		
		if (ifs.bad()) return false;
		ifs.close();
		
		adj.resize(numOfVertices);
		adj_weight.resize(numOfVertices);
		assert(adj.size() == numOfVertices && adj_weight.size() == numOfVertices && "No enough memory for adj lists!");
		if (DIRECTED_FLAG == true) {
			r_adj.resize(numOfVertices);
			r_adj_weight.resize(numOfVertices);
			assert(r_adj.size() == numOfVertices && r_adj_weight.size() == numOfVertices  && "No enough memory for adj lists!");
		}

		if(graph_type==1){
			for (size_t i = 0; i < es.size(); ++i) {
				NodeID u = es[i].first, v = es[i].second;
				EdgeWeight w = es_weight[i];
				adj[u].push_back(v);
				adj_weight[u].push_back(w);
				if (DIRECTED_FLAG == true) {
					r_adj[v].push_back(u);
					r_adj_weight[v].push_back(w);
				}
			}
		}else{
			for (size_t i = 0; i < es.size(); ++i) {
				NodeID u = es[i].first, v = es[i].second;
				EdgeWeight w = es_weight[i];
				adj[u].push_back(v);
				adj_weight[u].push_back(w);
				adj[v].push_back(u);
				adj_weight[v].push_back(w);
				if(DIRECTED_FLAG==true){
					r_adj[v].push_back(u);
					r_adj_weight[v].push_back(w);
					r_adj[u].push_back(v);
					r_adj_weight[u].push_back(w);
				}
			}
		}

		for (NodeID v = 0; v < numOfVertices; ++v) {
			vector<NodeID>& adj_v = adj[v];
			vector<EdgeWeight>& adj_weight_v = adj_weight[v];
			vector<pair<NodeID, EdgeWeight> > adj_vw_v(adj_v.size());
			for (size_t i = 0; i < adj_v.size(); ++i) {
				NodeID u = adj_v[i];
				EdgeWeight w = adj_weight_v[i];
				adj_vw_v[i] = make_pair(u, w);
			}
			sort(adj_vw_v.begin(), adj_vw_v.end());
			for (size_t i = 0; i < adj_v.size(); ++i) {
				adj_v[i] = adj_vw_v[i].first;
				adj_weight_v[i] = adj_vw_v[i].second;
			}
			adj_vw_v.clear();

			if (DIRECTED_FLAG == true) {
				vector<NodeID>& r_adj_v = r_adj[v];
				vector<EdgeWeight>& r_adj_weight_v = r_adj_weight[v];
				vector<pair<NodeID, EdgeWeight> > r_adj_vw_v(r_adj_v.size());
				for (size_t i = 0; i < r_adj_v.size(); ++i) {
					NodeID u = r_adj_v[i];
					EdgeWeight w = r_adj_weight_v[i];
					r_adj_vw_v[i] = make_pair(u, w);
				}
				sort(r_adj_vw_v.begin(), r_adj_vw_v.end());
				for (size_t i = 0; i < r_adj_v.size(); ++i) {
					r_adj_v[i] = r_adj_vw_v[i].first;
					r_adj_weight_v[i] = r_adj_vw_v[i].second;
				}
				r_adj_vw_v.clear();
			}
		}

		long long sum_of_edges = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			sum_of_edges += adj[v].size();
		}


		vertices.resize(numOfVertices + 1);
		edges.reserve(sum_of_edges);
		EdgeID duplicated = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			vertices[v] = edges.size();

			NodeID last_v = numOfVertices;
			EdgeWeight last_weight = -1;
			for (NodeID i = 0; i < adj[v].size(); ++i) {
				if (adj[v][i] == last_v) {
					if (last_weight <= adj_weight[v][i]) {
						last_v = adj[v][i];
						last_weight = adj_weight[v][i];
						duplicated++;
						continue;
					}
				}
				edges.push_back(NodeEdgeWeightPair(adj[v][i], adj_weight[v][i]));
				last_v = adj[v][i];
				last_weight = adj_weight[v][i];
			}
		}

		vertices[numOfVertices] = edges.size();
		cout << sum_of_edges << " edges in forward adj." << endl;
		cout << edges.size() << " edges inserted." << endl;
		cout << duplicated << " edges deleted." << endl;

		if (DIRECTED_FLAG == true) {

			sum_of_edges = 0;
			for (NodeID v = 0; v < numOfVertices; ++v) {
				sum_of_edges += r_adj[v].size();
			}
			r_vertices.resize(numOfVertices + 1);
			r_edges.reserve(sum_of_edges);
			EdgeID duplicated = 0;
			for (NodeID v = 0; v < numOfVertices; ++v) {
				r_vertices[v] = r_edges.size();
				NodeID last_v = numOfVertices;
				EdgeWeight last_weight = -1;
				for (NodeID i = 0; i < r_adj[v].size(); ++i) {
					if (r_adj[v][i] == last_v) {
						if (last_weight <= r_adj_weight[v][i]) {
							last_v = r_adj[v][i];
							last_weight = r_adj_weight[v][i];
							++duplicated;
							continue;
						}							
					}
					r_edges.push_back(NodeEdgeWeightPair(r_adj[v][i], r_adj_weight[v][i]));
					last_v = r_adj[v][i];
					last_weight = r_adj_weight[v][i];
				}
			}
			r_vertices[numOfVertices] = r_edges.size();

			cout << sum_of_edges << " edges in backward adj." << endl;
			cout << r_edges.size() << " edges inserted." << endl;
			cout << duplicated << " edges deleted." << endl;
		}
		numOfEdges=edges.size();
		std::cout<<"numOfEdges = "<<numOfEdges<<std::endl;

		adj.clear();
		r_adj.clear();
		adj.shrink_to_fit();
		r_adj.shrink_to_fit();

		adj_weight.clear();
		r_adj_weight.clear();
		adj_weight.shrink_to_fit();
		r_adj_weight.shrink_to_fit();

		return true;
	}

	void save_overlay_graph(const char* write_filename,const vector<int>& newToOriginal,int graph_type=0){
			//************output new graph file using original ids************
			int numOfOverlayvertices=newToOriginal.size();
			string output_file_original(write_filename);
			output_file_original.append("_overlay_graph.original");
			ofstream ofs1(output_file_original);
			if(!ofs1.is_open()) std::cerr<<output_file_original<<" cannot be opened!"<<std::endl;
			if(graph_type==0){
				for(NodeID u=0;u<numOfOverlayvertices;++u){
					for(size_t i=vertices[u];i<vertices[u+1];++i){
						if(newToOriginal[edges[i].first]<newToOriginal[u]) continue;
						ofs1<<newToOriginal[u]<<" "<<newToOriginal[edges[i].first]<<" "<<edges[i].second<<endl;
					}
				}
			}else if(graph_type==1){
				for(NodeID u=0;u<numOfOverlayvertices;++u){
					for(size_t i=vertices[u];i<vertices[u+1];++i){
						ofs1<<newToOriginal[u]<<" "<<newToOriginal[edges[i].first]<<" "<<edges[i].second<<endl;
					}
				}
			}
			ofs1.close();
			std::cout<<"write "<<output_file_original<<" successfully!"<<std::endl;
	}
    
};

#endif