/*
 * @Descripttion:  define the ordering class
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-10 21:24:40
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-10-16 11:07:30
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
#ifndef ORDERING_H
#define ORDERING_H

#include<algorithm>
#include<unordered_set>
//#include<unordered_map>
#include <unordered_map>
#include<time.h>
#ifdef _WIN32
	#include<google/sparsehash/sparseconfig.h>
#endif
    #include<google/dense_hash_map>
	#include<google/dense_hash_set>

#include "graph.h"
#include "labels.h"
#include "time_util.h" 
#include "heap.h"
#include "paras.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define CNT 16
#define cover_value_type long long

using namespace time_util;

class Ordering{
    public:
        vector<NodeID> inv; // Fetch the original vertex id by a given ranking.
	    vector<NodeID> rank; // Fetch the ranking of a given vertex id.

        /*
         *@description: relabel the graph by order
         *@author: wanjingyi
         *@date: 2021-01-18
        */
        void Relabel(Graph& graph) {
            for (NodeID v = 0; v < numOfVertices; ++v) rank[inv[v]] = v;
            vector<EdgeID> new_vertices(numOfVertices + 1);
            vector<NodeID> new_edges;
            new_edges.reserve(graph.edges.size());
            for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
                NodeID originalVertex = inv[ranking];
                for (EdgeID eid = graph.vertices[originalVertex]; eid < graph.vertices[originalVertex + 1]; ++eid)
                    new_edges.push_back(rank[graph.edges[eid]]);
                new_vertices[ranking + 1] = new_edges.size();
            }
            graph.vertices.swap(new_vertices);
            graph.edges.swap(new_edges);

            if (DIRECTED_FLAG == true) {
                vector<EdgeID> r_new_vertices(numOfVertices + 1);
                vector<NodeID> r_new_edges;
                r_new_edges.reserve(graph.r_edges.size());
                for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
                    NodeID originalVertex = inv[ranking];
                    for (EdgeID eid = graph.r_vertices[originalVertex]; eid < graph.r_vertices[originalVertex + 1]; ++eid)
                        r_new_edges.push_back(rank[graph.r_edges[eid]]);
                    r_new_vertices[ranking + 1] = r_new_edges.size();
                }
                graph.r_vertices.swap(r_new_vertices);
                graph.r_edges.swap(r_new_edges);
            }

        }

        void Relabel(WGraph& wgraph,bool isOrder) {
		    for (NodeID v = 0; v < numOfVertices; ++v) rank[inv[v]] = v;
            // Array Representation
            vector<EdgeID> new_vertices(numOfVertices + 1);
            vector<NodeEdgeWeightPair> new_edges;
            vector<NodeEdgeWeightPair> new_edges_v;
            new_edges.reserve(wgraph.edges.size());
            for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
                NodeID originalVertex = inv[ranking];
                for (EdgeID eid = wgraph.vertices[originalVertex]; eid < wgraph.vertices[originalVertex + 1]; ++eid){
                    new_edges_v.push_back(make_pair(rank[wgraph.edges[eid].first],wgraph.edges[eid].second));
                }
                sort(new_edges_v.rbegin(),new_edges_v.rend());
                new_edges.insert(new_edges.end(),new_edges_v.begin(),new_edges_v.end());
                new_vertices[ranking + 1] = new_edges.size();
                new_edges_v.clear();
            }
            wgraph.vertices.swap(new_vertices);
		    wgraph.edges.swap(new_edges);

            if (DIRECTED_FLAG == true) {
                vector<EdgeID> r_new_vertices(numOfVertices + 1);
                vector<NodeEdgeWeightPair> r_new_edges;
                vector<NodeEdgeWeightPair> r_new_edges_v;
                r_new_edges.reserve(wgraph.r_edges.size());
                for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
                    NodeID originalVertex = inv[ranking];
                    for (EdgeID eid = wgraph.r_vertices[originalVertex]; eid < wgraph.r_vertices[originalVertex + 1]; ++eid)
                    {
                         r_new_edges_v.push_back(make_pair(rank[wgraph.r_edges[eid].first],wgraph.r_edges[eid].second));
                    }
                    sort(r_new_edges_v.rbegin(),r_new_edges_v.rend());
                    r_new_edges.insert(r_new_edges.end(),r_new_edges_v.begin(),r_new_edges_v.end());
                    r_new_vertices[ranking + 1] = r_new_edges.size();
                    r_new_edges_v.clear();
                }
                wgraph.r_vertices.swap(r_new_vertices);
                wgraph.r_edges.swap(r_new_edges);
            }
        }

	    void Relabel(WGraph& wgraph) {
		    for (NodeID v = 0; v < numOfVertices; ++v) rank[inv[v]] = v;
            // Array Representation
            vector<EdgeID> new_vertices(numOfVertices + 1);
            vector<NodeEdgeWeightPair> new_edges;
            new_edges.reserve(wgraph.edges.size());
            for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
                NodeID originalVertex = inv[ranking];
                for (EdgeID eid = wgraph.vertices[originalVertex]; eid < wgraph.vertices[originalVertex + 1]; ++eid)
                    new_edges.push_back(make_pair(rank[wgraph.edges[eid].first],wgraph.edges[eid].second));
                new_vertices[ranking + 1] = new_edges.size();
            }
            wgraph.vertices.swap(new_vertices);
		    wgraph.edges.swap(new_edges);

            if (DIRECTED_FLAG == true) {
                vector<EdgeID> r_new_vertices(numOfVertices + 1);
                vector<NodeEdgeWeightPair> r_new_edges;
                r_new_edges.reserve(wgraph.r_edges.size());
                for (NodeID ranking = 0; ranking < numOfVertices; ++ranking) {
                    NodeID originalVertex = inv[ranking];
                    for (EdgeID eid = wgraph.r_vertices[originalVertex]; eid < wgraph.r_vertices[originalVertex + 1]; ++eid)
                        r_new_edges.push_back(make_pair(rank[wgraph.r_edges[eid].first],wgraph.r_edges[eid].second));
                    r_new_vertices[ranking + 1] = r_new_edges.size();
                }
                wgraph.r_vertices.swap(r_new_vertices);
                wgraph.r_edges.swap(r_new_edges);
            }
        }

	void ReswapLabel(Graph& graph) {

		vector<vector<NodeID> > new_adj(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < graph.adj[v].size(); ++i) {
				new_adj[v].push_back(inv[graph.adj[rank[v]][i]]);
			}
		}
		graph.adj.swap(new_adj);
	}

	void save_rank(char* order_file) {
        string orderFileName(order_file);
        orderFileName.append(".order");
		ofstream ofs(orderFileName);
        if(!ofs.is_open()) cout<<"Cannot open "<<orderFileName<<endl;
		for (int i = 0; i < numOfVertices; ++i) {
			ofs << inv[i] << endl;
		}
		ofs.close();
	}

	~Ordering() {
		inv.clear();
		rank.clear();
	}

};

class Betweenness_Ordering : public Ordering{

	typedef	vector<NodeID> large_tree; // A V-sized parent-pointer array representing sampled trees. -1 for those vertices which do not appear in the tree.
	//typedef	unordered_map<NodeID, NodeID> small_tree; // A smaller hashmap, mapping a tree node vertex to its parent vertex.
	typedef	google::dense_hash_map<NodeID, NodeID> small_tree; // A smaller hashmap, mapping a tree node vertex to its parent vertex.
	
public:
	vector<double> iteration_time;
	HFLabel labels;
	bool SECOND_LEVEL;

	//Add for directed need to be modefied
	DLabel dlabels;
	
	long children_size;
	long r_children_size;

	vector<large_tree> ltrees;
	vector<small_tree> strees;
	vector<NodeID> trees_pointers;
	double init_time;
	double updating_time;
	double labeling_time;
	double adding_time;
	double selecting_time;
	long long bounded_resources;
	int num_of_trees;
	long total_resources;
	long long alive_resources;
	double u1_time;
	double u2_time;
	double u3_time = 0;
	double u4_time = 0;
	double u5_time = 0;
	double u6_time = 0;
	double u7_time = 0;
	double u8_time = 0;
	double u9_time = 0;
	double u10_time = 0;
	double u11_time = 0;
	double u12_time = 0;
	double u81_time = 0;
	double u82_time = 0;
	double u91_time = 0;
	double u92_time = 0;

	small_tree empty_small_tree = small_tree();
	bool BUILD_SMALL_TREE_FLAG = false;
	double get_small_parent_time;
	bool switch_small;
	
	class nodeid_vector_hasher {
	public:
		std::size_t operator()(std::pair<NodeID, std::vector<NodeID> > const& pairvec) const {
		  std::size_t seed = pairvec.second.size();
		  seed ^= pairvec.first + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		  for(auto& i : pairvec.second) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		  }
		  return seed;
		}
	};

	//using boost::hash_combine
	template <class T>
	inline void hash_combine(std::size_t& seed, T const& v)
	{
		seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	}

	int build_large_tree_dij(NodeID source, large_tree& ltree, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, WGraph& wgraph, int accum_new_arcs, int labeled_arcs_bound) {
		NodeID new_tree_arcs = 0;

		descendants.clear();
		ltree.resize(numOfVertices, -1);
	//	ltree[source] = -1; // As the root, it has no parent.
		descendants_parents[source] = -1;

		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

		/*vector<vector<NodeID> >& adj = wgraph.adj;
		vector<vector<EdgeWeight> >& adj_w = wgraph.adj_weight;*/


		pqueue.update(source, 0);
		//++new_tree_arcs;

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		visited_que.push(source);
		while(!pqueue.empty()){
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				vis[v] = true;


				if (usd[v] == true) continue;
				//int adj_v_size = adj[v].size();

				if (DIRECTED_FLAG == false) {
					pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						NodeID w = tmp_idx_v.first[i];
						EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= v_d) {
							goto pruned;
						}
					}
				}
				else {
					pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

					// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
					for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
						NodeID w = r_tmp_idx_v.first[i];
						EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];

						if (td <= v_d) {
							goto pruned;
						}
					}
				}

				// num of hubs touched so far when building a new tree.
				++new_tree_arcs;
				descendants.push_back(v);
				
				if (accum_new_arcs + new_tree_arcs > labeled_arcs_bound) {
					break;
				}

				/*for (size_t i = 0; i < adj_v_size; ++i) {
					NodeID w = adj[v][i];
					EdgeWeight w_d = adj_w[v][i] + v_d;		*/

				 //Array Representations.
				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){
					
					NodeID w = wgraph.edges[eid].first;
					EdgeWeight w_d = wgraph.edges[eid].second + v_d;

					//if (wgraph.adj[v][eid - wgraph.vertices[v]] != w || wgraph.adj_weight[v][eid - wgraph.vertices[v]] + v_d != w_d) cout << "we got problems build large tree." << endl;

					if (!vis[w]) {
						if (distances[w] > w_d) {
							pqueue.update(w, w_d);
							distances[w] = w_d;
							//ltree[w] = v; // Set parent nodes. 
							descendants_parents[w] = v;
							visited_que.push(w);
						}
					}
				}
			pruned:
				{}
		}
		while (!visited_que.empty()) {
			NodeID vis_v = visited_que.front();
			visited_que.pop();
			vis[vis_v] = false;
			distances[vis_v] = INF_WEIGHT;
			pqueue.clear(vis_v);
		}
	/*	while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			vis[v] = false;
			distances[v] = INF_WEIGHT;
		}*/

		pqueue.clear_n();

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}

		for (int i = 0; i < descendants.size(); ++i)
			ltree[descendants[i]] = descendants_parents[descendants[i]];

		return new_tree_arcs;
	}

	void shrink_tree(int tid, NodeID root, large_tree& ltree_tid) {
		small_tree new_stree;
		//GOOGLE DENSE MAP
		new_stree.set_empty_key(numOfVertices+1);
		new_stree.set_deleted_key(numOfVertices+2);
		for (NodeID v = 0; v < ltree_tid.size(); ++v) {
			int pt = ltree_tid[v];
			if (pt != -1)
				new_stree.insert(make_pair(v, pt));
		}
		new_stree.insert(make_pair(root, -1));
		strees.push_back(new_stree);
		trees_pointers[tid] = strees.size() - 1;
		alive_resources += (long long)new_stree.size();

		//ltrees[tid].clear();
		// FIX
		ltree_tid.clear();
		ltree_tid.shrink_to_fit();
		alive_resources -= (long long)numOfVertices;
	}

	void enlarge_tree(int tid, NodeID root, small_tree& stree_tid, large_tree& ltree_tid) {
		ltree_tid.resize(numOfVertices,-1);
		
		for (small_tree::iterator it = stree_tid.begin(); it != stree_tid.end(); ++it) {			
			NodeID v = (*it).first;
			NodeID pt = (*it).second;
			ltree_tid[v] = pt;
		}

		ltree_tid[root] = -1;

		ltrees.push_back(ltree_tid);

		trees_pointers[tid] = -1;
		alive_resources += (long long)numOfVertices;

		alive_resources -= (long long)stree_tid.size();
		stree_tid.clear();
		stree_tid.resize(0);
	}

	int get_large_parent(NodeID v, large_tree& ltree) {
		return ltree[v];
	}

	int get_small_parent_withouttest(NodeID v, small_tree& stree) {

		return stree[v];

	}

	int get_small_parent(NodeID v, small_tree& stree) {
		small_tree::iterator it = stree.find(v);

		if (it != stree.end())
			return stree[v];
		else
			return -1;
	}

	void remove_large_node(NodeID v, large_tree& ltree) {
		ltree[v] = -1;
	}
	void remove_small_node(NodeID v, small_tree& stree) {
		stree.erase(v);
	//	stree.resize(0);
		//stree[v] = -1;
	}

	//GOOGLE
	void cnt16_2max_queue(vector<cover_value_type>& picked_value, vector<cover_value_type>& sum, vector<bool>& has_picked, /*unordered_set<NodeID>& affected_v*/ google::dense_hash_set<NodeID>& affected_v, vector<bool>& affected_list, const vector<vector<cover_value_type> >& cnt16, bool& switch_small, bool& switch_set, benchmark::heap<2, cover_value_type, NodeID>& value_heap) {

		// First Round, build the entire cnt16 one time.
		if (picked_value.size() == 0) {
			picked_value.resize(numOfVertices, 0);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<cover_value_type>& cnt16_v = cnt16[v];
				cover_value_type &picked_value_v = picked_value[v];
				cover_value_type& sum_v = sum[v];
				sum_v = 0;
				for (int i = 0; i < CNT; ++i) {
					sum_v += cnt16_v[i];
					if (cnt16_v[i] > max2)
						max2 = cnt16_v[i];
					if (max2 > max1) {
						NodeID tmp = max1;
						max1 = max2;
						max2 = tmp;
					}
				}
				picked_value_v = sum_v - max1 - max2;
				if (switch_small == true) {
					value_heap.update(v, -picked_value_v);
				}
			}
		}

		//for (unordered_set<NodeID>::iterator it = affected_v.begin(); it != affected_v.end(); ++it) {
		if (!switch_set) {
			for (int v = 0; v < numOfVertices; ++v) {
				//NodeID v = *it;
				if (affected_list[v] == false) continue;
				affected_list[v] = false;
				if (has_picked[v] == true) continue;
				//if (v == best_vid) continue;
				cover_value_type& sum_v = sum[v];

				/*	if (sum_v < best_so_far) {
				picked_value[v] = sum_v;
				if (switch_small == true)
				value_heap.update(v, sum_v);
				lazy_updates[v] = true;
				continue;
				}*/
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<cover_value_type>& cnt16_v = cnt16[v];
				cover_value_type& picked_value_v = picked_value[v];


				sum_v = 0;
				for (int i = 0; i < CNT; ++i) {
						sum_v += cnt16_v[i];
					if (cnt16_v[i] > max2)
						max2 = cnt16_v[i];
					if (max2 > max1) {
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
					}
				}
				picked_value_v = sum_v - max1 - max2;
				//if (picked_value_v > best_so_far) best_so_far = picked_value_v;
				if (switch_small == true) {
					value_heap.update(v, -picked_value_v);
				}
			}
		}
		else {
			//GOOGLE
			
			//for (unordered_set<NodeID>::iterator it = affected_v.begin(); it != affected_v.end(); ++it) {
			for (google::dense_hash_set<NodeID>::iterator it = affected_v.begin(); it != affected_v.end(); ++it) {

				NodeID v = *it;
				affected_list[v] = false;
				if (has_picked[v] == true) continue;
				//if (v == best_vid) continue;
				cover_value_type& sum_v = sum[v];

				//if (sum_v < best_so_far) {
				//	picked_value[v] = sum_v;
				//	if (switch_small == true)
				//		value_heap.update(v, sum_v);
				//	lazy_updates[v] = true;
				//	continue;
				//}
				NodeID max1 = 0;
				NodeID max2 = 0;
				const vector<cover_value_type>& cnt16_v = cnt16[v];
				cover_value_type& picked_value_v = picked_value[v];


				sum_v = 0;
				for (int i = 0; i < CNT; ++i) {
						sum_v += cnt16_v[i];
					if (cnt16_v[i] > max2)
						max2 = cnt16_v[i];
					if (max2 > max1) {
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
					}
				}
				picked_value_v = sum_v - max1 - max2;
				//if (picked_value_v > best_so_far) best_so_far = picked_value_v;
				if (switch_small == true) {
					value_heap.update(v, -picked_value_v);
				}
			}
			affected_v.clear();
		}
	}

	int labeling_source_dij(NodeID source, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, NodeID ranking, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, WGraph& wgraph) {
		NodeID labeled_arcs = 0;

		NodeID visited_arcs = 0;
	/*	vector<vector<NodeID> >& adj = wgraph.adj;
		vector<vector<EdgeWeight> >& adj_w = wgraph.adj_weight;*/

		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}

		pqueue.update(source, 0);


		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			vis[v] = true;
			visited_que.push(v);

			if (usd[v]) continue;


			if (DIRECTED_FLAG == false) {
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_forward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;

				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);
			}
			else {
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

				// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
				for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
					NodeID w = r_tmp_idx_v.first[i];
					EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_forward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;
				
				r_tmp_idx_v.first.back() = ranking;
				r_tmp_idx_v.second.back() = v_d;
				r_tmp_idx_v.first.push_back(numOfVertices);
				r_tmp_idx_v.second.push_back(INF_WEIGHT);
			}

			/*		for (size_t i = 0; i < adj[v].size(); ++i) {
			
				NodeID w = adj[v][i];
				EdgeWeight w_d = adj_w[v][i] + v_d;*/

			//Array Representation
			
			visited_arcs += (wgraph.vertices[v + 1] - wgraph.vertices[v]);
			
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){

				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;

				//		if (wgraph.adj[v][eid - wgraph.vertices[v]] != w || wgraph.adj_weight[v][eid - wgraph.vertices[v]] + v_d != w_d) cout << "we got problems labeling tree." << endl;
			
				if (!vis[w]) {
					if (distances[w] > w_d) {
						pqueue.update(w, w_d);
						distances[w] = w_d;
					}
				}
			}
			pruned_forward:
			{}
		}
		
		// Clear foward search structures.
		while (!visited_que.empty()) {
			NodeID vis_v = visited_que.front();
			visited_que.pop();
			vis[vis_v] = false;
			distances[vis_v] = INF_WEIGHT;
		}

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}


		// Backward search.
		if (DIRECTED_FLAG == true) {

			/*		vector<vector<NodeID> >& r_adj = wgraph.r_adj;
			vector<vector<EdgeWeight> >& r_adj_w = wgraph.r_adj_weight;*/

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}


			pqueue.update(source, 0);
			while (!pqueue.empty()) {
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				vis[v] = true;
				visited_que.push(v);


				if (usd[v]) continue;


				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];

				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned_backward;
					}
				}

				// arcs touched so far when building this labels
				labeled_arcs++;
				tmp_idx_v.first.back() = ranking;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

		/*		for (size_t i = 0; i < r_adj[v].size(); ++i) {
					NodeID w = r_adj[v][i];
					EdgeWeight w_d = r_adj_w[v][i] + v_d;*/
				
				////Array Representation
				
				visited_arcs += (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]);
				
				for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid){

					NodeID w = wgraph.r_edges[eid].first;
					EdgeWeight w_d = wgraph.r_edges[eid].second + v_d ;

					if (!vis[w]) {
						if (distances[w] > w_d) {
							pqueue.update(w, w_d);
							distances[w] = w_d;
						}
					}
				}
			pruned_backward:
				{}
			}

		while (!visited_que.empty()) {
			NodeID vis_v = visited_que.front();
			visited_que.pop();
			vis[vis_v] = false;
			distances[vis_v] = INF_WEIGHT;
			pqueue.clear(vis_v);
		}
		pqueue.clear_n();

		for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
			dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;
		}
	
		
		//cout << visited_arcs << endl;
	
		usd[source] = true;
		return labeled_arcs;
	}

	void get_descendants(int tree_id, NodeID picked_v, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, vector<NodeID>& que, vector<bool>& vis, WGraph& wgraph, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {

		descendants.clear();
		//	descendants.reserve(V / 8);

		descendants.push_back(picked_v);
		descendants_parents[picked_v] = -1;

		vis[picked_v] = true;

		/*	bool isLarge = false;
		int stree_idx = trees_pointers[tree_id];
		if (stree_idx == -1)
		isLarge = true;

		small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
		large_tree& ltree_tid = ltrees[tree_id];*/


		for (int i = 0; i < descendants.size(); ++i) {
			NodeID v = descendants[i];

		/*	vector<EdgeWeight>& adj_v = wgraph.adj[v];
			NodeID adj_v_size = adj_v.size();
			
			for (int j = 0; j < adj_v_size; ++j) {
				NodeID w = adj_v[j];*/

			//Array Representation

			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1];++eid){
				
				NodeID w = wgraph.edges[eid].first;

				//if(wgraph.adj[v][eid - wgraph.vertices[v]] != w) cout << "we got problems. get descent" << endl;

				NodeID parentID = 0;


				if (isLarge)
					parentID = get_large_parent(w, ltree_tid);
				else {
					//	double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(w, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}


				if (!vis[w] && parentID == v) {
					descendants.push_back(w);
					descendants_parents[w] = v;
					vis[w] = true;
				}
			}
		}
		int dsize = descendants.size();
		for (int i = 0; i <dsize; ++i) vis[descendants[i]] = false;
	}

	void get_ascendants(int tree_id, NodeID picked_v, vector<NodeID>& ascendants, vector<NodeID>& que, vector<bool>& vis, WGraph& wgraph, const bool isLarge, small_tree& stree_tid, large_tree& ltree_tid) {
		ascendants.clear();

		NodeID que_t0 = 0, que_t1 = 0, que_h = 0;
		//	vector<vector<NodeID> >& adj = graph.adj;
		//vector<vector<NodeID> >& r_adj = graph.r_adj;

		que[que_h++] = picked_v;
		vis[picked_v] = true;
		que_t1 = que_h;


		//bool isLarge = false;
		//int stree_idx = trees_pointers[tree_id];
		//if (stree_idx == -1) 
		//	isLarge = true;
		//
		//small_tree& stree_tid = (stree_idx != -1 )? strees[stree_idx] : empty_small_tree;
		//large_tree& ltree_tid = ltrees[tree_id];

		for (EdgeWeight d = 0; que_t0 < que_h; d = d + 1) {
			for (NodeID que_i = que_t0; que_i < que_t1; ++que_i) {
				NodeID v = que[que_i];
				ascendants.push_back(v);

				NodeID parentID = 0;
				if (isLarge)
					parentID = get_large_parent(v, ltree_tid);
				else {
					//double get_small_parent_incre_time = GetCurrentTimeSec();
					parentID = get_small_parent(v, stree_tid);
					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
				}
				if (DIRECTED_FLAG == false) {
				/*	vector<EdgeWeight>& adj_v = wgraph.adj[v];
					int adj_v_size = adj_v.size();
					for (int i = 0; i < adj_v_size; ++i) {
						NodeID w = adj_v[i];*/

					//Array Representation
					for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){

						NodeID w = wgraph.edges[eid].first;

						//if (!vis[w] && get_parent(tree_id, v) == w) {
						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
				else {
				/*	vector<EdgeWeight>& r_adj_v = wgraph.r_adj[v];
					int r_adj_v_size = r_adj_v.size();
					for (int i = 0; i < r_adj_v_size; ++i) {
						NodeID w = r_adj_v[i];*/

					//Array Representation
					for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid){

						NodeID w = wgraph.r_edges[eid].first;

					//	if (wgraph.r_adj[v][eid - wgraph.r_vertices[v]] != w) cout << "we got problems get ascendants." << endl;

						if (!vis[w] && parentID == w) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
				}
			}
			que_t0 = que_t1;
			que_t1 = que_h;
		}

		for (size_t i = 0; i < que_h; ++i) vis[que[i]] = false;
	}

	int build_small_tree_dij(NodeID source, small_tree& stree, vector<NodeID>& descendants, vector<NodeID>& descendants_parents, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, queue<NodeID> &visited_que, vector<EdgeWeight>& distances, vector<bool>& vis, vector<EdgeWeight>& dst_r, vector<bool>& usd, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& tmp_idx, vector<pair<vector<NodeID>, vector<EdgeWeight> > >& r_tmp_idx, WGraph& wgraph, int accum_new_arcs, int labeled_arcs_bound) {
		NodeID new_tree_arcs = 0;

		descendants.clear();
		stree.clear();
		//	ltree[source] = -1; // As the root, it has no parent.
		descendants_parents[source] = -1;

		const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[source];
		const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[source];

		/*vector<vector<NodeID> >& adj = wgraph.adj;
		vector<vector<EdgeWeight> >& adj_w = wgraph.adj_weight;*/


		pqueue.update(source, 0);
		//++new_tree_arcs;

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
		}
		visited_que.push(source);
		while (!pqueue.empty()) {
			NodeID v;
			EdgeWeight v_d;
			pqueue.extract_min(v, v_d);
			vis[v] = true;
			//visited_que.push(v);
			


			if (usd[v] == true) continue;

			//int adj_v_size = adj[v].size();

			if (DIRECTED_FLAG == false) {
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						goto pruned;
					}
				}
			}
			else {
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];

				// Pruned by the forward labels of r and backward labels of v in the forward search from r when reaching v.
				for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
					NodeID w = r_tmp_idx_v.first[i];
					EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];

					if (td <= v_d) {
						goto pruned;
					}
				}
			}

			// num of hubs touched so far when building a new tree.
			++new_tree_arcs;
			descendants.push_back(v);

			if (accum_new_arcs + new_tree_arcs > labeled_arcs_bound) {
				break;
			}
			
		/*	for (size_t i = 0; i < adj_v_size; ++i) {
				NodeID w = adj[v][i];
				EdgeWeight w_d = adj_w[v][i] + v_d;*/

			//Array Representation
			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid){
				NodeID w = wgraph.edges[eid].first;
				EdgeWeight w_d = wgraph.edges[eid].second + v_d;


				//if (wgraph.adj[v][eid - wgraph.vertices[v]] != w || wgraph.adj_weight[v][eid - wgraph.vertices[v]] + v_d != w_d ) cout << "we got problems build small tree." << endl;

				if (!vis[w]) {
					if (distances[w] > w_d) {
						pqueue.update(w, w_d);
						distances[w] = w_d;
						//ltree[w] = v; // Set parent nodes. 
						descendants_parents[w] = v;
						visited_que.push(w);
					}
				}
			}
		pruned:
			{}
		}
		while (!visited_que.empty()) {
			NodeID vis_v = visited_que.front();
			visited_que.pop();
			vis[vis_v] = false;
			distances[vis_v] = INF_WEIGHT;
			pqueue.clear(vis_v);
		}
	
		//while (!pqueue.empty()) {
		//	NodeID v;
		//	EdgeWeight v_d;
		//	pqueue.extract_min(v, v_d);
		//	//vis[v] = false;
		//	//distances[v] = INF_WEIGHT;
		//}
		
		pqueue.clear_n();
	//	pqueue.clear();

		for (NodeID i = 0; i < tmp_idx_r.first.size(); ++i) {
			dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
		}

		for (int i = 0; i < descendants.size(); ++i)
			stree.insert(make_pair(descendants[i], descendants_parents[descendants[i]]));

		return new_tree_arcs;
	}

    
	/**
	 * @description: used for generating the betweenness values for overlay graph and store to the
	 * @param {NodeID} k
	 * @param {double} beta
	 * @param {NodeID} stopv
	 * @param {unsigned} int
	 * @param {long long} bound_times
	 * @param {long long} trans_times
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	// Betweenness_Ordering(NodeID k, double beta, WGraph& wgraph, NodeID stopv,double& _ordering_time, double& _hub_time, vector<unsigned int>& betweenness_values,long long bound_times = 10, long long trans_times = 8) {
	// 	double start_time = GetCurrentTimeSec();
	// 	betweenness_values.resize(numOfVertices,0);
	// 	if (k < CNT) {
	// 		cout << " Should pick at least 16 seed!" << endl;
	// 		return;
	// 	}
	// 	iteration_time.resize(numOfVertices);

	// 	updating_time = 0;
	// 	labeling_time = 0;
	// 	adding_time = 0;
	// 	selecting_time = 0;
	// 	init_time = 0;
	// 	u1_time = 0;
	// 	u2_time = 0;
	// 	u3_time = 0;
	// 	num_of_trees = 0;
	// 	total_resources = 0;


	// 	inv.resize(numOfVertices, -1);
	// 	rank.resize(numOfVertices, -1);

	// 	labels.index_.resize(numOfVertices);

	// 	// Resource Bounding.
	// 	NodeID new_tree_arcs = 0;
	// 	NodeID labeled_arcs = 0;
	// 	alive_resources = 0;
	// 	bounded_resources = bound_times * (long long)k * (long long)numOfVertices;

	// 	vector<NodeID> st_roots; // The root vertices of the sampled trees.
	// 	vector<NodeID> st_sizes; // The size of the sampled trees.

	// 	NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / trans_times; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.

	// 	init_time = GetCurrentTimeSec();
	// 	// Building the initial k trees.
	// 	//	srand((unsigned)time(NULL));
	// 	srand(0x728f9d0);


	// 	ltrees.reserve(numOfVertices);
	// 	strees.reserve(numOfVertices);
	// 	trees_pointers.reserve(numOfVertices);
	// 	st_roots.reserve(numOfVertices);
	// 	st_sizes.reserve(numOfVertices);
	// 	st_roots.resize(k);
	// 	st_sizes.resize(k);

	// 	ltrees.resize(k);
	// 	//strees.clear();
	// 	trees_pointers.resize(k, -1); // Point the sampling trees to its tree representations. large tree [0, V) small tree [V, 2V).
	// 	//vector<unordered_map<int, int> > cover(k, unordered_map<int, int>());

	// 	vector<NodeID> ascendants;
	// 	ascendants.reserve(numOfVertices);
	// 	vector<NodeID> descendants;
	// 	descendants.reserve(numOfVertices);
	// 	vector<NodeID> que(numOfVertices);
	// 	vector<bool> vis(numOfVertices);
	// 	vector<bool> has_picked(numOfVertices);
	// 	vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
	// 	//		vector<bool> affected(V, true);
	// 	vector<bool> picked_as_root(numOfVertices, false);
	// 	vector<bool> trivial_root(numOfVertices, false);
	// 	int picked_count = 0;
	// 	int used_as_cover = 0;
	// 	vector<cover_value_type> sum(numOfVertices);
	// 	bool stop_flag = false;


	// 	// tmp structure used in the dijkstra searches. 

	// 	// Binary Heap
	// 	benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

	// 	queue<NodeID> visited_que;
	// 	vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
	// 	vector<bool> usd(numOfVertices, false);

	// 	//GOOGLE
	// 	//		unordered_set<NodeID> affected_v;
	// 	google::dense_hash_set<NodeID> affected_v;
	// 	affected_v.set_empty_key(numOfVertices + 1);
	// 	affected_v.set_deleted_key(numOfVertices + 2);

	// 	vector<cover_value_type> cover_tid(numOfVertices, 0);
	// 	//vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14


	// 	vector<bool> affected_list(numOfVertices, false);
	// 	int affected_count = 0;
	// 	bool switch_set = false;

	// 	//vector<unordered_set<int> > reverse_index(0);
	// 	//GOOGLE DENSE SET
	// 	vector<google::dense_hash_set<NodeID> > reverse_index;
	// 	//vector<NodeID> in_tree_count(numOfVertices);

	// 	// Preparing basic structure for pl algorithm.
	// 	vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
	// 	vector<pair<vector<NodeID>, vector<EdgeWeight> > >
	// 		tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
	// 			vector<EdgeWeight>(1, INF_WEIGHT)));

	// 	vector<pair<vector<NodeID>, vector<EdgeWeight> > >
	// 		r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
	// 			vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

	// 	vector<NodeID> descendants_parents(numOfVertices, -1);

	// 	// Building the initial k sampled trees. 
	// 	for (int i = 0; i < k; ++i) {
	// 		NodeID r = rand() % numOfVertices;
	// 		if (picked_as_root[r] == true) { --i; continue; }

	// 		st_roots[i] = r;
	// 		large_tree &ltree_i = ltrees[i];
	// 		picked_as_root[r] = true;
	// 		picked_count++;
	// 		used_as_cover++;

	// 		build_large_tree_dij(r, ltree_i, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
	// 		st_sizes[i] = descendants.size();
	// 		alive_resources += (long long)numOfVertices;

	// 		// Turn large tree to small tree.
	// 		// 1. Build new small tree with ltree_i and add it to strees;
	// 		// 2. Set tree_pointers[i] to the position of the strees;
	// 		// 3. Clear the large tree ltree_i by setting its size to 0.
	// 		if (st_sizes[i] < size_thr)
	// 			shrink_tree(i, st_roots[i], ltrees[i]);

	// 		bool isLarge = false;
	// 		int stree_idx = trees_pointers[i];
	// 		if (stree_idx == -1)
	// 			isLarge = true;

	// 		small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

	// 		// Calculate the subtree size of each nodes.
	// 		//	unordered_map<int, int>& cover_i = cover[i];
	// 		//vector<int> cover_i(numOfVertices, 0);

	// 		for (int di = descendants.size() - 1; di > -1; --di) {
	// 			NodeID q = descendants[di];

	// 			++cover_tid[q]; // Increment because it hits the r-q path.
	// 			++cnt16[q][i%CNT];
	// 			//++sum[q];

	// 			int parent_node;

	// 			if (isLarge)
	// 				parent_node = get_large_parent(q, ltree_i);
	// 			else
	// 				parent_node = get_small_parent(q, stree_i);

	// 			if (parent_node > -1) {
	// 				cover_tid[parent_node] += cover_tid[q];
	// 				cnt16[parent_node][i%CNT] += cover_tid[q];
	// 				//	sum[parent_node] += cover_i[q];
	// 			}
	// 		}
	// 		for (int di = descendants.size() - 1; di > -1; --di) {
	// 			NodeID q = descendants[di];
	// 			cover_tid[q] = 0;
	// 		}
	// 	}



	// 	//cnt16_aggreate(cnt16, has_picked, cover);

	// 	vector<cover_value_type> picked_value;
	// 	benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
	// 	//bool switch_small = false;
	// 	switch_small = false;

	// 	//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
    //    //		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);
	// 	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

	// 	init_time = GetCurrentTimeSec() - init_time;

	// 	cout << "Finish initial sampling " << init_time << "secs." << endl;

	// 	//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

	// 	double max_value = -1;


	// 	//output betweenness
	// 	//ofstream ofs("../dataset/Manhattan/order_values/original.betweenness");
	// 	//if(!ofs.is_open()) std::cout<<"Betweenness.values file can't be opened!"<<endl;
	// 	// Perform V iterations to rank all vertices.
	// 	NodeID ranking = 0;
	// 	for (; ranking < numOfVertices; ++ranking) {
	// 		if (ranking % (numOfVertices / 10) == 0)
	// 			cout << "picking the " << ranking << "-th orders." << endl;

	// 		double time_this_iteration = GetCurrentTimeSec();

	// 		NodeID pv = 0;
	// 		long long pv_cover = -1;

	// 		if (ranking == stopv) {
	// 			stop_flag = true;
	// 			cout << "topv labeling time:" << labeling_time << endl;
	// 			break;
	// 		}

	// 		double selecting_acc_time = GetCurrentTimeSec();
	// 		// In the beginning, pick the best cover vertex by a linear scan.
	// 		if (switch_small == false) {
	// 			double u1_acc_time = GetCurrentTimeSec();

	// 			for (NodeID v = 0; v < numOfVertices; ++v) {
	// 				if (has_picked[v] == true) continue;
	// 				if (pv_cover == picked_value[v]) {
	// 					if (sum[pv] < sum[v]) {
	// 						pv_cover = picked_value[v];
	// 						pv = v;
	// 					}
	// 				}
	// 				else if (pv_cover < picked_value[v]) {
	// 					pv_cover = picked_value[v];
	// 					pv = v;
	// 				}
	// 			}

	// 			if (sum[pv] < numOfVertices / 8 && switch_small == false) {
	// 				//cout << "u3 time:" << u3_time << endl;
	// 				cout << ranking << endl;
	// 				cout << pv_cover << endl;
	// 				cout << sum[pv] << endl;
	// 				cout << picked_count << endl;
	// 				cout << "Switch to max-heap representations." << endl;
	// 				switch_small = true;
	// 			}
	// 			u1_time += (GetCurrentTimeSec() - u1_acc_time);
	// 		}
	// 		else {// Later, utilizing the max-heap directly
	// 			double u2_acc_time = GetCurrentTimeSec();

	// 			value_heap.extract_min(pv, pv_cover);
	// 			pv_cover = -pv_cover;
	// 			u2_time += (GetCurrentTimeSec() - u2_acc_time);
	// 			/*
	// 			int max1 = 0;
	// 			int max2 = 0;
	// 			int max1v = -1;
	// 			int max2v = -1;
	// 			for(int  i = 0; i < CNT; ++i){
	// 				if(max2 <= cnt16[pv][i]){
	// 					max2 = cnt16[pv][i];
	// 					max2v = i;
	// 				}
	// 				if(max1 <= max2){
	// 					int tmp = max1;
	// 					max1 = max2;
	// 					max2 = tmp;
	// 					int tmpv = max1v;
	// 					max1v = max2v;
	// 					max2v = max1v;
	// 				}
	// 			}
	// 			long long sumstd = 0;
	// 			long long summean = 0;
	// 			for(int  i = 0; i < CNT; ++i){
	// 				if(i != max1v && i != max2v){
	// 					sumstd += (long long)cnt16[pv][i] * (long long)cnt16[pv][i];
	// 					summean += (long long)cnt16[pv][i];
	// 				}
	// 			}
	// 			sumstd = sumstd/(long long)(CNT-2);
	// 			summean = summean/(long long)(CNT-2);
	// 			cout << sqrt(sumstd - summean * summean) << endl;*/
	// 		}
	// 		//output betweenness
	// 		//ofs<<pv<<" "<<pv_cover<<endl;
	// 		betweenness_values[pv]=pv_cover;
	// 		selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

	// 		//// FIX
	// 		//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size() * 10) {//sum[pv] <= in_tree_count[pv] * 10){//
	// 		//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
	// 		//	break;
	// 		//}

	// 		//if (has_picked[pv] == true)
	// 		//	cout << "wrong1" << endl;

	// 		has_picked[pv] = true;
	// 		double labeling_acc_time = GetCurrentTimeSec();
	// 		labeled_arcs = labeling_source_dij(pv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
	// 		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
	// 		inv[ranking] = pv;
	// 		rank[pv] = ranking;

	// 		int current_trees_num = st_roots.size();

	// 		// 4 steps in updating:
	// 		double updating_acc_time = GetCurrentTimeSec();
	// 		if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
	// 			//	double u3_acc_time = GetCurrentTimeSec();
	// 			for (int tid = 0; tid < current_trees_num; ++tid) {

	// 				// Get the descendants and ascendants of pv in tree tid.
	// 				if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

	// 				bool isLarge = false;
	// 				int stree_idx = trees_pointers[tid];
	// 				if (stree_idx == -1)
	// 					isLarge = true;

	// 				small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
	// 				large_tree& ltree_tid = ltrees[tid];

	// 				// pv has already been covered in tree tid.
	// 				int pvpid = 0;
	// 				if (isLarge)
	// 					pvpid = get_large_parent(pv, ltree_tid);
	// 				else {
	// 					//double get_small_parent_incre_time = GetCurrentTimeSec();
	// 					pvpid = get_small_parent(pv, stree_tid);
	// 					//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
	// 				}
	// 				if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
	// 					continue;


	// 				double u3_acc_time = GetCurrentTimeSec();

	// 				// //1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
	// 				//get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph.adj, isLarge, stree_tid, ltree_tid);
	// 				//get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

	// 				//Array Representation
	// 				get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
	// 				get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

	// 				u3_time += (GetCurrentTimeSec() - u3_acc_time);


	// 				double u4_acc_time = GetCurrentTimeSec();
	// 				NodeID subtract_size = descendants.size();
	// 				for (int ai = 0; ai < ascendants.size(); ++ai) {
	// 					NodeID av = ascendants[ai];
	// 					// source - av - pv: subtract the size of subtree rooted at pv.
	// 					if (pv == av) continue; // It will be subtracted in the descending search.
	// 				//	cover[tid][av] -= subtract_size;
	// 					cnt16[av][tid%CNT] -= subtract_size;
	// 					//sum[av] -= subtract_size;
	// 					//affected[av] = true;

	// 					if (affected_list[av] == false) {
	// 						affected_list[av] = true;
	// 						affected_count++;
	// 						if (switch_set)
	// 							affected_v.insert(av);
	// 					}
	// 				}
	// 				u4_time += (GetCurrentTimeSec() - u4_acc_time);


	// 				double u5_acc_time = GetCurrentTimeSec();
	// 				// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
	// 				//vector<int> cover_tid(numOfVertices, 0);
	// 				for (int di = descendants.size() - 1; di > -1; --di) {
	// 					NodeID dv = descendants[di];
	// 					// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
	// 					cover_tid[dv]++;
	// 					cnt16[dv][tid%CNT] -= cover_tid[dv];
	// 					//sum[dv] -= cover_tid[dv];

	// 					if (descendants_parents[dv] != -1)
	// 						cover_tid[descendants_parents[dv]] += cover_tid[dv];
	// 					if (cnt16[dv][tid%CNT] < 0)
	// 						cout << "a6 " << cnt16[dv][tid%CNT] << endl;

	// 					/*if (isLarge)
	// 						remove_large_node(dv, ltree_tid);
	// 					else
	// 						remove_small_node(dv, stree_tid);*/

	// 						//						remove_node(tid, dv);
	// 												//affected[dv] = true;

	// 					if (affected_list[dv] == false) {
	// 						affected_count++;
	// 						affected_list[dv] = true;
	// 						if (switch_set)
	// 							affected_v.insert(dv);
	// 					}
	// 				}
	// 				for (int di = descendants.size() - 1; di > -1; --di) {
	// 					NodeID dv = descendants[di];
	// 					cover_tid[dv] = 0;

	// 					if (isLarge)
	// 						remove_large_node(dv, ltree_tid);
	// 					else 
	// 						remove_small_node(dv, stree_tid);
	// 				}
	// 				u5_time += (GetCurrentTimeSec() - u5_acc_time);

	// 				// Drop the size of subtree.
	// 				st_sizes[tid] -= subtract_size;
	// 				if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;

	// 				// 3. Convert shrink trees to small trees representataions.	
	// 				if (st_sizes[tid] > 0) {
	// 					if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
	// 						shrink_tree(tid, st_roots[tid], ltrees[tid]);
	// 					}
	// 				}

	// 				if (st_roots[tid] == pv)
	// 					if (isLarge) {
	// 						ltree_tid.clear();
	// 						ltree_tid.shrink_to_fit();
	// 						//ltrees.shrink_to_fit();
	// 					}
	// 					else {
	// 						stree_tid.clear();
	// 						stree_tid.resize(0);
	// 						//strees.shrink_to_fit();
	// 					}

	// 			}
	// 			//	u3_time += (GetCurrentTimeSec() - u3_acc_time);
	// 		}
	// 		else { // We use the reverse index to update trees.
	// 		   // Build the reverse_index first (one time).
	// 			if (reverse_index.size() == 0) {
	// 				reverse_index.resize(numOfVertices);
	// 				for (NodeID v = 0; v < numOfVertices; ++v) {
	// 					reverse_index[v].set_empty_key(numOfVertices + 1);
	// 					reverse_index[v].set_deleted_key(numOfVertices + 2);
	// 				}
	// 				for (int tid = 0; tid < current_trees_num; ++tid) {
	// 					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

	// 					bool isLarge = false;
	// 					int stree_idx = trees_pointers[tid];
	// 					if (stree_idx == -1)
	// 						isLarge = true;

	// 					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
	// 					large_tree& ltree_tid = ltrees[tid];

	// 					for (int v = 0; v < numOfVertices; ++v) {

	// 						int vpid = 0;
	// 						if (isLarge)
	// 							vpid = get_large_parent(v, ltree_tid);
	// 						else {
	// 							//double get_small_parent_incre_time = GetCurrentTimeSec();
	// 							vpid = get_small_parent(v, stree_tid);
	// 							//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
	// 						}
	// 						if ((vpid != -1 || st_roots[tid] == v)) {
	// 							reverse_index[v].insert(tid);
	// 							//in_tree_count[v]++;
	// 						}
	// 					}
	// 				}
	// 			}

	// 			//cout << reverse_index[pv].size() << "," << sum[pv] << "," <<  (double)sum[pv] / (double)reverse_index[pv].size() << endl;
	// 			// Access the alive trees directly by the reverse index.
	// 			//for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
	// 			for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
	// 				int tid = *it;
	// 				if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

	// 				bool isLarge = false;
	// 				int stree_idx = trees_pointers[tid];
	// 				if (stree_idx == -1)
	// 					isLarge = true;

	// 				small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
	// 				large_tree& ltree_tid = ltrees[tid];

	// 				double u6_acc_time = GetCurrentTimeSec();
	// 				/*
	// 						int parent_node;
	// 						if (isLarge)
	// 							parent_node = get_large_parent(pv, ltree_tid);
	// 						else {
	// 							//double get_small_parent_incre_time = GetCurrentTimeSec();
	// 							parent_node = get_small_parent(pv, stree_tid);
	// 							//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
	// 						}
	// 						if(parent_node == -1 && pv != st_roots[tid])
	// 							continue;
	// 					*/
	// 				double u7_acc_time = GetCurrentTimeSec();

	// 				//	// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
	// 				//	//double u4_acc_time = getcurrenttimesec();
	// 				//	get_descendants(tid, pv, descendants, descendants_parents,que, vis, wgraph.adj, islarge, stree_tid, ltree_tid);
	// 				////	u4_time += (getcurrenttimesec() - u4_acc_time);

	// 				//	get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

	// 					// Array Representation
	// 				get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
	// 				get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

	// 				NodeID subtract_size = descendants.size();
	// 				for (int ai = 0; ai < ascendants.size(); ++ai) {
	// 					NodeID av = ascendants[ai];
	// 					// source - av - pv: subtract the size of subtree rooted at pv.
	// 					if (pv == av) continue; // It will be subtracted in the descending search.
	// 					//cover[tid][av] -= subtract_size;
	// 					cnt16[av][tid%CNT] -= subtract_size;
	// 					//affected[av] = true;
	// 				//	sum[av] -= subtract_size;

	// 					if (affected_list[av] == false) {
	// 						affected_count++;
	// 						if (switch_set)
	// 							affected_v.insert(av);
	// 						affected_list[av] = true;
	// 					}
	// 				}
	// 				u7_time += (GetCurrentTimeSec() - u7_acc_time);

	// 				//double u2_acc_time = GetCurrentTimeSec();
	// 				// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
	// 				double u8_acc_time = GetCurrentTimeSec();
	// 				for (int di = descendants.size() - 1; di > -1; --di) {
	// 					NodeID dv = descendants[di];
	// 					// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
	// 					cover_tid[dv]++;
	// 					cnt16[dv][tid%CNT] -= cover_tid[dv];
	// 					//	sum[dv] -= cover_tid[dv];

	// 					if (descendants_parents[dv] != -1)
	// 						cover_tid[descendants_parents[dv]] += cover_tid[dv];
	// 					//	if (cnt16[dv][tid%CNT] < 0)
	// 					//		cout << "a6 " << cnt16[dv][tid%CNT] << endl;

	// 						//if(dv == pv){
	// 							/*if (isLarge)
	// 								remove_large_node(dv, ltree_tid);
	// 							else
	// 								remove_small_node(dv, stree_tid);*/
	// 								//}
	// 								//affected[dv] = true;

	// 								//double u82_acc_time = GetCurrentTimeSec();
	// 					if (dv != pv)
	// 						reverse_index[dv].erase(tid);
	// 					//if (dv != pv)
	// 					//	in_tree_count[dv]--;

	// 				//	double u81_acc_time = GetCurrentTimeSec();
	// 					if (affected_list[dv] == false) {
	// 						affected_count++;
	// 						if (switch_set)
	// 							affected_v.insert(dv);
	// 						affected_list[dv] = true;
	// 					}

	// 					//u81_time += (GetCurrentTimeSec() - u81_acc_time);

	// 				//	u82_time += (GetCurrentTimeSec() - u82_acc_time);
	// 				}
	// 				for (int di = descendants.size() - 1; di > -1; --di) {
	// 					NodeID dv = descendants[di];
	// 					cover_tid[dv] = 0;
	// 					if (isLarge)
	// 						remove_large_node(dv, ltree_tid);
	// 					else 
	// 						remove_small_node(dv, stree_tid);
	// 				}
	// 				u8_time += (GetCurrentTimeSec() - u8_acc_time);

	// 				if (st_roots[tid] == pv)
	// 					if (isLarge) {
	// 						ltree_tid.clear();
	// 						ltree_tid.shrink_to_fit();
	// 					}
	// 					else{
	// 						stree_tid.clear();
	// 						stree_tid.resize(0);
	// 					}


	// 				st_sizes[tid] -= subtract_size;
	// 				if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;
	// 				//	if (st_sizes[tid] < 0) {
	// 				//		cout << "a5" << st_sizes[tid] << endl;
	// 				//	}

	// 					// 3. Convert shrink trees to small trees representataions.				
	// 				if (st_sizes[tid] > 0) {
	// 					if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
	// 						shrink_tree(tid, st_roots[tid], ltrees[tid]);
	// 					}
	// 				}

	// 				//u2_time += (GetCurrentTimeSec() - u2_acc_time);
	// 			}
	// 			reverse_index[pv].clear();
	// 			reverse_index[pv].resize(0);

	// 		}
	// 		updating_time += (GetCurrentTimeSec() - updating_acc_time);

	// 		// 4. Sample new trees based on.
	// 		// If we sample a new tree, we need to dynamically increment the following data structures:
	// 		// 1) ltrees
	// 		// 2) trees_pointers (-1 or pointer to strees)
	// 		// 3) cover
	// 		// 4) st_roots
	// 		// 5) st_sizes 
	// 		// 6) strees (maybe)
	// 		// two constraints: cl > ct and 10kn total vertices
	// 		double adding_acc_time = GetCurrentTimeSec();


	// 		/*cout << "ranking:" << ranking << endl;
	// 		cout << "pick cover:" << pv_cover << endl;
	// 		cout << "pick count:" << picked_count << endl;
	// 		cout << "labeled_arcs:" << labeled_arcs << endl;
	// 		cout << "new_tree_arcs:" << new_tree_arcs << endl;
	// 		cout << "alive_resources:" << alive_resources << endl;
	// 		cout << "bounded_resources:" << bounded_resources << endl;
	// 		cout << "#tree:" << cover.size() << endl;
	// 		cout << "u5/adding" << u5_time / adding_time << endl;
	// 		*/
	// 		new_tree_arcs = 0;
	// 		int how_many = 0;
	// 		while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {

	// 			if (picked_count == numOfVertices) break;

	// 			NodeID r = rand() % numOfVertices;
	// 			if (picked_as_root[r] == true) continue;
	// 			if (trivial_root[r] == true) continue;

	// 			large_tree new_ltree;
	// 			small_tree new_stree;
	// 			//GOOGLE DENSE MAP
	// 			new_stree.set_empty_key(numOfVertices + 1);
	// 			new_stree.set_deleted_key(numOfVertices + 2);
	// 			bool isLarge = false;
	// 			int tmp_tree_arcs = 0;

	// 			if (!BUILD_SMALL_TREE_FLAG) {

	// 				//double u5_acc_time = GetCurrentTimeSec();
	// 				double u9_acc_time = GetCurrentTimeSec();
	// 				tmp_tree_arcs = build_large_tree_dij(r, new_ltree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
	// 				u9_time += (GetCurrentTimeSec() - u9_acc_time);
	// 				//u5_time += (GetCurrentTimeSec() - u5_acc_time);
	// 				new_tree_arcs += tmp_tree_arcs;

	// 				if (descendants.size() < size_thr)
	// 					BUILD_SMALL_TREE_FLAG = true;

	// 				if (descendants.size() < 100) {
	// 					//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
	// 					trivial_root[r] = true;
	// 					picked_count++;
	// 					continue; // we have sampled a trivial tree.
	// 				}

	// 				//	double u91_acc_time = GetCurrentTimeSec();

	// 				isLarge = true;
	// 				how_many++;
	// 				//cout << descendants.size() << endl;
	// 				picked_as_root[r] = true;
	// 				picked_count++;
	// 				used_as_cover++;
	// 				st_roots.push_back(r);
	// 				trees_pointers.push_back(-1);
	// 				st_sizes.push_back(descendants.size());

	// 				//new_tree_arcs += descendants.size();
	// 				alive_resources += (long long)numOfVertices;

	// 				//	u91_time += (GetCurrentTimeSec() - u91_acc_time);
	// 					// Turn large tree to small tree.
	// 					// 1. Build new small tree with ltree_i and add it to strees;
	// 					// 2. Set tree_pointers[i] to the position of the strees;
	// 					// 3. Clear the large tree ltree_i by setting its size to 0.

	// 				//	double u92_acc_time = GetCurrentTimeSec();
	// 				int i = st_roots.size() - 1;
	// 				if (st_sizes[i] < size_thr) {
	// 					shrink_tree(i, st_roots[i], new_ltree);
	// 					// FIX
	// 					isLarge = false;
	// 					/*for (int di = descendants.size() - 1; di > -1; --di) {
	// 						NodeID q = descendants[di];
	// 						NodeID v_parent = get_large_parent(q, new_ltree);
	// 						NodeID s_parent = get_small_parent(q, strees[strees.size() - 1]);
	// 						if (v_parent != s_parent) cout << v_parent << " vs " << s_parent << endl;
	// 					}*/
	// 				}
	// 				ltrees.push_back(new_ltree);
	// 				//	u92_time += (GetCurrentTimeSec() - u92_acc_time);
	// 			}
	// 			else {//Build Small Tree Directly
	// 				//double u5_acc_time = GetCurrentTimeSec();
	// 				double u9_acc_time = GetCurrentTimeSec();
	// 				tmp_tree_arcs = build_small_tree_dij(r, new_stree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
	// 				//u5_time += (GetCurrentTimeSec() - u5_acc_time);
	// 				new_tree_arcs += tmp_tree_arcs;
	// 				u9_time += (GetCurrentTimeSec() - u9_acc_time);

	// 				if (descendants.size() >= size_thr)
	// 					BUILD_SMALL_TREE_FLAG = false;

	// 				if (descendants.size() < 100) {
	// 					//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
	// 					trivial_root[r] = true;
	// 					picked_count++;
	// 					continue; // we have sampled a trivial tree.
	// 				}

	// 				isLarge = false;

	// 				how_many++;
	// 				used_as_cover++;
	// 				//cout << descendants.size() << endl;
	// 				picked_as_root[r] = true;
	// 				picked_count++;
	// 				st_roots.push_back(r);
	// 				trees_pointers.push_back(strees.size()); //we assume it is inserted.
	// 				st_sizes.push_back(descendants.size());

	// 				alive_resources += (long long)new_stree.size();
	// 				//strees.push_back(new_stree);
	// 				//// new_ltree is empty.
	// 				//ltrees.push_back(new_ltree);

	// 				int i = st_roots.size() - 1;
	// 				if (st_sizes[i] >= size_thr) {
	// 					enlarge_tree(i, st_roots[i], new_stree, new_ltree);
	// 					new_stree.clear();
	// 					new_stree.resize(0);
	// 					isLarge = true;
	// 				}
	// 				else {
	// 					strees.push_back(new_stree);
	// 					ltrees.push_back(new_ltree);
	// 				}

	// 			}

	// 			int i = st_roots.size() - 1;

	// 			int numOfStrees = strees.size();
	// 			small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

	// 			// Calculate the subtree size of each nodes.

	// 			double u10_acc_time = GetCurrentTimeSec();
	// 			for (int di = descendants.size() - 1; di > -1; --di) {
	// 				NodeID q = descendants[di];

	// 				++cover_tid[q]; // Increment because it hitsthe r-q path.
	// 				++cnt16[q][i%CNT];
	// 				//++sum[q];
	// 				//affected[q] = true;

	// 				if (affected_list[q] == false) {
	// 					affected_count++;
	// 					if (switch_set)
	// 						affected_v.insert(q);
	// 					affected_list[q] = true;
	// 				}

	// 				if (switch_small == true) {
	// 					reverse_index[q].insert(i);
	// 					//in_tree_count[q]++;
	// 				}


	// 				int parent_node;
	// 				if (isLarge)
	// 					parent_node = get_large_parent(q, new_ltree);
	// 				else {
	// 					//double get_small_parent_incre_time = GetCurrentTimeSec();
	// 					parent_node = get_small_parent(q, new_stree_ref);
	// 					//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
	// 				}
	// 				if (parent_node > -1) {
	// 					cover_tid[parent_node] += cover_tid[q];
	// 					cnt16[parent_node][i%CNT] += cover_tid[q];
	// 					//sum[parent_node] += cover_new[q];
	// 				}
	// 			}
	// 			u10_time += (GetCurrentTimeSec() - u10_acc_time);
	// 			for (int di = descendants.size() - 1; di > -1; --di) {
	// 				NodeID q = descendants[di];
	// 				cover_tid[q] = 0;
	// 			}
	// 			//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;

	// 		}
	// 		adding_time += (GetCurrentTimeSec() - adding_acc_time);


	// 		//double u3_acc_time = GetCurrentTimeSec();

	// 		updating_acc_time = GetCurrentTimeSec();

	// 		// Build the cover value heap for one time.
	// 		if (switch_small == true && value_heap.empty())
	// 			for (int nid = 0; nid < numOfVertices; ++nid)
	// 				if (has_picked[nid] == false)
	// 					value_heap.update(nid, -picked_value[nid]);

	// 		//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);
	// 		//c	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);

	// 		double u11_acc_time = GetCurrentTimeSec();
	// 		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
	// 		u11_time += (GetCurrentTimeSec() - u11_acc_time);
	// 		//	iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
	// 		updating_time += (GetCurrentTimeSec() - updating_acc_time);

	// 		if (affected_count < numOfVertices / 8)
	// 			switch_set = true;
	// 		else
	// 			switch_set = false;
	// 		affected_count = 0;

	// 		time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
	// 		iteration_time[ranking] = time_this_iteration;

	// 	}

	// 	cout << "picked as tree covers:" << used_as_cover << endl;
	// 	cout << "picked count:" << picked_count << endl;
		
	// 	//// FIX
	// 	//for (int v = 0; v < numOfVertices; ++v){
	// 	//	for (int i = 0; i < CNT; ++i) {
	// 	//		if (cnt16[v][i] != 0)
	// 	//			cout << "mi?" << endl;
	// 	//	}
	// 	//}
	// 	if (stop_flag == true) return;

	// 	double labeling_acc_time = GetCurrentTimeSec();
	// 	// Order the remaining vertices with degree, if any.
	// 	if (ranking != numOfVertices) {
	// 		vector<NodeID> degree(numOfVertices);
	// 		vector<NodeID> r_degree(numOfVertices);
	// 		for (NodeID v = 0; v < numOfVertices; ++v) {
	// 			if (has_picked[v] == true) continue;
	// 			/*vector<NodeID>& adj_v = graph.adj[v];
	// 			for (NodeID i = 0; i < adj_v.size(); ++i) {
	// 			if (has_picked[adj_v[i]] != true)
	// 			degree[v]++;
	// 			}*/

	// 			for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
	// 				NodeID w = wgraph.edges[eid].first;
	// 				if (has_picked[w] != true)
	// 					degree[v]++;
	// 			}

	// 		}
			
	// 		vector<pair<NodeID, NodeID> > degree_v;
	// 		degree_v.reserve(numOfVertices - ranking);
	// 		for (NodeID v = 0; v < numOfVertices; ++v) {
	// 			if (has_picked[v] == true) continue;
	// 			degree_v.push_back(make_pair(degree[v], v));
	// 		}
	// 		sort(degree_v.rbegin(), degree_v.rend());

	// 		int count = 0;
	// 		for (; ranking < numOfVertices; ++ranking) {
	// 			NodeID tpv = degree_v[count++].second;
	// 			inv[ranking] = tpv;
	// 			rank[tpv] = ranking;
	// 			//output betweenness
	// 			//ofs<<" degree ranking ="<<ranking <<" pv="<<tpv<<" pv_cover="<< degree_v[count++].first<<" picked_count="<<picked_count<<endl;
	// 			labeling_source_dij(tpv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
	// 		}
	// 	}

	// 	labeling_time += (GetCurrentTimeSec() - labeling_acc_time);

	// 	for (NodeID v = 0; v < numOfVertices; ++v) {

    //         NodeID k = tmp_idx[v].first.size();
    //         labels.index_[v].spt_v.resize(k);
    //         labels.index_[v].spt_d.resize(k);
    //         for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
    //         for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
    //         tmp_idx[v].first.clear();
    //         tmp_idx[v].second.clear();
    //         tmp_idx[v].first.shrink_to_fit();
    //         tmp_idx[v].second.shrink_to_fit();

	// 	}

	// 	num_of_trees = st_roots.size();

	// 	double small_tree_num = 0;

	// 	for (int i = 0; i < num_of_trees; ++i) {
	// 		if (trees_pointers[i] != -1) {
	// 			small_tree_num++;
	// 			//small_tree_num += strees[trees_pointers[i]].size();
	// 		}
	// 	}
		
	// 	total_resources = small_tree_num;
	// 	//output time
	// 	start_time=GetCurrentTimeSec()-start_time;
	// 	_hub_time=labeling_time;
	// 	_ordering_time=start_time-_hub_time;
	// }

	Betweenness_Ordering(NodeID k, double beta, WGraph& wgraph, NodeID stopv,double& _ordering_time, double& _hub_time, vector<unsigned int>& betweenness_values,long long bound_times = 10, long long trans_times = 8) {
		double start_time = GetCurrentTimeSec();

		betweenness_values.resize(numOfVertices,0);

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;

		bool isOutputBet=false; //indict whether to output betwenness modified by wanwan
		ofstream bet_ofs; 

		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		if (DIRECTED_FLAG == false) {
			labels.index_.resize(numOfVertices);
		}
		else { //if (DIRECTED_FLAG == true) {
			dlabels.index_.resize(numOfVertices);
			dlabels.bindex_.resize(numOfVertices);
		}

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = bound_times * (long long)k * (long long)numOfVertices;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / trans_times; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.

		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
		//	srand((unsigned)time(NULL));
		srand(0x728f9d0);


		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		//strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its tree representations. large tree [0, V) small tree [V, 2V).
		//vector<unordered_map<int, int> > cover(k, unordered_map<int, int>());

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		bool stop_flag = false;


		// tmp structure used in the dijkstra searches. 

		// Binary Heap
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> usd(numOfVertices, false);

		//GOOGLE
		//		unordered_set<NodeID> affected_v;
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<cover_value_type> cover_tid(numOfVertices, 0);
		//vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14


		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<unordered_set<int> > reverse_index(0);
		//GOOGLE DENSE SET
		vector<google::dense_hash_set<NodeID> > reverse_index;
		//vector<NodeID> in_tree_count(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

		vector<NodeID> descendants_parents(numOfVertices, -1);

		// Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			build_large_tree_dij(r, ltree_i, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			// Calculate the subtree size of each nodes.
			//	unordered_map<int, int>& cover_i = cover[i];
			//vector<int> cover_i(numOfVertices, 0);

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hits the r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];

				int parent_node;

				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);

				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_i[q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
		}



		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		//bool switch_small = false;
		switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
       //		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;

		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;

		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
			if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;

			double time_this_iteration = GetCurrentTimeSec();

			NodeID pv = 0;
			long long pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}

			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();

				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					//cout << "u3 time:" << u3_time << endl;
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);
			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();

				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;
				u2_time += (GetCurrentTimeSec() - u2_acc_time);
				/*
				int max1 = 0;
				int max2 = 0;
				int max1v = -1;
				int max2v = -1;
				for(int  i = 0; i < CNT; ++i){
					if(max2 <= cnt16[pv][i]){
						max2 = cnt16[pv][i];
						max2v = i;
					}
					if(max1 <= max2){
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
						int tmpv = max1v;
						max1v = max2v;
						max2v = max1v;
					}
				}
				long long sumstd = 0;
				long long summean = 0;
				for(int  i = 0; i < CNT; ++i){
					if(i != max1v && i != max2v){
						sumstd += (long long)cnt16[pv][i] * (long long)cnt16[pv][i];
						summean += (long long)cnt16[pv][i];
					}
				}
				sumstd = sumstd/(long long)(CNT-2);
				summean = summean/(long long)(CNT-2);
				cout << sqrt(sumstd - summean * summean) << endl;*/
			}
			//output betweenness
			//if(isOutputBet) bet_ofs<<pv<<" "<<pv_cover<<endl;
			betweenness_values[pv]=pv_cover;
			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size() * 10) {//sum[pv] <= in_tree_count[pv] * 10){//
			//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
			//	break;
			//}

			//if (has_picked[pv] == true)
			//	cout << "wrong1" << endl;

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_dij(pv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
				//	double u3_acc_time = GetCurrentTimeSec();
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// Get the descendants and ascendants of pv in tree tid.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					double u3_acc_time = GetCurrentTimeSec();

					// //1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph.adj, isLarge, stree_tid, ltree_tid);
					//get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					//Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					u3_time += (GetCurrentTimeSec() - u3_acc_time);


					double u4_acc_time = GetCurrentTimeSec();
					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
					//	cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//affected[av] = true;

						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}
					u4_time += (GetCurrentTimeSec() - u4_acc_time);


					double u5_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						/*if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);*/

							//						remove_node(tid, dv);
													//affected[dv] = true;

						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
					}
					u5_time += (GetCurrentTimeSec() - u5_acc_time);

					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;

					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
				//	u3_time += (GetCurrentTimeSec() - u3_acc_time);
			}
			else { // We use the reverse index to update trees.
			   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {
						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								vpid = get_small_parent(v, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if ((vpid != -1 || st_roots[tid] == v)) {
								reverse_index[v].insert(tid);
								//in_tree_count[v]++;
							}
						}
					}
				}

				//cout << reverse_index[pv].size() << "," << sum[pv] << "," <<  (double)sum[pv] / (double)reverse_index[pv].size() << endl;
				// Access the alive trees directly by the reverse index.
				//for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					double u6_acc_time = GetCurrentTimeSec();
					/*
							int parent_node;
							if (isLarge)
								parent_node = get_large_parent(pv, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								parent_node = get_small_parent(pv, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if(parent_node == -1 && pv != st_roots[tid])
								continue;
						*/
					double u7_acc_time = GetCurrentTimeSec();

					//	// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//	//double u4_acc_time = getcurrenttimesec();
					//	get_descendants(tid, pv, descendants, descendants_parents,que, vis, wgraph.adj, islarge, stree_tid, ltree_tid);
					////	u4_time += (getcurrenttimesec() - u4_acc_time);

					//	get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

						// Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
						//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//affected[av] = true;
					//	sum[av] -= subtract_size;

						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}
					u7_time += (GetCurrentTimeSec() - u7_acc_time);

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					double u8_acc_time = GetCurrentTimeSec();
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						//	if (cnt16[dv][tid%CNT] < 0)
						//		cout << "a6 " << cnt16[dv][tid%CNT] << endl;

							//if(dv == pv){
								/*if (isLarge)
									remove_large_node(dv, ltree_tid);
								else
									remove_small_node(dv, stree_tid);*/
									//}
									//affected[dv] = true;

									//double u82_acc_time = GetCurrentTimeSec();
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//if (dv != pv)
						//	in_tree_count[dv]--;

					//	double u81_acc_time = GetCurrentTimeSec();
						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}

						//u81_time += (GetCurrentTimeSec() - u81_acc_time);

					//	u82_time += (GetCurrentTimeSec() - u82_acc_time);
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
					}
					u8_time += (GetCurrentTimeSec() - u8_acc_time);

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else{
							stree_tid.clear();
							stree_tid.resize(0);
						}


					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;
					//	if (st_sizes[tid] < 0) {
					//		cout << "a5" << st_sizes[tid] << endl;
					//	}

						// 3. Convert shrink trees to small trees representataions.				
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					//u2_time += (GetCurrentTimeSec() - u2_acc_time);
				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);

			}
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();


			/*cout << "ranking:" << ranking << endl;
			cout << "pick cover:" << pv_cover << endl;
			cout << "pick count:" << picked_count << endl;
			cout << "labeled_arcs:" << labeled_arcs << endl;
			cout << "new_tree_arcs:" << new_tree_arcs << endl;
			cout << "alive_resources:" << alive_resources << endl;
			cout << "bounded_resources:" << bounded_resources << endl;
			cout << "#tree:" << cover.size() << endl;
			cout << "u5/adding" << u5_time / adding_time << endl;
			*/
			new_tree_arcs = 0;
			int how_many = 0;
			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {

				if (picked_count == numOfVertices) break;

				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);
				bool isLarge = false;
				int tmp_tree_arcs = 0;

				if (!BUILD_SMALL_TREE_FLAG) {

					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_large_tree_dij(r, new_ltree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					u9_time += (GetCurrentTimeSec() - u9_acc_time);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() < size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					//	double u91_acc_time = GetCurrentTimeSec();

					isLarge = true;
					how_many++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					used_as_cover++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += (long long)numOfVertices;

					//	u91_time += (GetCurrentTimeSec() - u91_acc_time);
						// Turn large tree to small tree.
						// 1. Build new small tree with ltree_i and add it to strees;
						// 2. Set tree_pointers[i] to the position of the strees;
						// 3. Clear the large tree ltree_i by setting its size to 0.

					//	double u92_acc_time = GetCurrentTimeSec();
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						// FIX
						isLarge = false;
						/*for (int di = descendants.size() - 1; di > -1; --di) {
							NodeID q = descendants[di];
							NodeID v_parent = get_large_parent(q, new_ltree);
							NodeID s_parent = get_small_parent(q, strees[strees.size() - 1]);
							if (v_parent != s_parent) cout << v_parent << " vs " << s_parent << endl;
						}*/
					}
					ltrees.push_back(new_ltree);
					//	u92_time += (GetCurrentTimeSec() - u92_acc_time);
				}
				else {//Build Small Tree Directly
					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_small_tree_dij(r, new_stree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;
					u9_time += (GetCurrentTimeSec() - u9_acc_time);

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += (long long)new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}

				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.

				double u10_acc_time = GetCurrentTimeSec();
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path.
					++cnt16[q][i%CNT];
					//++sum[q];
					//affected[q] = true;

					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
						//in_tree_count[q]++;
					}


					int parent_node;
					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						parent_node = get_small_parent(q, new_stree_ref);
						//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//sum[parent_node] += cover_new[q];
					}
				}
				u10_time += (GetCurrentTimeSec() - u10_acc_time);
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}
				//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;

			}
			adding_time += (GetCurrentTimeSec() - adding_acc_time);


			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);
			//c	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);

			double u11_acc_time = GetCurrentTimeSec();
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			u11_time += (GetCurrentTimeSec() - u11_acc_time);
			//	iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;

		}

		cout << "picked as tree covers:" << used_as_cover << endl;
		cout << "picked count:" << picked_count << endl;
		
		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}
		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
					NodeID w = wgraph.edges[eid].first;
					if (has_picked[w] != true)
						degree[v]++;
				}

				if (DIRECTED_FLAG == true) {
					/*vector<NodeID>& r_adj_v = graph.r_adj[v];
					for (NodeID i = 0; i < r_adj_v.size(); ++i) {
					if (has_picked[r_adj_v[i]] != true)
					r_degree[v]++;
					}*/
					for (EdgeID eid = wgraph.r_vertices[v]; eid < wgraph.r_vertices[v + 1]; ++eid) {
						NodeID w = wgraph.r_edges[eid].first;
						if (has_picked[w] != true)
							r_degree[v]++;
					}
				}
			}
			
			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				if (DIRECTED_FLAG == true)
					degree_v.push_back(make_pair(degree[v] * r_degree[v], v));
				else
					degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			std::cout<<"==================push by degree==================="<<std::endl;
			std::cout<<"Num of vertices ordered by degree is "<<numOfVertices - ranking<<std::endl;//debug by wanjingyi
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				//output betweenness
				//ofs<<" degree ranking ="<<ranking <<" pv="<<tpv<<" pv_cover="<< degree_v[count++].first<<" picked_count="<<picked_count<<endl;
				labeling_source_dij(tpv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			}
		}
		//output time
		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);

		// Finish the constructions.
		for (NodeID v = 0; v < numOfVertices; ++v) {
			if (DIRECTED_FLAG == false) {
				NodeID k = tmp_idx[v].first.size();
				labels.index_[v].spt_v.resize(k);
				labels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

			}
			else {
				NodeID k = tmp_idx[v].first.size();
				dlabels.index_[v].spt_v.resize(k);
				dlabels.index_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_v[i] = tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.index_[v].spt_d[i] = tmp_idx[v].second[i];
				tmp_idx[v].first.clear();
				tmp_idx[v].second.clear();
				tmp_idx[v].first.shrink_to_fit();
				tmp_idx[v].second.shrink_to_fit();

				k = r_tmp_idx[v].first.size();
				dlabels.bindex_[v].spt_v.resize(k);
				dlabels.bindex_[v].spt_d.resize(k);
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_v[i] = r_tmp_idx[v].first[i];
				for (NodeID i = 0; i < k; ++i) dlabels.bindex_[v].spt_d[i] = r_tmp_idx[v].second[i];
				r_tmp_idx[v].first.clear();
				r_tmp_idx[v].second.clear();
				r_tmp_idx[v].first.shrink_to_fit();
				r_tmp_idx[v].second.shrink_to_fit();
			}
		}
		

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}
		
		total_resources = small_tree_num;
		//output time
		_hub_time=labeling_time;
		start_time=GetCurrentTimeSec()-start_time;
		_ordering_time=start_time-_hub_time;
	}

	Betweenness_Ordering(NodeID k, double beta, WGraph& wgraph, NodeID stopv,double& _ordering_time, double& _hub_time, long long bound_times = 10, long long trans_times = 8) {
		double start_time = GetCurrentTimeSec();

		if (k < CNT) {
			cout << " Should pick at least 16 seed!" << endl;
			return;
		}
		iteration_time.resize(numOfVertices);

		updating_time = 0;
		labeling_time = 0;
		adding_time = 0;
		selecting_time = 0;
		init_time = 0;
		u1_time = 0;
		u2_time = 0;
		u3_time = 0;
		num_of_trees = 0;
		total_resources = 0;


		inv.resize(numOfVertices, -1);
		rank.resize(numOfVertices, -1);

		labels.index_.resize(numOfVertices);

		// Resource Bounding.
		NodeID new_tree_arcs = 0;
		NodeID labeled_arcs = 0;
		alive_resources = 0;
		bounded_resources = bound_times * (long long)k * (long long)numOfVertices;

		vector<NodeID> st_roots; // The root vertices of the sampled trees.
		vector<NodeID> st_sizes; // The size of the sampled trees.

		NodeID size_thr = numOfVertices < 8 ? numOfVertices : numOfVertices / trans_times; // When a tree shrink to a size smaller than the size_thr, it should be converted into a small_tree representation.

		init_time = GetCurrentTimeSec();
		// Building the initial k trees.
		//	srand((unsigned)time(NULL));
		srand(0x728f9d0);


		ltrees.reserve(numOfVertices);
		strees.reserve(numOfVertices);
		trees_pointers.reserve(numOfVertices);
		st_roots.reserve(numOfVertices);
		st_sizes.reserve(numOfVertices);
		st_roots.resize(k);
		st_sizes.resize(k);

		ltrees.resize(k);
		//strees.clear();
		trees_pointers.resize(k, -1); // Point the sampling trees to its tree representations. large tree [0, V) small tree [V, 2V).
		//vector<unordered_map<int, int> > cover(k, unordered_map<int, int>());

		vector<NodeID> ascendants;
		ascendants.reserve(numOfVertices);
		vector<NodeID> descendants;
		descendants.reserve(numOfVertices);
		vector<NodeID> que(numOfVertices);
		vector<bool> vis(numOfVertices);
		vector<bool> has_picked(numOfVertices);
		vector<vector<cover_value_type> > cnt16(numOfVertices, vector<cover_value_type>(CNT));
		//		vector<bool> affected(V, true);
		vector<bool> picked_as_root(numOfVertices, false);
		vector<bool> trivial_root(numOfVertices, false);
		int picked_count = 0;
		int used_as_cover = 0;
		vector<cover_value_type> sum(numOfVertices);
		bool stop_flag = false;


		// tmp structure used in the dijkstra searches. 

		// Binary Heap
		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
		vector<bool> usd(numOfVertices, false);

		//GOOGLE
		//		unordered_set<NodeID> affected_v;
		google::dense_hash_set<NodeID> affected_v;
		affected_v.set_empty_key(numOfVertices + 1);
		affected_v.set_deleted_key(numOfVertices + 2);

		vector<cover_value_type> cover_tid(numOfVertices, 0);
		//vector<cover_value_type> cover_new(numOfVertices, 0); // modified-5-14


		vector<bool> affected_list(numOfVertices, false);
		int affected_count = 0;
		bool switch_set = false;

		//vector<unordered_set<int> > reverse_index(0);
		//GOOGLE DENSE SET
		vector<google::dense_hash_set<NodeID> > reverse_index;
		//vector<NodeID> in_tree_count(numOfVertices);

		// Preparing basic structure for pl algorithm.
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT))); // Backward labels.

		vector<NodeID> descendants_parents(numOfVertices, -1);

		// Building the initial k sampled trees. 
		for (int i = 0; i < k; ++i) {
			NodeID r = rand() % numOfVertices;
			if (picked_as_root[r] == true) { --i; continue; }

			st_roots[i] = r;
			large_tree &ltree_i = ltrees[i];
			picked_as_root[r] = true;
			picked_count++;
			used_as_cover++;

			build_large_tree_dij(r, ltree_i, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, 0, numOfVertices);
			st_sizes[i] = descendants.size();
			alive_resources += (long long)numOfVertices;

			// Turn large tree to small tree.
			// 1. Build new small tree with ltree_i and add it to strees;
			// 2. Set tree_pointers[i] to the position of the strees;
			// 3. Clear the large tree ltree_i by setting its size to 0.
			if (st_sizes[i] < size_thr)
				shrink_tree(i, st_roots[i], ltrees[i]);

			bool isLarge = false;
			int stree_idx = trees_pointers[i];
			if (stree_idx == -1)
				isLarge = true;

			small_tree& stree_i = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;

			// Calculate the subtree size of each nodes.
			//	unordered_map<int, int>& cover_i = cover[i];
			//vector<int> cover_i(numOfVertices, 0);

			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];

				++cover_tid[q]; // Increment because it hits the r-q path.
				++cnt16[q][i%CNT];
				//++sum[q];

				int parent_node;

				if (isLarge)
					parent_node = get_large_parent(q, ltree_i);
				else
					parent_node = get_small_parent(q, stree_i);

				if (parent_node > -1) {
					cover_tid[parent_node] += cover_tid[q];
					cnt16[parent_node][i%CNT] += cover_tid[q];
					//	sum[parent_node] += cover_i[q];
				}
			}
			for (int di = descendants.size() - 1; di > -1; --di) {
				NodeID q = descendants[di];
				cover_tid[q] = 0;
			}
		}



		//cnt16_aggreate(cnt16, has_picked, cover);

		vector<cover_value_type> picked_value;
		benchmark::heap<2, cover_value_type, NodeID> value_heap(numOfVertices);
		//bool switch_small = false;
		switch_small = false;

		//cnt16_2max(picked_value, sum, has_picked,affected, cnt16, switch_small, value_heap);
       //		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);
		cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);

		init_time = GetCurrentTimeSec() - init_time;

		cout << "Finish initial sampling " << init_time << "secs." << endl;

		//for (NodeID v = 0; v < V; ++v)	value_heap.update(v, -picked_value[v]); // Min-heap.

		double max_value = -1;


		//output betweenness
		//ofstream ofs("../dataset/Manhattan/order_values/original.betweenness");
		//if(!ofs.is_open()) std::cout<<"Betweenness.values file can't be opened!"<<endl;
		// Perform V iterations to rank all vertices.
		NodeID ranking = 0;
		for (; ranking < numOfVertices; ++ranking) {
			if (ranking % (numOfVertices / 10) == 0)
				cout << "picking the " << ranking << "-th orders." << endl;

			double time_this_iteration = GetCurrentTimeSec();

			NodeID pv = 0;
			long long pv_cover = -1;

			if (ranking == stopv) {
				stop_flag = true;
				cout << "topv labeling time:" << labeling_time << endl;
				break;
			}

			double selecting_acc_time = GetCurrentTimeSec();
			// In the beginning, pick the best cover vertex by a linear scan.
			if (switch_small == false) {
				double u1_acc_time = GetCurrentTimeSec();

				for (NodeID v = 0; v < numOfVertices; ++v) {
					if (has_picked[v] == true) continue;
					if (pv_cover == picked_value[v]) {
						if (sum[pv] < sum[v]) {
							pv_cover = picked_value[v];
							pv = v;
						}
					}
					else if (pv_cover < picked_value[v]) {
						pv_cover = picked_value[v];
						pv = v;
					}
				}

				if (sum[pv] < numOfVertices / 8 && switch_small == false) {
					//cout << "u3 time:" << u3_time << endl;
					cout << ranking << endl;
					cout << pv_cover << endl;
					cout << sum[pv] << endl;
					cout << picked_count << endl;
					cout << "Switch to max-heap representations." << endl;
					switch_small = true;
				}
				u1_time += (GetCurrentTimeSec() - u1_acc_time);
			}
			else {// Later, utilizing the max-heap directly
				double u2_acc_time = GetCurrentTimeSec();

				value_heap.extract_min(pv, pv_cover);
				pv_cover = -pv_cover;
				u2_time += (GetCurrentTimeSec() - u2_acc_time);
				/*
				int max1 = 0;
				int max2 = 0;
				int max1v = -1;
				int max2v = -1;
				for(int  i = 0; i < CNT; ++i){
					if(max2 <= cnt16[pv][i]){
						max2 = cnt16[pv][i];
						max2v = i;
					}
					if(max1 <= max2){
						int tmp = max1;
						max1 = max2;
						max2 = tmp;
						int tmpv = max1v;
						max1v = max2v;
						max2v = max1v;
					}
				}
				long long sumstd = 0;
				long long summean = 0;
				for(int  i = 0; i < CNT; ++i){
					if(i != max1v && i != max2v){
						sumstd += (long long)cnt16[pv][i] * (long long)cnt16[pv][i];
						summean += (long long)cnt16[pv][i];
					}
				}
				sumstd = sumstd/(long long)(CNT-2);
				summean = summean/(long long)(CNT-2);
				cout << sqrt(sumstd - summean * summean) << endl;*/
			}
			//output betweenness
			//ofs<<pv<<" "<<pv_cover<<endl;
			selecting_time += (GetCurrentTimeSec() - selecting_acc_time);

			//// FIX
			//if (switch_small == true && reverse_index.size() > 0 && sum[pv] <= reverse_index[pv].size() * 10) {//sum[pv] <= in_tree_count[pv] * 10){//
			//	cout << "I think we can stop here " << GetCurrentTimeSec() - start_time << " ranking:" << ranking << endl;
			//	break;
			//}

			//if (has_picked[pv] == true)
			//	cout << "wrong1" << endl;

			has_picked[pv] = true;
			double labeling_acc_time = GetCurrentTimeSec();
			labeled_arcs = labeling_source_dij(pv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			labeling_time += (GetCurrentTimeSec() - labeling_acc_time);
			inv[ranking] = pv;
			rank[pv] = ranking;

			int current_trees_num = st_roots.size();

			// 4 steps in updating:
			double updating_acc_time = GetCurrentTimeSec();
			if (switch_small == false) { // We still iterate all alive trees to calculate the cover.	
				//	double u3_acc_time = GetCurrentTimeSec();
				for (int tid = 0; tid < current_trees_num; ++tid) {

					// Get the descendants and ascendants of pv in tree tid.
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					// pv has already been covered in tree tid.
					int pvpid = 0;
					if (isLarge)
						pvpid = get_large_parent(pv, ltree_tid);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						pvpid = get_small_parent(pv, stree_tid);
						//	get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (pvpid == -1 && st_roots[tid] != pv) // pv has been covered in tree tid 
						continue;


					double u3_acc_time = GetCurrentTimeSec();

					// //1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;				 
					//get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph.adj, isLarge, stree_tid, ltree_tid);
					//get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

					//Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					u3_time += (GetCurrentTimeSec() - u3_acc_time);


					double u4_acc_time = GetCurrentTimeSec();
					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
					//	cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//sum[av] -= subtract_size;
						//affected[av] = true;

						if (affected_list[av] == false) {
							affected_list[av] = true;
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
						}
					}
					u4_time += (GetCurrentTimeSec() - u4_acc_time);


					double u5_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					//vector<int> cover_tid(numOfVertices, 0);
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						if (cnt16[dv][tid%CNT] < 0)
							cout << "a6 " << cnt16[dv][tid%CNT] << endl;

						/*if (isLarge)
							remove_large_node(dv, ltree_tid);
						else
							remove_small_node(dv, stree_tid);*/

							//						remove_node(tid, dv);
													//affected[dv] = true;

						if (affected_list[dv] == false) {
							affected_count++;
							affected_list[dv] = true;
							if (switch_set)
								affected_v.insert(dv);
						}
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;

						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
					}
					u5_time += (GetCurrentTimeSec() - u5_acc_time);

					// Drop the size of subtree.
					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;

					// 3. Convert shrink trees to small trees representataions.	
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
							//ltrees.shrink_to_fit();
						}
						else {
							stree_tid.clear();
							stree_tid.resize(0);
							//strees.shrink_to_fit();
						}

				}
				//	u3_time += (GetCurrentTimeSec() - u3_acc_time);
			}
			else { // We use the reverse index to update trees.
			   // Build the reverse_index first (one time).
				if (reverse_index.size() == 0) {
					reverse_index.resize(numOfVertices);
					for (NodeID v = 0; v < numOfVertices; ++v) {
						reverse_index[v].set_empty_key(numOfVertices + 1);
						reverse_index[v].set_deleted_key(numOfVertices + 2);
					}
					for (int tid = 0; tid < current_trees_num; ++tid) {
						if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue; //This whole tree has already been covered. 

						bool isLarge = false;
						int stree_idx = trees_pointers[tid];
						if (stree_idx == -1)
							isLarge = true;

						small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
						large_tree& ltree_tid = ltrees[tid];

						for (int v = 0; v < numOfVertices; ++v) {

							int vpid = 0;
							if (isLarge)
								vpid = get_large_parent(v, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								vpid = get_small_parent(v, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if ((vpid != -1 || st_roots[tid] == v)) {
								reverse_index[v].insert(tid);
								//in_tree_count[v]++;
							}
						}
					}
				}

				//cout << reverse_index[pv].size() << "," << sum[pv] << "," <<  (double)sum[pv] / (double)reverse_index[pv].size() << endl;
				// Access the alive trees directly by the reverse index.
				//for (unordered_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
				for (google::dense_hash_set<int>::iterator it = reverse_index[pv].begin(); it != reverse_index[pv].end(); ++it) {
					int tid = *it;
					if (has_picked[st_roots[tid]] == true && st_roots[tid] != pv) continue;

					bool isLarge = false;
					int stree_idx = trees_pointers[tid];
					if (stree_idx == -1)
						isLarge = true;

					small_tree& stree_tid = (stree_idx != -1) ? strees[stree_idx] : empty_small_tree;
					large_tree& ltree_tid = ltrees[tid];

					double u6_acc_time = GetCurrentTimeSec();
					/*
							int parent_node;
							if (isLarge)
								parent_node = get_large_parent(pv, ltree_tid);
							else {
								//double get_small_parent_incre_time = GetCurrentTimeSec();
								parent_node = get_small_parent(pv, stree_tid);
								//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
							}
							if(parent_node == -1 && pv != st_roots[tid])
								continue;
						*/
					double u7_acc_time = GetCurrentTimeSec();

					//	// 1. Subtract the size of descendants of pv in tid for all ascendants of pv in the same tree;
					//	//double u4_acc_time = getcurrenttimesec();
					//	get_descendants(tid, pv, descendants, descendants_parents,que, vis, wgraph.adj, islarge, stree_tid, ltree_tid);
					////	u4_time += (getcurrenttimesec() - u4_acc_time);

					//	get_ascendants(tid, pv, ascendants, que, vis, wgraph.adj, wgraph.r_adj, isLarge, stree_tid, ltree_tid);

						// Array Representation
					get_descendants(tid, pv, descendants, descendants_parents, que, vis, wgraph, isLarge, stree_tid, ltree_tid);
					get_ascendants(tid, pv, ascendants, que, vis, wgraph, isLarge, stree_tid, ltree_tid);

					NodeID subtract_size = descendants.size();
					for (int ai = 0; ai < ascendants.size(); ++ai) {
						NodeID av = ascendants[ai];
						// source - av - pv: subtract the size of subtree rooted at pv.
						if (pv == av) continue; // It will be subtracted in the descending search.
						//cover[tid][av] -= subtract_size;
						cnt16[av][tid%CNT] -= subtract_size;
						//affected[av] = true;
					//	sum[av] -= subtract_size;

						if (affected_list[av] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(av);
							affected_list[av] = true;
						}
					}
					u7_time += (GetCurrentTimeSec() - u7_acc_time);

					//double u2_acc_time = GetCurrentTimeSec();
					// 2. Set all parent pointers of dv to -1(NULL), cover of dv to 0, and remove the subtree rooted at pv.
					double u8_acc_time = GetCurrentTimeSec();
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						// source - pv - dv: set cover of dv to 0(because they have been covered by pv).
						cover_tid[dv]++;
						cnt16[dv][tid%CNT] -= cover_tid[dv];
						//	sum[dv] -= cover_tid[dv];

						if (descendants_parents[dv] != -1)
							cover_tid[descendants_parents[dv]] += cover_tid[dv];
						//	if (cnt16[dv][tid%CNT] < 0)
						//		cout << "a6 " << cnt16[dv][tid%CNT] << endl;

							//if(dv == pv){
								/*if (isLarge)
									remove_large_node(dv, ltree_tid);
								else
									remove_small_node(dv, stree_tid);*/
									//}
									//affected[dv] = true;

									//double u82_acc_time = GetCurrentTimeSec();
						if (dv != pv)
							reverse_index[dv].erase(tid);
						//if (dv != pv)
						//	in_tree_count[dv]--;

					//	double u81_acc_time = GetCurrentTimeSec();
						if (affected_list[dv] == false) {
							affected_count++;
							if (switch_set)
								affected_v.insert(dv);
							affected_list[dv] = true;
						}

						//u81_time += (GetCurrentTimeSec() - u81_acc_time);

					//	u82_time += (GetCurrentTimeSec() - u82_acc_time);
					}
					for (int di = descendants.size() - 1; di > -1; --di) {
						NodeID dv = descendants[di];
						cover_tid[dv] = 0;
						if (isLarge)
							remove_large_node(dv, ltree_tid);
						else 
							remove_small_node(dv, stree_tid);
					}
					u8_time += (GetCurrentTimeSec() - u8_acc_time);

					if (st_roots[tid] == pv)
						if (isLarge) {
							ltree_tid.clear();
							ltree_tid.shrink_to_fit();
						}
						else{
							stree_tid.clear();
							stree_tid.resize(0);
						}


					st_sizes[tid] -= subtract_size;
					if (trees_pointers[tid] != -1) alive_resources -= (long long)subtract_size;
					//	if (st_sizes[tid] < 0) {
					//		cout << "a5" << st_sizes[tid] << endl;
					//	}

						// 3. Convert shrink trees to small trees representataions.				
					if (st_sizes[tid] > 0) {
						if (trees_pointers[tid] == -1 && st_sizes[tid] < size_thr) {
							shrink_tree(tid, st_roots[tid], ltrees[tid]);
						}
					}

					//u2_time += (GetCurrentTimeSec() - u2_acc_time);
				}
				reverse_index[pv].clear();
				reverse_index[pv].resize(0);

			}
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			// 4. Sample new trees based on.
			// If we sample a new tree, we need to dynamically increment the following data structures:
			// 1) ltrees
			// 2) trees_pointers (-1 or pointer to strees)
			// 3) cover
			// 4) st_roots
			// 5) st_sizes 
			// 6) strees (maybe)
			// two constraints: cl > ct and 10kn total vertices
			double adding_acc_time = GetCurrentTimeSec();


			/*cout << "ranking:" << ranking << endl;
			cout << "pick cover:" << pv_cover << endl;
			cout << "pick count:" << picked_count << endl;
			cout << "labeled_arcs:" << labeled_arcs << endl;
			cout << "new_tree_arcs:" << new_tree_arcs << endl;
			cout << "alive_resources:" << alive_resources << endl;
			cout << "bounded_resources:" << bounded_resources << endl;
			cout << "#tree:" << cover.size() << endl;
			cout << "u5/adding" << u5_time / adding_time << endl;
			*/
			new_tree_arcs = 0;
			int how_many = 0;
			while (beta * labeled_arcs > new_tree_arcs && alive_resources < bounded_resources) {

				if (picked_count == numOfVertices) break;

				NodeID r = rand() % numOfVertices;
				if (picked_as_root[r] == true) continue;
				if (trivial_root[r] == true) continue;

				large_tree new_ltree;
				small_tree new_stree;
				//GOOGLE DENSE MAP
				new_stree.set_empty_key(numOfVertices + 1);
				new_stree.set_deleted_key(numOfVertices + 2);
				bool isLarge = false;
				int tmp_tree_arcs = 0;

				if (!BUILD_SMALL_TREE_FLAG) {

					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_large_tree_dij(r, new_ltree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					u9_time += (GetCurrentTimeSec() - u9_acc_time);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;

					if (descendants.size() < size_thr)
						BUILD_SMALL_TREE_FLAG = true;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					//	double u91_acc_time = GetCurrentTimeSec();

					isLarge = true;
					how_many++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					used_as_cover++;
					st_roots.push_back(r);
					trees_pointers.push_back(-1);
					st_sizes.push_back(descendants.size());

					//new_tree_arcs += descendants.size();
					alive_resources += (long long)numOfVertices;

					//	u91_time += (GetCurrentTimeSec() - u91_acc_time);
						// Turn large tree to small tree.
						// 1. Build new small tree with ltree_i and add it to strees;
						// 2. Set tree_pointers[i] to the position of the strees;
						// 3. Clear the large tree ltree_i by setting its size to 0.

					//	double u92_acc_time = GetCurrentTimeSec();
					int i = st_roots.size() - 1;
					if (st_sizes[i] < size_thr) {
						shrink_tree(i, st_roots[i], new_ltree);
						// FIX
						isLarge = false;
						/*for (int di = descendants.size() - 1; di > -1; --di) {
							NodeID q = descendants[di];
							NodeID v_parent = get_large_parent(q, new_ltree);
							NodeID s_parent = get_small_parent(q, strees[strees.size() - 1]);
							if (v_parent != s_parent) cout << v_parent << " vs " << s_parent << endl;
						}*/
					}
					ltrees.push_back(new_ltree);
					//	u92_time += (GetCurrentTimeSec() - u92_acc_time);
				}
				else {//Build Small Tree Directly
					//double u5_acc_time = GetCurrentTimeSec();
					double u9_acc_time = GetCurrentTimeSec();
					tmp_tree_arcs = build_small_tree_dij(r, new_stree, descendants, descendants_parents, pqueue, visited_que, distances, vis, dst_r, usd, tmp_idx, r_tmp_idx, wgraph, new_tree_arcs, beta * labeled_arcs);
					//u5_time += (GetCurrentTimeSec() - u5_acc_time);
					new_tree_arcs += tmp_tree_arcs;
					u9_time += (GetCurrentTimeSec() - u9_acc_time);

					if (descendants.size() >= size_thr)
						BUILD_SMALL_TREE_FLAG = false;

					if (descendants.size() < 100) {
						//	cout << ranking << ":" << labeled_arcs << " vs. " << new_tree_arcs << " current picked:" << picked_count << " r:" << r << " alive resrouces ratio:" << (double)alive_resources/(double)bounded_resources <<   endl;
						trivial_root[r] = true;
						picked_count++;
						continue; // we have sampled a trivial tree.
					}

					isLarge = false;

					how_many++;
					used_as_cover++;
					//cout << descendants.size() << endl;
					picked_as_root[r] = true;
					picked_count++;
					st_roots.push_back(r);
					trees_pointers.push_back(strees.size()); //we assume it is inserted.
					st_sizes.push_back(descendants.size());

					alive_resources += (long long)new_stree.size();
					//strees.push_back(new_stree);
					//// new_ltree is empty.
					//ltrees.push_back(new_ltree);

					int i = st_roots.size() - 1;
					if (st_sizes[i] >= size_thr) {
						enlarge_tree(i, st_roots[i], new_stree, new_ltree);
						new_stree.clear();
						new_stree.resize(0);
						isLarge = true;
					}
					else {
						strees.push_back(new_stree);
						ltrees.push_back(new_ltree);
					}

				}

				int i = st_roots.size() - 1;

				int numOfStrees = strees.size();
				small_tree &new_stree_ref = numOfStrees != 0 ? strees[strees.size() - 1] : empty_small_tree;

				// Calculate the subtree size of each nodes.

				double u10_acc_time = GetCurrentTimeSec();
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];

					++cover_tid[q]; // Increment because it hitsthe r-q path.
					++cnt16[q][i%CNT];
					//++sum[q];
					//affected[q] = true;

					if (affected_list[q] == false) {
						affected_count++;
						if (switch_set)
							affected_v.insert(q);
						affected_list[q] = true;
					}

					if (switch_small == true) {
						reverse_index[q].insert(i);
						//in_tree_count[q]++;
					}


					int parent_node;
					if (isLarge)
						parent_node = get_large_parent(q, new_ltree);
					else {
						//double get_small_parent_incre_time = GetCurrentTimeSec();
						parent_node = get_small_parent(q, new_stree_ref);
						//get_small_parent_time += (GetCurrentTimeSec() - get_small_parent_incre_time);
					}
					if (parent_node > -1) {
						cover_tid[parent_node] += cover_tid[q];
						cnt16[parent_node][i%CNT] += cover_tid[q];
						//sum[parent_node] += cover_new[q];
					}
				}
				u10_time += (GetCurrentTimeSec() - u10_acc_time);
				for (int di = descendants.size() - 1; di > -1; --di) {
					NodeID q = descendants[di];
					cover_tid[q] = 0;
				}
				//cout << "beta*labeld_arcs: " << beta * labeled_arcs << " new tree _arcs: " << new_tree_arcs << " alive " << alive_resources << " bounded " << bounded_resources << endl;

			}
			adding_time += (GetCurrentTimeSec() - adding_acc_time);


			//double u3_acc_time = GetCurrentTimeSec();

			updating_acc_time = GetCurrentTimeSec();

			// Build the cover value heap for one time.
			if (switch_small == true && value_heap.empty())
				for (int nid = 0; nid < numOfVertices; ++nid)
					if (has_picked[nid] == false)
						value_heap.update(nid, -picked_value[nid]);

			//	cnt16_2max(picked_value, sum, has_picked, affected, cnt16, switch_small, value_heap);
			//c	cnt16_2max_queue(picked_value, sum, has_picked, affected_v, cnt16, switch_small, value_heap);

			double u11_acc_time = GetCurrentTimeSec();
			cnt16_2max_queue(picked_value, sum, has_picked, affected_v, affected_list, cnt16, switch_small, switch_set, value_heap);
			u11_time += (GetCurrentTimeSec() - u11_acc_time);
			//	iteration_time[ranking] = (GetCurrentTimeSec() - updating_acc_time);
			updating_time += (GetCurrentTimeSec() - updating_acc_time);

			if (affected_count < numOfVertices / 8)
				switch_set = true;
			else
				switch_set = false;
			affected_count = 0;

			time_this_iteration = GetCurrentTimeSec() - time_this_iteration;
			iteration_time[ranking] = time_this_iteration;

		}

		cout << "picked as tree covers:" << used_as_cover << endl;
		cout << "picked count:" << picked_count << endl;
		
		//// FIX
		//for (int v = 0; v < numOfVertices; ++v){
		//	for (int i = 0; i < CNT; ++i) {
		//		if (cnt16[v][i] != 0)
		//			cout << "mi?" << endl;
		//	}
		//}
		if (stop_flag == true) return;

		double labeling_acc_time = GetCurrentTimeSec();
		// Order the remaining vertices with degree, if any.
		if (ranking != numOfVertices) {
			vector<NodeID> degree(numOfVertices);
			vector<NodeID> r_degree(numOfVertices);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				/*vector<NodeID>& adj_v = graph.adj[v];
				for (NodeID i = 0; i < adj_v.size(); ++i) {
				if (has_picked[adj_v[i]] != true)
				degree[v]++;
				}*/

				for (EdgeID eid = wgraph.vertices[v]; eid < wgraph.vertices[v + 1]; ++eid) {
					NodeID w = wgraph.edges[eid].first;
					if (has_picked[w] != true)
						degree[v]++;
				}

			}
			
			vector<pair<NodeID, NodeID> > degree_v;
			degree_v.reserve(numOfVertices - ranking);
			for (NodeID v = 0; v < numOfVertices; ++v) {
				if (has_picked[v] == true) continue;
				degree_v.push_back(make_pair(degree[v], v));
			}
			sort(degree_v.rbegin(), degree_v.rend());

			int count = 0;
			for (; ranking < numOfVertices; ++ranking) {
				NodeID tpv = degree_v[count++].second;
				inv[ranking] = tpv;
				rank[tpv] = ranking;
				//output betweenness
				//ofs<<" degree ranking ="<<ranking <<" pv="<<tpv<<" pv_cover="<< degree_v[count++].first<<" picked_count="<<picked_count<<endl;
				labeling_source_dij(tpv, pqueue, visited_que, distances, vis, dst_r, usd, ranking, tmp_idx, r_tmp_idx, wgraph);
			}
		}

		labeling_time += (GetCurrentTimeSec() - labeling_acc_time);

		for (NodeID v = 0; v < numOfVertices; ++v) {

            NodeID k = tmp_idx[v].first.size();
            labels.index_[v].spt_v.resize(k);
            labels.index_[v].spt_d.resize(k);
            for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_v[i] = tmp_idx[v].first[i];
            for (NodeID i = 0; i < k; ++i) labels.index_[v].spt_d[i] = tmp_idx[v].second[i];
            tmp_idx[v].first.clear();
            tmp_idx[v].second.clear();
            tmp_idx[v].first.shrink_to_fit();
            tmp_idx[v].second.shrink_to_fit();

		}

		num_of_trees = st_roots.size();

		double small_tree_num = 0;

		for (int i = 0; i < num_of_trees; ++i) {
			if (trees_pointers[i] != -1) {
				small_tree_num++;
				//small_tree_num += strees[trees_pointers[i]].size();
			}
		}
		
		total_resources = small_tree_num;
		//output time
		start_time=GetCurrentTimeSec()-start_time;
		_hub_time=labeling_time;
		_ordering_time=start_time-_hub_time;
	}
};


class Degree_Ordering : public Ordering {
    public:
        Degree_Ordering(Graph& graph) {
            inv.resize(numOfVertices);
            rank.resize(numOfVertices);

            vector<pair<float, NodeID> > deg(numOfVertices);
            //srand((unsigned)time(NULL));
            srand(100);
            for (size_t v = 0; v < numOfVertices; ++v) {
                if (DIRECTED_FLAG == true)
                    deg[v] = make_pair((graph.vertices[v+1] - graph.vertices[v]) * (graph.r_vertices[v + 1] - graph.r_vertices[v]) + float(rand()) / RAND_MAX, v); // In-degree + Out-degree.
                else
                    deg[v] = make_pair((graph.vertices[v + 1] - graph.vertices[v]) + float(rand()) / RAND_MAX, v);//
            }
            sort(deg.rbegin(), deg.rend());
            for (size_t v = 0; v < numOfVertices; ++v) inv[v] = deg[v].second;
            Relabel(graph);	
        }	

        Degree_Ordering(WGraph& wgraph) {
            inv.resize(numOfVertices);
            rank.resize(numOfVertices);

            vector<pair<float, NodeID> > deg(numOfVertices);
            srand(100);
            for (size_t v = 0; v < numOfVertices; ++v) {
                if (DIRECTED_FLAG == true)
                    deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]) + float(rand()) / RAND_MAX, v);
                else
                    deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) + float(rand()) / RAND_MAX, v);
            }
            sort(deg.rbegin(), deg.rend());
            for (size_t v = 0; v < numOfVertices; ++v) inv[v] = deg[v].second;
            Relabel(wgraph);
        }

        Degree_Ordering(WGraph& wgraph,bool isOrder) {
            inv.resize(numOfVertices);
            rank.resize(numOfVertices);

            vector<pair<float, NodeID> > deg(numOfVertices);
            srand(100);
            for (size_t v = 0; v < numOfVertices; ++v) {
                if (DIRECTED_FLAG == true)
                    deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]) + float(rand()) / RAND_MAX, v);
                else
                    deg[v] = make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) + float(rand()) / RAND_MAX, v);
            }
            sort(deg.rbegin(), deg.rend());
            for (size_t v = 0; v < numOfVertices; ++v) inv[v] = deg[v].second;
            Relabel(wgraph,isOrder);
        }
};

/**
 * @Author: wanjingyi
 * @description: given order to construct labels
 */
class Given_Ordering : public Ordering{
public:
	Given_Ordering(char* order_file, Graph& graph) {
		
		inv.resize(numOfVertices);
		rank.resize(numOfVertices);
		
		ifstream ifs(order_file);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ifs >> inv[v];
		}
		ifs.close();
		Relabel(graph);
	}

	Given_Ordering(char* order_file, WGraph& wgraph) {

		inv.resize(numOfVertices);
		rank.resize(numOfVertices);

		ifstream ifs(order_file);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ifs >> inv[v];
		}
		ifs.close();
		for (NodeID v = 0; v < inv.size(); ++v) rank[inv[v]] = v;
		Relabel(wgraph);
	}

	Given_Ordering(char* order_file, WGraph& wgraph,const vector<int>& originalToNew) {

		inv.resize(numOfVertices);
		rank.resize(numOfVertices);

		ifstream ifs(order_file);
		NodeID tmp;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ifs >> tmp;
			inv[v]=originalToNew[tmp];
		}
		ifs.close();
		for (NodeID v = 0; v < inv.size(); ++v) rank[inv[v]] = v;
		Relabel(wgraph);
	}

	Given_Ordering(char* order_file, WGraph& wgraph,bool isOrder) {

		inv.resize(numOfVertices);
		rank.resize(numOfVertices);

		ifstream ifs(order_file);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			ifs >> inv[v];
		}
		ifs.close();
		for (NodeID v = 0; v < inv.size(); ++v) rank[inv[v]] = v;
		Relabel(wgraph,isOrder);
	}

	Given_Ordering(char* order_file) {

		inv.resize(0);
		ifstream ifs(order_file);
		NodeID inv_id;
		while (ifs >> inv_id) {
			inv.push_back(inv_id);
		}
		ifs.close();
		rank.resize(inv.size());
		for (NodeID v = 0; v < inv.size(); ++v) rank[inv[v]] = v;
	}

	Given_Ordering(WGraph& wgraph,vector<NodeID> inv1,vector<NodeID> rank1) {
		inv.reserve(numOfVertices);
		rank.reserve(numOfVertices);
		inv.assign(inv1.begin(),inv1.end());
		rank.assign(rank1.begin(),rank1.end());
		Relabel(wgraph);
	}

};

#endif