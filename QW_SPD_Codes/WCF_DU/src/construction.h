#ifndef CONSTRUCTION_H
#define CONSTRUCTION_H

#include <omp.h>
#include <unordered_map>
#include<math.h>
#include "./graph.h"
#include "./paras.h"
#include "./labels.h"
#include "./ordering.h"
#include "./heap.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define INF_WEIGHT SP_Constants::INF_WEIGHT

class construction {
    public:
        HFDLabel dlabels;
        HFLabel labels;
        Ordering orders;
};

class PL_W : public construction {
    public:
        vector<double> iteration_generated;
	    vector<double> pruning_power;
        PL_W(WGraph &wgraph, Ordering &orders) {
            iteration_generated.resize(numOfVertices);
            pruning_power.resize(numOfVertices);
            vector<index_t1>& index_ = labels.index_;
            vector<NodeID> &inv = orders.inv;
            vector<NodeID> &rank = orders.rank;
            vector<EdgeID>& vertices = wgraph.vertices;
		    vector<NodeEdgeWeightPair>& edges = wgraph.edges;
            vector<bool> usd(numOfVertices, false);
            vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			    tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),vector<EdgeWeight>(1, INF_WEIGHT)));
            vector<bool> vis(numOfVertices);
            vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
            queue<NodeID> visited_que;
            vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);
            benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);
            long pop = 0;
		    double hsize = 0;

            for (size_t r = 0; r < numOfVertices; ++r) {
                if (usd[r]) continue;
                const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];
                for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				    dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			    }
                pqueue.update(r, 0);
                long max_heap_size = 0;
                long heap_size = 1;

                while (!pqueue.empty()) {
                    pop++;
                    heap_size--;
                    NodeID v;
                    EdgeWeight v_d;
                    pqueue.extract_min(v, v_d);
                    pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
                    vis[v] = true;
                    visited_que.push(v);
                    if (usd[v]) continue;
                        for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
                        NodeID w = tmp_idx_v.first[i];
                        EdgeWeight td = tmp_idx_v.second[i] + dst_r[w];
                        if (td <= v_d) {
                            pruning_power[w]++;
                            goto pruned;
                        }
                    }
                    // Traverse
                    tmp_idx_v.first.back() = r;
                    tmp_idx_v.second.back() = v_d;
                    tmp_idx_v.first.push_back(numOfVertices);
                    tmp_idx_v.second.push_back(INF_WEIGHT);
                    iteration_generated[r]++;
                    for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
                        NodeID w = edges[eid].first;
                        EdgeWeight w_d = edges[eid].second + v_d;
                        if (!vis[w]&&r<w) {
                            if (distances[w] == INF_WEIGHT) {
                                heap_size++;
                                if (max_heap_size < heap_size)
                                    max_heap_size = heap_size;
                            }

                            if( distances[w] > w_d ){
                                pqueue.update(w, w_d);
                                distances[w] = w_d;
                            }
                        }
                    }
                    pruned:
                        {}
                }

                hsize = hsize + max_heap_size;//count the hsize
                while (!visited_que.empty()) {
                    NodeID vis_v = visited_que.front();
                    visited_que.pop();
                    vis[vis_v] = false;
                    distances[vis_v] = INF_WEIGHT;
                    pqueue.clear(vis_v);
                }
                pqueue.clear_n();
                for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				    dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;
                usd[r] = true;
            }

            std::cout << "total pop = " << pop << endl;
            std::cout << "toal heap size = " << hsize<<endl;
		    std::cout << "average heap size = " << (double)hsize / (double)numOfVertices << endl;

            double count = 0;
            for (size_t v = 0; v < numOfVertices; ++v) {
                NodeID k = tmp_idx[v].first.size();
                count = count + k - 1;
                index_[inv[v]].spt_v.resize(k);
                index_[inv[v]].spt_d.resize(k);
                for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
                for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
                tmp_idx[v].first.clear(); 
                tmp_idx[v].second.clear();
                tmp_idx[v].first.shrink_to_fit();
                tmp_idx[v].second.shrink_to_fit();
            }
            cout << "Average Label Size:" << count / numOfVertices << endl;

        }

    //Add fordirected
	/**
     * @description: consturct pll labels
     * @param {WGraph} &wgraph
     * @param {Ordering} &orders
     * @param {bool} D_FLAGS
     * @return {*}
     * @author: Wan Jingyi
     */ 
    PL_W(WGraph &wgraph, Ordering &orders, bool D_FLAGS) {

		iteration_generated.resize(numOfVertices);
		pruning_power.resize(numOfVertices);

		vector<index_t1>& index_ = dlabels.index_;

		vector<NodeID> &inv = orders.inv;
		vector<NodeID> &rank = orders.rank;
		/*vector<vector<NodeID> > &adj = wgraph.adj;
		vector<vector<EdgeWeight> > &adj_weight = wgraph.adj_weight;*/

		vector<index_t1>& bindex_ = dlabels.bindex_;/*
		vector<vector<NodeID> > &r_adj = wgraph.r_adj;
		vector<vector<EdgeWeight> > &r_adj_weight = wgraph.r_adj_weight;*/

		//Array Representation
		vector<EdgeID>& vertices = wgraph.vertices;
		vector<EdgeID>& r_vertices = wgraph.r_vertices;
		vector<NodeEdgeWeightPair>& edges = wgraph.edges;
		vector<NodeEdgeWeightPair>& r_edges = wgraph.r_edges;

		vector<bool> usd(numOfVertices, false);

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<pair<vector<NodeID>, vector<EdgeWeight> > >
			r_tmp_idx(numOfVertices, make_pair(vector<NodeID>(1, numOfVertices),
				vector<EdgeWeight>(1, INF_WEIGHT)));

		vector<bool> vis(numOfVertices);
		vector<EdgeWeight> dst_r(numOfVertices + 1, INF_WEIGHT);
		vector<EdgeWeight> r_dst_r(numOfVertices + 1, INF_WEIGHT);
		queue<NodeID> visited_que;
		vector<EdgeWeight> distances(numOfVertices, INF_WEIGHT);


		benchmark::heap<2, EdgeWeight, NodeID> pqueue(numOfVertices);

		// Forward search from r.
		for (size_t r = 0; r < numOfVertices; ++r) {
			if (usd[r]) continue;

			const pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_r = tmp_idx[r];

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			pqueue.update(r, 0);
			//vis[r] = true;

			while (!pqueue.empty()) {
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_v = r_tmp_idx[v];
				vis[v] = true;
				visited_que.push(v);

				if (usd[v]) continue;
				for (size_t i = 0; i < r_tmp_idx_v.first.size(); ++i) {
					NodeID w = r_tmp_idx_v.first[i];
					EdgeWeight td = r_tmp_idx_v.second[i] + dst_r[w];
					if (td <= v_d) {
						pruning_power[w]++;
						goto pruned_forward;
					}
				}

				// Traverse
				r_tmp_idx_v.first.back() = r;
				r_tmp_idx_v.second.back() = v_d;
				r_tmp_idx_v.first.push_back(numOfVertices);
				r_tmp_idx_v.second.push_back(INF_WEIGHT);
				iteration_generated[r]++;
				
			
			/*	for (size_t i = 0; i < adj[v].size(); ++i) {
					NodeID w = adj[v][i];
					EdgeWeight w_d = adj_weight[v][i] + v_d;*/

					//Array Representation

			for (EdgeID eid = vertices[v]; eid < vertices[v + 1]; ++eid) {
				NodeID w = edges[eid].first;
				EdgeWeight w_d = edges[eid].second + v_d;
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

			while (!visited_que.empty()) {
				NodeID vis_v = visited_que.front();
				visited_que.pop();
				vis[vis_v] = false;
				distances[vis_v] = INF_WEIGHT;
				//pqueue.clear(vis_v);
			}

			//pqueue.clear_n();

			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i)
				dst_r[tmp_idx_r.first[i]] = INF_WEIGHT;

			
			// Backward search from r.
			const pair<vector<NodeID>, vector<EdgeWeight> > &r_tmp_idx_r = r_tmp_idx[r];

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i) {
				r_dst_r[r_tmp_idx_r.first[i]] = r_tmp_idx_r.second[i];
			}

			pqueue.update(r, 0);

			while (!pqueue.empty()) {
				NodeID v;
				EdgeWeight v_d;
				pqueue.extract_min(v, v_d);
				pair<vector<NodeID>, vector<EdgeWeight> > &tmp_idx_v = tmp_idx[v];
				vis[v] = true;
				visited_que.push(v);

				if (usd[v]) continue;
				for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
					NodeID w = tmp_idx_v.first[i];
					EdgeWeight td = tmp_idx_v.second[i] + r_dst_r[w];
					if (td <= v_d) {
						pruning_power[w]++;
						goto pruned_backward;
					}
				}

				// Traverse
				tmp_idx_v.first.back() = r;
				tmp_idx_v.second.back() = v_d;
				tmp_idx_v.first.push_back(numOfVertices);
				tmp_idx_v.second.push_back(INF_WEIGHT);

				iteration_generated[r]++;

				/*for (size_t i = 0; i < r_adj[v].size(); ++i) {
					NodeID w = r_adj[v][i];
					EdgeWeight w_d = r_adj_weight[v][i] + v_d;*/

				//Array Representation
				for (EdgeID eid = r_vertices[v]; eid < r_vertices[v + 1]; ++eid) {
					NodeID w = r_edges[eid].first;
					EdgeWeight w_d = r_edges[eid].second + v_d;

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
				//pqueue.clear(vis_v);
			}

		//	pqueue.clear_n();

			for (size_t i = 0; i < r_tmp_idx_r.first.size(); ++i)
				r_dst_r[r_tmp_idx_r.first[i]] = INF_WEIGHT;

			usd[r] = true;			
		}

		for (size_t v = 0; v < numOfVertices; ++v) {
			NodeID k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v.resize(k);
			index_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
			tmp_idx[v].first.shrink_to_fit();
			tmp_idx[v].second.shrink_to_fit();

			k = r_tmp_idx[v].first.size();
			bindex_[inv[v]].spt_v.resize(k);
			bindex_[inv[v]].spt_d.resize(k);
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_v[i] = r_tmp_idx[v].first[i];
			for (NodeID i = 0; i < k; ++i) bindex_[inv[v]].spt_d[i] = r_tmp_idx[v].second[i];
			r_tmp_idx[v].first.clear();
			r_tmp_idx[v].second.clear();		
			r_tmp_idx[v].first.shrink_to_fit();
			r_tmp_idx[v].second.shrink_to_fit();	
		}

	}

};

#endif