/*
 * @Descripttion: vertex-centric padradigm to construct labels parallelly
 * @version: 
 * @Author: wanjingyi
 * @Date: 2021-02-22 15:19:41
 * @LastEditors: Please set LastEditors
 * @LastEditTime: 2021-03-01 20:27:47
 */
#ifndef _PARALLEL_CONSTRUCTION_H
#define _PARALLEL_CONSTRUCTION_H

#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <climits>
#include <xmmintrin.h>
#include <immintrin.h>
#include <bitset>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "./graph.h"
#include "./paras.h"
#include "./time_util.h"
#include "./utils.h"

using std::vector;
using std::unordered_map;
using std::map;
using std::bitset;
using std::stable_sort;
using std::min;
using std::fill;
using namespace time_util;

namespace PADO {
// Structure for the type of label
struct IndexType {
    vector<NodeID> vertices; // Vertices in the label, preresented as temperory ID
    vector<EdgeWeight> distances;
    NodeID size() {
    return vertices.size();
    }
}; //__attribute__((aligned(64)));

// Structure for the type of temporary label
struct ShortIndex {
    // Use a queue to store candidates
    vector<NodeID> candidates_que;
    NodeID end_candidates_que = 0;
    // Use a array to store distances of candidates; length of roots_size
    // record the distances to candidates. If candidates_dists[c] = INF, then c is NOT a candidate
    // The candidates_dists is also used for distance query.
    vector<EdgeWeight> candidates_dists; 
    // Use a queue to store temporary labels in this batch; use it so don't need to traverse vertices_dists.
    // Elements in vertices_que are roots_id, not real id
    vector<NodeID> vertices_que; 
    NodeID end_vertices_que = 0;
    // Use an array to store distances to roots; length of roots_size
    vector<EdgeWeight> vertices_dists; // labels_table
    // Usa a queue to record which roots have reached this vertex in this batch.
    // It is used for reset the dists_table
    vector<NodeID> reached_roots_que;
    NodeID end_reached_roots_que = 0;
    // Use a queue to store last inserted labels (IDs); distances are stored in vertices_dists.
    vector<NodeID> last_new_roots;
    NodeID end_last_new_roots = 0;
    ShortIndex(){}
    ShortIndex (int num_roots){
        candidates_que.resize(num_roots);
        candidates_dists.resize(num_roots, INF_WEIGHT);
        last_new_roots.resize(num_roots);
        vertices_que.resize(num_roots);
        vertices_dists.resize(num_roots, INF_WEIGHT);
        reached_roots_que.resize(num_roots);
    }
};//__attribute__((aligned(64)));

// Structure of the public ordered index for distance queries.
struct index_t {
	vector<NodeID> spt_v;
	vector<EdgeWeight> spt_d;
	NodeID size() {
		return spt_v.size();
	}
};

class WeightedVCPLL{
public:
    //***************variables*************
    vector<IndexType> L;
	vector<index_t> index_; // Ordered labels for original vertex ID
    NodeID BATCH_SIZE;
    //***************functions*************
    void reorder_labels(const vector<NodeID> &rank){
        int labels_count = 0;
        index_.resize(numOfVertices);
        for (NodeID v_id = 0; v_id < numOfVertices; ++v_id) {
            NodeID v_rank = rank[v_id];
            const IndexType &Lv = L[v_rank];
            NodeID size_labels = Lv.vertices.size();
            labels_count += size_labels;
            vector< std::pair<NodeID, EdgeWeight> > ordered_labels(size_labels);
            // Traverse v_id's all labels
            for (NodeID l_i = 0; l_i < size_labels; ++l_i) {
                ordered_labels[l_i].first = Lv.vertices[l_i];
                ordered_labels[l_i].second = Lv.distances[l_i];
            }
            // Sort
 		    sort(ordered_labels.begin(), ordered_labels.end());
            index_[v_id].spt_v.resize(size_labels+1);
            index_[v_id].spt_d.resize(size_labels+1);
            for (NodeID l_i = 0; l_i < size_labels; ++l_i) {
                index_[v_id].spt_v[l_i]=ordered_labels[l_i].first;
                index_[v_id].spt_d[l_i] = ordered_labels[l_i].second;
            }
            index_[v_id].spt_v[size_labels]=numOfVertices;
            index_[v_id].spt_d[size_labels]=INF_WEIGHT;
        }
        cout<<"Total label size="<<labels_count<<" Average label size="<<static_cast<double>(labels_count)/numOfVertices<<endl;
    }

    void save_labels(const char *save_filename)
    {
        string label_filename(save_filename);
        label_filename.append(".label");
        ofstream ofs(label_filename, ios::binary | ios::out);
        if(!ofs.is_open()){
            cout<<"Cannot open "<<label_filename<<endl;
            exit(EXIT_FAILURE);
        }
        // Store into file the number of vertices and the number of bit-parallel roots.
        ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
        for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID isize = index_[v].size();
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
				ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
			}
		}
        ofs.close();
    }

	void save_label_size(char* label_size_file,const vector<NodeID>& inv) {
		string labelSizefile(label_size_file);
		labelSizefile.append(".size");
		ofstream ofs(labelSizefile);
		if(!ofs.is_open()) {cerr<<"Cannot open "<<labelSizefile<<endl;}
		for (int i = 0; i < numOfVertices; ++i) {
			ofs << index_[i].size()-1 << endl;
		}
		ofs.close();
		//output the size by descending order,format:label_size
		string labelSizefile_prefix1(label_size_file);
		string labelSizefile1=labelSizefile_prefix1.append("_byOrder.size");
		ofstream out(labelSizefile1.c_str());
		if(!out.is_open()) {cerr<<"Cannot open "<<labelSizefile1<<endl;}
		for (int i = 0; i < numOfVertices; ++i) {
			out<<index_[inv[i]].size()-1 << endl;
		}		
		out.close();
	}

	void write_labels(char* write_filename,const vector<NodeID>& inv,bool isOrder=false)
	{
			string write_filename_prefix(write_filename);
			string write_filename1=write_filename_prefix.append(".list");
			ofstream ofs(write_filename1.c_str());
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			NodeID isize = index_[v].size();
			ofs <<isize<<" "<<v;
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				ofs<<" "<<'('<<index_[v].spt_v[i]<<","<<index_[v].spt_d[i]<<")";
			}
			ofs<<endl;
		}
		ofs.close();

		//write labels with original NodeId in graph
		string write_filename2=write_filename_prefix.append("_original");
		ofstream out(write_filename2.c_str());
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			NodeID isize = index_[v].size();
			out<<isize-1<<" "<<v;
			for (NodeID i = 0; i < index_[v].size()-1; ++i) {
				out<<" "<<'('<<inv[index_[v].spt_v[i]]<<","<<index_[v].spt_d[i]<<")";
			}
			out<<endl;
		}			
		out.close();

		if(isOrder){
			string write_filename_prefix1(write_filename);
			string write_filename3=write_filename_prefix1.append(".list_order");
			ofstream out1(write_filename3.c_str());
			for (NodeID r= 0; r < numOfVertices; ++r) 
			{
				NodeID v=inv[r];
				NodeID isize = index_[v].size();
				out1<<isize-1<<" "<<v;
				for (NodeID i = 0; i < index_[v].size()-1; ++i) {
					out1<<" "<<'('<<inv[index_[v].spt_v[i]]<<","<<index_[v].spt_d[i]<<")";
				}
				out1<<endl;
			}			
			out1.close();
			}
	}
};

class WeightedVertexCentricPLL:public WeightedVCPLL{
public:
    //****************variables***********
    //*********************constructions and deconstructions***********
    WeightedVertexCentricPLL(){
    }

    WeightedVertexCentricPLL(WGraph& wgraph,int batchSize){
        BATCH_SIZE=batchSize;
        construct(wgraph);
    }

    ~WeightedVertexCentricPLL(){ }

protected:
    void construct(WGraph& wgraph){
        L.resize(numOfVertices);
        NodeID remainer = numOfVertices % BATCH_SIZE;
	    NodeID b_i_bound = numOfVertices / BATCH_SIZE;
        vector<NodeID> active_queue(numOfVertices);//Active queue
        NodeID end_active_queue = 0;
        // is_active[v] is true means vertex v is in the active queue.
        vector<bool> is_active(numOfVertices, false);	
        // Queue of vertices having candidates
        vector<NodeID> has_cand_queue(numOfVertices);
        NodeID end_has_cand_queue = 0;
        // has_candidates[v] is true means vertex v is in the queue has_cand_queue.
        vector<bool> has_candidates(numOfVertices, false);
        // Distance table of shortest distance from roots to other vertices.
        // The distance table is roots_sizes by N. 
        // 1. record the shortest distance so far from every root to every vertex; (from v to root r)
        vector< vector<EdgeWeight> > dists_table(numOfVertices, vector<EdgeWeight>(BATCH_SIZE, INF_WEIGHT));
        // r's old label distance from r to v
        //The distance buffer, recording old label distances of every root. It needs to be initialized every batch by labels of roots. (dists_table[r][l], r > l)
        vector< vector<EdgeWeight> > roots_labels_buffer(BATCH_SIZE, vector<EdgeWeight>(numOfVertices, INF_WEIGHT)); 
        // Here the size of short_index actually is fixed because it is static.
        vector<ShortIndex> short_index(numOfVertices, ShortIndex(BATCH_SIZE));
        // A queue to store vertices which have got new labels in this batch. This queue is used for update the index.
        vector<NodeID> has_new_labels_queue(numOfVertices);
        NodeID end_has_new_labels_queue = 0;
        vector<bool> has_new_labels(numOfVertices, false);
        // A queue to store vertices which have got candidates in this batch. This queue is used fore reset dists_table, accompanying with reached_roots_que in the ShortIndex.
        vector<NodeID> once_had_cand_queue(numOfVertices);
        NodeID end_once_had_cand_queue = 0;
        vector<bool> once_had_cand(numOfVertices, false);
        double hub_time=GetCurrentTimeSec();
        for (NodeID b_i = 0; b_i < b_i_bound; ++b_i) {
            vertex_centric_labeling_in_batches(
                    wgraph,
                    b_i,
                    b_i * BATCH_SIZE,
                    BATCH_SIZE,
                    L,
                    active_queue,
                    end_active_queue,
                    is_active,
                    has_cand_queue,
                    end_has_cand_queue,
                    has_candidates,
                    dists_table,
                    roots_labels_buffer,
                    short_index,
                    has_new_labels_queue,
                    end_has_new_labels_queue,
                    has_new_labels,
                    once_had_cand_queue,
                    end_once_had_cand_queue,
                    once_had_cand);
        }
        if (0 != remainer) {
            vertex_centric_labeling_in_batches(
                    wgraph,
                    b_i_bound,
                    b_i_bound * BATCH_SIZE,
                    remainer,
                    L,
                    active_queue,
                    end_active_queue,
                    is_active,
                    has_cand_queue,
                    end_has_cand_queue,
                    has_candidates,
                    dists_table,
                    roots_labels_buffer,
                    short_index,
                    has_new_labels_queue,
                    end_has_new_labels_queue,
                    has_new_labels,
                    once_had_cand_queue,
                    end_once_had_cand_queue,
                    once_had_cand);
        }
        hub_time=GetCurrentTimeSec()-hub_time;
        cout << "labeling time:" << hub_time *1e6 <<  " microseconds" << endl;
        double ave_hub_time=hub_time/(double) numOfVertices;
        cout<<"average labeling time:"<<ave_hub_time*1e6 <<  " microseconds" << endl;
    }

    void initialize_tables(
			vector<ShortIndex> &short_index,
			vector< vector<EdgeWeight> > &roots_labels_buffer,
			vector<NodeID> &active_queue,
			NodeID &end_active_queue,
			NodeID roots_start,
			NodeID roots_size,
			vector<IndexType> &L,
			vector<NodeID> &has_new_labels_queue,
			NodeID &end_has_new_labels_queue,
			vector<bool> &has_new_labels)
    {
            NodeID roots_bound = roots_start + roots_size;
            // **********Distance Tables**********
            for (NodeID r_root_id = 0; r_root_id < roots_size; ++r_root_id) {// Traverse all roots
                NodeID r_real_id = r_root_id + roots_start;
                const IndexType &Lr = L[r_real_id];
                NodeID l_i_bound = Lr.vertices.size();
                _mm_prefetch(&Lr.vertices[0], _MM_HINT_T0);
                _mm_prefetch(&Lr.distances[0], _MM_HINT_T0);
                for (NodeID l_i = 0; l_i < l_i_bound; ++l_i) {
				roots_labels_buffer[r_root_id][Lr.vertices[l_i]] = Lr.distances[l_i];
			    }
                roots_labels_buffer[r_root_id][r_real_id] = 0;
            }
            // **********Distance Tables**********

            // **********Short Index**********
            for (NodeID r_real_id = roots_start; r_real_id < roots_bound; ++r_real_id) {
                ShortIndex &SI_r = short_index[r_real_id];
                NodeID r_root_id = r_real_id - roots_start;
                // Record itself as last inserted label distance
                SI_r.vertices_que[SI_r.end_vertices_que++] = r_root_id;
                SI_r.vertices_dists[r_root_id] = 0;
                SI_r.last_new_roots[SI_r.end_last_new_roots++] = r_root_id;
            }

            // **********Active queue**********
            for (NodeID r_real_id = roots_start; r_real_id < roots_bound; ++r_real_id) {
                active_queue[end_active_queue++] = r_real_id;
            }
            // **********Active queue**********

            // **********Has new labels queue**********
            // has_new_labels_queue: Put all roots into the has_new_labels_queue
            for (NodeID r_real_id = roots_start; r_real_id < roots_bound; ++r_real_id) {
                has_new_labels_queue[end_has_new_labels_queue++] = r_real_id;
                has_new_labels[r_real_id] = true;
            }
            // **********Has new labels queue**********
    }

    void send_messages(
        NodeID v_head,
        NodeID roots_start,
        WGraph &wgraph,
        vector< vector<EdgeWeight> > &dists_table,
        vector<ShortIndex> &short_index,
        vector<NodeID> &has_cand_queue,
        NodeID &end_has_cand_queue,
        vector<bool> &has_candidates,
        vector<NodeID> &once_had_cand_queue,
        NodeID &end_once_had_cand_queue,
        vector<bool> &once_had_cand
    )
    {
        //graph
        vector<EdgeID>& wgraph_vertices = wgraph.vertices;
        vector<NodeEdgeWeightPair>& wgraph_edges = wgraph.edges;
        ShortIndex &SI_v_head = short_index[v_head];
        // Traverse v_head's every neighbor v_tail
        for (EdgeID eid = wgraph_vertices[v_head]; eid < wgraph_vertices[v_head + 1]; ++eid) {
            NodeID v_tail = wgraph_edges[eid].first;
            if(v_tail<=roots_start){ // v_tail has higher rank than any roots, then no roots can push new labels to it.
                break;
            }
            EdgeWeight weight_h_t = wgraph_edges[eid].second;
            ShortIndex &SI_v_tail = short_index[v_tail];
            bool got_candidates = false; // A flag indicates if v_tail got new candidates.
            // Traverse v_head's last inserted labels
            NodeID bound_r_i = SI_v_head.end_last_new_roots;
            for (NodeID r_i = 0; r_i < bound_r_i; ++r_i) {
                NodeID r_root_id = SI_v_head.last_new_roots[r_i]; // last inserted label r_root_id of v_head
                if (v_tail <= r_root_id + roots_start) {
				    continue;
                }
                EdgeWeight tmp_dist_r_t = SI_v_head.vertices_dists[r_root_id] + weight_h_t;
                if (tmp_dist_r_t < dists_table[v_tail][r_root_id] && tmp_dist_r_t < SI_v_tail.candidates_dists[r_root_id]) {
                    // Mark r_root_id as a candidate of v_tail
                    if (INF_WEIGHT == SI_v_tail.candidates_dists[r_root_id]) {
                        // Add r_root_id into v_tail's candidates_que
                        SI_v_tail.candidates_que[SI_v_tail.end_candidates_que++] = r_root_id;
                    }
                    SI_v_tail.candidates_dists[r_root_id] = tmp_dist_r_t;
                    if (INF_WEIGHT == dists_table[v_tail][r_root_id]) {
                        // Add r_root_id into v_tail's reached_roots_que so to reset dists_table[r_root_id][v_tail] at the end of this batch
                        SI_v_tail.reached_roots_que[SI_v_tail.end_reached_roots_que++] = r_root_id;
                    }
                    dists_table[v_tail][r_root_id] = tmp_dist_r_t; // For filtering out longer distances coming later
                    got_candidates = true;
                }
            }

            if (got_candidates && !has_candidates[v_tail]) {
                // Put v_tail into has_cand_queue
                has_candidates[v_tail] = true;
                has_cand_queue[end_has_cand_queue++] = v_tail;
                if (!once_had_cand[v_tail]) {
                    once_had_cand[v_tail] = true;
                    once_had_cand_queue[end_once_had_cand_queue++] = v_tail;
                }
            }
        }
        SI_v_head.end_last_new_roots = 0; // clear v_head's last_new_roots
    }


    //Function: to see if v_id and cand_root_id can have other path already cover the candidate distance
    // If there was other path, return the shortest distance. If there was not, return INF
    EdgeWeight distance_query(
        NodeID v_id,
        NodeID cand_root_id,
        const vector< vector<EdgeWeight> > &roots_labels_buffer,
        const vector<ShortIndex> &short_index,
        const vector<IndexType> &L,
        NodeID roots_start,
        EdgeWeight tmp_dist_v_c
    )
    {
        // Traverse all available hops of v, to see if they reach c
	    // 1. Labels in L[v]
        NodeID cand_real_id = cand_root_id + roots_start;
        const IndexType &Lv = L[v_id];
        NodeID bound_i_l = Lv.vertices.size();
        for (NodeID i_l = 0; i_l < bound_i_l; ++i_l) {
            NodeID r = Lv.vertices[i_l];
            if (cand_real_id <= r || INF_WEIGHT == roots_labels_buffer[cand_root_id][r]) {
                continue;
            }
            EdgeWeight label_dist_v_c = Lv.distances[i_l] + roots_labels_buffer[cand_root_id][r];
            if (label_dist_v_c <= tmp_dist_v_c) {
                return label_dist_v_c;
            }
        }
        // 2. Labels in short_index[v_id].vertices_que
        const ShortIndex &SI_v = short_index[v_id];
        const ShortIndex &SI_c = short_index[cand_root_id];
        NodeID bound_i_que = SI_v.end_vertices_que;
        for (NodeID i_que = 0; i_que < bound_i_que; ++i_que) {
            NodeID r_root_id = SI_v.vertices_que[i_que];
            if (cand_real_id <= r_root_id + roots_start) {
                continue;
            }
            // Check r_root_id in cand_root_id's labels inserted in this batch
            if (INF_WEIGHT != SI_c.vertices_dists[r_root_id]) {
                EdgeWeight label_dist_v_c = SI_v.vertices_dists[r_root_id] + SI_c.vertices_dists[r_root_id];
                if (label_dist_v_c <= tmp_dist_v_c) {
                    return label_dist_v_c;
                }
            }
        }
        return INF_WEIGHT;
    }

    // Function: at the end of every batch, filter out wrong and redundant labels
    inline void filter_out_labels(
        const vector<NodeID> &has_new_labels_queue,
		NodeID end_has_new_labels_queue,
		vector<ShortIndex> &short_index,
		NodeID roots_start
    )
    {
        // Traverse has_new_labels_queue for every vertex who got new labels in this batch
        for (NodeID i_q = 0; i_q < end_has_new_labels_queue; ++i_q) {
            NodeID v_id = has_new_labels_queue[i_q];
            IndexType &Lv = L[v_id];
		    ShortIndex &SI_v = short_index[v_id];
            NodeID bound_label_i = SI_v.end_vertices_que;
            // Traverse v_id's every new label
            for (NodeID i_c = 0; i_c < bound_label_i; ++i_c) {
                NodeID c_root_id = SI_v.vertices_que[i_c];
                EdgeWeight dist_v_c = SI_v.vertices_dists[c_root_id];
                ShortIndex &SI_c = short_index[c_root_id + roots_start];
                // Traverse v_id's other new label to see if has shorter distance
                for (NodeID i_r = 0; i_r < bound_label_i; ++i_r) {
                    // r_root_id is the middle hop.
                    NodeID r_root_id = SI_v.vertices_que[i_r];
                    if (c_root_id <= r_root_id) {
                        continue;
                    }
                    EdgeWeight dist_v_r = SI_v.vertices_dists[r_root_id];
                    if (dist_v_r > dist_v_c) {
                        continue;
                    }
                    EdgeWeight dist_c_r = SI_c.vertices_dists[r_root_id];
                    if (dist_c_r > dist_v_c) {
                        continue;
                    }
                    if (dist_v_r + dist_c_r <= dist_v_c) {
                        // Shorter distance exists, then label c is not necessary
                        SI_v.vertices_dists[c_root_id] = INF_WEIGHT; // Filter out c_root_id from v_id's labels
                    }
                }
            }
        }
    }

    void reset_tables(
        vector<ShortIndex> &short_index,
		NodeID roots_start,
		NodeID roots_size,
		const vector<IndexType> &L,
		vector< vector<EdgeWeight> > &dists_table,
		vector< vector<EdgeWeight> > &roots_labels_buffer,
		vector<NodeID> &once_had_cand_queue,
		NodeID &end_once_had_cand_queue,
		vector<bool> &once_had_cand
    )
    {
        // Reset roots_labels_buffer according to L (old labels)
        for (NodeID r_roots_id = 0; r_roots_id < roots_size; ++r_roots_id) {
            NodeID r_real_id = r_roots_id + roots_start;
            // Traverse labels of r
            const IndexType &Lr = L[r_real_id];
            NodeID bound_i_l = Lr.vertices.size();
            for (NodeID i_l = 0; i_l < bound_i_l; ++i_l) {
                roots_labels_buffer[r_roots_id][Lr.vertices[i_l]] = INF_WEIGHT;
            }
            roots_labels_buffer[r_roots_id][r_real_id] = INF_WEIGHT;
        }
        
        // Reset dists_table according to short_index[v].reached_roots_que and once_had_cand_queue
        for (NodeID i_q = 0; i_q < end_once_had_cand_queue; ++i_q) {
            NodeID v_id = once_had_cand_queue[i_q];
            once_had_cand[v_id] = false; // Reset once_had_cand
            ShortIndex &SI_v = short_index[v_id];
            // Traverse roots which have reached v_id
            NodeID bound_i_r = SI_v.end_reached_roots_que;
            for (NodeID i_r = 0; i_r < bound_i_r; ++i_r) {
                dists_table[v_id][SI_v.reached_roots_que[i_r]] = INF_WEIGHT;
            }
            SI_v.end_reached_roots_que = 0; // Clear v_id's reached_roots_que
        }
        end_once_had_cand_queue = 0; // Clear once_had_cand_queue.

    }

    void update_index(
        vector<IndexType> &L,
		vector<ShortIndex> &short_index,
		NodeID roots_start,
		vector<NodeID> &has_new_labels_queue,
		NodeID &end_has_new_labels_queue,
		vector<bool> &has_new_labels
    )
    {
        for (NodeID i_q = 0; i_q < end_has_new_labels_queue; ++i_q) {
            NodeID v_id = has_new_labels_queue[i_q];
            has_new_labels[v_id] = false; // Reset has_new_labels
            IndexType &Lv = L[v_id];
		    ShortIndex &SI_v = short_index[v_id];
            NodeID bound_i_r = SI_v.end_vertices_que;
            for (NodeID i_r = 0; i_r < bound_i_r; ++i_r) {
                NodeID r_root_id = SI_v.vertices_que[i_r];
                EdgeWeight dist = SI_v.vertices_dists[r_root_id];
                if (INF_WEIGHT == dist) {
                    continue;
                }
                SI_v.vertices_dists[r_root_id] = INF_WEIGHT; // Reset v_id's vertices_dists
                Lv.vertices.push_back(r_root_id + roots_start);
                Lv.distances.push_back(dist);
            }
            SI_v.end_vertices_que = 0; // Clear v_id's vertices_que
        }
        end_has_new_labels_queue = 0; // Clear has_new_labels_queue
    }

    void vertex_centric_labeling_in_batches(
        WGraph& wgraph,
        NodeID b_id,
        NodeID roots_start,// start id of roots
        NodeID roots_size, // how many roots in the batch
        vector<IndexType> &L,
        vector<NodeID> &active_queue,
        NodeID end_active_queue,
		vector<bool> &is_active,
        vector<NodeID> &has_cand_queue,
        NodeID end_has_cand_queue,
        vector<bool> &has_candidates,
        vector< vector<EdgeWeight> > &dists_table,
        vector< vector<EdgeWeight> > &roots_labels_buffer,
        vector<ShortIndex> &short_index,
        vector<NodeID> &has_new_labels_queue,
        NodeID end_has_new_labels_queue,
        vector<bool> &has_new_labels,
        vector<NodeID> &once_had_cand_queue,
        NodeID end_once_had_cand_queue,
        vector<bool> &once_had_cand
    )
    {
        // At the beginning of a batch, initialize the labels L and distance buffer dist_matrix;
        initialize_tables(
                short_index,
                roots_labels_buffer,
                active_queue,
                end_active_queue,
                roots_start,
                roots_size,
                L,
                has_new_labels_queue,
                end_has_new_labels_queue,
                has_new_labels);
        
        while (0 != end_active_queue) {
            // ************First stage, sending distances.************
            // Traverse the active queue, every active vertex sends distances to its neighbors
            //cout<<"************************************"<<endl;//to be deleted
            for (NodeID i_queue = 0; i_queue < end_active_queue; ++i_queue) {
                NodeID v_head = active_queue[i_queue];
                //cout<<v_head<<":"<<endl;//to be deleted
                is_active[v_head] = false; // reset is_active
                send_messages(
					v_head,
					roots_start,
					wgraph,
					dists_table,
					short_index,
					has_cand_queue,
					end_has_cand_queue,
					has_candidates,
					once_had_cand_queue,
					end_once_had_cand_queue,
					once_had_cand);
                // //to be deleted
                // cout<<"has_cand_queue:end_has_cand_queue="<<end_has_cand_queue;
                // for(int i=0;i<end_has_cand_queue;++i) cout<<" "<<has_cand_queue[i];
                // cout<<endl;
                // //to be deleted
                // cout<<"once_had_cand_queue:end_has_cand_queue="<<end_once_had_cand_queue;
                // for(int i=0;i<end_once_had_cand_queue;++i) cout<<" "<<once_had_cand_queue[i];
                // cout<<endl;
            }

            end_active_queue = 0; // Set the active_queue empty
            // Traverse vertices in the has_cand_queue to insert labels
            for (NodeID i_queue = 0; i_queue < end_has_cand_queue; ++i_queue) {
                NodeID v_id = has_cand_queue[i_queue];
                bool need_activate = false;
                ShortIndex &SI_v = short_index[v_id];
                // Traverse v_id's all candidates
                NodeID bound_cand_i = SI_v.end_candidates_que;
                for (NodeID cand_i = 0; cand_i < bound_cand_i; ++cand_i) {
                    NodeID cand_root_id = SI_v.candidates_que[cand_i];
                    EdgeWeight tmp_dist_v_c = SI_v.candidates_dists[cand_root_id];
                    // Distance check for pruning
                    EdgeWeight query_dist_v_c;
                    if(INF_WEIGHT==
                        (query_dist_v_c = distance_query(
                                        v_id,
                                        cand_root_id,
                                        roots_labels_buffer,
                                        short_index,
                                        L,
                                        roots_start,
                                        tmp_dist_v_c)))
                    {
                        if (INF_WEIGHT == SI_v.vertices_dists[cand_root_id]) {
                            // Record cand_root_id as v_id's label
                            SI_v.vertices_que[SI_v.end_vertices_que++] = cand_root_id;
                        }
                        // Record the new distance in the label table
                        SI_v.vertices_dists[cand_root_id] = tmp_dist_v_c;
                        SI_v.last_new_roots[SI_v.end_last_new_roots++] = cand_root_id;
                        need_activate = true;
                    }else if(query_dist_v_c < tmp_dist_v_c){
                        // Update the dists_table
                        dists_table[v_id][cand_root_id] = query_dist_v_c;
                    }
                }
                if (need_activate) {
                    if (!is_active[v_id]) {
                        is_active[v_id] = true;
                        active_queue[end_active_queue++] = v_id;
                    }
                    if (!has_new_labels[v_id]) {
                        has_new_labels[v_id] = true;
                        has_new_labels_queue[end_has_new_labels_queue++] = v_id;
                    }
                }
            }
            // Reset vertices' candidates_que and candidates_dists
            for (NodeID i_queue = 0; i_queue < end_has_cand_queue; ++i_queue) {
                NodeID v_id = has_cand_queue[i_queue];
                has_candidates[v_id] = false; // reset has_candidates
                ShortIndex &SI_v = short_index[v_id];
                NodeID bound_cand_i = SI_v.end_candidates_que;
                for (NodeID cand_i = 0; cand_i < bound_cand_i; ++cand_i) {
                    NodeID cand_root_id = SI_v.candidates_que[cand_i];
                    SI_v.candidates_dists[cand_root_id] = INF_WEIGHT; // Reset candidates_dists after using in distance_query.
                }
                SI_v.end_candidates_que = 0; // Clear v_id's candidates_que
            }
            end_has_cand_queue = 0; // Set the has_cand_queue empty
        }

		// Filter out wrong and redundant labels in this batch
        filter_out_labels(
                has_new_labels_queue,
                end_has_new_labels_queue,
                short_index,
                roots_start);

         // Reset dists_table and short_index
        reset_tables(
			short_index,
			roots_start,
			roots_size,
			L,
			dists_table,
			roots_labels_buffer,
			once_had_cand_queue,
			end_once_had_cand_queue,
			once_had_cand);
        
        // Update the index according to labels_table
        update_index(
			L,
			short_index,
			roots_start,
			has_new_labels_queue,
			end_has_new_labels_queue,
			has_new_labels);

    }

};

/**
 * @description: 
 * @Author: wanjingyi
 * @Date: 2021-02-23 12:57:34
 * @param {*}
 * @return {*}
 */
class WeightedParaVertexCentricPLL:public WeightedVCPLL{
public:
    //****************variables**************
    int num_of_threads{5};
    //****************constructions**************
    WeightedParaVertexCentricPLL() {}

    WeightedParaVertexCentricPLL(WGraph& wgraph,int batchSize,int numThreads) {
        BATCH_SIZE=batchSize;
        num_of_threads=numThreads;
        construct(wgraph);
    }

    ~WeightedParaVertexCentricPLL() {}

protected:
    void construct(WGraph& wgraph){
        L.resize(numOfVertices);
        NodeID remainer = numOfVertices % BATCH_SIZE;
	    NodeID b_i_bound = numOfVertices / BATCH_SIZE;
        vector<NodeID> active_queue(numOfVertices);//Active queue
        NodeID end_active_queue = 0;
        // is_active[v] is true means vertex v is in the active queue.
        vector<uint8_t> is_active(numOfVertices, 0);	
        // Queue of vertices having candidates
        vector<NodeID> has_cand_queue(numOfVertices);
        NodeID end_has_cand_queue = 0;
        // has_candidates[v] is true means vertex v is in the queue has_cand_queue.
        vector<uint8_t> has_candidates(numOfVertices, 0);
        // Distance table of shortest distance from roots to other vertices.
        // The distance table is roots_sizes by N. 
        // 1. record the shortest distance so far from every root to every vertex; (from v to root r)
        vector< vector<EdgeWeight> > dists_table(numOfVertices, vector<EdgeWeight>(BATCH_SIZE, INF_WEIGHT));
        // r's old label distance from r to v
        //The distance buffer, recording old label distances of every root. It needs to be initialized every batch by labels of roots. (dists_table[r][l], r > l)
        vector< vector<EdgeWeight> > roots_labels_buffer(BATCH_SIZE, vector<EdgeWeight>(numOfVertices, INF_WEIGHT)); 
        // Here the size of short_index actually is fixed because it is static.
        vector<ShortIndex> short_index(numOfVertices, ShortIndex(BATCH_SIZE));
        // A queue to store vertices which have got new labels in this batch. This queue is used for update the index.
        vector<NodeID> has_new_labels_queue(numOfVertices);
        NodeID end_has_new_labels_queue = 0;
        vector<uint8_t> has_new_labels(numOfVertices, 0);
        // A queue to store vertices which have got candidates in this batch. This queue is used fore reset dists_table, accompanying with reached_roots_que in the ShortIndex.
        vector<NodeID> once_had_cand_queue(numOfVertices);
        NodeID end_once_had_cand_queue = 0;
        vector<uint8_t> once_had_cand(numOfVertices, 0);
        double hub_time=GetCurrentTimeSec();
        for (NodeID b_i = 0; b_i < b_i_bound; ++b_i) {
            vertex_centric_labeling_in_batches(
                    wgraph,
                    b_i,
                    b_i * BATCH_SIZE,
                    BATCH_SIZE,
                    L,
                    active_queue,
                    end_active_queue,
                    is_active,
                    has_cand_queue,
                    end_has_cand_queue,
                    has_candidates,
                    dists_table,
                    roots_labels_buffer,
                    short_index,
                    has_new_labels_queue,
                    end_has_new_labels_queue,
                    has_new_labels,
                    once_had_cand_queue,
                    end_once_had_cand_queue,
                    once_had_cand);
        }
        if (0 != remainer) {
            vertex_centric_labeling_in_batches(
                    wgraph,
                    b_i_bound,
                    b_i_bound * BATCH_SIZE,
                    remainer,
                    L,
                    active_queue,
                    end_active_queue,
                    is_active,
                    has_cand_queue,
                    end_has_cand_queue,
                    has_candidates,
                    dists_table,
                    roots_labels_buffer,
                    short_index,
                    has_new_labels_queue,
                    end_has_new_labels_queue,
                    has_new_labels,
                    once_had_cand_queue,
                    end_once_had_cand_queue,
                    once_had_cand);
        }
        hub_time=GetCurrentTimeSec()-hub_time;
        cout << "labeling time:" << hub_time *1e6 <<  " microseconds" << endl;
        double ave_hub_time=hub_time/(double) numOfVertices;
        cout<<"average hub time:"<<ave_hub_time*1e6 <<  " microseconds" << endl;
    }

    void initialize_tables(
			vector<ShortIndex> &short_index,
			vector< vector<EdgeWeight> > &roots_labels_buffer,
			vector<NodeID> &active_queue,
			NodeID &end_active_queue,
			NodeID roots_start,
			NodeID roots_size,
			vector<IndexType> &L,
			vector<NodeID> &has_new_labels_queue,
			NodeID &end_has_new_labels_queue,
			vector<uint8_t> &has_new_labels)
    {
            NodeID roots_bound = roots_start + roots_size;
            // **********Distance Tables**********
            {
            #pragma omp parallel for
            for (NodeID r_root_id = 0; r_root_id < roots_size; ++r_root_id) {// Traverse all roots
                NodeID r_real_id = r_root_id + roots_start;
                const IndexType &Lr = L[r_real_id];
                NodeID l_i_bound = Lr.vertices.size();
                _mm_prefetch(&Lr.vertices[0], _MM_HINT_T0);
                _mm_prefetch(&Lr.distances[0], _MM_HINT_T0);
                for (NodeID l_i = 0; l_i < l_i_bound; ++l_i) {
				roots_labels_buffer[r_root_id][Lr.vertices[l_i]] = Lr.distances[l_i];
			    }
                roots_labels_buffer[r_root_id][r_real_id] = 0;
            }
            }
            // **********Distance Tables**********

            // **********Short Index**********
            {
            #pragma omp parallel for
            for (NodeID r_real_id = roots_start; r_real_id < roots_bound; ++r_real_id) {
                ShortIndex &SI_r = short_index[r_real_id];
                NodeID r_root_id = r_real_id - roots_start;
                // Record itself as last inserted label distance
                SI_r.vertices_que[SI_r.end_vertices_que++] = r_root_id;
                SI_r.vertices_dists[r_root_id] = 0;
                SI_r.last_new_roots[SI_r.end_last_new_roots++] = r_root_id;
            }
            }

            // **********Active queue**********
            for (NodeID r_real_id = roots_start; r_real_id < roots_bound; ++r_real_id) {
                active_queue[end_active_queue++] = r_real_id;
            }
            // **********Active queue**********

            // **********Has new labels queue**********
            // has_new_labels_queue: Put all roots into the has_new_labels_queue
            for (NodeID r_real_id = roots_start; r_real_id < roots_bound; ++r_real_id) {
                has_new_labels_queue[end_has_new_labels_queue++] = r_real_id;
                has_new_labels[r_real_id] = true;
            }
            // **********Has new labels queue**********
    }

    // Function that pushes v_head's labels to v_head's every neighbor
    void send_messages(
        NodeID v_head,
        NodeID roots_start,
        WGraph &wgraph,
        vector< vector<EdgeWeight> > &dists_table,
        vector<ShortIndex> &short_index,
        vector<NodeID> &tmp_has_cand_queue,
        NodeID &size_tmp_has_cand_queue,
        NodeID offset_tmp_queue,
        vector<uint8_t> &has_candidates,
        vector<NodeID> &tmp_once_had_cand_queue,
        NodeID &size_tmp_once_had_cand_queue,
        vector<uint8_t> &once_had_cand
    ){
        ShortIndex &SI_v_head = short_index[v_head];
        vector<EdgeID>& wgraph_vertices = wgraph.vertices;
        vector<NodeEdgeWeightPair>& wgraph_edges = wgraph.edges;
        // Traverse v_head's every neighbor v_tail
        for (EdgeID eid = wgraph_vertices[v_head]; eid < wgraph_vertices[v_head + 1]; ++eid) {
            NodeID v_tail = wgraph_edges[eid].first;
            if (v_tail <= roots_start) { // v_tail has higher rank than any roots, then no roots can push new labels to it.
                break;
            }
            EdgeWeight weight_h_t = wgraph_edges[eid].second;
            ShortIndex &SI_v_tail = short_index[v_tail];
            bool got_candidates = false; // A flag indicates if v_tail got new candidates.
            // Traverse v_head's last inserted labels
            // Sequential Version
            NodeID bound_r_i = SI_v_head.end_last_new_roots;
            for (NodeID r_i = 0; r_i < bound_r_i; ++r_i) {
                NodeID r_root_id = SI_v_head.last_new_roots[r_i]; // last inserted label r_root_id of v_head
                if (v_tail <= r_root_id + roots_start) {
                    continue;
                } 
                EdgeWeight tmp_dist_r_t = SI_v_head.vertices_dists[r_root_id] + weight_h_t;
                if (tmp_dist_r_t < dists_table[v_tail][r_root_id] && tmp_dist_r_t < SI_v_tail.candidates_dists[r_root_id]) {
                    // Mark r_root_id as a candidate of v_tail
                    if (CAS<EdgeWeight>(&SI_v_tail.candidates_dists[r_root_id], (EdgeWeight) INF_WEIGHT, tmp_dist_r_t)) {
                        TS_enqueue<NodeID,NodeID>(SI_v_tail.candidates_que, SI_v_tail.end_candidates_que, r_root_id);
                    }
                    volatile EdgeWeight old_v = SI_v_tail.candidates_dists[r_root_id];
                    while (tmp_dist_r_t < old_v	&& !CAS<EdgeWeight>(&SI_v_tail.candidates_dists[r_root_id], old_v, tmp_dist_r_t)) {
                        old_v = SI_v_tail.candidates_dists[r_root_id];
                    }
                    if (CAS<EdgeWeight>(&dists_table[v_tail][r_root_id], (EdgeWeight) INF_WEIGHT, tmp_dist_r_t)) {
                        TS_enqueue<NodeID,NodeID>(SI_v_tail.reached_roots_que, SI_v_tail.end_reached_roots_que, r_root_id);
                    }
                    old_v = dists_table[v_tail][r_root_id];
                    while (tmp_dist_r_t < old_v && !CAS<EdgeWeight>(&dists_table[v_tail][r_root_id], old_v, tmp_dist_r_t)) {
                        old_v = dists_table[v_tail][r_root_id];
                    }
                    got_candidates = true;
                }
            }
            if (got_candidates) {
                if (CAS(&has_candidates[v_tail], (uint8_t) 0, (uint8_t) 1)) {
                    tmp_has_cand_queue[offset_tmp_queue + size_tmp_has_cand_queue++] = v_tail;
                }

                if (CAS(&once_had_cand[v_tail], (uint8_t) 0, (uint8_t) 1)) {
                    tmp_once_had_cand_queue[offset_tmp_queue + size_tmp_once_had_cand_queue++] = v_tail;
                }
            }
        }
        SI_v_head.end_last_new_roots = 0; // clear v_head's last_new_roots
    }

    // Function: thread-save enqueue. The queue has enough size already. An index points the end of the queue.
    template <typename T, typename Int>
    inline void TS_enqueue(
            vector<T> &queue,
            Int &end_queue,
            const T &e)
    {
        volatile Int old_i = end_queue;
        volatile Int new_i = old_i + 1;
        while (!CAS(&end_queue, old_i, new_i)) {
            old_i = end_queue;
            new_i = old_i + 1;
        }
        queue[old_i] = e;
    }
    template <typename T>
    inline void collect_into_queue(					
        vector<T> &tmp_queue,
        vector<NodeID> &offsets_tmp_queue, // the locations in tmp_queue for writing from tmp_queue
        vector<NodeID> &offsets_queue, // the locations in queue for writing into queue.
        NodeID num_elements, // total number of elements which need to be added from tmp_queue to queue
        vector<T> &queue,
        NodeID &end_queue)
    {
        if (0 == num_elements) {
            return;
        }
        NodeID i_bound = offsets_tmp_queue.size();
    #pragma omp parallel for
        for (NodeID i = 0; i < i_bound; ++i) {
            NodeID i_q_start = end_queue + offsets_queue[i];
            NodeID i_q_bound;
            if (i_bound - 1 != i) {
                i_q_bound = end_queue + offsets_queue[i + 1];
            } else {
                i_q_bound = end_queue + num_elements;
            }
            if (i_q_start == i_q_bound) {
                // If the group has no elements to be added, then continue to the next group
                continue;
            }
            NodeID end_tmp = offsets_tmp_queue[i];
            for (NodeID i_q = i_q_start; i_q < i_q_bound; ++i_q) {
                queue[i_q] = tmp_queue[end_tmp++];
            }
        }
        end_queue += num_elements;
    }

    // Function to get the prefix sum of elements in offsets
    NodeID prefix_sum_for_offsets(vector<NodeID> &offsets){
        NodeID size_offsets = offsets.size();
        if (1 == size_offsets){
            NodeID tmp = offsets[0];
            offsets[0] = 0;
            return tmp;
        }else if(size_offsets < 2048){
            NodeID offset_sum = 0;
            NodeID size = size_offsets;
            for (NodeID i = 0; i < size; ++i) {
                NodeID tmp = offsets[i];
                offsets[i] = offset_sum;
                offset_sum += tmp;
            }
            return offset_sum;
        }else{
            // Parallel Prefix Sum, based on Guy E. Blelloch's Prefix Sums and Their Applications
            NodeID last_element = offsets[size_offsets - 1];
            NodeID size = 1 << ((NodeID) log2(size_offsets));
            NodeID tmp_element = offsets[size - 1];
            // Up-Sweep (Reduce) Phase
            NodeID log2size = log2(size);
            for (NodeID d = 0; d < log2size; ++d) {
                NodeID by = 1 << (d + 1);
                #pragma omp parallel for
                for (NodeID k = 0; k < size; k += by) {
                    offsets[k + (1 << (d + 1)) - 1] += offsets[k + (1 << d) - 1];
                }
            }
            // Down-Sweep Phase
            offsets[size - 1] = 0;
            for (NodeID d = log2(size) - 1; d != (NodeID) -1 ; --d) {
                NodeID by = 1 << (d + 1);
                #pragma omp parallel for
                for (NodeID k = 0; k < size; k += by) {
                    NodeID t = offsets[k + (1 << d) - 1];
                    offsets[k + (1 << d) - 1] = offsets[k + (1 << (d + 1)) - 1];
                    offsets[k + (1 << (d + 1)) - 1] += t;
                }
            }
            if (size != size_offsets) {
                NodeID tmp_sum = offsets[size - 1] + tmp_element;
                for (NodeID i = size; i < size_offsets; ++i) {
                    NodeID t = offsets[i];
                    offsets[i] = tmp_sum;
                    tmp_sum += t;
                }
            }
            return offsets[size_offsets - 1] + last_element;     
        }
    }

    inline void filter_out_labels(
        const vector<NodeID> &has_new_labels_queue,
		NodeID end_has_new_labels_queue,
		vector<ShortIndex> &short_index,
		NodeID roots_start)
    {
        // Traverse has_new_labels_queue for every vertex who got new labels in this batch
        #pragma omp parallel for
        for (NodeID i_q = 0; i_q < end_has_new_labels_queue; ++i_q) {
            NodeID v_id = has_new_labels_queue[i_q];
            IndexType &Lv = L[v_id];
            ShortIndex &SI_v = short_index[v_id];
            NodeID bound_label_i = SI_v.end_vertices_que;
            // Traverse v_id's every new label
            for (NodeID i_c = 0; i_c < bound_label_i; ++i_c) {
                NodeID c_root_id = SI_v.vertices_que[i_c];
                EdgeWeight dist_v_c = SI_v.vertices_dists[c_root_id];
                ShortIndex &SI_c = short_index[c_root_id + roots_start];
                // Traverse v_id's other new label to see if has shorter distance
                for (NodeID i_r = 0; i_r < bound_label_i; ++i_r) {
                    NodeID r_root_id = SI_v.vertices_que[i_r];
                    if (c_root_id <= r_root_id) {
                        continue;
                    }
                    EdgeWeight dist_v_r = SI_v.vertices_dists[r_root_id];
                    if (dist_v_r > dist_v_c) {
                        continue;
                    }
                    EdgeWeight dist_c_r = SI_c.vertices_dists[r_root_id];
                    if (dist_c_r > dist_v_c) {
                        continue;
                    }
                    if (dist_v_r + dist_c_r <= dist_v_c) {
                        // Shorter distance exists, then label c is not necessary
                        SI_v.vertices_dists[c_root_id] = INF_WEIGHT; // Filter out c_root_id from v_id's labels
                    }
                }
            }
        }
    }    

    inline void update_index(
        vector<IndexType> &L,
		vector<ShortIndex> &short_index,
		NodeID roots_start,
		vector<NodeID> &has_new_labels_queue,
		NodeID &end_has_new_labels_queue,
		vector<uint8_t> &has_new_labels)
    {
        // Traverse has_new_labels_queue for every vertex who got new labels in this batch
        #pragma omp parallel for
        for (NodeID i_q = 0; i_q < end_has_new_labels_queue; ++i_q) {
            NodeID v_id = has_new_labels_queue[i_q];
            has_new_labels[v_id] = 0; // Reset has_new_labels
            IndexType &Lv = L[v_id];
            ShortIndex &SI_v = short_index[v_id];
            NodeID bound_i_r = SI_v.end_vertices_que;
            for (NodeID i_r = 0; i_r < bound_i_r; ++i_r) {
                NodeID r_root_id = SI_v.vertices_que[i_r];
                EdgeWeight dist = SI_v.vertices_dists[r_root_id];
                if (INF_WEIGHT == dist) {
                    continue;
                }
                SI_v.vertices_dists[r_root_id] = INF_WEIGHT; // Reset v_id's vertices_dists
                Lv.vertices.push_back(r_root_id + roots_start);
                Lv.distances.push_back(dist);
            }
            SI_v.end_vertices_que = 0; // Clear v_id's vertices_que
        }
        end_has_new_labels_queue = 0; // Clear has_new_labels_queue
    }

    inline EdgeWeight distance_query(
        NodeID v_id,
        NodeID cand_root_id,
        const vector< vector<EdgeWeight> > &roots_labels_buffer,
        const vector<ShortIndex> &short_index,
        const vector<IndexType> &L,
        NodeID roots_start,
        EdgeWeight tmp_dist_v_c)
    {
        // 1. Labels in L[v]
        NodeID cand_real_id = cand_root_id + roots_start;
	    const IndexType &Lv = L[v_id];
        NodeID bound_i_l = Lv.vertices.size();
        for (NodeID i_l = 0; i_l < bound_i_l; ++i_l) {
            NodeID r = Lv.vertices[i_l];
            if (cand_real_id <= r || INF_WEIGHT == roots_labels_buffer[cand_root_id][r]) {
                continue;
            }
            EdgeWeight label_dist_v_c = Lv.distances[i_l] + roots_labels_buffer[cand_root_id][r];
            if (label_dist_v_c <= tmp_dist_v_c) {
                //			distance_query_time += WallTimer::get_time_mark();
                //++l_l_hit_count;
                return label_dist_v_c;
            }
        }

        // 2. Labels in short_index[v_id].vertices_que
        const ShortIndex &SI_v = short_index[v_id];
        const ShortIndex &SI_c = short_index[cand_root_id];
        NodeID bound_i_que = SI_v.end_vertices_que;
        for (NodeID i_que = 0; i_que < bound_i_que; ++i_que) {
            NodeID r_root_id = SI_v.vertices_que[i_que];
            if (cand_real_id <= r_root_id + roots_start) {
                continue;
            }
            // Check r_root_id in cand_root_id's labels inserted in this batch
            if (INF_WEIGHT != SI_c.vertices_dists[r_root_id]) {
                EdgeWeight label_dist_v_c = SI_v.vertices_dists[r_root_id] + SI_c.vertices_dists[r_root_id];
                if (label_dist_v_c <= tmp_dist_v_c) {
                    return label_dist_v_c;
                }
            }
        }
        return INF_WEIGHT;
    }

    inline void reset_tables(
        vector<ShortIndex> &short_index,
		NodeID roots_start,
		NodeID roots_size,
		const vector<IndexType> &L,
		vector< vector<EdgeWeight> > &dists_table,
		vector< vector<EdgeWeight> > &roots_labels_buffer,
		vector<NodeID> &once_had_cand_queue,
		NodeID &end_once_had_cand_queue,
		vector<uint8_t> &once_had_cand
    )
    {
	    // Reset roots_labels_buffer according to L (old labels)
        #pragma omp parallel for
        for (NodeID r_roots_id = 0; r_roots_id < roots_size; ++r_roots_id) {
            NodeID r_real_id = r_roots_id + roots_start;
            // Traverse labels of r
            const IndexType &Lr = L[r_real_id];
            NodeID bound_i_l = Lr.vertices.size();
            for (NodeID i_l = 0; i_l < bound_i_l; ++i_l) {
                roots_labels_buffer[r_roots_id][Lr.vertices[i_l]] = INF_WEIGHT;
            }
            roots_labels_buffer[r_roots_id][r_real_id] = INF_WEIGHT;
        }

        // Reset dists_table according to short_index[v].reached_roots_que
        for (NodeID i_q = 0; i_q < end_once_had_cand_queue; ++i_q) {
            NodeID v_id = once_had_cand_queue[i_q];
            once_had_cand[v_id] = false; // Reset once_had_cand
            ShortIndex &SI_v = short_index[v_id];
            // Traverse roots which have reached v_id
            NodeID bound_i_r = SI_v.end_reached_roots_que;
            for (NodeID i_r = 0; i_r < bound_i_r; ++i_r) {
                dists_table[v_id][SI_v.reached_roots_que[i_r]] = INF_WEIGHT;
            }
            SI_v.end_reached_roots_que = 0; // Clear v_id's reached_roots_que
        }
        end_once_had_cand_queue = 0; // Clear once_had_cand_queue.
    }

    void vertex_centric_labeling_in_batches(
        WGraph& wgraph,
        NodeID b_id,
        NodeID roots_start,// start id of roots
        NodeID roots_size, // how many roots in the batch
        vector<IndexType> &L,
        vector<NodeID> &active_queue,
        NodeID end_active_queue,
		vector<uint8_t> &is_active,
        vector<NodeID> &has_cand_queue,
        NodeID end_has_cand_queue,
        vector<uint8_t> &has_candidates,
        vector< vector<EdgeWeight> > &dists_table,
        vector< vector<EdgeWeight> > &roots_labels_buffer,
        vector<ShortIndex> &short_index,
        vector<NodeID> &has_new_labels_queue,
        NodeID end_has_new_labels_queue,
        vector<uint8_t> &has_new_labels,
        vector<NodeID> &once_had_cand_queue,
        NodeID end_once_had_cand_queue,
        vector<uint8_t> &once_had_cand
    )
    {
        omp_set_num_threads(num_of_threads);
        // At the beginning of a batch, initialize the labels L and distance buffer dist_matrix;
        initialize_tables(
                short_index,
                roots_labels_buffer,
                active_queue,
                end_active_queue,
                roots_start,
                roots_size,
                L,
                has_new_labels_queue,
                end_has_new_labels_queue,
                has_new_labels);

        while (0 != end_active_queue) {
            //*****************first stage sending distances***************
            // Prepare for parallel processing the active_queue and adding to candidate_queue.
            // Every vertex's offset location in tmp_candidate_queue
            // It's used for every thread to write into tmp_candidate_queue and tmp_once_candidated_queue
            vector<NodeID> offsets_tmp_queue(end_active_queue);
            #pragma omp parallel for
           	for (NodeID i_queue = 0; i_queue < end_active_queue; ++i_queue) {
				// Traverse all active vertices, get their out degrees.
                NodeID v=active_queue[i_queue];
				offsets_tmp_queue[i_queue] = wgraph.vertices[v+1]-wgraph.vertices[v];
			} 
            NodeID num_neighbors = prefix_sum_for_offsets(offsets_tmp_queue);
            // every thread writes to tmp_candidate_queue at its offset location
			vector<NodeID> tmp_has_cand_queue(num_neighbors);
			// A vector to store the true number of pushed neighbors of every active vertex.
			vector<NodeID> sizes_tmp_has_cand_queue(end_active_queue, 0);
            vector<NodeID> tmp_once_had_cand_queue(num_neighbors);
			vector<NodeID> sizes_tmp_once_had_cand_queue(end_active_queue, 0);
            // Traverse the active queue, every active vertex sends distances to its neighbors
            //omp_set_num_threads(num_of_threads);
            #pragma omp parallel for
            for (NodeID i_queue = 0; i_queue < end_active_queue; ++i_queue) {
				NodeID v_head = active_queue[i_queue];
				is_active[v_head] = 0; // reset is_active
				send_messages(
						v_head,
						roots_start,
						wgraph,
						dists_table,
						short_index,
						tmp_has_cand_queue,
						sizes_tmp_has_cand_queue[i_queue],
						offsets_tmp_queue[i_queue],
						has_candidates,
						tmp_once_had_cand_queue,
						sizes_tmp_once_had_cand_queue[i_queue],
						once_had_cand);
			}
			// According to sizes_tmp_candidate_queue, get the offset for inserting to the real queue
			NodeID total_new = prefix_sum_for_offsets(sizes_tmp_has_cand_queue);
            // Collect all candidate vertices from tmp_candidate_queue into candidate_queue.
            collect_into_queue<NodeID>(
                tmp_has_cand_queue,
                offsets_tmp_queue, // the locations in tmp_queue for writing from tmp_queue
                sizes_tmp_has_cand_queue, // the locations in queue for writing into queue.
                total_new, // total number of elements which need to be added from tmp_queue to queue
                has_cand_queue,
                end_has_cand_queue);
            total_new = prefix_sum_for_offsets(sizes_tmp_once_had_cand_queue);
            collect_into_queue<NodeID>(
					tmp_once_had_cand_queue,
					offsets_tmp_queue,
					sizes_tmp_once_had_cand_queue,
					total_new,
					once_had_cand_queue,
					end_once_had_cand_queue);
			end_active_queue = 0; // Set the active_queue empty
            //*****************second stage-adding***************
            offsets_tmp_queue.clear();
            offsets_tmp_queue.resize(end_has_cand_queue);
            //vector<NodeID> offsets_tmp_queue(end_has_cand_queue);
            #pragma omp parallel for
            for (NodeID i_queue = 0; i_queue < end_has_cand_queue; ++i_queue) {
				// Traverse all active vertices, get their out degrees.
				offsets_tmp_queue[i_queue] = i_queue;
			}
            // every thread writes to tmp_active_queue at its offset location
			vector<NodeID> tmp_active_queue(end_has_cand_queue);
            // A vector to store the true number of pushed neighbors of every active vertex.
			vector<NodeID> sizes_tmp_active_queue(end_has_cand_queue, 0);
            // every thread writes to its offset location.
			vector<NodeID> tmp_has_new_labels_queue(end_has_cand_queue);
            // Store the size of every group recorded by every thread.
			vector<NodeID> sizes_tmp_has_new_labels_queue(end_has_cand_queue, 0);
			// Traverse vertices in the has_cand_queue to insert labels
            #pragma omp parallel for
            for (NodeID i_queue = 0; i_queue < end_has_cand_queue; ++i_queue) {
				NodeID v_id = has_cand_queue[i_queue];
				bool need_activate = 0;
				ShortIndex &SI_v = short_index[v_id];
				// Traverse v_id's all candidates
				NodeID bound_cand_i = SI_v.end_candidates_que;
				for (NodeID cand_i = 0; cand_i < bound_cand_i; ++cand_i) {
					NodeID cand_root_id = SI_v.candidates_que[cand_i];
					EdgeWeight tmp_dist_v_c = SI_v.candidates_dists[cand_root_id];
					// Distance check for pruning
					EdgeWeight query_dist_v_c;
					if (INF_WEIGHT ==
							(query_dist_v_c = distance_query(
									v_id,
									cand_root_id,
									roots_labels_buffer,
									short_index,
									L,
									roots_start,
									tmp_dist_v_c))) {
						if (INF_WEIGHT == SI_v.vertices_dists[cand_root_id]) {
							// Record cand_root_id as v_id's label
							SI_v.vertices_que[SI_v.end_vertices_que++] = cand_root_id;

						}
							// Record the new distance in the label table
						SI_v.vertices_dists[cand_root_id] = tmp_dist_v_c;
						SI_v.last_new_roots[SI_v.end_last_new_roots++] = cand_root_id;
						need_activate = 1;
					} else if (query_dist_v_c < tmp_dist_v_c){
						// Update the dists_table
						dists_table[v_id][cand_root_id] = query_dist_v_c;
					}
				}
				if (need_activate) {
					if (!is_active[v_id]) {
						is_active[v_id] = 1;
						tmp_active_queue[offsets_tmp_queue[i_queue] + sizes_tmp_active_queue[i_queue]++] = v_id;
						//						active_queue[end_active_queue++] = v_id;
					}
					if (!has_new_labels[v_id]) {
						has_new_labels[v_id] = 1;
						tmp_has_new_labels_queue[offsets_tmp_queue[i_queue] + sizes_tmp_has_new_labels_queue[i_queue]++] = v_id;
						//						has_new_labels_queue[end_has_new_labels_queue++] = v_id;
					}
				}
			}
            // According to sizes_tmp_active_queue, get the offset for inserting to the real queue
			NodeID total_new_1 = prefix_sum_for_offsets(sizes_tmp_active_queue);
            // Collect all candidate vertices from tmp_active_queue into active_queue.
			collect_into_queue<NodeID>(
						tmp_active_queue,
						offsets_tmp_queue, // the locations in tmp_queue for writing from tmp_queue
						sizes_tmp_active_queue, // the locations in queue for writing into queue.
						total_new_1, // total number of elements which need to be added from tmp_queue to queue
						active_queue,
						end_active_queue);
            // According to sizes_tmp_active_queue, get the offset for inserting to the real queue
			total_new_1 = prefix_sum_for_offsets(sizes_tmp_has_new_labels_queue);
            // Collect all candidate vertices from tmp_candidate_queue into candidate_queue.
			collect_into_queue<NodeID>(
						tmp_has_new_labels_queue,
						offsets_tmp_queue, // the locations in tmp_queue for writing from tmp_queue
						sizes_tmp_has_new_labels_queue, // the locations in queue for writing into queue.
						total_new_1, // total number of elements which need to be added from tmp_queue to queue
						has_new_labels_queue,
						end_has_new_labels_queue);
            // Reset vertices' candidates_que and candidates_dists
            #pragma omp parallel for
            for (NodeID i_queue = 0; i_queue < end_has_cand_queue; ++i_queue) {
				NodeID v_id = has_cand_queue[i_queue];
				has_candidates[v_id] = 0; // reset has_candidates
				ShortIndex &SI_v = short_index[v_id];
                // Sequential Version
				NodeID bound_cand_i = SI_v.end_candidates_que;
				for (NodeID cand_i = 0; cand_i < bound_cand_i; ++cand_i) {
					NodeID cand_root_id = SI_v.candidates_que[cand_i];
					SI_v.candidates_dists[cand_root_id] = INF_WEIGHT; // Reset candidates_dists after using in distance_query.
				}
				SI_v.end_candidates_que = 0; // Clear v_id's candidates_que
            }
            end_has_cand_queue = 0; // Set the has_cand_queue empty
        }
        // Filter out wrong and redundant labels in this batch
        filter_out_labels(
			has_new_labels_queue,
			end_has_new_labels_queue,
			short_index,
			roots_start);
        
        // Reset dists_table and short_index
        reset_tables(
			short_index,
			roots_start,
			roots_size,
			L,
			dists_table,
			roots_labels_buffer,
			once_had_cand_queue,
			end_once_had_cand_queue,
			once_had_cand);

        // Update the index according to labels_table
	    update_index(
			L,
			short_index,
			roots_start,
			has_new_labels_queue,
			end_has_new_labels_queue,
			has_new_labels);
    }

};

}

#endif // !_PARALLEL_CONSTRUCTION_H