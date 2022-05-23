/*
 * @Descripttion:  define the label class and its util functions
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-10 21:13:18
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-02-14 10:24:30
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

#ifndef LABELS_H
#define LABELS_H


#include <limits>
#include <climits>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include "graph.h"
#include "paras.h" 
#include "./utils.h"
#include <malloc.h>
#include <xmmintrin.h>
//typedef unsigned __int64 BPSeed;
#include <omp.h>
#include<bitset>
#include <sstream>
#include <string>
#include <math.h>
#include <iomanip>
#ifdef _WIN32
	#include<google/sparsehash/sparseconfig.h>
#endif
    #include<google/dense_hash_map>
	#include<google/dense_hash_set>

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define INF_WEIGHT SP_Constants::INF_WEIGHT
#define AVG_LABEL_SIZE 200
#define DIVISION_FACTOR 100

struct index_pivot {
	int* pointNo;
	int pointNum;
};//每个pivot对应的所有结点编号的数据结构 modified by wanjingyi

struct index_t {
	vector<NodeID> spt_v;
	vector<EdgeWeight> spt_d;

	NodeID size() {
		return spt_v.size();
	}

};

struct index_t_p {
	NodeID* spt_v;
	EdgeWeight* spt_d;
}__attribute__((aligned(64)));  // Aligned for cache lines;

//label structure for cache defined by wanjingyi
struct index_cache_t_p {
	NodeID* spt_v;
	EdgeWeight* spt_d;
	google::dense_hash_map<NodeID, EdgeWeight> dis_cache_m;
}__attribute__((aligned(64)));  // Aligned for cache lines;

//label structure for parallel defined by wanjingyi
struct index_parallel_t_p {
	NodeID* spt_v;
	EdgeWeight* spt_d;
	NodeID lsiz;
	NodeID* index;//store the parallel index for sections
}__attribute__((aligned(64)));  // Aligned for cache lines;
// struct index_parallel_t_p {
// 	NodeID* spt_v_l;
// 	EdgeWeight* spt_d_l;
// 	NodeID* spt_v_h;
// 	EdgeWeight* spt_d_h;
// };  // Aligned for cache lines;


struct two_index_t_p {
	NodeID* spt_v;
	EdgeWeight* spt_d;
	uint8_t* spt_lv;
	EdgeWeight* spt_ld;
}__attribute__((aligned(64)));  // Aligned for cache lines;

struct index_t_path {
	vector<NodeID> spt_v;
	vector<NodeID> spt_p;//parent nodes
	vector<EdgeWeight> spt_d;

	NodeID size() {
		return spt_v.size();
	}

};

struct index_t_path_p {
	NodeID* spt_v;
	NodeID* spt_p;
	EdgeWeight* spt_d;
};

struct query_info {
	NodeID meet_node;
	NodeID search_len;
	double time_cost;
	EdgeWeight distance;
};

template<int kNumBitParallelRoots = 50>
struct index_t_bp {
	NodeID* spt_v;
	EdgeWeight* spt_d;
	EdgeWeight bpspt_d[kNumBitParallelRoots];
	uint64_t bpspt_s[kNumBitParallelRoots][2];
}__attribute__((aligned(64)));  // Aligned for cache lines;


struct token_t {
	NodeID* sptc_v; // sptc_v[0] is the root
	EdgeWeight* sptc_d;	 // |*| = k + 1, sptc_d[0] is the number of children - k
	unsigned char* sptc_fbv; // first-level bit vector
	unsigned char* sptc_sbv; // second-level bit vector
	NodeID* sptc_pathv; // intermediate point for a path
}__attribute__((aligned(64)));

//*************written by wanjingyi***************
struct queryPair{
	NodeID s,t;
	int freq;
	queryPair(NodeID s1,NodeID t1,int f) : s(s1), t(t1), freq(f){ }
	bool operator < (const queryPair& qp){
		if(qp.freq!=freq) return freq<qp.freq;
		else if(qp.s!=s) return s<qp.s;
		else return t<qp.t;
	}
};

bool cmp(pair<NodeID, int> a, pair<NodeID, int> b)
{
	if(a.second!=b.second) return a.second<b.second;//根据second的值升序排序
	else return a.first<b.first;
}

class Label {

public:
	vector<index_t> index_;	
	index_t_p* index_p;
	two_index_t_p* two_index_p;


	double GetCurrentTimeSec() {
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec + tv.tv_usec * 1e-6;
	}


	Label() {
		index_.resize(numOfVertices);
	}

	~Label() {
		Free();
	}

	EdgeWeight query_p(NodeID s, NodeID t) {
		//
		//EdgeWeight distance = INF_WEIGHT;

		//NodeID *vs = index_p[s].spt_v;
		//NodeID *vt = index_p[t].spt_v;
		//EdgeWeight* ws = index_p[s].spt_d;
		//EdgeWeight* wt = index_p[t].spt_d;

		//_mm_prefetch(vs, _MM_HINT_T0);
		//_mm_prefetch(vt, _MM_HINT_T0);
		//_mm_prefetch(ws, _MM_HINT_T0);
		//_mm_prefetch(wt, _MM_HINT_T0);

		//for (unsigned i = 0, j = 0; ; ) {
		//	if (*(vs + i) == *(vt + j)) {
		//		if (*(vs + i) == numOfVertices) break;  // Sentinel
		//		EdgeWeight td = *(ws + i) + *(wt + j);
		//		if (td < distance) distance = td;
		//		++i;
		//		++j;
		//	}
		//	else {
		//		i += *(vs + i) < *(vt + j) ? 1 : 0;
		//		j += *(vs + i) > *(vt + j) ? 1 : 0;
		//	}
		//}
		//return distance;

		EdgeWeight distance = INF_WEIGHT;

		const index_t_p &idx_s = index_p[s];
		const index_t_p &idx_t = index_p[t];

		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) distance = td;
				++i;
				++j;
			} 
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}
		return distance;
	}

	EdgeWeight two_query_p_sequential(NodeID s, NodeID t) {
		
		EdgeWeight distance = INF_WEIGHT;
		EdgeWeight ldistance = INF_WEIGHT;

		const two_index_t_p &idx_s = two_index_p[s];
		const two_index_t_p &idx_t = two_index_p[t];

		_mm_prefetch(&idx_s.spt_lv[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_lv[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_ld[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_ld[0], _MM_HINT_T0);

		for (uint8_t i = 0, j = 0; ; ) {
			uint8_t uv8_1 = idx_s.spt_lv[i], uv8_2 = idx_t.spt_lv[j];

			if (uv8_1 == UCHAR_MAX) break;  // Sentinel

			if (uv8_1 == uv8_2) {
				EdgeWeight td = idx_s.spt_ld[i] + idx_t.spt_ld[j];
				if (td < ldistance) ldistance = td;
				++i;
				++j;
			}
			else {
				i += uv8_1 < uv8_2 ? 1 : 0;
				j += uv8_1 > uv8_2 ? 1 : 0;
			}
		}
	
		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) distance = td;
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}
		
		if(distance < ldistance) 
			return distance;
		else
			return ldistance;		
	}

	
	EdgeWeight two_query_p_parallel(NodeID s, NodeID t) {
		
		EdgeWeight distance = INF_WEIGHT;
		EdgeWeight ldistance = INF_WEIGHT;

		const two_index_t_p &idx_s = two_index_p[s];
		const two_index_t_p &idx_t = two_index_p[t];

		
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

				for (int i = 0, j = 0; ; ) {
					NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

					if (v1 == numOfVertices) break;  // Sentinel

					if (v1 == v2) {
						EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
						if (td < distance) distance = td;
						++i;
						++j;
					}
					else {
						i += v1 < v2 ? 1 : 0;
						j += v1 > v2 ? 1 : 0;
					}
				}
			}
			
			#pragma omp section
			{
				_mm_prefetch(&idx_s.spt_lv[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_lv[0], _MM_HINT_T0);
				_mm_prefetch(&idx_s.spt_ld[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_ld[0], _MM_HINT_T0);

				for (uint8_t i = 0, j = 0; ; ) {
					uint8_t uv8_1 = idx_s.spt_lv[i], uv8_2 = idx_t.spt_lv[j];

					if (uv8_1 == UCHAR_MAX) break;  // Sentinel

					if (uv8_1 == uv8_2) {
						EdgeWeight td = idx_s.spt_ld[i] + idx_t.spt_ld[j];
						if (td < ldistance) ldistance = td;
						++i;
						++j;
					}
					else {
						i += uv8_1 < uv8_2 ? 1 : 0;
						j += uv8_1 > uv8_2 ? 1 : 0;
					}
				}
			}
		}
		if(distance < ldistance) 
			return distance;
		else
			return ldistance;		
	}
	
	EdgeWeight query_p_with_nums(NodeID s, NodeID t, int k) {
		//
		//EdgeWeight distance = INF_WEIGHT;

		//NodeID *vs = index_p[s].spt_v;
		//NodeID *vt = index_p[t].spt_v;
		//EdgeWeight* ws = index_p[s].spt_d;
		//EdgeWeight* wt = index_p[t].spt_d;

		//_mm_prefetch(vs, _MM_HINT_T0);
		//_mm_prefetch(vt, _MM_HINT_T0);
		//_mm_prefetch(ws, _MM_HINT_T0);
		//_mm_prefetch(wt, _MM_HINT_T0);

		//for (unsigned i = 0, j = 0; ; ) {
		//	if (*(vs + i) == *(vt + j)) {
		//		if (*(vs + i) == numOfVertices) break;  // Sentinel
		//		EdgeWeight td = *(ws + i) + *(wt + j);
		//		if (td < distance) distance = td;
		//		++i;
		//		++j;
		//	}
		//	else {
		//		i += *(vs + i) < *(vt + j) ? 1 : 0;
		//		j += *(vs + i) > *(vt + j) ? 1 : 0;
		//	}
		//}
		//return distance;


		EdgeWeight distance = INF_WEIGHT;

		const index_t_p &idx_s = index_p[s];
		const index_t_p &idx_t = index_p[t];

		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
		int k1 = k, k2 = k;
		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) distance = td;
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}

			if (i > k1 || j > k2) break;
		}
		return distance;
	}

	//used to query distances at first stage
	EdgeWeight query(NodeID s, NodeID t) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;

		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j]) 
				distance = min(distance, (EdgeWeight)(index_s_d[i++] + index_t_d[j++]));
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		}
		return distance;
	}
	
	

	EdgeWeight query(NodeID s, NodeID t, NodeID& meet, EdgeWeight& dis1, EdgeWeight& dis2) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;
		meet = numeric_limits<NodeID>::max();
		dis1 = numeric_limits<EdgeWeight>::max();
		dis2 = numeric_limits<EdgeWeight>::max();
		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j]) {
				if (distance > (EdgeWeight)(index_s_d[i] + index_t_d[j])) {
					distance = (EdgeWeight)(index_s_d[i] + index_t_d[j]);
					meet = index_s[i];
					dis1 = index_s_d[i];
					dis2 = index_t_d[j];
				}
				++i; ++j;
			}
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		}
		return distance;
	}

	/*EdgeWeight query_new(NodeID s, NodeID t, Ordering& ordering) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;



		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j])
				distance = min(distance, (EdgeWeight)(index_s_d[i++] + index_t_d[j++]));
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		}
		return distance;
	}
	*/
	double avg_size() {
		
			double total = 0;
		if(index_.size()!=0){
			for (int i = 0; i < numOfVertices; ++i) total += index_[i].spt_v.size();

			double avg = total / numOfVertices - 1; // We do not count the trivial label (V, INF_WEIGHT).
			return avg;
		}
		
		total = 0;
		for (int i = 0; i < numOfVertices; ++i) {
			int unit_count = 0;
			const index_t_p &idx_s = index_p[i];
			for(int j = 0; ;){
				NodeID v = idx_s.spt_v[j++];
				++unit_count;
				if( v == numOfVertices) break;
			}
			total += unit_count;
		}

		double avg = total / numOfVertices - 1; // We do not count the trivial label (V, INF_WEIGHT).

		return avg;
	}
	/*
	NodeID max_size() {
		NodeID maxsize = numeric_limits<NodeID>::min();
		for (int i = 0; i < V; ++i) maxsize = max(maxsize, index_[i].spt_v.size());
		return maxsize;
	}*/

	void append(NodeID v, NodeID root, EdgeWeight distance) {
		index_[v].spt_v.push_back(root);
		index_[v].spt_d.push_back(distance);
	}

	void print_stat() {
		cout << "Average Label Size: " << avg_size() << endl;
		//cout << "Maximum Label Size: " << max_size() << endl;
	}

	void Free() {
		if (index_.size() == 0) return;
		for (int v = 0; v < numOfVertices; ++v) {
			index_[v].spt_v.clear();
			index_[v].spt_d.clear();
		}
		index_.clear();
	}

		 /**
		  * function used to serialize vertice's labels to file,format as follows:
		  * label_size v (w0,d(v,w0)) (w1,d(v,w1))...(wn,d(v,wn))
		  * written by wanjingyi
		  * */
		 void write_labels(const char* write_filename)
		 {
			 ofstream ofs(write_filename);
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
		 }

	void save_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

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
	
	void load_labels(const char* load_filename) {
		if(index_p)
		{
			for (NodeID v = 0; v < numOfVertices; ++v) {
					free(index_p[v].spt_v);
					free(index_p[v].spt_d);
				} 
			free(index_p);
		}	
		index_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		cout<<"numOfVertices = "<<numOfVertices<<endl;

		index_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));
		


		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));
			
			idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));

		//	index_[v].spt_v.resize(isize);
		//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_v[i] = hub;
				idx.spt_d[i] = hub_weight;

			}
		}
		ifs.close();

		/*
		index_.clear();
		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		index_.resize(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {

			ifs.read((char*)&isize, sizeof(isize));
			index_[v].spt_v.resize(isize);
			index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < index_[v].size(); ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				index_[v].spt_v[i] = hub;
				index_[v].spt_d[i] = hub_weight;
			}
		}
		ifs.close();
		*/
	}

	void load_labels(const char* load_filename,const vector<bool>& flag) {
		index_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		//cout<<"numOfVertices = "<<numOfVertices<<endl;
		index_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));
	
		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));
			if(flag[v]){
				idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
				idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
			}

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;
				if(flag[v]){
					idx.spt_v[i] = hub;
					idx.spt_d[i] = hub_weight;
				}
			}
		}
		ifs.close();
	}

	void convert_to_fewerbit(){
		
		two_index_p = NULL;
		two_index_p = (two_index_t_p*)memalign(64, numOfVertices * sizeof(two_index_t_p));

		double compressed_size = 0;
		double total_size = 0;
		for (NodeID v = 0; v < numOfVertices; ++v) {
			two_index_t_p &idx = two_index_p[v];
			
			index_t_p &idx_original = index_p[v];
			
			NodeID isize = 0;
			for(NodeID i = 0; idx_original.spt_v[i] < UCHAR_MAX; ++i){
				++isize;
			}
			

			idx.spt_lv = (uint8_t*)memalign(64, (isize + 1) * sizeof(uint8_t));
			idx.spt_ld = (EdgeWeight*)memalign(64, (isize + 1) * sizeof(EdgeWeight));

		//	index_[v].spt_v.resize(isize);
		//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_lv[i] = idx_original.spt_v[i];
				idx.spt_ld[i] = idx_original.spt_d[i];

			}
			
			compressed_size += 4 * (isize - 1)  - isize;
			
			idx.spt_lv[isize] = UCHAR_MAX;
			idx.spt_ld[isize] = INF_WEIGHT;
			
			NodeID larger_size = 0;
			for(NodeID i = isize; idx_original.spt_v[i] != numOfVertices; ++i){
				++larger_size;
			}
			
			larger_size++;
			
			idx.spt_v = (NodeID*)memalign(64, larger_size * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, larger_size * sizeof(EdgeWeight));
			
			for (NodeID i = 0; i < larger_size; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_v[i] = idx_original.spt_v[i + isize];
				idx.spt_d[i] = idx_original.spt_d[i + isize];
			}			
			
			total_size += 4 * (isize - 1 + larger_size) * 2;
			
		}		
		cout << "reduce size :" << compressed_size << " out of " << total_size << " saving " << int(compressed_size * 100 / total_size) << "%" << endl;

	}
	
	void load_labels_with_k(const char* load_filename, int k) {
		/*	for (NodeID v = 0; v < numOfVertices; ++v) {
		free(index_p[v].spt_v);
		free(index_p[v].spt_d);
		}
		*/
		//free(index_p);

		long total_amount = 0;
		long actual_amount = 0;

		index_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		index_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));



		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));
			int actual_isize = k;
			if (isize > k) actual_isize = k;
			else actual_isize = isize;

			total_amount += isize;
			actual_amount += actual_isize;

			idx.spt_v = (NodeID*)memalign(64, actual_isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, actual_isize * sizeof(EdgeWeight));

			//	index_[v].spt_v.resize(isize);
			//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				if (i > actual_isize) continue;
				if (i == actual_isize - 1) {
					idx.spt_v[i] = numOfVertices;
					idx.spt_d[i] = INF_WEIGHT;
				}else {
					idx.spt_v[i] = hub;
					idx.spt_d[i] = hub_weight;
				}
			}
		}
		ifs.close();

		cout << "Total Labels:" << total_amount << endl;
		cout << "Actual Labels:" << actual_amount << endl;

		/*
		index_.clear();
		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		index_.resize(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {

		ifs.read((char*)&isize, sizeof(isize));
		index_[v].spt_v.resize(isize);
		index_[v].spt_d.resize(isize);

		for (NodeID i = 0; i < index_[v].size(); ++i) {
		NodeID hub;
		EdgeWeight hub_weight;
		ifs.read((char*)&hub, sizeof(hub));
		ifs.read((char*)&hub_weight, sizeof(hub_weight));
		index_[v].spt_v[i] = hub;
		index_[v].spt_d[i] = hub_weight;
		}
		}
		ifs.close();
		*/
	}


	void save_labels_iteration_stats(const char* save_filename) {

		vector<NodeID> stat(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < index_[v].size(); ++i)
				stat[index_[v].spt_v[i]]++;
		}

		ofstream ofs(save_filename);

		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs << stat[v] << endl;
		}
		ofs.close();
	}


	EdgeWeight query_with_info(NodeID s, NodeID t, query_info& q_info) {

		double stime = GetCurrentTimeSec();

		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& index_t = index_[t].spt_v;
		vector<EdgeWeight>& index_t_d = index_[t].spt_d;		

		q_info.meet_node = numOfVertices;
		double meet_distance;

		for (int i = 0, j = 0; i < index_s.size(), j < index_t.size(); ) {
			if (index_s[i] == index_t[j]) {
				meet_distance = (EdgeWeight)(index_s_d[i++] + index_t_d[j++]);
				if ( distance >  meet_distance) {
					distance = meet_distance;
					q_info.meet_node = index_s[i];
				}
			}
			else {
				if (index_s[i] < index_t[j])
					++i;
				else
					++j;
			}
		};

		stime = GetCurrentTimeSec() - stime;

		q_info.time_cost = stime;

		if (index_s.size() < index_t.size())
			q_info.search_len = index_s.size();
		else
			q_info.search_len = index_t.size();

		return distance;
	}

	/**
	 * function used to save label size to file
	 * written by wanjingyi
	 * */
	void save_label_size(const char* label_size_file) {
		string labelSizefile_prefix(label_size_file);
		string labelSizefile=labelSizefile_prefix.append(".size");
		ofstream ofs(labelSizefile);
		if(!ofs.is_open()) {cerr<<"Cannot open "<<labelSizefile<<endl;}
		for (int i = 0; i < numOfVertices; ++i) {
			ofs << index_[i].size()-1 << endl;
		}
		ofs.close();
	}

	/**
	 * function used to save label size to file
	 * written by wanjingyi
	 * */
	void save_label_size(const char* label_size_file,const vector<NodeID>& inv) {
		string labelSizefile_prefix(label_size_file);
		string labelSizefile=labelSizefile_prefix.append(".size");
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

	/**
	 * function used to serialize vertice's labels to file,format as follows:
	 * label_size v (w0,d(v,w0)) (w1,d(v,w1))...(wn,d(v,wn))
	 * written by wanjingyi
	 * */
	void write_labels(char* write_filename)
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
	}

	/**
	 * function used to serialize vertice's labels to file,format as follows:
	 * label_size v (w0,d(v,w0)) (w1,d(v,w1))...(wn,d(v,wn))
	 * written by wanjingyi
	 * */
	void write_labels(const char* write_filename,const vector<NodeID>& inv,bool isOrder=false)
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
		string write_filename2(write_filename);
		write_filename2.append("_original.list");
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
			string write_filename3(write_filename);
			write_filename3.append("_original_ordered.list");
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

/***
 * class HFLabel inherented from Label 
 * written by wanjingyi
 * **/
class HFLabel : public Label{
	public:
		vector<NodeID> HFPoint; //store the high frequency point referred orders
		vector<bool> HFinGraphIndex;//HFinGraphIndex is index of high frequency point in original graph, 0 represents low frequency point，1 represents high frequency point
		vector<NodeID> HFOripoint; //store the high frequency points'ids
		vector<NodeID> LFOripoint; //store the low frequency points'ids
		vector<NodeID> spt_v_num; //store the size of each vertice's labels
		vector<NodeID> rank;//store the verices' rank
		vector<int> queryTime;//store the query time of each point by id
		vector< pair<unsigned int,NodeID> > queryTime_pair; //store the query time of each point which is not null
		int numOfHFpoint=0; //number of High frequency points
		long long totalQueryTime=0;//num of total query time
		long long totalQueryTime_hf=0;//num of total query time
		int maxQueryTime=0;//num of max query time
		int minQueryTime=INT_MAX;//num of min query time
		double hf_rate=0;
		vector< vector<double> > hf_chache;//use to cache the h-h point distance
		vector< vector<int> > queryPairTime;//store the pair query frequency read from file
		vector<NodeID> HFPoint_inv; //order the node by descending freq
		vector<NodeID> HFPoint_rank; //fetch the rank by node index
		index_cache_t_p* index_cache_p;//cache_laebl
		index_parallel_t_p* index_parallel_p;//parallel_laebl
		typedef google::dense_hash_map<NodeID, EdgeWeight> dis_map;//store the cache distance
		int cache_size;
	
		HFLabel(){}
		~HFLabel(){
			Free();
			//delete pointer
			// delete index_cache_p;
			// index_cache_p=NULL;
			// delete index_parallel_p;
			// index_parallel_p=NULL;
		}

		double get_query_cost(char* load_filename){
			//initialize iterms
			vector<int> queryTime_tmp(numOfVertices,0);
			double performance_result=0;
			std::ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int s,t,qt,i=0;
			char line[24];
			//read each line representing HFpoint to vector 
			while (in.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>s>>t>>qt;
				queryTime_tmp[s]+=qt;
				queryTime_tmp[t]+=qt;
			}
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				NodeID isize = spt_v_num[v];
				double ratio=(double)queryTime_tmp[v]/(double)DIVISION_FACTOR;
				performance_result+=ratio*(double)isize;
			}
			//std::cout<<"get_query_cost successfully!"<<std::endl;
			return performance_result;
		}

		void append_experiment_result_pll(char* write_filename){
			long long total_sum_size=0,hf_sum_size=0;
			double total_ave_size=0,hf_ave_size=0;
			double performance_result=0;//total performance function 
			ofstream ofs(write_filename,ios::app|ios::out);//append way
			if(!ofs.is_open()) cout<<"Cannot open "<<write_filename<<endl;
			ofs.setf(ios::fixed);
			ofs.precision(4);
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				NodeID isize = spt_v_num[v];
				total_sum_size+=isize;
				if(HFinGraphIndex[v]) hf_sum_size+=isize;
				//compute the query performance function
				double ratio=(double)queryTime[v]/((double)(maxQueryTime*DIVISION_FACTOR));
				performance_result+=ratio*(double)isize;
			}
			total_ave_size= (double) total_sum_size/(double) numOfVertices;
			hf_ave_size= (double) hf_sum_size/(double) numOfHFpoint;
			double standard_performance_result=performance_result*((double)maxQueryTime);
			std::cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			std::cout<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			std::cout<<"nomalization performance_result = "<<performance_result<<endl;
			std::cout<<"standard performance_result = "<<standard_performance_result<<endl;
			ofs<<total_sum_size<<" "<<total_ave_size<<" "<<hf_sum_size<<" "<<hf_ave_size<<" "<<standard_performance_result<<" ";
			ofs.close();
		}

		void append_experiment_result(char* write_filename,double _labeling_time,double _ordering_time=0){
			long long total_sum_size=0,hf_sum_size=0;
			double total_ave_size=0,hf_ave_size=0;
			double performance_result=0;//total performance function 
			ofstream ofs(write_filename,ios::app|ios::out);//append way
			if(!ofs.is_open()) cout<<"Cannot open "<<write_filename<<endl;
			ofs.setf(ios::fixed);
			ofs.precision(4);
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				NodeID isize = index_[v].size()-1;
				total_sum_size+=isize;
				if(HFinGraphIndex[v]) hf_sum_size+=isize;
				//compute the query performance function
				double ratio=(double)queryTime[v]/((double)(maxQueryTime*DIVISION_FACTOR));
				performance_result+=ratio*(double)isize;
			}
            total_ave_size= (double) total_sum_size/(double) numOfVertices;
			hf_ave_size= (double) hf_sum_size/(double) numOfHFpoint;
			double standard_performance_result=performance_result*((double)maxQueryTime);
			// std::cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			// std::cout<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			// std::cout<<"nomalization performance_result = "<<performance_result<<endl;
			// std::cout<<"standard performance_result = "<<standard_performance_result<<endl;
			ofs<<_ordering_time* 1e6 <<" "<<_labeling_time* 1e6 <<" "<<total_sum_size<<" "<<total_ave_size<<" "<<hf_sum_size<<" "<<hf_ave_size<<" "<<standard_performance_result<<" ";
			ofs.close();
		}

		//*************functions for no optimazation pll*************
		// void load_hfpoint_and_qt(char* load_filename,int hf_rate=0){//5%%
		// 	//initialize iterms
		// 	queryTime.resize(numOfVertices,0);
		// 	HFinGraphIndex.resize(numOfVertices,false);
		// 	if(hf_rate==0) numOfHFpoint=numOfVertices;
		// 	else numOfHFpoint=static_cast<int>((double)(numOfVertices*hf_rate)/(double)HF_DIVIDION);
		// 	std::cout<<"initial numOfHfpoint = "<<numOfHFpoint;
		// 	std::ifstream in(load_filename);//input HFPoint file to ifstream
		// 	if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
		// 	int id;int t,i=0;
		// 	char line[24];
		// 	//read each line representing HFpoint to vector 
		// 	while (in.getline(line,sizeof(line))){
		// 		std::stringstream hp(line);
		// 		hp>>id>>t;		
		// 		totalQueryTime+=t;
		// 		if(i<numOfHFpoint){
		// 			totalQueryTime_hf+=t;
		// 			if(t>maxQueryTime) maxQueryTime=t;
		// 			if(t<minQueryTime) minQueryTime=t;
		// 			queryTime[id]=t;
		// 			HFPoint.push_back(id);
		// 			HFinGraphIndex[id]=true;
		// 			i++;
		// 		}
		// 	}
		// 	if(i<numOfHFpoint){
		// 		numOfHFpoint=i;
		// 	}
		// 	//compute real hf_ratio
		// 	hf_rate=((double)numOfHFpoint/(double)numOfVertices)*HF_DIVIDION;
		// 	std::cout<<" real numOfHfpoint = "<<numOfHFpoint<<" maxQueryTime = "<<maxQueryTime<<" totalQueryTime = "<<totalQueryTime<<std::endl;
		// }

		void load_hfpoint_and_qt(char* load_filename,int hf_rate=0){//5%%
			//initialize iterms
			queryTime.resize(numOfVertices,0);
			HFinGraphIndex.resize(numOfVertices,false);
			//if(hf_rate==0) numOfHFpoint=numOfVertices;
			//else numOfHFpoint=static_cast<int>((double)(numOfVertices*hf_rate)/(double)HF_DIVIDION);
			std::ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int s,t,qt,i=0;
			char line[24];
			//read each line representing HFpoint to vector 
			while (in.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>s>>t>>qt;
				queryTime[s]+=qt;
				queryTime[t]+=qt;
				HFinGraphIndex[s]=true;
				HFinGraphIndex[t]=true;
			}
			for(size_t v=0;v<numOfVertices;++v){
				int qt=queryTime[v];
				if(qt>maxQueryTime) maxQueryTime=qt;
				if(qt<minQueryTime) minQueryTime=qt;
			}
			in.close();
			std::cout<<"maxQueryTime = "<<maxQueryTime<<std::endl;
			std::cout<<"minQueryTime = "<<minQueryTime<<std::endl;
			std::cout<<"Load_hfpoint_and_qt successfully!"<<std::endl;
			return;
		}

        void save_anaylysis_size(const char* write_filename){
			//write size analysis to file
			long long total_sum_size=0,hf_sum_size=0;
			double total_ave_size=0,hf_ave_size=0;
			double performance_result=0;//total performance function 
            string asize_filename(write_filename);
			asize_filename.append("_analysis.size");
			ofstream ofs_size(asize_filename.c_str());
			if(!ofs_size.is_open()) {cerr<<"Cannot open "<<asize_filename<<endl;}
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				NodeID isize = index_[v].size()-1;
				total_sum_size+=isize;
				if(HFinGraphIndex[v]) hf_sum_size+=isize;
				//compute the query performance function
				double ratio=(double)(queryTime[v])/((double)(maxQueryTime*DIVISION_FACTOR));
				performance_result+=ratio*(double)isize;
				//cout<<v<<": qtime="<<queryTime[v]<<", size="<<isize<<endl;
			}
            total_ave_size= (double) total_sum_size/(double) numOfVertices;
			hf_ave_size= (double) hf_sum_size/(double) numOfHFpoint;
			std::cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			std::cout<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			ofs_size<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			ofs_size<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			//write performance function to file
			double standard_performance_result=performance_result*((double)maxQueryTime);
			std::cout<<"nomalization performance_result = "<<performance_result<<endl;
			std::cout<<"standard performance_result = "<<standard_performance_result<<endl;
			ofs_size<<"nomalization performance_result = "<<performance_result<<endl;
			ofs_size<<"standard performance_result = "<<standard_performance_result<<endl;	
			ofs_size.close();	
        }
		//*************functions for no optimazation pll*************

		/*
		 *@description: load_query_pair_frequency and build cache
		 *@author: wanjingyi
		 *@parameter:k-cache size
		 *@date: 2021-01-15
		*/
		// int refresh_cache(char* queryDataFileName,char* pointFreqFileName,int hfRate=5){
		// 	numOfHFpoint = 0;//first line is the number of HFpoints
		// 	numOfHFpoint = static_cast<int> ( (double)numOfVertices*hfRate/(double)1000);
		// 	if(numOfHFpoint<=0) cout<<"error:numOfHFpoint<=0"<<endl;
		// 	cout<<"initial numOfHFpoint  = "<<numOfHFpoint <<endl;
		// 	std::cout<<"refresh cache numOfVertices="<<numOfVertices<<endl;//to be deleted
		// 	//clear cache dense_hash_map
		// 	for(int v=0;v<numOfVertices;++v){
		// 		//dis_map& v_dis_cache=index_cache_p[v].dis_cache_m;
		// 		if(!index_cache_p[v].dis_cache_m.empty()){
		// 			index_cache_p[v].dis_cache_m.clear();
		// 		}
		// 		index_cache_p[v].dis_cache_m.set_empty_key(numOfVertices+1);//size=0
		// 		//index_cache_p[v].dis_cache_m.set_deleted_key(numOfVertices+2);
		// 		HFinGraphIndex[v]=false;
		// 	}
		// 	//to be deleted
		// 	std::cout<<"dense_hash_map initialize finished!"<<endl;
		// 	//std::cout<<"0 dense_hash_map_size()="<<index_cache_p[0].dis_cache_m.size()<<endl;
		// 	//load hfpoint that needed to be cached
		// 	load_cache_hfPoint(pointFreqFileName);
		// 	std::cout<<"real numOfHFpoint = "<<numOfHFpoint<<endl;
		// 	//cache construction
		// 	int cnt=0,both_cnt; 
		// 	NodeID s,t;
		// 	bool s_hf,t_hf;
		// 	std::cout<<"start to read "<<queryDataFileName<<endl;
		// 	vector<pair<int,int> > cache_pair;
		// 	int num_cache_pair=load_query_pair_time(queryDataFileName,cache_pair);
		// 	std::cout<<"num_cache_pair="<<num_cache_pair<<endl;
		// 	for(int i=0;i<num_cache_pair;++i)
		// 	{
		// 		s_hf=false;t_hf=false;
		// 		s=cache_pair[i].first;
		// 		t=cache_pair[i].second;
		// 		dis_map& s_dis_cache=index_cache_p[s].dis_cache_m;
		// 		dis_map& t_dis_cache=index_cache_p[t].dis_cache_m;
		// 		s_dis_cache.insert(make_pair(0,0));
		// 		std::cout<<s<<"-"<<t<<endl;
		// 		//undirected graph
		// 		if(HFinGraphIndex[s]==true&&s_dis_cache.size()<cache_size){
		// 			s_dis_cache.insert(make_pair(t,INF_WEIGHT));
		// 			//s_dis_cache[t]=INF_WEIGHT;
		// 			cnt++;
		// 			s_hf=true;
		// 			std::cout<<s<<" isCached"<<endl;//to be deleted
		// 		}
		// 		if(HFinGraphIndex[t]&&t_dis_cache.size()<cache_size){
		// 			t_dis_cache.insert(make_pair(s,INF_WEIGHT));
		// 			cnt++;
		// 			t_hf=true;
		// 			std::cout<<t<<" isCached"<<endl;//to be deleted
		// 		}
		// 		if(s_hf&t_hf) both_cnt++;
		// 	}
		// 	std::cout<<"cache s-t pair size = "<<cnt<<endl;
		// 	std::cout<<"cache both hf size = "<<both_cnt<<endl;
		// }

		/*
		 *@description: load hfpoint for cache construction
		 *@author: wanjingyi
		 *@date: 2021-01-17
		*/
		// void load_cache_hfPoint(char* pointFreqFileName){
		// 	ifstream in(pointFreqFileName);
		// 	if(!in.is_open()) {cerr<<"Cannot open "<<pointFreqFileName<<endl;}
		// 	NodeID id;int t,i=0;
		// 	char line[24];
		// 	//read each line representing HFpoint to vector 
		// 	while (in.getline(line,sizeof(line))){
		// 		stringstream hp(line);
		// 		hp>>id>>t;		
		// 		if(i>=numOfHFpoint)	break;
		// 		HFinGraphIndex[id]=true;
		// 		i++;
		// 	}
		// 	if(i<numOfHFpoint){
		// 		numOfHFpoint=i;
		// 	}
		// 	in.close();
		// }

		/*
		 *@description: load cache labels
		 *@author: wanjingyi
		 *@date: 2021-01-17
		*/
		// void load_cache_labels(const char* load_filename,int k){
		// 	ifstream ifs(load_filename);
		// 	NodeID isize = 0;
		// 	ifs.read((char*)&isize, sizeof(isize));
		// 	numOfVertices = isize;
		// 	cout<<"numOfVertices = "<<numOfVertices<<endl;
		// 	if(index_cache_p)
		// 	{
		// 		for (NodeID v = 0; v < numOfVertices; ++v) {
		// 				free(index_cache_p[v].spt_v);
		// 				free(index_cache_p[v].spt_d);
		// 			} 
		// 		free(index_cache_p);
		// 	}else{
		// 		HFinGraphIndex.resize(numOfVertices,false);
		// 		std::cout<<"initialize HFinGraphIndex to false!"<<endl;
		// 	}
		// 	index_cache_p = NULL;
		// 	index_cache_p = (index_cache_t_p*)memalign(64, numOfVertices * sizeof(index_cache_t_p));
		// 	for (NodeID v = 0; v < numOfVertices; ++v) {
		// 		index_cache_t_p &idx = index_cache_p[v];
		// 		ifs.read((char*)&isize, sizeof(isize));
		// 		idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
		// 		idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
		// 		for (NodeID i = 0; i < isize; ++i) {
		// 			NodeID hub;
		// 			EdgeWeight hub_weight;
		// 			ifs.read((char*)&hub, sizeof(hub));
		// 			ifs.read((char*)&hub_weight, sizeof(hub_weight));
		// 			idx.spt_v[i] = hub;
		// 			idx.spt_d[i] = hub_weight;
		// 		}
		// 	}
		// 	ifs.close();
		// 	//initialize iterms
		// }

		/*
		 *@description: load hf and parallel labels for parallel
		 *@author: wanjingyi
		 *@date: 2021-01-15
		*/
		// void load_cache_parallel_labels(const char* load_filename,int k,int hfRate){
		// 	//clear label
		// 	if(index_cache_p)
		// 	{
		// 		for(NodeID v=0;v<numOfVertices;++v){
		// 			free(index_cache_p[v].spt_d_h);
		// 			free(index_cache_p[v].spt_d_l);
		// 			free(index_cache_p[v].spt_v_h);
		// 			free(index_cache_p[v].spt_v_l);
		// 		}
		// 		free(index_cache_p);
		// 		vector<bool> ().swap(HFinGraphIndex);
		// 		HFinGraphIndex.clear();
		// 	}
		// 	index_cache_p=NULL;
		// 	cache_size=k;
		// 	HFinGraphIndex.resize(numOfVertices);
		// 	//read labels from binary file
		// 	ifstream ifs(load_filename);
		// 	NodeID isize = 0;
		// 	ifs.read((char*)&isize, sizeof(isize));
		// 	numOfVertices = isize;
		// 	cout<<"numOfVertices = "<<numOfVertices<<endl;

		// 	index_cache_p=(index_cache_t_p*)memalign(64,numOfVertices*sizeof(index_cache_t_p));
		// 	for (NodeID v = 0; v < numOfVertices; ++v) {
		// 		index_cache_t_p &idx = index_cache_p[v];
		// 		ifs.read((char*)&isize, sizeof(isize));
		// 		NodeID isize_h=(NodeID)isize/2;
		// 		NodeID isize_l=isize-1-isize_h;
		// 		idx.spt_v_h = (NodeID*)memalign(64, (isize_h+1) * sizeof(NodeID));
		// 		idx.spt_v_l = (NodeID*)memalign(64, (isize_l+1) * sizeof(NodeID));
		// 		idx.spt_d_h = (EdgeWeight*)memalign(64, (isize_h+1) * sizeof(EdgeWeight));
		// 		idx.spt_d_l = (EdgeWeight*)memalign(64, (isize_l+1) * sizeof(EdgeWeight));
		// 		NodeID hub;
		// 		EdgeWeight hub_weight;
		// 		NodeID i;
		// 		for(i=0;i<isize_h;++i){
		// 			ifs.read((char*)&hub, sizeof(hub));
		// 			ifs.read((char*)&hub_weight, sizeof(hub_weight));
		// 			idx.spt_v_h[i] = hub;
		// 			idx.spt_d_h[i] = hub_weight;
		// 		}
		// 		idx.spt_v_h[i] = numOfVertices;
		// 		idx.spt_d_h[i] = INF_WEIGHT;
		// 		for(;i<isize;++i){
		// 			ifs.read((char*)&hub, sizeof(hub));
		// 			ifs.read((char*)&hub_weight, sizeof(hub_weight));
		// 			idx.spt_v_l[i] = hub;
		// 			idx.spt_d_l[i] = hub_weight;			
		// 		}
		// 	}
		// 	ifs.close();
		// }

		/*
		 *@description: load parallel labels for parallel
		 *@author: wanjingyi
		 *@date: 2021-01-15
		*/
		// void load_parallel_labels(const char* load_filename,int threads_num){
		// 	//clear label
		// 	if(index_parallel_p)
		// 	{
		// 		for(NodeID v=0;v<numOfVertices;++v){
		// 			free(index_parallel_p[v].spt_v);
		// 			free(index_parallel_p[v].spt_d);
		// 			free(index_parallel_p[v].index);
		// 		}
		// 		free(index_parallel_p);
		// 	}
		// 	index_parallel_p=NULL;

		// 	//read labels from binary file
		// 	ifstream ifs(load_filename);
		// 	NodeID isize = 0;
		// 	ifs.read((char*)&isize, sizeof(isize));
		// 	numOfVertices = isize;
		// 	cout<<"numOfVertices = "<<numOfVertices<<endl;

		// 	index_parallel_p=(index_parallel_t_p*)memalign(64,numOfVertices*sizeof(index_parallel_t_p));
		// 	//index_parallel_p=(index_parallel_t_p*)malloc(numOfVertices*sizeof(index_parallel_t_p));
		// 	std::cout<<"memalign index_parallel_p finished!"<<endl;
		// 	for (NodeID v = 0; v < numOfVertices; ++v) {
		// 		index_parallel_t_p &idx = index_parallel_p[v];
		// 		ifs.read((char*)&isize, sizeof(isize));
		// 		//store th label size
		// 		idx.lsiz=isize-1;
		// 		//compute section index
		// 		idx.index=new int[threads_num+1];
		// 		int j,k=0,p_siz=(isize-1)/threads_num;
		// 		for(j=0;j<threads_num;++j){
		// 			idx.index[j]=k;
		// 			k+=p_siz;
		// 		}
		// 		idx.index[j]=isize-1;
		// 		//malloc hub and edgeWeight pointer
		// 		idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
		// 		idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
		// 		//read and store in list
		// 		for (NodeID i = 0; i < isize; ++i) {
		// 			NodeID hub;
		// 			EdgeWeight hub_weight;
		// 			ifs.read((char*)&hub, sizeof(hub));
		// 			ifs.read((char*)&hub_weight, sizeof(hub_weight));
		// 			idx.spt_v[i] = hub;
		// 			idx.spt_d[i] = hub_weight;
		// 		}
		// 	}
		// 	std::cout<<"load_parallel_labels finished! "<<endl;
		// 	ifs.close();
		// }

		/*
		 *@description: warmUp distance for cache labels
		 *@author: wanjingyi
		 *@date: 2021-01-17
		*/
		// EdgeWeight cache_warmUp(NodeID s,NodeID t){
		// 	EdgeWeight distance = INF_WEIGHT;
		// 	index_cache_t_p &idx_s = index_cache_p[s];
		// 	index_cache_t_p &idx_t = index_cache_p[t];
		// 	_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
			
		// 	for (int i = 0, j = 0; ; ) {
		// 		NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
		// 		if (v1 == numOfVertices) break;  // Sentinel
		// 		if (v1 == v2) {
		// 			EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
		// 			if (td < distance) distance = td;
		// 			++i;
		// 			++j;
		// 		} 
		// 		else {
		// 			i += v1 < v2 ? 1 : 0;
		// 			j += v1 > v2 ? 1 : 0;
		// 		}
		// 	}
			
		// 	return distance;
		// }

		/*
		 *@description: query distance for cache labels
		 *@author: wanjingyi
		 *@date: 2021-01-17
		*/
		// EdgeWeight cache_query(NodeID s,NodeID t){
		// 	EdgeWeight distance = INF_WEIGHT;
		// 	index_cache_t_p &idx_s = index_cache_p[s];
		// 	index_cache_t_p &idx_t = index_cache_p[t];
		// 	_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		// 	//find s-cache
		// 	if(idx_s.dis_cache_m.find(t)!=idx_s.dis_cache_m.end()){
		// 		distance=idx_s.dis_cache_m[t];
		// 		if(distance!=INF_WEIGHT) return distance;
		// 		else{
		// 			//*******merge browse*********
		// 			for (int i = 0, j = 0; ; ) {
		// 				NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
		// 				if (v1 == numOfVertices) break;  // Sentinel
		// 				if (v1 == v2) {
		// 					EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
		// 					if (td < distance) distance = td;
		// 					++i;
		// 					++j;
		// 				} 
		// 				else {
		// 					i += v1 < v2 ? 1 : 0;
		// 					j += v1 > v2 ? 1 : 0;
		// 				}
		// 			}
		// 			//*******merge browse*********
		// 			idx_s.dis_cache_m[t]=distance;
		// 			return distance;
		// 		}
		// 	}

		// 	//find t-cache
		// 	if(idx_t.dis_cache_m.find(s)!=idx_t.dis_cache_m.end()){
		// 		distance=idx_s.dis_cache_m[s];
		// 		if(distance!=INF_WEIGHT) return distance;
		// 		else{
		// 			//*******merge browse*********
		// 			for (int i = 0, j = 0; ; ) {
		// 				NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
		// 				if (v1 == numOfVertices) break;  // Sentinel
		// 				if (v1 == v2) {
		// 					EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
		// 					if (td < distance) distance = td;
		// 					++i;
		// 					++j;
		// 				} 
		// 				else {
		// 					i += v1 < v2 ? 1 : 0;
		// 					j += v1 > v2 ? 1 : 0;
		// 				}
		// 			}
		// 			//*******merge browse*********
		// 			idx_t.dis_cache_m[s]=distance;
		// 			return distance;
		// 		}
		// 	}

		// 	//*******cache miss*********
		// 	for (int i = 0, j = 0; ; ) {
		// 		NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
		// 		if (v1 == numOfVertices) break;  // Sentinel
		// 		if (v1 == v2) {
		// 			EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
		// 			if (td < distance) distance = td;
		// 			++i;
		// 			++j;
		// 		} 
		// 		else {
		// 			i += v1 < v2 ? 1 : 0;
		// 			j += v1 > v2 ? 1 : 0;
		// 		}
		// 	}
		// 	//*******merge browse*********
		// 	return distance;
		// }

		// EdgeWeight thread_query(const index_parallel_t_p& idx_s,const index_parallel_t_p& idx_t,int l,int h){
		// 	EdgeWeight distance=INF_WEIGHT;
		// 	for (int i = l, j = 0;i<h; ) {
		// 		NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
		// 		if (v1 == numOfVertices) break;  // Sentinel

		// 		if (v1 == v2) {
		// 			EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
		// 			//if(td==0) std::cout<<" ("<<v1<<","<<idx_s.spt_d[i]<<") "<<"("<<v2<<","<<idx_t.spt_d[j]<<")"<<endl;
		// 			if (td < distance) {
		// 				distance = td;
		// 			}
		// 			++i;
		// 			++j;
		// 		}
		// 		else {
		// 			i += v1 < v2 ? 1 : 0;
		// 			j += v1 > v2 ? 1 : 0;
		// 		}
		// 	}
		// 	return distance;
		// }

		/*
		 *@description: query distances parallel by sections
		 *@author: wanjingyi
		 *@date: 2021-01-19
		*/
		// EdgeWeight parallel_query(NodeID s, NodeID t){
		// 	EdgeWeight l_distance=INF_WEIGHT,h_distance=INF_WEIGHT;
		// 	index_parallel_t_p& idx_s=index_parallel_p[s];
		// 	index_parallel_t_p& idx_t=index_parallel_p[t];
		// 	_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
		// 	omp_set_num_threads(2);
			
		// 	// if(idx_s.lsiz<idx_t.lsiz){
		// 	// 	index_parallel_t_p& tmp=index_parallel_p[s];
		// 	// 	idx_s=idx_t;
		// 	// 	idx_t=tmp;
		// 	// }			
		// 	#pragma omp parallel sections
		// 	{
		// 		#pragma omp section
		// 		{
		// 			l_distance=thread_query(idx_s,idx_t,idx_s.index[0],idx_s.index[1]);
		// 		}

		// 		#pragma omp section
		// 		{
		// 		h_distance=thread_query(idx_s,idx_t,idx_s.index[1],idx_s.index[2]);
		// 		}
		// 	}
		// 	return l_distance<h_distance?l_distance:h_distance;
		// }

		/*
		*@description: query the pivot in parallel model
		*@author: wanjingyi
		*@date: 2021-01-15
		*/
		// EdgeWeight parallel_query(NodeID s, NodeID t,int threads_num){
		// 	EdgeWeight shared_distance=INF_WEIGHT;
		// 	index_parallel_t_p& idx_s=index_parallel_p[s];
		// 	index_parallel_t_p& idx_t=index_parallel_p[t];
		// 	int k=0;
		// 	// if(idx_s.lsiz<idx_t.lsiz){
		// 	// 	idx_s=index_parallel_p[t];
		// 	// 	idx_t=index_parallel_p[s];
		// 	// }
		// 	_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		// 	_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		// 	//parallel
		// 	omp_set_num_threads(threads_num);
		// 	#pragma omp parallel
		// 	{
		// 		EdgeWeight distance=INF_WEIGHT;
		// 		#pragma omp for nowait
		// 		// #pragma omp for
		// 		for(k=0;k<threads_num;++k){
		// 			for(int i=idx_s.index[k],j=0;i<idx_s.index[k+1];){
		// 				NodeID v1=idx_s.spt_v[i],v2=idx_t.spt_v[j];
		// 				if (v2 == numOfVertices) break;  // Sentinel
		// 				if(v1==v2){
		// 					EdgeWeight td=idx_s.spt_d[i]+idx_t.spt_d[j];
		// 					if(td<distance) distance=td;
		// 					++i;
		// 					++j;
		// 				}else{
		// 					i += v1 < v2 ? 1 : 0;
		// 					j += v1 > v2 ? 1 : 0;
		// 				}
		// 			}
		// 			#pragma omp critical 
		// 			{
		// 				shared_distance = std::min(shared_distance, distance);
		// 			}
		// 		}
		// 	}
		// 	//parallel sections finished!
		// 	return shared_distance;
		// }

		/*
		*@description: query the pivot in parallel model
		*@author: wanjingyi
		*@date: 2021-01-15
		*/
		// EdgeWeight parallel_query(NodeID s,NodeID t){
		// 	EdgeWeight distance_h = INF_WEIGHT,distance_l=INF_WEIGHT;
		// 	index_parallel_t_p& idx_s=index_parallel_p[s];
		// 	index_parallel_t_p& idx_t=index_parallel_p[t];
		// 	#pragma omp parallel sections
		// 	{
		// 		#pragma omp section
		// 		{
		// 			_mm_prefetch(&idx_s.spt_v_h[0], _MM_HINT_T0);
		// 			_mm_prefetch(&idx_t.spt_v_h[0], _MM_HINT_T0);	
		// 			_mm_prefetch(&idx_s.spt_d_h[0], _MM_HINT_T0);
		// 			_mm_prefetch(&idx_t.spt_d_h[0], _MM_HINT_T0);
		// 			for(int i=0,j=0;;){
		// 				NodeID v1=idx_s.spt_v_h[i],v2=idx_t.spt_v_h[j];
		//				if (v1 == numOfVertices) break;  // Sentinel
		// 				if(v1==v2){
		// 					EdgeWeight td=idx_s.spt_d_h[i]+idx_t.spt_d_h[j];
		// 					if(td<distance_h) distance_h=td;
		// 					++i;
		// 					++j;
		// 				}else{
		// 					i += v1 < v2 ? 1 : 0;
		// 					j += v1 > v2 ? 1 : 0;
		// 				}
		// 			}
		// 		}

		// 		#pragma omp section
		// 		{
		// 			_mm_prefetch(&idx_s.spt_v_l[0], _MM_HINT_T0);
		// 			_mm_prefetch(&idx_t.spt_v_l[0], _MM_HINT_T0);
		// 			_mm_prefetch(&idx_s.spt_d_l[0], _MM_HINT_T0);
		// 			_mm_prefetch(&idx_t.spt_d_l[0], _MM_HINT_T0);
		//				if (v1 == numOfVertices) break;  // Sentinel
		// 			for(int i=0,j=0;;){
		// 				NodeID v1=idx_s.spt_v_l[i],v2=idx_t.spt_v_l[j];
		// 				if(v1==v2){
		// 					EdgeWeight td=idx_s.spt_d_l[i]+idx_t.spt_d_l[j];
		// 					if(td<distance_l) distance_l=td;
		// 					++i;
		// 					++j;
		// 				}else{
		// 					i += v1 < v2 ? 1 : 0;
		// 					j += v1 > v2 ? 1 : 0;
		// 				}
		// 			}
		// 		}

		// 	}
		// 	return distance_h<distance_l?distance_h:distance_l;
		// }

		/**
		 * function update load_query_time
		 * return total query_time
		 * written by wanjingyi
		 * */
		pair<int,long long> load_deordered_query_time(char* load_filename){
			ifstream in(load_filename);//input query file to ifstream
			if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
			NodeID q_id;
			unsigned int query_time=0; //each line representing the label size of each vertice
			char line[24];
			long long  query_cnt=0;
			int max_query_time=0;
			 while (in.getline(line,sizeof(line)))
			 {
				 stringstream ls(line);
				 ls>>q_id>>query_time;
				 queryTime[q_id]=query_time;
				 query_cnt+=query_time;
				queryTime_pair.push_back(make_pair(query_time,q_id));
			 }
			 in.close();
			 //sort
			//sort(queryTime_pair.rbegin(),queryTime_pair.rend());

			for(NodeID i=0;i<numOfVertices;++i){
				NodeID v=queryTime_pair[i].second;
				HFPoint_inv[i]=v;
				HFPoint_rank[v]=i;
				if(i<numOfHFpoint) HFinGraphIndex[v]=true;
			}
			std::cout<<"total query time = "<<query_cnt<<endl;
			double ave_total_query_time=(double)query_cnt/(double)numOfVertices;
			max_query_time=queryTime[HFPoint_inv[0]];
			if(max_query_time!=queryTime_pair[0].first) std::cout<<"max_query_time compute error!"<<endl;
			std::cout<<"average query time = "<<ave_total_query_time<<endl;
			std::cout<<"max_query_time = "<<max_query_time<<endl;
			return make_pair(max_query_time,query_cnt);
		}


		/**
		 * function used to read  query time from file
		 * return total query_time
		 * written by wanjingyi
		 * */
		pair<int,long long> load_query_time(char* load_filename){
			ifstream in(load_filename);//input query file to ifstream
			if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
			NodeID q_id;
			unsigned int query_time=0; //each line representing the label size of each vertice
			char line[24];
			long long  query_cnt=0;
			int max_query_time=0;
			 while (in.getline(line,sizeof(line)))
			 {
				 stringstream ls(line);
				 ls>>q_id>>query_time;
				 queryTime[q_id]=query_time;
				 query_cnt+=query_time;
				queryTime_pair.push_back(make_pair(query_time,q_id));
			 }
			 in.close();
			 //sort
			//sort(queryTime_pair.rbegin(),queryTime_pair.rend());
			long long cover_cnt=0;
			unsigned int hfpoint_cnt=0;
			for(NodeID i=0;i<queryTime_pair.size();++i){
				NodeID v=queryTime_pair[i].second;
				HFPoint_inv[i]=v;
				HFPoint_rank[v]=i;
				if(i<numOfHFpoint) 
				{
					HFinGraphIndex[v]=true;
					++hfpoint_cnt;
					cover_cnt+=queryTime_pair[i].first;
				}
			}
			NodeID rank=queryTime_pair.size();
			for(NodeID v=0;v<numOfVertices;++v){
				if(queryTime[v]==0){
					HFPoint_inv[rank]=v;
					HFPoint_rank[v]=rank;
					++rank;
				}		
			}
			std::cout<<"total query time = "<<query_cnt<<endl;
			double ave_total_query_time=(double)query_cnt/(double)numOfVertices;
			max_query_time=queryTime[HFPoint_inv[0]];
			if(max_query_time!=queryTime_pair[0].first) std::cout<<"max_query_time compute error!"<<endl;
			std::cout<<"average query time = "<<ave_total_query_time<<endl;
			std::cout<<"max_query_time = "<<max_query_time<<endl;
			double hfpoint_cover_ratio=(double)hfpoint_cnt/(double)numOfVertices;
			double freq_cover_ratio=(double)cover_cnt/(double)query_cnt;
			std::cout<<"real_hfpoint_cover_ratio = "<<hfpoint_cover_ratio<<endl;
			std::cout<<"real_freq_cover_ratio = "<<freq_cover_ratio<<endl;
			return make_pair(max_query_time,query_cnt);
		}

		/**
		 * function used to read  spt_v_num from file
		 * written by wanjingyi
		 * */
		void load_label_size(char* load_filename){
			//spt_v_num.resize(numOfVertices,0);
			ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
			NodeID label_size=0; //each line representing the label size of each vertice
			int i=0;
			 //read each line representing HFpoint to vector 
			 while (in>>label_size)
			{
				spt_v_num.push_back(label_size);
			}
			in.close();
			if(numOfVertices == 0){
				numOfVertices=spt_v_num.size();
				//cout<<"load_label_size: numOfVertices="<<spt_v_num.size()<<endl;
			}
			return;
		}

		/**
		 * function used to read orders from file
		 * written by wanjingyi
		 * */
		void load_order(char* load_filename)
		{
			rank.resize(numOfVertices,0);
			ifstream in(load_filename);//input orderfile to ifstream
			if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
			NodeID id,r=0;
			char line[24];
			 //read each line representing HFpoint to vector 
			 while (in.getline(line,sizeof(line)))
			 {
				 stringstream id_ss(line);
				 id_ss>>id;
				 rank[id]=r++;
			 }
			in.close();
			if(r!=numOfVertices) cout<<"r!=numOfVertices"<<endl;
			
		}

		/**
		 * function used to read HFpoint from file
		 * written by wanjingyi
		 * */
		 void load_HFpoint(char* load_filename,int hfRate=50){
			 numOfHFpoint = 0;//first line is the number of HFpoints
			 numOfHFpoint = static_cast<int> ( (double)numOfVertices*hfRate/(double)1000);
			 if(numOfHFpoint<=0) cout<<"error:numOfHFpoint<=0"<<endl;
			 cout<<"numOfHFpoint  = "<<numOfHFpoint <<endl;
			 ifstream in(load_filename);//input HFPoint file to ifstream
			 if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
			 NodeID t,id;int i=0;
			 char line[24];
			 //read each line representing HFpoint to vector 
			 while (in.getline(line,sizeof(line)))
			 {
				 stringstream hp(line);
				 hp>>t;
				 if(i<numOfHFpoint) HFOripoint.push_back(t);
				 else LFOripoint.push_back(t);
				 ++i;
			 }
			 in.close();
			 if(i!=numOfVertices) cout<<"i!=numOfVertices"<<endl;
			 cout<<"HFOripoint.size()="<<HFOripoint.size()<<endl;
			 cout<<"LFOripoint.size()="<<LFOripoint.size()<<endl;
			cout<<"numOfHFpoint="<<numOfHFpoint<<endl;

			//update HFPoint and HFinGraphIndex by order instead of node id
			HFinGraphIndex.resize(numOfVertices,0);//initial 0
			HFPoint.resize(numOfHFpoint,0);
			for(i=0;i<numOfHFpoint;i++)
			{
				id=HFOripoint[i];
				HFPoint[i]=rank[id]; //fetch the order rank by NodeId
				HFinGraphIndex[id]=true;
			}
		 }

		 /**
		  * function used to find  the minist-index label of current point
		  * modified by wanjingyi
		  * */
		 int findminIndex(vector<pair<NodeID,NodeID>>& index)
		 {
			NodeID min{ INT_MAX };//min代表当前找到的label中结点的最小编号
			int result{ -1 };//表示当前找到最小值的点对应的数组位置
			NodeID currentPoint;//表示当前比较的label结点
			//cout<<"findminIndex:numOfVertices="<<numOfVertices<<endl;
			for (int i = 0; i < numOfHFpoint; i++) {
				//if(i%100==0) cout<<"i="<<i<<" spt_v_num[i]="<<spt_v_num[i]<<";";
				if (index[i].second < spt_v_num[i]) {
					//cout<<"index=("<<index[i].first<<","<<index[i].second<<")"<<";";
					currentPoint = index_p[index[i].first].spt_v[index[i].second];
					//cout<<" currentPoint="<<currentPoint<<endl;
					if (currentPoint < min) {
						min = currentPoint;
						result = i;
					}
				}
			}
			for (int i = 0; i < HFPoint.size(); i++) {
				if (index_p[index[i].first].spt_v[index[i].second] == min) {
					index[i].second++;
				}	
			}
		if(result!=-1)
		index[result].second--;
			if (INT_MAX == min)
				result == -1;
			return result;
		 }

		 /**
		  * function used to serialize vertice's labels to file,format as follows:
		  * label_size v (w0,d(v,w0)) (w1,d(v,w1))...(wn,d(v,wn))
		  * written by wanjingyi
		  * */
		 void write_update_labels(const char* write_filename)
		 {
			 cout<<"write_update_labels begins!"<<endl;
			 string write_filename_prefix(write_filename);
			 string write_filename1=write_filename_prefix.append("_update.list");
			 ofstream ofs(write_filename1.c_str());
			 if(!ofs.is_open()) {cerr<<"Cannot open "<<write_filename1<<endl;}
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				NodeID isize = spt_v_num[v]+1;
				ofs <<isize<<" "<<v;
				for (NodeID i = 0; i < isize; ++i) {
					ofs<<" "<<'('<<index_p[v].spt_v[i]<<","<<index_p[v].spt_d[i]<<")";
				}
				ofs<<endl;
			}
			ofs.close();
			cout<<"write_update_labels finished!"<<endl;
		 }

		/*
		 *@description: get the query cost fo one time slicing need to load label size
		 *@return:double query cost
		 *@author: wanjingyi
		 *@date: 2021-01-06
		*/
		double get_time_slicing_query_cost(string freqid_filename)
		{
		}
		/*
		 *@description: get the query cost fo one time slicing
		 *@return:double query cost
		 *@author: wanjingyi
		 *@date: 2020-12-29
		*/
		double get_time_slicing_query_cost(string freqid_filename,bool hasLabel)
		{
			vector<int> idQueryTime(numOfVertices,0);
			ifstream in(freqid_filename);//input query file to ifstream
			if(!in.is_open()) {cerr<<"Cannot open "<<freqid_filename<<endl;}
			NodeID q_id;
			unsigned int query_time=0; //each line representing the label size of each vertice
			char line[24];
			long long  total_query_time=0;
			 while (in.getline(line,sizeof(line)))
			 {
				 stringstream ls(line);
				 ls>>q_id>>query_time;
				 idQueryTime[q_id]=query_time;
				 total_query_time+=query_time;
			 }
			 in.close();
			double performance_result=0;//total performance function 
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				NodeID isize = index_[v].size()-1 ;
				//compute the query performance function
				double ratio=(double)idQueryTime[v]/(double)DIVISION_FACTOR;
				performance_result+=ratio*(double)isize;
			}
			//write performance function to file
			//std::cout<<"standard performance_result = "<<performance_result<<endl;
			return performance_result;
		}

		 /**
		  * function used to output the total and ave size comparation
		  * written by wanjingyi
		  * */
		 void save_anaylysis_size(char* freqid_filename,char* write_filename,int hfRate=5)
		 {
			 //load hfPoint
			if(hfRate==0) numOfHFpoint=numOfVertices;
			else numOfHFpoint = static_cast<int> ( (double)numOfVertices*hfRate/(double)HF_DIVIDION);
			if(numOfHFpoint<=0) cout<<"error:numOfHFpoint<=0"<<endl;
			cout<<"numOfHFpoint  = "<<numOfHFpoint <<endl;
			//initialize vector
			HFOripoint.resize(numOfHFpoint);
			HFPoint_rank.resize(numOfVertices);
			HFPoint_inv.resize(numOfVertices);
			queryTime.resize(numOfVertices,0);
			queryTime_pair.reserve(numOfVertices);
			HFinGraphIndex.resize(numOfVertices,false);
			pair<int,long long> query_pair=load_query_time(freqid_filename);
			int max_query_time=query_pair.first;
			long long total_query_time=query_pair.second;
			//write size analysis to file
			long long total_sum_size=0,hf_sum_size=0;
			double total_ave_size=0,hf_ave_size=0;
			double performance_result=0;//total performance function 
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				NodeID isize = index_[v].size()-1;
				total_sum_size+=isize;
				if(HFinGraphIndex[v]) hf_sum_size+=isize;
				//compute the query performance function
				double ratio=(double)queryTime[v]/((double)(max_query_time*DIVISION_FACTOR));
				performance_result+=ratio*(double)isize;
			}
			total_ave_size= (double) total_sum_size/(double) numOfVertices;
			hf_ave_size= (double) hf_sum_size/(double) numOfHFpoint;
			cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			cout<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			string write_filename_prefix(write_filename);
			string asize_filename=write_filename_prefix.append("_analysis.size");
			ofstream ofs(asize_filename.c_str());
			if(!ofs.is_open()) {cerr<<"Cannot open "<<asize_filename<<endl;}
			ofs<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			ofs<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			ofs.close();
			//write performance function to file
			std::cout<<"nomalization performance_result = "<<performance_result<<endl;
			performance_result=performance_result*((double)max_query_time);
			std::cout<<"standard performance_result = "<<performance_result<<endl;
		 }

		 /**
		  * function used to serialize vertice's labels to file,format as follows:
		  * label_size v (w0,d(v,w0)) (w1,d(v,w1))...(wn,d(v,wn))
		  * written by wanjingyi
		  * */
		 void write_labels(const char* write_filename,const vector<NodeID>& inv)
		 {
			 string write_filename_prefix(write_filename);
			 string write_filename1=write_filename_prefix.append("_label.list");
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
			string write_filename2=write_filename_prefix.append("o");
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
		 }

		/**
		 * function used to quickly updtae all labels
		 * written by wanjingyi
		 * */
		void update_all_labels()
		{

			//数据量小时，注释掉所有omp部分，数据量在50万以上时，恢复所有omp部分，打开VS的openMP并行支持

			//以下部分代码生成pivot的label，即每个pivot对应的结点编号，作为预处理，不计入Update时间，其结果pivot_p作为quickUpdateLabel函数的参数
				vector<vector<int>> pivotLabel(numOfVertices);//存放每个pivot对应的所有结点编号
				for (int i = 0; i < numOfVertices; i++) {
					pivotLabel[i].reserve(AVG_LABEL_SIZE);
				}
				for (int i = 0; i < numOfVertices; i++) {//对于每一个结点
					omp_set_num_threads(16);
			#pragma omp parallel for
					for (int j = 0; j < spt_v_num[i]; j++) {//对于结点中的每一个pivot
						pivotLabel[index_p[i].spt_v[j]].push_back(i);//将当前结点号存入该pivot，遍历过程确保了存储结果有序
					}
				}

				/****************quickUpdate***********/
				//time_t start_quick = clock();
				cout<<"update all point begins!"<<endl;
				double update_time=GetCurrentTimeSec();
				//unsigned __int64 numOfHFpoint{=HFPoint.size() };
				vector<int> HFhopPoint;//存放所有HFPoint中无重复的hop点
				vector<vector<int>> tempLabel(numOfVertices);//用来临时存放更新的所有结点的label
				//为tempLabel中每一个结点预留存储空间,设定预留大小为结点总数的1%
				omp_set_num_threads(16);
			#pragma omp parallel for
				for (int i = 0; i < numOfVertices; i++) {
					tempLabel[i].reserve(AVG_LABEL_SIZE+1);
				}																																			
				vector<bool> HFPivotFlag(numOfVertices, 0);//临时存放高频点pivto标记
				omp_set_num_threads(16);
			#pragma omp parallel for
				for (int i = 0; i < HFOripoint.size(); i++) {//对于每一个高频点
					for (int j = 0; j < spt_v_num[HFOripoint[i]]; j++) {//对于高频点中的每一个pivot
						HFPivotFlag[index_p[HFOripoint[i]].spt_v[j]] = 1;//将HFPivotFlag赋值为1
					}
				}
				for (int i = 0; i < numOfVertices; i++) {//对于HFPivotFlag中的每一个pivot
					if (1 == HFPivotFlag[i]) {//如果标记为1
						HFhopPoint.push_back(i);//将pivot号放入HFhopPoint
					}
				}
				HFPivotFlag.clear();
				//将所有高频点中的pivot对应结点的label按照pivot大小顺序存放在低频点的tempLabel中
				for (int i = 0; i < HFhopPoint.size(); i++) {//对于所有高频点中的每一个pivot
					omp_set_num_threads(16);
			#pragma omp parallel for
					for (int j = 0; j < pivotLabel[HFhopPoint[i]].size(); j++) {//对于pivot对应的每一个结点
						if (HFinGraphIndex[pivotLabel[HFhopPoint[i]][j]] == 0) {//如果当前点是低频点
							tempLabel[pivotLabel[HFhopPoint[i]][j]].push_back(HFhopPoint[i]);//将该pivot存入对应结点的临时label中
						}	
					}
				}
				HFhopPoint.clear();
				//以下代码用于查找HFPoint中共有的pivot，并将计数大于2的pivot按编号由小到大放入HFPivot中，作为更新HFPoint时调用updateHFLabel()的参数
				vector<int> pivotCount(numOfVertices, 0);//用来记录所有高频点中重复的pivot，其值大于2则代表此pivot应排在高频点label前部
				vector<int> HFPivot;//存放所有HFPoint中计数大于2的pivot
				for (int i = 0; i < numOfHFpoint; i++) {//对于每一个高频点
					for (int j = 0; j < spt_v_num[HFOripoint[i]]; j++)//对于高频点中的每个pivot
						pivotCount[index_p[HFOripoint[i]].spt_v[j]]++;//使得pivotCount中对应的pivot计数器+1
				}
				for (int i = 0; i < numOfVertices; i++) {
					if (pivotCount[i] > 1) //如果pivot计数大于1，则将该pivot放入HFpivot中
						HFPivot.push_back(i);
				}
				//将所有高频点共有的pivot放入所有高频点的tempLabel中
				for (int i = 0; i < HFPivot.size(); i++) {//对于所有高频点中的每一个pivot
					omp_set_num_threads(16);
			#pragma omp parallel for
					for (int j = 0; j < pivotLabel[HFPivot[i]].size(); j++) {//对于pivot对应的每一个结点
						if (HFinGraphIndex[pivotLabel[HFPivot[i]][j]] == 1)//如果当前点是高频点
							tempLabel[pivotLabel[HFPivot[i]][j]].push_back(HFPivot[i]);//将该pivot存入对应结点的临时label中
					}
				}
				HFPivot.clear();
				pivotCount.clear();

				vector<int> tempJoinCur(numOfVertices,0);//用于表示在tempLabel中，当前比较的位置
				vector<int> CurPreMax(numOfVertices,0);//用来表示每个结点在合并前tempLabel中包含pivot的最大值
				//将tempLabel中所有结果与原label合并，使得label中前半部分为tempLabel中的pivot，后半部分为其它pivot，前后都按降序排列
				omp_set_num_threads(16);
			#pragma omp parallel for
				for (int i = 0; i < numOfVertices; i++) {//对于每一个结点
					tempJoinCur[i] = 0;
					CurPreMax[i] = tempLabel[i].size();
					index_p[i].spt_v[spt_v_num[i]] = CurPreMax[i];//将后半部分pivot的索引号放在无穷大的位置
					for (int j = 0; j < spt_v_num[i]; j++) {//对于结点中的每一个pivot
						if (tempJoinCur[i] < CurPreMax[i]) {//如果tempLabel中前半部分pivot未比较结束
							if (index_p[i].spt_v[j] != tempLabel[i][tempJoinCur[i]])
								tempLabel[i].push_back(index_p[i].spt_v[j]);
							else
								tempJoinCur[i]++;
						}
						else
							tempLabel[i].push_back(index_p[i].spt_v[j]);//如果前半部分已比较结束，则只需将剩下的结点全部存入tempLabel
					}
					for (int j = 0; j < spt_v_num[i]; j++) {//将更新完后的Label赋值给原label，实现label更新
						index_p[i].spt_v[j] = tempLabel[i][j];
					}
				}
				//time_t end_quick = clock();

				update_time = GetCurrentTimeSec()-update_time;
				cout << " Total update time is : " << update_time* 1e6 << endl;
				double utime=update_time/numOfVertices;
				cout << " Each point update time is : " << utime* 1e6 << endl;
				cout<<"update all point successffully!"<<endl;
				//清除缓存
				tempLabel.clear();
				tempJoinCur.clear();
				CurPreMax.clear();

		}

		/**
		 * function used to updtae labels by low frequency points
		 * written by wanjingyi
		 * */
		 void update_low_labels()
		 {
			 vector<NodeID> HFhopPoint;//存放所有HFPoint中无重复的hop点
			 //建立高频点与遍历到的点的对应关系,第一个点为HFPoint序号，第二个值为HFPoint中当前搜索到的标签位置
			 vector<pair<NodeID,NodeID> > indexHFhopPoint(numOfHFpoint ,make_pair(0,0));
			 int i = numOfHFpoint;//i表示HFPoint各标签中剩余未取出的hop点
			 //用于临时存放更新label时选择的高频点对应的hop以及最终更新后的label
			 pair<vector<NodeID>, vector<NodeID> > tempHOP;
			 tempHOP.first.resize(numOfVertices);//初始化临时存放pair,预留10000储存区
			tempHOP.second.resize(numOfVertices);//初始化临时存放pair,预留10000储存区
			cout<<"tempHOP. initialize size is numOfVertices="<<numOfVertices<<endl;
			int tempcur;//表示tempHOP中当前存储到第几个label
			int hIndex, lIndex;//分别存储更新label时高频点和低频点对应的当前检索的label位置
			//初始化indexHFhopPoint为label对应点	
			for (int i = 0; i < HFOripoint.size(); i++) {
				indexHFhopPoint[i].first = HFOripoint[i];
			}
			cout<<"update_labels initialize successffully!"<<endl;

			//将各HFPoint中的label由小到大无重复地放入HFhopPoint中
			for (int i = 0; i != -1;) {
				i = findminIndex(indexHFhopPoint);
				//cout<<"i = "<<i<<endl;
				if (i != -1) {
					//将查找到的最小Label对应的Hop号放入HFhopPoint
					HFhopPoint.push_back(index_p[indexHFhopPoint[i].first].spt_v[indexHFhopPoint[i].second]);
					indexHFhopPoint[i].second++;
				}
			}
			cout<<"findminIndex successffully!"<<endl;

			double update_time=GetCurrentTimeSec();
			//更新图中所有低频点的label
			for (int i = 0; i < HFinGraphIndex.size(); i++) {
				if (!HFinGraphIndex[i]) {//如果该点不是高频点，则更新其label
					hIndex = 0;//初始化hIndex
					lIndex = 0;//初始化lIndex
					tempcur = 0;//初始化tempcur
					//对于每一个低频点，比对HFhopPoint和低频点中的label，将HFhopPoint中的部分放在低频点label的最前面，更新后的label仍按序排列
					while (hIndex < HFhopPoint.size() && lIndex < spt_v_num[i] ) {
						if (HFhopPoint[hIndex] == index_p[i].spt_v[lIndex]) {
							tempHOP.first[tempcur] = index_p[i].spt_v[lIndex];//将选中的共有hop存入tempHOP
							tempHOP.second[tempcur] = index_p[i].spt_d[lIndex];
							index_p[i].spt_v[lIndex] = -1;//将选中的共有hop标志为-1
							hIndex++;
							lIndex++;
							tempcur++;
						}
						else if (HFhopPoint[hIndex] < index_p[i].spt_v[lIndex])
							hIndex++;
						else
							lIndex++;
					}
					//将当前低频点中剩余的非共有hop点放入tempHOP,未比较部分不用处理
					if (tempcur > 0) {
						index_p[i].spt_v[spt_v_num[i]] = tempcur;//将Label后半部分的起始位置索引放在无穷大位置
						for (int j = 0; j < spt_v_num[i]; j++) {
							if (index_p[i].spt_v[j] != -1) {
								tempHOP.first[tempcur] = index_p[i].spt_v[j];
								tempHOP.second[tempcur++] = index_p[i].spt_d[j];
							}
						}
						//将更新后的tempHOP中的内容赋值给原label，实现原label更新
						for (int j = 0; j < tempcur; j++) {
							index_p[i].spt_v[j] = tempHOP.first[j];
							index_p[i].spt_d[j] = tempHOP.second[j];
						}
					}
				}
			}
			update_time = GetCurrentTimeSec()-update_time;
			cout << " Total update time is : " << update_time* 1e6 << endl;
			double utime=update_time/(numOfVertices-numOfHFpoint);
			cout << " Each LFPoint update time is : " << utime* 1e6 << endl;
			cout<<"update low frequency point successffully!"<<endl;

		 }

		/**
		 * function used to updtae labels by high frequency points
		 * written by wanjingyi
		 * */		 
		void update_high_labels()
		{
			 vector<NodeID> HFhopPoint;//存放所有HFPoint中无重复的hop点
			 //建立高频点与遍历到的点的对应关系,第一个点为HFPoint序号，第二个值为HFPoint中当前搜索到的标签位置
			 vector<pair<NodeID,NodeID> > indexHFhopPoint(numOfHFpoint ,make_pair(0,0));
			 int i = numOfHFpoint;//i表示HFPoint各标签中剩余未取出的hop点
			 //用于临时存放更新label时选择的高频点对应的hop以及最终更新后的label
			 pair<vector<NodeID>, vector<NodeID> > tempHOP;
			 tempHOP.first.resize(numOfVertices);//初始化临时存放pair,预留10000储存区
			tempHOP.second.resize(numOfVertices);//初始化临时存放pair,预留10000储存区
			cout<<"tempHOP. initialize size is numOfVertices="<<numOfVertices<<endl;
			int tempcur;//表示tempHOP中当前存储到第几个label
			int hIndex, lIndex;//分别存储更新label时高频点和低频点对应的当前检索的label位置
			//初始化indexHFhopPoint为label对应点	
			for (int i = 0; i < HFOripoint.size(); i++) {
				indexHFhopPoint[i].first = HFOripoint[i];
			}
			cout<<"update_labels initialize successffully!"<<endl;

			//将各HFPoint中的label由小到大无重复地放入HFhopPoint中
			for (int i = 0; i != -1;) {
				i = findminIndex(indexHFhopPoint);
				//cout<<"i = "<<i<<endl;
				if (i != -1) {
					//将查找到的最小Label对应的Hop号放入HFhopPoint
					HFhopPoint.push_back(index_p[indexHFhopPoint[i].first].spt_v[indexHFhopPoint[i].second]);
					indexHFhopPoint[i].second++;
				}
			}
			cout<<"findminIndex successffully!"<<endl;
			vector<int> pivotCount(numOfVertices, 0);//用来记录所有高频点中重复的pivot，其值大于2则代表此pivot应排在高频点label前部
			vector<int> HFPivot;//存放所有HFPoint中计数大于2的pivot
			for (int i = 0; i < numOfHFpoint; i++) {//对于每一个高频点
				for (int j = 0; j < spt_v_num[HFOripoint[i]]; j++)//对于高频点中的每个pivot
					pivotCount[index_p[HFOripoint[i]].spt_v[j]]++;//使得pivotCount中对应的pivot计数器+1
			}
			for (int i = 0; i < numOfVertices; i++) {
				if (pivotCount[i] > 1) //如果pivot计数大于1，则将该pivot放入HFpivot中
					HFPivot.push_back(i);
			}

			//更新图中所有低频点的label
			for (int i = 0; i < HFinGraphIndex.size(); i++) {
				if (!HFinGraphIndex[i]) {//如果该点不是高频点，则更新其label
					hIndex = 0;//初始化hIndex
					lIndex = 0;//初始化lIndex
					tempcur = 0;//初始化tempcur
					//对于每一个低频点，比对HFhopPoint和低频点中的label，将HFhopPoint中的部分放在低频点label的最前面，更新后的label仍按序排列
					while (hIndex < HFhopPoint.size() && lIndex < spt_v_num[i] ) {
						if (HFhopPoint[hIndex] == index_p[i].spt_v[lIndex]) {
							tempHOP.first[tempcur] = index_p[i].spt_v[lIndex];//将选中的共有hop存入tempHOP
							tempHOP.second[tempcur] = index_p[i].spt_d[lIndex];
							index_p[i].spt_v[lIndex] = -1;//将选中的共有hop标志为-1
							hIndex++;
							lIndex++;
							tempcur++;
						}
						else if (HFhopPoint[hIndex] < index_p[i].spt_v[lIndex])
							hIndex++;
						else
							lIndex++;
					}
					//将当前低频点中剩余的非共有hop点放入tempHOP,未比较部分不用处理
					if (tempcur > 0) {
						for (int j = 0; j < spt_v_num[i]; j++) {
							if (index_p[i].spt_v[j] != -1) {
								tempHOP.first[tempcur] = index_p[i].spt_v[j];
								tempHOP.second[tempcur++] = index_p[i].spt_d[j];
							}
						}
						//将更新后的tempHOP中的内容赋值给原label，实现原label更新
						for (int j = 0; j < tempcur; j++) {
							index_p[i].spt_v[j] = tempHOP.first[j];
							index_p[i].spt_d[j] = tempHOP.second[j];
						}
					}
				}else{//更新高频点的label，使得所有高频点的共有pivot在label前按序排列，原label大小不变
					updateHFLabel(HFPivot);
				}
			}
		}

//更新高频点的label，使得所有高频点的共有pivot在label前部按序排列，原label大小不变
void updateHFLabel(const vector<int> mutualPivot) {
	int hIndex{ 0 };//高频点label索引
	int pivotIndex{ 0 };//mutualPivot索引
	vector<pair<int, int>> tempLabel(numOfVertices, make_pair(0, 0));//临时存放更新过程中的高频点label
	int tempcur{ 0 };//表示tempLabel中当前存储到第几个label
	//对于每一个高频点，比对mutualPivot中的元素和高频点中的label，将包含在pivot中的部分放在高频点label的最前面，更新后的label仍按序排列
	for (int i = 0; i < HFOripoint.size(); i++) {
		int hIndex = 0;//初始化hIndex
		int pivotIndex = 0;//初始化pivotIndex
		int tempcur = 0;//初始化tempcur
		while (hIndex < spt_v_num[HFOripoint[i]] && pivotIndex < mutualPivot.size()) {
			if (mutualPivot[pivotIndex] == index_p[HFOripoint[i]].spt_v[hIndex]) {
				tempLabel[tempcur].first = index_p[HFOripoint[i]].spt_v[hIndex];//将选中的pivot存入tempLabel
				tempLabel[tempcur].second = index_p[HFOripoint[i]].spt_d[hIndex];
				index_p[HFOripoint[i]].spt_v[hIndex] = -1;//将选中的共有hop标志为-1
				hIndex++;
				pivotIndex++;
				tempcur++;
			}
			else if (mutualPivot[pivotIndex] < index_p[HFOripoint[i]].spt_v[hIndex])
				pivotIndex++;
			else
				hIndex++;
		}
		//将当前高频点中剩余的非共有pivot放入tempLabel,未比较部分不用处理
		if (tempcur > 0) {
			for (int j = 0; j < spt_v_num[HFOripoint[i]]; j++) {
				if (index_p[HFOripoint[i]].spt_v[j] != -1) {
					tempLabel[tempcur].first = index_p[HFOripoint[i]].spt_v[j];
					tempLabel[tempcur++].second = index_p[HFOripoint[i]].spt_d[j];
				}
			}
			//将更新后的tempLabel中的内容赋值给原label，实现原label更新
			for (int j = 0; j < tempcur; j++) {
				index_p[HFOripoint[i]].spt_v[j] = tempLabel[j].first;
				index_p[HFOripoint[i]].spt_d[j] = tempLabel[j].second;
			}
		}
	}
	
}

		/*
		 *@description: count the label size
		 *@return:pair<long long,double>(<total_label_size,avg_label_size>)
		 *@author: wanjingyi
		 *@date: 2021-01-07
		*/
		void count_label_size(long long& total_label_size,double& avg_label_size){
			long long cnt_total=0;
			double avg_size=0;
			for (int i = 0; i < numOfVertices; ++i) {
				cnt_total+=index_[i].size()-1;
			}
			avg_size=(double)cnt_total/(double)numOfVertices;
			avg_label_size=avg_size;
			total_label_size=cnt_total;
		}

		/*
		 *@description: save and output label size to file spefified for time slices
		 *@author: wanjingyi
		 *@date: 2021-01-07
		*/
		void save_label_size(const char* label_size_file){//vertice index
			ofstream ofs(label_size_file);
			if(!ofs.is_open()) {cerr<<"Cannot open "<<label_size_file<<endl;}
			for (int i = 0; i < numOfVertices; ++i) {
				ofs << index_[i].size()-1 << endl;
			}
			ofs.close();
		}

		/**
		 * function used to save label size to file
		 * written by wanjingyi
		 * */
		void save_label_size(const char* label_size_file,const vector<NodeID>& inv) {
			// string labelSizefile_prefix(label_size_file);
			// string labelSizefile=labelSizefile_prefix.append(".size");
			ofstream ofs(label_size_file);
			if(!ofs.is_open()) {cerr<<"Cannot open "<<label_size_file<<endl;}
			for (int i = 0; i < numOfVertices; ++i) {
				ofs << index_[i].size()-1 << endl;
			}
			ofs.close();
			//output the size by descending order,format:label_size
			string labelSizefile_prefix1(label_size_file);
			string labelSizefile1=labelSizefile_prefix1.append("_byOrder");
			ofstream out(labelSizefile1.c_str());
			if(!out.is_open()) {cerr<<"Cannot open "<<labelSizefile1<<endl;}
			for (int i = 0; i < numOfVertices; ++i) {
				out<<index_[inv[i]].size()-1 << endl;
			}		
			out.close();
		}

		/**
		 * function used to save query distance to file
		 * written by wanjingyi
		 * */
		void save_query_distance(vector<pair<int, int> > queries,const char* query_distance_file,const vector<EdgeWeight>& distance_p, const vector<EdgeWeight>& distance_hf,int numQuery,int warmup,int queryModel=0,int updateModel=0) {
			string query_distance_file_prefix(query_distance_file);
			string appendix="";
			appendix+="_u"+to_string(updateModel)+"_q"+to_string(queryModel)+"_"+to_string(numQuery)+".dis";
			string query_distance_filename=query_distance_file_prefix.append(appendix);
			ofstream ofs(query_distance_filename);
			if(!ofs.is_open()) {cerr<<"Cannot open "<<query_distance_filename<<endl;}
			for (int i = warmup; i < warmup + numQuery; ++i) {
				ofs <<"("<<queries[i].first<<","<<queries[i].second<<") "<<distance_p[i-warmup];
				if(updateModel!=0) ofs<<" "<<distance_hf[i-warmup];
				ofs<<endl;
			}
			ofs.close();
		}

		/**
		 * function used to save query distance to file
		 * written by wanjingyi
		 * */
		void save_query_distance_data(vector<pair<int, int> > queries,const char* query_distance_file,const vector<EdgeWeight>& distance_p, const vector<EdgeWeight>& distance_hf,int timeQuery,int queryModel=0,int updateModel=0) {
			string query_distance_file_prefix(query_distance_file);
			string appendix="";
			appendix+="_u"+to_string(updateModel)+"_q"+to_string(queryModel)+"_"+to_string(timeQuery)+".dis";
			string query_distance_filename=query_distance_file_prefix.append(appendix);
			ofstream ofs(query_distance_filename);
			if(!ofs.is_open()) {cerr<<"Cannot open "<<query_distance_filename<<endl;}
			for (int i = 0; i < timeQuery; ++i) {
				ofs <<"("<<queries[i].first<<","<<queries[i].second<<") "<<distance_p[i];
				if(updateModel!=0) ofs<<" "<<distance_hf[i];
				ofs<<endl;
			}
			ofs.close();
		}

		EdgeWeight query_lf(NodeID s, NodeID t)
		{
			EdgeWeight distance = INF_WEIGHT;

			const index_t_p &idx_s = index_p[s];
			const index_t_p &idx_t = index_p[t];

			_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
			bool isBreak=false; //judge break or not
			for(int i=0,j=0; ; )
				{
					//if(i==0&&j==0) cout<<" s-t:h-l";
					if(j==idx_t.spt_v[spt_v_num[t]] ) break;
					NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
					if (v1 == v2) {
						EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
						if (td < distance) distance = td;
						++i;
						++j;
					} 
					else {
						i += v1 < v2 ? 1 : 0;
						j += v1 > v2 ? 1 : 0;
					}
					
				}
				return distance;
		}
		
		/**
		 * function used to query s-t case:h-h including cache and comparison between no cache
		 * written by wanjingyi
		 * */
		EdgeWeight query_h_h(NodeID s,NodeID t)
		{
			EdgeWeight distance = INF_WEIGHT;
			const index_t_p &idx_s = index_p[s];
			const index_t_p &idx_t = index_p[t];

			_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

			//所有label都Update后，高频点之间的查询,无穷大位置存放后半部分起始点
			for (int i = 0, j = 0; ; ) {//i,j从label前半部分开始搜索
				if (i == idx_s.spt_v[spt_v_num[s]] || j == idx_t.spt_v[spt_v_num[t]]) break;//前半部分搜索到分界点位置时结束
				NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
				if (v1 == v2) {
					EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
					if (td < distance) distance = td;
					++i;
					++j;
				}
				else {
					i += v1 < v2 ? 1 : 0;
					j += v1 > v2 ? 1 : 0;
				}
			}
			return distance;
		}

		/**
		 * function used to query s-t case:h-l
		 * written by wanjingyi
		 * */
		EdgeWeight query_h_l(NodeID s,NodeID t)
		{
			EdgeWeight distance = INF_WEIGHT;
			const index_t_p &idx_s = index_p[s];
			const index_t_p &idx_t = index_p[t];

			_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
			//所有label都update之后，高频点和低频点之间的查询,其中idx_s为高频点，idx_t为低频点
			for (int i = 0, imiddle = idx_s.spt_v[spt_v_num[s]], j = 0; ; )
			 {//imiddle表示高频点后半部分的起始点(该索引号存放在无穷大的位置)
				if (j==idx_t.spt_v[spt_v_num[t]] || (i == idx_s.spt_v[spt_v_num[s]] && imiddle == spt_v_num[s])) break;//低频点当前pivot大于下一个pivot或高频点中前半部分和后半部分都查询结束时，整个查询结束，低频点的无穷大保证了正确性
				NodeID v1 = idx_s.spt_v[i], v1middle = idx_s.spt_v[imiddle], v2 = idx_t.spt_v[j];
				if (v1 == v2) {
					EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
					if (td < distance) distance = td;
					++i;
					++j;
				}else{
					if (v1middle == v2) {
						EdgeWeight td = idx_s.spt_d[imiddle] + idx_t.spt_d[j];
						if (td < distance) distance = td;
						++imiddle;
						++j;
					}else {//高频点前后查询中都未找到与低频点共有的pivot,则pivot编号最小的+1
						if (v1 > v2) {
							if (v2 > v1middle)
								++imiddle;//v1middle最小
							else
								++j;//v2最小
						}
						else {
							if (v1 > v1middle) {
								++imiddle;//v1middle最小
							}
							else
								++i;//v1最小
						}
					}
				}
			 }
			return distance;

		}

		/**
		 * function used to query s-t case:l-h
		 * written by wanjingyi
		 * */
		EdgeWeight query_l_h(NodeID s,NodeID t)
		{
			NodeID tmp=s;s=t;t=tmp; //swap s and t
			EdgeWeight distance = INF_WEIGHT;
			const index_t_p &idx_s = index_p[s];
			const index_t_p &idx_t = index_p[t];

			_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

			//所有label都update之后，高频点和低频点之间的查询,其中idx_s为高频点，idx_t为低频点
			for (int i = 0, imiddle = idx_s.spt_v[spt_v_num[s]], j = 0; ; )
			 {//imiddle表示高频点后半部分的起始点(该索引号存放在无穷大的位置)
				if (j==idx_t.spt_v[spt_v_num[t]] || (i == idx_s.spt_v[spt_v_num[s]] && imiddle == spt_v_num[s])) break;//低频点当前pivot大于下一个pivot或高频点中前半部分和后半部分都查询结束时，整个查询结束，低频点的无穷大保证了正确性
				NodeID v1 = idx_s.spt_v[i], v1middle = idx_s.spt_v[imiddle], v2 = idx_t.spt_v[j];
				if (v1 == v2) {
					EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
					if (td < distance) distance = td;
					++i;
					++j;
				}else{
					if (v1middle == v2) {
						EdgeWeight td = idx_s.spt_d[imiddle] + idx_t.spt_d[j];
						if (td < distance) distance = td;
						++imiddle;
						++j;
					}else {//高频点前后查询中都未找到与低频点共有的pivot,则pivot编号最小的+1
						if (v1 > v2) {
							if (v2 > v1middle)
								++imiddle;//v1middle最小
							else
								++j;//v2最小
						}
						else {
							if (v1 > v1middle) {
								++imiddle;//v1middle最小
							}
							else
								++i;//v1最小
						}
					}
				}
			 }
			return distance;
		}

		/**
		 * function used to query s-t case:l-l
		 * written by wanjingyi
		 * */
		EdgeWeight query_l_l(NodeID s,NodeID t)
		{
			EdgeWeight distance = INF_WEIGHT;
			const index_t_p &idx_s = index_p[s];
			const index_t_p &idx_t = index_p[t];

			_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
			//所有label都Update后两点都为低频点的查询,无穷大位置存放后半部分pivot的起始位置
			for (int i = 0, j = 0; ; ) 
			{//i,j从label前半部分开始搜索
				if (i == idx_s.spt_v[spt_v_num[s]] || j == idx_t.spt_v[spt_v_num[t]]) break;//前半部分搜索到分界点位置时结束
				NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
				if (v1 == v2) {
					EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
					if (td < distance) distance = td;
					++i;
					++j;
				}
				else {
					i += v1 < v2 ? 1 : 0;
					j += v1 > v2 ? 1 : 0;
				}
			}
			for (int iMiddle = idx_s.spt_v[spt_v_num[s]], jMiddle = idx_t.spt_v[spt_v_num[t]]; ; )
			 {//imiddle和jmiddle从label后半部分开始搜索
				if (iMiddle == spt_v_num[s] || jMiddle == spt_v_num[t]) break;//后半部分搜索到pivot最后一个位置之后结果
				NodeID v1Middle = idx_s.spt_v[iMiddle], v2Middle = idx_t.spt_v[jMiddle];
				if (v1Middle == v2Middle) {
					EdgeWeight tdMiddle = idx_s.spt_d[iMiddle] + idx_t.spt_d[jMiddle];
					if (tdMiddle < distance) distance = tdMiddle;
					++iMiddle;
					++jMiddle;
				}
				else {
					iMiddle += v1Middle < v2Middle ? 1 : 0;
					jMiddle += v1Middle > v2Middle ? 1 : 0;
				}
			 }	
			return distance;
		}

		/**
		 * function used to query all s-t cases
		 * written by wanjingyi
		 * */
		EdgeWeight query_all(NodeID s,NodeID t)
		{
			//0-random(default) 1-(s-t):(h-l),2-(s-t):(l:h),3-(s-t):(h-h),4-(s-t):(l-l)
			if(HFinGraphIndex[s]&&!HFinGraphIndex[t]) return query_h_l(s,t);
			else if(!HFinGraphIndex[s]&&HFinGraphIndex[t]) return query_l_h(s,t);
			else if(HFinGraphIndex[s]&&!HFinGraphIndex[t]) return query_h_h(s,t);
			else return query_l_l(s,t);
		}

		/**
		 * function used to query all s-t cases use cache
		 * written by wanjingyi
		 * */
		EdgeWeight query_all_cached(NodeID s,NodeID t)
		{
			//0-random(default) 1-(s-t):(h-l),2-(s-t):(l:h),3-(s-t):(h-h),4-(s-t):(l-l)
			if(HFinGraphIndex[s]&&!HFinGraphIndex[t]) return query_h_l(s,t);
			else if(!HFinGraphIndex[s]&&HFinGraphIndex[t]) return query_l_h(s,t);
			else if(HFinGraphIndex[s]&&!HFinGraphIndex[t]){
				if(hf_chache[s][t]!=INF_WEIGHT) return hf_chache[s][t];
				double result_distance=query_h_h(s,t);
				hf_chache[s][t]=result_distance;
				return result_distance;
			}
			else return query_l_l(s,t);
		}

		/**
		 * function used to query distinguished by high frequency and low frequency
		 * written by wanjingyi
		 * */
		EdgeWeight query_hf(NodeID s, NodeID t) {
			EdgeWeight distance = INF_WEIGHT;
			if(!HFinGraphIndex[s]&&HFinGraphIndex[t]){ //s-t:l-h
				NodeID tmp=s;
				s=t;
				t=tmp;
			}
			const index_t_p &idx_s = index_p[s];
			const index_t_p &idx_t = index_p[t];

			_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

			if(HFinGraphIndex[s]&&HFinGraphIndex[t]){ // (s,t)both high 
				for(int i=0,j=0; ; )
				{
					NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
					if (v1 == numOfVertices|| v2== numOfVertices) break;  // Sentinel
					if(v1>idx_s.spt_v[i+1] || v2>idx_t.spt_v[j+1]){
						break;
					}
					if (v1 == v2) {
						EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
						if (td < distance) distance = td;
						++i;
						++j;
					} 
					else {
						i += v1 < v2 ? 1 : 0;
						j += v1 > v2 ? 1 : 0;
					}

				}
			}else if(HFinGraphIndex[s]&&!HFinGraphIndex[t]){ //s-high,t-low
			//imiddle表示高频点后半部分的起始点(该索引号存放在无穷大的位置)
				for (int i = 0, imiddle = idx_s.spt_v[spt_v_num[s]], j = 0; ; ) {
					//低频点当前pivot大于下一个pivot或高频点中前半部分和后半部分都查询结束时，整个查询结束，低频点的无穷大保证了正确性
					if (idx_t.spt_v[j] > idx_t.spt_v[j + 1] || (i == idx_s.spt_v[spt_v_num[s]] && imiddle == spt_v_num[s])) break;
					NodeID v1 = idx_s.spt_v[i], v1middle = idx_s.spt_v[imiddle], v2 = idx_t.spt_v[j];
					if (v1 == numOfVertices) break;  // Sentinel
					if (i < idx_s.spt_v[spt_v_num[s]]) {//i只查询高频点前半部分
						if (v1 == v2) {
							EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
							if (td < distance) distance = td;
							++i;
							++j;
						}else{//若前半部分当前未查到相同的pivot，则要进入后半部分进行查询
							while (v1middle != v2 && imiddle < spt_v_num[s]) {//imiddle只查询高频点后半部分
								++imiddle;
								v1middle == idx_s.spt_v[imiddle];//无穷大位置用来存放后半部分起始点标记，保证了数组不会越界
							}
							if (v1middle == v2) {
								EdgeWeight td = idx_s.spt_d[imiddle] + idx_t.spt_d[j];
								if (td < distance) distance = td;
								++imiddle;
								++j;
							}else{
								if (v1 > v2) {
									if (v2 > v1middle)
										++imiddle;//v1middle最小
									else
										++j;//v2最小
								}
								else {
									if (v1 > v1middle) {
										++imiddle;//v1middle最小
									}
									else ++i;//v1最小
								}
							}
						}

					}
					else {//前半部分都查询结束了，那么只剩后半部分要查询
						if (v1middle == v2) {
							EdgeWeight td = idx_s.spt_d[imiddle] + idx_t.spt_d[j];
							if (td < distance) distance = td;
							++imiddle;
							++j;
						}
						else {
							imiddle += v1middle < v2 ? 1 : 0;
							j += v1middle > v2 ? 1 : 0;
						}
				 	}
				}
			}else if(!HFinGraphIndex[s]&&HFinGraphIndex[t]){ //s-low,t-high
				cout<<" s-t:l-h"<<endl;
			}else{ //low-low frequency point query
				// for (int i = 0, j = 0; ; ) {
				// NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

				// if (v1 == numOfVertices) break;  // Sentinel

				// 	if (v1 == v2) {
				// 	EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				// 	if (td < distance) distance = td;
				// 	++i;
				// 	++j;
				// 	} 
				// 	else {
				// 		i += v1 < v2 ? 1 : 0;
				// 		j += v1 > v2 ? 1 : 0;
				// 	}
				// }
					int flag{ 1 };
					for (int i = 0, j = 0, iBack = spt_v_num[s]-1, jBack = spt_v_num[t]; ; ) {//i,j从label前向搜索，iBack和jBack从label后向搜索
						if (i == iBack || j == jBack|| 0==flag) break;
						flag = 0;//flag用于排除两个搜索悬空无法相等的特殊情况
						NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j], v1Back = idx_s.spt_v[iBack], v2Back = idx_t.spt_v[jBack];
						if (v1 == numOfVertices || v2 == numOfVertices) break;  // Sentinel
						if (v1 < idx_s.spt_v[i + 1] && v2 < idx_t.spt_v[j + 1]) {//进入前向搜索的条件是label中当前点序号小于下一点序号
							flag = 1;
							if (v1 == v2) {
								EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
								if (td < distance) distance = td;
								++i;
								++j;
							}
							else {
								i += v1 < v2 ? 1 : 0;
								j += v1 > v2 ? 1 : 0;
							}
						}
						if (v1Back > idx_s.spt_v[iBack - 1] && v2Back > idx_t.spt_v[jBack - 1]) {//进入后向搜索的条件是label中当前点序号大于下一点序号
							flag = 1;
							if (v1Back == v2Back) {
								EdgeWeight tdBack = idx_s.spt_d[iBack] + idx_t.spt_d[jBack];
								if (tdBack < distance) distance = tdBack;
								--iBack;
								--jBack;
							}
							else {
								iBack -= v1 > v2 ? 1 : 0;
								jBack -= v1 < v2 ? 1 : 0;
							}
						}
					}

			}
			
			return distance;
	}


};

class DLabel : public Label {
	public:
	vector<index_t> bindex_; // Backward labels.

	index_t_p* bindex_p;
	
	
	two_index_t_p* b_two_index_p;


	DLabel() {
		index_.resize(numOfVertices);
		bindex_.resize(numOfVertices);
	}
	
	~DLabel() {  
		Free();
	}

	EdgeWeight query(NodeID s, NodeID t) {

		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& bindex_t = bindex_[t].spt_v;
		vector<EdgeWeight>& bindex_t_d = bindex_[t].spt_d;

		for (int i = 0, j = 0; i < index_s.size(), j < bindex_t.size(); ) {
			if (index_s[i] == bindex_t[j]) {
				distance = min(distance, (EdgeWeight)(index_s_d[i] + bindex_t_d[j]));
				++i;
				++j;
			}
			else {
				if (index_s[i] < bindex_t[j])
					++i;
				else
					++j;
			}
		}

		return distance;
	}

	EdgeWeight query(NodeID s, NodeID t, NodeID& meet, EdgeWeight& dis1, EdgeWeight& dis2) {
		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

		vector<NodeID>& bindex_t = bindex_[t].spt_v;
		vector<EdgeWeight>& bindex_t_d = bindex_[t].spt_d;

		meet = numeric_limits<NodeID>::max();
		dis1 = numeric_limits<EdgeWeight>::max();
		dis2 = numeric_limits<EdgeWeight>::max();
		for (int i = 0, j = 0; i < index_s.size(), j < bindex_t.size(); ) {
			if (index_s[i] == bindex_t[j]) {
				if (distance > (EdgeWeight)(index_s_d[i] + bindex_t_d[j])) {
					distance = (EdgeWeight)(index_s_d[i] + bindex_t_d[j]);
					meet = index_s[i];
					dis1 = index_s_d[i];
					dis2 = bindex_t_d[j];
				}
				++i;
				++j;
			}
			else {
				if (index_s[i] < bindex_t[j])
					++i;
				else
					++j;
			}
		}

		return distance;
	}
	
	inline EdgeWeight query_p(NodeID s, NodeID t) {

		//EdgeWeight distance = INF_WEIGHT;
		//
		////const index_t_p &idx_s = index_p[s];
		////const index_t_p &idx_t = bindex_p[t];

		//NodeID *vs = index_p[s].spt_v;
		//NodeID *vt = bindex_p[t].spt_v;
		//EdgeWeight* ws = index_p[s].spt_d;
		//EdgeWeight* wt = bindex_p[t].spt_d;


		//_mm_prefetch(vs, _MM_HINT_T0);
		//_mm_prefetch(vt, _MM_HINT_T0);
		//_mm_prefetch(ws, _MM_HINT_T0); 
		//_mm_prefetch(wt, _MM_HINT_T0);

		//for (unsigned i = 0, j = 0; ; ) {
		//	if (*(vs + i) == *(vt + j)) {
		//		if (*(vs + i) == numOfVertices) break;  // Sentinel
		//		EdgeWeight td = *(ws + i) + *(wt + j);
		//		if (td < distance) distance = td;
		//		++i;
		//		++j;
		//	}
		//	else {
		//		i += *(vs + i) < *(vt + j) ? 1 : 0;
		//		j += *(vs + i) > *(vt + j) ? 1 : 0;
		//	}
		//}
		//return distance;

		
		EdgeWeight distance = INF_WEIGHT;

		const index_t_p &idx_s = index_p[s];
		const index_t_p &idx_t = bindex_p[t];


		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == v2) {
				if (v1 == numOfVertices) break;  //two_index_t_p Sentinel
				EdgeWeight td = idx_s.spt_d[i] + idx_t.spt_d[j];
				if (td < distance) distance = td;
				++i;
				++j;
			}
			else {
				i += v1 < v2 ? 1 : 0;
				j += v1 > v2 ? 1 : 0;
			}
		}
		return distance; 
	}

	EdgeWeight query_with_info(NodeID s, NodeID t, query_info& q_info) {

		double stime = GetCurrentTimeSec();

		EdgeWeight distance = INF_WEIGHT;
		vector<NodeID>& index_s = index_[s].spt_v;
		vector<EdgeWeight>& index_s_d = index_[s].spt_d;

	//	vector<NodeID>& index_t = index_[t].spt_v;
	//	vector<EdgeWeight>& index_t_d = index_[t].spt_d;	
		vector<NodeID>& bindex_t = bindex_[t].spt_v;
		vector<EdgeWeight>& bindex_t_d = bindex_[t].spt_d;

		q_info.meet_node = numOfVertices;
		double meet_distance;

		for (int i = 0, j = 0; i < index_s.size(), j < bindex_t.size(); ) {
			if (index_s[i] == bindex_t[j]) {
				meet_distance = (EdgeWeight)(index_s_d[i++] + bindex_t[j++]);
				if (distance >  meet_distance) {
					distance = meet_distance;
					q_info.meet_node = index_s[i];
				}
			}
			else {
				if (index_s[i] < bindex_t[j])
					++i;
				else
					++j;
			}
		};

		stime = GetCurrentTimeSec() - stime;

		q_info.time_cost = stime;

		if (index_s.size() < bindex_t.size())
			q_info.search_len = index_s.size();
		else
			q_info.search_len = bindex_t.size();

		return distance;
	}


	void append(NodeID v, NodeID root, EdgeWeight distance, bool forward) { // forward(backward) search from root to vertex v.
		if (forward) { // forward search from root to vertex v, hence append (root, distance) to backward index of vertex v.
			bindex_[v].spt_v.push_back(root);
			bindex_[v].spt_d.push_back(distance);
		}
		else { // backward search from root to vertex v, hence append (root, distance) to forward index of vertex v.
			index_[v].spt_v.push_back(root);
			index_[v].spt_d.push_back(distance);
		}
	}

	void Free() {
		if (index_.size() == 0 || bindex_.size() == 0) return;
		for (int v = 0; v < numOfVertices; ++v) {
			index_[v].spt_v.clear();
			index_[v].spt_d.clear();
			if (DIRECTED_FLAG == true) {
				bindex_[v].spt_v.clear();
				bindex_[v].spt_d.clear();
			}
		}
		index_.clear();
		bindex_.clear();
	}

	double avg_size() {
		double total = 0;
		for (int i = 0; i < numOfVertices; ++i) {
			total += index_[i].spt_v.size() ;
			total += bindex_[i].spt_v.size();
		}

		double avg = total / numOfVertices / 2 - 1; // We do not count the trivial labels (V, INF_WEIGHT).

		return avg;
	}

	void print_stat() {
		cout << "Average Label Size: " << avg_size() << endl;
		//cout << "Maximum Label Size: " << max_size() << endl;
	}

	void save_labels(const char* save_filename) {
		ofstream ofs(save_filename, ios::binary | ios::out);

		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			int isize = index_[v].size();
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
				ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
			}
			int bisize = bindex_[v].size();
			ofs.write((const char*)&bisize, sizeof(bisize));
			for (NodeID i = 0; i < bindex_[v].size(); ++i) {
				ofs.write((const char*)&bindex_[v].spt_v[i], sizeof(bindex_[v].spt_v[i]));
				ofs.write((const char*)&bindex_[v].spt_d[i], sizeof(bindex_[v].spt_d[i]));
			}
		}
		ofs.close();
	}

	void load_labels(const char* load_filename) {
		cout << "Loading Labels" << endl;
	/*
		for (NodeID v = 0; v < numOfVertices; ++v) {
			free(index_p[v].spt_v);
			free(index_p[v].spt_d);
		}*/

		//free(index_p);
		index_p = NULL;
		bindex_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		index_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));
		bindex_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));

		cout << numOfVertices << " vertices." << endl;

		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));

			idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_v[i] = hub;
				idx.spt_d[i] = hub_weight;

			}

			//	index_[v].spt_v.resize(isize);
			//	index_[v].spt_d.resize(isize);


			index_t_p &bidx = bindex_p[v];
			ifs.read((char*)&isize, sizeof(isize));
			bidx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
			bidx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));

			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;
				bidx.spt_v[i] = hub;
				bidx.spt_d[i] = hub_weight;

			}

		}
		ifs.close();

		/*
		index_.clear();
		bindex_.clear();
		ifs.open(load_filename, ios::binary | ios::in);
		
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		index_.resize(numOfVertices);
		bindex_.resize(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {

			ifs.read((char*)&isize, sizeof(isize));
			index_[v].spt_v.resize(isize);
			index_[v].spt_d.resize(isize);
			for (NodeID i = 0; i < index_[v].size(); ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				index_[v].spt_v[i] = hub;
				index_[v].spt_d[i] = hub_weight;
			}

			ifs.read((char*)&isize, sizeof(isize));
			bindex_[v].spt_v.resize(isize);
			bindex_[v].spt_d.resize(isize); 			
			for (NodeID i = 0; i < bindex_[v].size(); ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				bindex_[v].spt_v[i] = hub;
				bindex_[v].spt_d[i] = hub_weight;
			}
		}
		ifs.close();
		*/
	/*	for (int i = 0; i < numOfVertices; ++i) {
			for (int j = 0; j < index_[i].size(); ++j)
				if (index_[i].spt_v[j] != index_p[i].spt_v[j])
					cout << "warning." << endl;
		}*/
		
	}

	void load_labels(const char* load_filename,const vector<bool>& flag) {
		// for (NodeID v = 0; v < numOfVertices; ++v) {
		// 	free(index_p[v].spt_v);
		// 	free(index_p[v].spt_d);
		// }
		// free(index_p);
		index_p = NULL;
		bindex_p = NULL;

		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;

		index_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));
		bindex_p = (index_t_p*)memalign(64, numOfVertices * sizeof(index_t_p));

		//cout << numOfVertices << " vertices." << endl;

		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t_p &idx = index_p[v];
			ifs.read((char*)&isize, sizeof(isize));
			if(flag[v]){
				idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
				idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
			}
			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;
				if(flag[v]){
					idx.spt_v[i] = hub;
					idx.spt_d[i] = hub_weight;
				}
			}

			//	index_[v].spt_v.resize(isize);
			//	index_[v].spt_d.resize(isize);


			index_t_p &bidx = bindex_p[v];
			ifs.read((char*)&isize, sizeof(isize));
			if(flag[v]){
				bidx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
				bidx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
			}
			for (NodeID i = 0; i < isize; ++i) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;
				if(flag[v]){
					bidx.spt_v[i] = hub;
					bidx.spt_d[i] = hub_weight;
				}
			}
		}
		ifs.close();
	}
	
	void convert_to_fewerbit(){
		
		two_index_p = NULL;
		b_two_index_p = NULL;
		two_index_p = (two_index_t_p*)memalign(64, numOfVertices * sizeof(two_index_t_p));
		b_two_index_p = (two_index_t_p*)memalign(64, numOfVertices * sizeof(two_index_t_p));

		for (NodeID v = 0; v < numOfVertices; ++v) {
			two_index_t_p &idx = two_index_p[v];
			
			index_t_p &idx_original = index_p[v];
			
			NodeID isize = 0;
			for(NodeID i = 0; idx_original.spt_v[i] < UCHAR_MAX; ++i){
				++isize;
			}
			

			idx.spt_lv = (uint8_t*)memalign(64, (isize + 1) * sizeof(uint8_t));
			idx.spt_ld = (EdgeWeight*)memalign(64, (isize + 1) * sizeof(EdgeWeight));

		//	index_[v].spt_v.resize(isize);
		//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_lv[i] = idx_original.spt_v[i];
				idx.spt_ld[i] = idx_original.spt_d[i];

			}
			 
			idx.spt_lv[isize] = UCHAR_MAX;
			idx.spt_ld[isize] = INF_WEIGHT;
			
			NodeID larger_size = 0;
			for(NodeID i = isize; idx_original.spt_v[i] != numOfVertices; ++i){
				++larger_size;
			}
			
			
			idx.spt_v = (NodeID*)memalign(64, larger_size * sizeof(NodeID));
			idx.spt_d = (EdgeWeight*)memalign(64, larger_size * sizeof(EdgeWeight));
			
			for (NodeID i = 0; i < larger_size; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				idx.spt_v[i] = idx_original.spt_v[i + isize];
				idx.spt_d[i] = idx_original.spt_d[i + isize];

			}			
			
			two_index_t_p &b_idx = b_two_index_p[v];
			
			index_t_p &b_idx_original = bindex_p[v];
			
			isize = 0;
			for(NodeID i = 0; b_idx_original.spt_v[i] < UCHAR_MAX; ++i){
				++isize;
			}
			

			b_idx.spt_lv = (uint8_t*)memalign(64, (isize + 1) * sizeof(uint8_t));
			b_idx.spt_ld = (EdgeWeight*)memalign(64, (isize + 1) * sizeof(EdgeWeight));

		//	index_[v].spt_v.resize(isize);
		//	index_[v].spt_d.resize(isize);

			for (NodeID i = 0; i < isize; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;

				b_idx.spt_lv[i] = b_idx_original.spt_v[i];
				b_idx.spt_ld[i] = b_idx_original.spt_d[i];

			}
			 
			b_idx.spt_lv[isize] = UCHAR_MAX;
			b_idx.spt_ld[isize] = INF_WEIGHT;
			
			larger_size = 0;
			for(NodeID i = isize; b_idx_original.spt_v[i] != numOfVertices; ++i){
				++larger_size;
			}
			
			
			b_idx.spt_v = (NodeID*)memalign(64, larger_size * sizeof(NodeID));
			b_idx.spt_d = (EdgeWeight*)memalign(64, larger_size * sizeof(EdgeWeight));
			
			for (NodeID i = 0; i < larger_size; ++i) {
				uint8_t hub;
				EdgeWeight hub_weight;
				//index_[v].spt_v[i] = hub;
				//index_[v].spt_d[i] = hub_weight;
 
				b_idx.spt_v[i] = b_idx_original.spt_v[i + isize];
				b_idx.spt_d[i] = b_idx_original.spt_d[i + isize];

			}			
		}		
	}
	
	
	void save_labels_iteration_stats(const char* save_filename) {

		vector<NodeID> stat(numOfVertices);
		for (NodeID v = 0; v < numOfVertices; ++v) {
			for (NodeID i = 0; i < index_[v].size(); ++i)
				stat[index_[v].spt_v[i]]++;
			for (NodeID i = 0; i < bindex_[v].size(); ++i)
				stat[bindex_[v].spt_v[i]]++;
		}

		ofstream ofs(save_filename);

		for (NodeID v = 0; v < numOfVertices; ++v) {
			ofs << stat[v] << endl;
		}
		ofs.close();
	}

};

/**
 * @description: directed hf label defined by wanwan
 * @param {*}
 * @return {*}
 * @author: Wan Jingyi
 */
class HFDLabel : public DLabel {
	public:
		vector<int> queryTime;//store each vertex's query frequency 
		vector<int> queryTime_s;//store each vertex's query frequency as start
		vector<int> queryTime_t;//store each vertex's query frequency as destination
		vector<NodeID> spt_v_num; //store the size of each vertice's forward labels
		vector<NodeID> spt_v_num_f; //store the size of each vertice's forward labels
		vector<NodeID> spt_v_num_r; //store the size of each vertice's backward labels
		vector<NodeID> rank;//store the verices' rank
		int numOfHFpoint=0; //number of High frequency points
		vector<NodeID> HFPoint; //store the high frequency point referred orders
		vector<bool> HFinGraphIndex;//HFinGraphIndex is index of high frequency point in original graph, 0 represents low frequency point，1 represents high frequency point
		int maxQueryTime=0;//num of max query time total
		int minQueryTime=INT_MAX;//num of min query time total
		int maxQueryTime_s=0;//num of max query time as start
		int minQueryTime_s=INT_MAX;//num of min query time as start
		int maxQueryTime_t=0;//num of max query time as destination
		int minQueryTime_t=INT_MAX;//num of min query time as destination
		/**************constructions and deconstructions*************/
		HFDLabel(){}
		~HFDLabel(){
			Free();
		}
	
	//*******************util functions*****************
		/*
		 *@description: save and output label size to file
		 *format:label_size r_label_size 
		 *@author: wanjingyi
		 *@date: 2021-01-07
		*/
		void save_label_size(const char* label_size_file){//vertice index
			ofstream ofs(label_size_file);
			if(!ofs.is_open()) {cerr<<"Cannot open "<<label_size_file<<endl;}
			for (int i = 0; i < numOfVertices; ++i) {
				ofs << index_[i].size()-1 <<" "<<bindex_[i].size()-1<<endl;
			}
			ofs.close();
		}

		/**
		 * function used to save label size to file
		 * written by wanjingyi
		 * label_size r_label_size 
		 * */
		void save_label_size(const char* label_size_file,const vector<NodeID>& inv) {
			// string labelSizefile_prefix(label_size_file);
			// string labelSizefile=labelSizefile_prefix.append(".size");
			ofstream ofs(label_size_file);
			if(!ofs.is_open()) {cerr<<"Cannot open "<<label_size_file<<endl;}
			for (int i = 0; i < numOfVertices; ++i) {
				ofs << index_[i].size()-1<<" "<<bindex_[i].size()-1 << endl;
			}
			ofs.close();
			//output the size by descending order,format:label_size
			string labelSizefile_prefix1(label_size_file);
			string labelSizefile1=labelSizefile_prefix1.append("_byOrder");
			ofstream out(labelSizefile1.c_str());
			if(!out.is_open()) {cerr<<"Cannot open "<<labelSizefile1<<endl;}
			for (int i = 0; i < numOfVertices; ++i) {
				out<<index_[inv[i]].size()-1<<" "<<bindex_[inv[i]].size()-1 << endl;
			}		
			out.close();
		}


	/**
	 * function used to serialize vertice's labels to file,format as follows:
	 *	format:
	 * v:
	 *  label_size (w0,d(v,w0)) (w1,d(v,w1))...(wn,d(v,wn))
	 *  r_label_size (w0,d(v,w0)) (w1,d(v,w1))...(wn,d(v,wn))
	 * written by wanjingyi
	 * */
	void write_labels(const char* write_filename,const vector<NodeID>& inv){
		string write_filename_prefix(write_filename);
		string write_filename1=write_filename_prefix.append("_label.list");
		ofstream ofs(write_filename1.c_str());
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			ofs<<v<<":"<<endl;
			//backward
			NodeID isize = index_[v].size();
			ofs <<isize<<" ";
			for (NodeID i = 0; i < isize; ++i) {
				ofs<<" "<<'('<<index_[v].spt_v[i]<<","<<index_[v].spt_d[i]<<")";
			}
			ofs<<endl;
			//forward
			isize = bindex_[v].size();
			ofs <<isize<<" ";
			for (NodeID i = 0; i < isize; ++i) {
				ofs<<" "<<'('<<bindex_[v].spt_v[i]<<","<<bindex_[v].spt_d[i]<<")";
			}
			ofs<<endl;	
		}
		ofs.close();
		
		//write labels with original NodeId in graph
		string write_filename2=write_filename_prefix.append("o");
		ofstream out(write_filename2.c_str());
		NodeID isize;
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			out<<v<<":"<<endl;
			//backward
			isize = index_[v].size();
			out<<isize-1<<" ";
			for (NodeID i = 0; i < isize-1; ++i) {
				out<<" "<<'('<<inv[index_[v].spt_v[i]]<<","<<index_[v].spt_d[i]<<")";
			}
			out<<endl;
			//forward
			isize = bindex_[v].size();
			out<<isize-1<<" ";
			for (NodeID i = 0; i < isize-1; ++i) {
				out<<" "<<'('<<inv[bindex_[v].spt_v[i]]<<","<<bindex_[v].spt_d[i]<<")";
			}
			out<<endl;
		}			
		out.close();
		return;
	}

	/**
	 * @description: load hfpoint and its queryTime_s and queryTime_t
	 * @param {int} hf_rate
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void load_hfpoint_and_qt(char* load_filename,int hf_rate=0){
		//initialize iterms
		//queryTime.resize(numOfVertices,0);
		queryTime_s.resize(numOfVertices,0);
		queryTime_t.resize(numOfVertices,0);
		HFinGraphIndex.resize(numOfVertices,false);
		numOfHFpoint=0;
		//if(hf_rate==0) numOfHFpoint=numOfVertices;
		//else numOfHFpoint=static_cast<int>((double)(numOfVertices*hf_rate)/(double)HF_DIVIDION);
		std::ifstream in(load_filename);//input HFPoint file to ifstream
		if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
		int s,t,qt,i=0;
		char line[24];
		//read each line representing HFpoint to vector 
		while (in.getline(line,sizeof(line))){
			std::stringstream hp(line);
			hp>>s>>t>>qt;
			queryTime_s[s]+=qt;
			queryTime_t[t]+=qt;
			HFinGraphIndex[s]=true;
			HFinGraphIndex[t]=true;

		}
		for(size_t v=0;v<numOfVertices;++v){
			int qt_s=queryTime_s[v];
			int qt_t=queryTime_t[v];
			if(qt_s>maxQueryTime_s) maxQueryTime_s=qt_s;
			if(qt_s<minQueryTime_s) minQueryTime_s=qt_s;
			if(qt_t>maxQueryTime_t) maxQueryTime_t=qt_t;
			if(qt_t<minQueryTime_t) minQueryTime_t=qt_t;
		}
		in.close();
		std::cout<<"maxQueryTime_s = "<<maxQueryTime_s<<std::endl;
		std::cout<<"maxQueryTime_t = "<<maxQueryTime_t<<std::endl;
		std::cout<<"Load_hfpoint_and_qt successfully!"<<std::endl;
		return;
	}

	/**
  * @description: compute the query cost for directed graph
  * @param {const} char
  * @return {*}
  * @author: Wan Jingyi
  */ 
 	void save_anaylysis_size(const char* write_filename){
		double performance_result=0;//total performance function
		//total
		long long total_sum_size=0,hf_sum_size=0;
		double total_ave_size=0,hf_ave_size=0;
		//open file
		string asize_filename(write_filename);
		asize_filename.append("_analysis.size");
		ofstream ofs_size(asize_filename.c_str());
		if(!ofs_size.is_open()) {cerr<<"Cannot open "<<asize_filename<<endl;}

		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			NodeID isize_s = index_[v].size()-1;
			NodeID isize_t = bindex_[v].size()-1;
			total_sum_size+=isize_s+isize_t;
			if(HFinGraphIndex[v]){
				hf_sum_size+=isize_s+isize_t;
			}
			//compute the query performance function
			double ratio_s=(double)(queryTime_s[v])/(double)DIVISION_FACTOR;
			double ratio_t=(double)(queryTime_t[v])/(double)DIVISION_FACTOR;
			performance_result+=ratio_s*(double)isize_s+ratio_t*(double)isize_t;
			//cout<<v<<": qtime="<<queryTime[v]<<", size="<<isize<<endl;
		}
		total_ave_size= (double) total_sum_size/(double) numOfVertices;
		hf_ave_size= (double) hf_sum_size/(double) numOfHFpoint;
		std::cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
		std::cout<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
		ofs_size<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
		ofs_size<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
		//write performance function to file
		std::cout<<"standard performance_result = "<<performance_result<<endl;
		ofs_size<<"nomalization performance_result = "<<performance_result<<endl;
		ofs_size.close();	
		return;
	}

	/** experiment model output for PLL
	 * @description: 
	 * @param {*}
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void append_experiment_result_pll(char* write_filename){
		long long total_sum_size=0,hf_sum_size=0;
		double total_ave_size=0,hf_ave_size=0;
		double performance_result=0;//total performance function 
		ofstream ofs(write_filename,ios::app|ios::out);//append way
		if(!ofs.is_open()) cout<<"Cannot open "<<write_filename<<endl;
		ofs.setf(ios::fixed);
		ofs.precision(4);
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			NodeID isize_s = spt_v_num_f[v];
			NodeID isize_t = spt_v_num_r[v];
			NodeID isize =  isize_s + isize_t;
			total_sum_size+=isize;
			if(HFinGraphIndex[v]) hf_sum_size+=isize;
			//compute the query performance function
			double ratio_s=(double)queryTime_s[v]/(double)DIVISION_FACTOR;
			double ratio_t=(double)queryTime_t[v]/(double)DIVISION_FACTOR;
			performance_result+=ratio_s*(double)isize_s+ratio_t*(double)isize_t;
		}
		total_ave_size= (double) total_sum_size/(double) numOfVertices;
		hf_ave_size= (double) hf_sum_size/(double) numOfHFpoint;
		std::cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
		std::cout<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
		std::cout<<"standard performance_result = "<<performance_result<<endl;
		ofs<<total_sum_size<<" "<<total_ave_size<<" "<<hf_sum_size<<" "<<hf_ave_size<<" "<<performance_result<<" ";
		ofs.close();
	}

	/** experiment model output for WHP
	 * @description: 
	 * @param {double} _labeling_time
	 * @param {double} _ordering_time
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void append_experiment_result(char* write_filename,double _labeling_time,double _ordering_time=0){
		long long total_sum_size=0,hf_sum_size=0;
		double total_ave_size=0,hf_ave_size=0;
		double performance_result=0;//total performance function 
		ofstream ofs(write_filename,ios::app|ios::out);//append way
		if(!ofs.is_open()) cout<<"Cannot open "<<write_filename<<endl;
		ofs.setf(ios::fixed);
		ofs.precision(4);
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			NodeID isize_s = spt_v_num_f[v];
			NodeID isize_t = spt_v_num_r[v];
			NodeID isize = spt_v_num[v];
			total_sum_size+=isize;
			if(HFinGraphIndex[v]) hf_sum_size+=isize;
			//compute the query performance function
			double ratio_s=(double)queryTime_s[v]/(double)DIVISION_FACTOR;
			double ratio_t=(double)queryTime_t[v]/(double)DIVISION_FACTOR;
			performance_result+=ratio_s*(double)isize_s+ratio_t*(double)isize_t;
		}
		total_ave_size= (double) total_sum_size/(double) numOfVertices;
		hf_ave_size= (double) hf_sum_size/(double) numOfHFpoint;
		ofs<<_ordering_time* 1e6 <<" "<<_labeling_time* 1e6 <<" "<<total_sum_size<<" "<<total_ave_size<<" "<<hf_sum_size<<" "<<hf_ave_size<<" "<<performance_result<<" ";
		ofs.close();
		return;
	}

 	double get_query_cost(char* load_filename){
		//initialize iterms
		vector<int> queryTime_s_tmp(numOfVertices,0);
		vector<int> queryTime_t_tmp(numOfVertices,0);
		double performance_result=0;
		std::ifstream in(load_filename);//input HFPoint file to ifstream
		if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
		int s,t,qt,i=0;
		char line[24];
		//read each line representing HFpoint to vector 
		while (in.getline(line,sizeof(line))){
			std::stringstream hp(line);
			hp>>s>>t>>qt;
			queryTime_s_tmp[s]+=qt;
			queryTime_t_tmp[t]+=qt;
		}
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			NodeID isize_s = spt_v_num_f[v];
			NodeID isize_t = spt_v_num_r[v];
			//compute the query performance function
			double ratio_s=(double)queryTime_s_tmp[v]/(double)DIVISION_FACTOR;
			double ratio_t=(double)queryTime_t_tmp[v]/(double)DIVISION_FACTOR;
			performance_result+=ratio_s*(double)isize_s+ratio_t*(double)isize_t;
		}
		//std::cout<<"get_query_cost successfully!"<<std::endl;
		return performance_result;
	}

	/**
	 * @description: load each vertice's label size 
	 * @param {*}
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void load_label_size(char* load_filename){
		// spt_v_num.resize(numOfVertices,0);
		// spt_v_num_f.resize(numOfVertices,0);
		// spt_v_num_r.resize(numOfVertices,0);
		ifstream in(load_filename);//input HFPoint file to ifstream
		if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
		NodeID isize_f=0,isize_r=0; //each line representing the label size of each vertice
		//read each line representing HFpoint to vector 
		while (in>>isize_f>>isize_r)
		{
			spt_v_num_f.push_back(isize_f);
			spt_v_num_r.push_back(isize_r);
			spt_v_num.push_back(isize_f+isize_r);
		}
		in.close();
		if(numOfVertices==0){
			numOfVertices=spt_v_num.size();
			cout<<"load_label_size: numOfVertices="<<spt_v_num.size()<<endl;
		}
		return;
	}

};


#endif
