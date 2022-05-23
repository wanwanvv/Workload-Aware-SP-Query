/*
 * @Descripttion:  define the label class and its util functions
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-10 21:13:18
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-09-13 10:14:14
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
#include <malloc.h>
#include <xmmintrin.h>
//typedef unsigned __int64 BPSeed;
#include <omp.h>
#include <bitset>
#include <sstream>
#include <math.h>
#ifdef _WIN32
	#include<google/sparsehash/sparseconfig.h>
#endif
    #include<google/dense_hash_map>
	#include<google/dense_hash_set>

#include "./graph.h"
#include "./paras.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define INF_WEIGHT SP_Constants::INF_WEIGHT
#define AVG_LABEL_SIZE 200

typedef google::dense_hash_map<NodeID, EdgeWeight> dis_map;//store the cache distance

struct index_t1 {
	vector<NodeID> spt_v;
	vector<EdgeWeight> spt_d;

	NodeID size() {
		return spt_v.size();
	}
};

struct index_t1_p {
	NodeID* spt_v;
	EdgeWeight* spt_d;
}__attribute__((aligned(64)));  // Aligned for cache lines;

struct two_index_t1_p {
	NodeID* spt_v;
	EdgeWeight* spt_d;
	uint8_t* spt_lv;
	EdgeWeight* spt_ld;
}__attribute__((aligned(64)));  // Aligned for cache lines;

//label structure for cache defined by wanjingyi
struct index_hash_t_p {
	NodeID* spt_v; 
	EdgeWeight* spt_d;
	dis_map* dis_hash_m;
}__attribute__((aligned(64)));  // Aligned for cache lines;

//label structure for parallel defined by wanjingyi
// struct index_parallel_t_p {
// 	NodeID* spt_v_l;
// 	EdgeWeight* spt_d_l;
// 	NodeID* spt_v_h;
// 	EdgeWeight* spt_d_h;
// }__attribute__((aligned(64)));  // Aligned for cache lines;
struct index_parallel_t_p {
	NodeID* spt_v_l;
	EdgeWeight* spt_d_l;
	NodeID* spt_v_h;
	EdgeWeight* spt_d_h;
};  // Aligned for cache lines;

struct query_info {
	NodeID meet_node;
	NodeID search_len;
	double time_cost;
	EdgeWeight distance;
};

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
	vector<index_t1> index_;	
	index_t1_p* index_p;

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
		EdgeWeight distance = INF_WEIGHT;

		const index_t1_p &idx_s = index_p[s];
		const index_t1_p &idx_t = index_p[t];

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
	
	EdgeWeight query_p_with_nums(NodeID s, NodeID t, int k) {
		EdgeWeight distance = INF_WEIGHT;

		const index_t1_p &idx_s = index_p[s];
		const index_t1_p &idx_t = index_p[t];

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
			const index_t1_p &idx_s = index_p[i];
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
		if(!ofs.is_open()) cout<<"Cannot open"<<save_filename<<endl;
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

		index_p = (index_t1_p*)memalign(64, numOfVertices * sizeof(index_t1_p));
		


		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t1_p &idx = index_p[v];
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

		index_p = (index_t1_p*)memalign(64, numOfVertices * sizeof(index_t1_p));



		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t1_p &idx = index_p[v];
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

	void save_label_size(char* label_size_file){
		string labelSizefile(label_size_file);
		labelSizefile.append(".size");
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
	void save_label_size(char* label_size_file,const vector<NodeID>& inv) {
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

	void write_labels(char* write_filename){
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

/***
 * class HFLabel inherented from Label 
 * written by wanjingyi
 * **/
class HFLabel : public Label{
	public:
		vector<bool> HFinGraphIndex;//HFinGraphIndex is index of high frequency point in original graph, 0 represents low frequency point，1 represents high frequency point
		vector<NodeID> HFOripoint; //store the high frequency points'ids
		vector<int> queryTime;//store the query time of each point by id
		vector< pair<unsigned int,NodeID> > queryTime_pair; //store the query time of each point which is not null
		int numOfHFpoint; //number of High frequency points
		vector< vector<double> > hf_chache;//use to cache the h-h point distance
		vector< vector<int> > queryPairTime;//store the pair query frequency read from file
		vector<NodeID> HFPoint_inv; //order the node by descending freq
		vector<NodeID> HFPoint_rank; //fetch the rank by node index
		index_hash_t_p* index_hash_p;//cache_laebl
		vector<dis_map> dis_maps;//vecor of dis hash table
		index_parallel_t_p* index_parallel_p;//parallel_laebl
		int cache_size;
		bool* isDeleted;
		NodeID numOfOverlayVertices=0;
		NodeID numOfOriginalVertices=0;
		NodeID numOfDeletedVertices=0;

		void load_labels_overlay_hash(const char* load_filename){
			ifstream ifs(load_filename);
			if(!ifs.is_open()) cout<<"Cannot open "<<load_filename<<endl;
			NodeID isize = 0,i,dis_map_index=0;
			ifs.read((char*)&isize, sizeof(isize));
			numOfVertices = isize;
			cout<<"numOfVertices = "<<numOfVertices<<" numOfDeletedVertices="<<numOfDeletedVertices<<endl;
			//initialize vector dis map
			dis_maps.resize(numOfDeletedVertices);
			if(index_hash_p)
			{
				for (NodeID v = 0; v < numOfVertices; ++v) {
						free(index_hash_p[v].spt_v);
						free(index_hash_p[v].spt_d);
						free(index_hash_p[v].dis_hash_m);
					} 
				free(index_hash_p);
			}
			index_hash_p = NULL;
			index_hash_p = (index_hash_t_p*)memalign(64, numOfVertices * sizeof(index_hash_t_p));
			NodeID hub;
			EdgeWeight hub_weight;
			for (NodeID v = 0; v < numOfVertices; ++v) {
					index_hash_t_p &idx = index_hash_p[v];
					ifs.read((char*)&isize, sizeof(isize));
				if(!isDeleted[v]){
					idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
					idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
					for (i = 0; i < isize; ++i) {
						ifs.read((char*)&hub, sizeof(hub));
						ifs.read((char*)&hub_weight, sizeof(hub_weight));
						idx.spt_v[i] = hub;
						idx.spt_d[i] = hub_weight;
					}
				}else{
					NodeID lf_start,k;
					EdgeWeight lf_end;
					ifs.read((char*)&lf_start, sizeof(lf_start));
					ifs.read((char*)&lf_end, sizeof(lf_end));
					idx.spt_v = (NodeID*)memalign(64, lf_start * sizeof(NodeID));
					idx.spt_d = (EdgeWeight*)memalign(64, lf_start * sizeof(EdgeWeight));
					for (i = 1,k=0; i < lf_start; ++i,++k) {
						ifs.read((char*)&hub, sizeof(hub));
						ifs.read((char*)&hub_weight, sizeof(hub_weight));
						idx.spt_v[k] = hub;
						idx.spt_d[k] = hub_weight;
					}		
					idx.spt_v[k] = numOfVertices;
					idx.spt_d[k] = INF_WEIGHT;
					//store to dis hash map
					idx.dis_hash_m=&dis_maps[dis_map_index++];
					idx.dis_hash_m->set_empty_key(numOfVertices+1);
					for(;i<lf_end;++i){
						ifs.read((char*)&hub, sizeof(hub));
						ifs.read((char*)&hub_weight, sizeof(hub_weight));
						idx.dis_hash_m->insert(make_pair(hub,hub_weight));
					}
				}
			}
			ifs.close();
			//to be deleted
			// for(NodeID v=0;v<numOfVertices;++v){
			// 	cout<<v<<":"<<endl;
			// 	const index_hash_t_p& idx=index_hash_p[v];
			// 	for(size_t i=0;idx.spt_v[i]!=numOfVertices;++i){
			// 		cout<<"("<<idx.spt_v[i]<<","<<idx.spt_d[i]<<") ";
			// 	}
			// 	if(isDeleted[v]){
			// 		cout<<"hash_map:";
			// 		for(dis_map::iterator iter=idx.dis_hash_m->begin();iter!=idx.dis_hash_m->end();++iter) cout<<"("<<iter->first<<","<<iter->second<<") ";
			// 	}
			// 	cout<<endl;
			// }
		}

		/**
		 * @description: query distances for hierarchy labels
		 * @Author: wanjingyi
		 * @Date: 2021-02-05 09:08:38
		 * @param {NodeID} s
		 * @param {NodeID} t
		 * @return {*}
		 */
		EdgeWeight query_overlay(NodeID s, NodeID t) {
			EdgeWeight distance = INF_WEIGHT;
			
			if(!isDeleted[s]&&!isDeleted[t]){//h-h
				return query_p(s,t);
			}else if(!isDeleted[s]&&isDeleted[t]){//h-l
				EdgeWeight distance_tmp = INF_WEIGHT;
				const index_t1_p &idx_l = index_p[t];
				_mm_prefetch(&idx_l.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_l.spt_d[0], _MM_HINT_T0);
				for(int k=1;k<idx_l.spt_v[0];++k){
					NodeID l=idx_l.spt_v[k];
					EdgeWeight l_w=idx_l.spt_d[k];
					distance_tmp=query_p(l,s);
					if(distance_tmp+l_w<distance) distance=distance_tmp+l_w;
				}
				return distance;
			}else if(isDeleted[s]&&!isDeleted[t]){//l-h
				EdgeWeight distance_tmp = INF_WEIGHT;
				const index_t1_p &idx_l = index_p[s];
				_mm_prefetch(&idx_l.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_l.spt_d[0], _MM_HINT_T0);
				for(int k=1;k<idx_l.spt_v[0];++k){		
					NodeID l=idx_l.spt_v[k];
					EdgeWeight l_w=idx_l.spt_d[k];
					distance_tmp=query_p(l,t);
					if(distance_tmp+l_w<distance) distance=distance_tmp+l_w;
				}
				return distance;
			}else{//l-l
				EdgeWeight distance_tmp = INF_WEIGHT;
				const index_t1_p &idx_s = index_p[s];
				const index_t1_p &idx_t = index_p[t];
				_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
				//lf segmentation traverse
				int end=idx_s.spt_d[0];
				for(int k=idx_s.spt_v[0];k<idx_s.spt_d[0];++k){
					if(idx_s.spt_v[k]==t){
						distance=idx_s.spt_d[k];
						break;
						//return idx_s.spt_d[k];
					}
				}
				for(int k=idx_t.spt_v[0];k<idx_t.spt_d[0];++k){
					if(idx_t.spt_v[k]==s){
						distance=min(distance,idx_t.spt_d[k]);
						break;
						//return idx_t.spt_d[k];
					}
				}
				if(distance!=INF_WEIGHT) return distance;
				//hf segmentation traverse
				for(int i=1;i<idx_s.spt_v[0];++i){
					NodeID s_l=idx_s.spt_v[i];
					EdgeWeight s_w=idx_s.spt_d[i];
					for(int j=1;j<idx_t.spt_v[0];++j){
						NodeID t_l=idx_t.spt_v[j];
						EdgeWeight t_w=idx_t.spt_d[j];
						distance_tmp=query_p(s_l,t_l);
						distance_tmp=distance_tmp+s_w+t_w;
						if(distance_tmp<distance) distance=distance_tmp;
					}
				}
				return distance;
			}
			return distance;
		}

		/**
		 * @Descripttion: 
		 * @version: 
		 * @Author: Wan Jingyi
		 * @Date: 2021-04-12 11:39:44
		 * @LastEditors: Wan Jingyi
		 * @LastEditTime: Do not Edit
		 * @param {NodeID} s
		 * @param {NodeID} t
		 */  
  		EdgeWeight query_p_hash(NodeID s, NodeID t) {
			EdgeWeight distance = INF_WEIGHT;

			const index_hash_t_p &idx_s = index_hash_p[s];
			const index_hash_t_p &idx_t = index_hash_p[t];

			// _mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
			// _mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
			// _mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
			// _mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

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

		EdgeWeight query_overlay_hash(NodeID s, NodeID t) {
			EdgeWeight distance = INF_WEIGHT;
			if(!isDeleted[s]&&!isDeleted[t]){//h-h
				return query_p_hash(s,t);
			}else if(!isDeleted[s]&&isDeleted[t]){//h-l
				EdgeWeight distance_tmp = INF_WEIGHT;
				const index_hash_t_p &idx_l = index_hash_p[t];
				_mm_prefetch(&idx_l.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_l.spt_d[0], _MM_HINT_T0);
				for(int k=0;idx_l.spt_v[k]!=numOfVertices;++k){
					NodeID l=idx_l.spt_v[k];
					EdgeWeight l_w=idx_l.spt_d[k];
					distance_tmp=query_p_hash(l,s);
					if(distance_tmp+l_w<distance) distance=distance_tmp+l_w;
				}
				return distance;
			}else if(isDeleted[s]&&!isDeleted[t]){//l-h
				EdgeWeight distance_tmp = INF_WEIGHT;
				const index_hash_t_p &idx_l = index_hash_p[s];
				_mm_prefetch(&idx_l.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_l.spt_d[0], _MM_HINT_T0);
				for(int k=0;idx_l.spt_v[k]!=numOfVertices;++k){		
					NodeID l=idx_l.spt_v[k];
					EdgeWeight l_w=idx_l.spt_d[k];
					distance_tmp=query_p_hash(l,t);
					if(distance_tmp+l_w<distance) distance=distance_tmp+l_w;
				}
				return distance;
			}else{//l-l
				EdgeWeight distance_tmp = INF_WEIGHT;
				const index_hash_t_p &idx_s = index_hash_p[s];
				const index_hash_t_p &idx_t = index_hash_p[t];
				_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
				_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
				_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
				//lf hash traverse
				if(idx_s.dis_hash_m->find(t)!=idx_s.dis_hash_m->end()){
					return (*idx_s.dis_hash_m)[t];
				}else if(idx_t.dis_hash_m->find(s)!=idx_t.dis_hash_m->end()){
					return (*idx_t.dis_hash_m)[s];
				}else{
					for(int i=0;idx_s.spt_v[i]!=numOfVertices;++i){
						NodeID s_l=idx_s.spt_v[i];
						EdgeWeight s_w=idx_s.spt_d[i];
						for(int j=0;idx_t.spt_v[j]!=numOfVertices;++j){
							NodeID t_l=idx_t.spt_v[j];
							EdgeWeight t_w=idx_t.spt_d[j];
							distance_tmp=query_p_hash(s_l,t_l);
							distance_tmp=distance_tmp+s_w+t_w;
							if(distance_tmp<distance) distance=distance_tmp;
						}
					}
					return distance;
				}
			}
			return distance;
		}

		/**
	 * @description: load whether the vertex is deleted
	 * @Author: wanjingyi
	 * @Date: 2021-02-03 15:58:03
	 * @param {*}
	 * @return {*}
	 */
		void load_is_deleted(char* isDeletedFIleName){
			ifstream in(isDeletedFIleName);
			if(!in.is_open()) cerr<<isDeletedFIleName<<" cannot be opened!"<<endl;
			in>>numOfOriginalVertices>>numOfOverlayVertices;
			numOfDeletedVertices=numOfOriginalVertices-numOfOverlayVertices;
			isDeleted=new bool[numOfOriginalVertices];
			int i,u,flag;
			//build the map relationship 
			for (i = 0; in >> u >> flag;i++) {
				isDeleted[u]=flag;
			}
			if(i!=numOfOriginalVertices) cerr<<"error:i!=numOfOriginalVertices!"<<endl;	
			in.close();
		}

		void save_analysis_overlay(char* write_filename,const vector<unsigned int>& _queryTime,const vector<bool>& _HFinGraphIndex,int _max_query_time,int _numOfHFpoint){
			long long total_sum_size=0,hf_sum_size=0;
			double total_ave_size=0,hf_ave_size=0;
			double performance_result=0;//total performance function 
			NodeID isize;
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				isize = index_[v].size()-1;
				total_sum_size+=isize;
				if(_HFinGraphIndex[v]) hf_sum_size+=isize;
				//compute the query performance function
				double ratio=(double)_queryTime[v]/((double)(_max_query_time*DIVISION_FACTOR));
				performance_result+=ratio*(double)isize;
			}
			total_ave_size= (double) total_sum_size/(double) numOfVertices;
			hf_ave_size= (double) hf_sum_size/(double) _numOfHFpoint;
			cout<<"total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			cout<<"hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			//append way
			ofstream ofs(write_filename,ios::app|ios::out);
			if(!ofs.is_open()) cout<<"Cannot open "<<write_filename<<endl;
			ofs.setf(ios::fixed);
			ofs.precision(4);
			performance_result=performance_result*((double)_max_query_time);
			ofs<<performance_result<<" ";
			ofs<<total_sum_size<<" "<<total_ave_size<<" ";
			ofs<<hf_sum_size<<" "<<hf_ave_size<<" ";
			ofs.close();
			cout<<"***************************save_anaylysis_overlay finished!*****************************"<<endl;
		}

		/**
		  * function used to output the analysis information
		  * hfpoint by assigned vector
		  * written by wanjingyi
		  * */
		void save_anaylysis_size(char* write_filename,vector<unsigned int>& _queryTime,vector<bool>& _HFinGraphIndex,int max_query_time,long long total_query_time,int _numOfHfpoint){
			cout<<"***************************save_anaylysis_size start!*****************************"<<endl;
			numOfHFpoint=_numOfHfpoint;
			long long total_sum_size=0,hf_sum_size=0;
			double total_ave_size=0,hf_ave_size=0;
			double performance_result=0;//total performance function 
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				NodeID isize = index_[v].size()-1;
				total_sum_size+=isize;
				if(_HFinGraphIndex[v]) hf_sum_size+=isize;
				//compute the query performance function
				double ratio=(double)_queryTime[v]/((double)(max_query_time*DIVISION_FACTOR));
				performance_result+=ratio*(double)isize;
			}
			total_ave_size= (double) total_sum_size/(double) numOfVertices;
			hf_ave_size= (double) hf_sum_size/(double) numOfHFpoint;
			cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			cout<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			string output_asizeFileName(write_filename);
        	output_asizeFileName.append(".analysis");
			ofstream ofs(output_asizeFileName);
			if(!ofs.is_open()) {cerr<<"Cannot open "<<output_asizeFileName<<endl;}
			ofs<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			ofs<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			//write performance function to file
			std::cout<<"nomalization performance_result = "<<performance_result<<endl;
			ofs<<"nomalization performance_result = "<<performance_result<<endl;
			performance_result=performance_result*((double)max_query_time);
			std::cout<<"standard performance_result = "<<performance_result<<endl;
			ofs<<"standard performance_result = "<<performance_result<<endl;
			ofs.close();
			cout<<"***************************save_anaylysis_size finished!*****************************"<<endl;
		}

		void append_experiment_overlay(char* write_filename,vector<unsigned int>& _queryTime,vector<bool>& _HFinGraphIndex,int max_query_time,long long total_query_time,int _numOfHfpoint,double _ordering_time,double _labeling_time){
			numOfHFpoint=_numOfHfpoint;
			long long total_sum_size=0,hf_sum_size=0;
			double total_ave_size=0,hf_ave_size=0;
			double performance_result=0;//total performance function 
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				NodeID isize = index_[v].size()-1;
				total_sum_size+=isize;
				if(_HFinGraphIndex[v]) hf_sum_size+=isize;
				//compute the query performance function
				double ratio=(double)_queryTime[v]/((double)(max_query_time*DIVISION_FACTOR));
				performance_result+=ratio*(double)isize;
			}
			//compute label_num
			total_ave_size= (double) total_sum_size/(double) numOfVertices;
			hf_ave_size= (double) hf_sum_size/(double) numOfHFpoint;
			cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
			cout<<"numOfHFpoint = "<<numOfHFpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
			//compute query cost
			std::cout<<"nomalization performance_result = "<<performance_result<<endl;
			performance_result=performance_result*((double)max_query_time);
			std::cout<<"standard performance_result = "<<performance_result<<endl;
			//append to file
			ofstream ofs(write_filename,ios::app|ios::out);//append way
			if(!ofs.is_open()) cout<<"Cannot open "<<write_filename<<endl;
			ofs.setf(ios::fixed);
			ofs.precision(4);
			ofs<<_ordering_time* 1e6 <<" "<<_labeling_time* 1e6<<" "<<hf_sum_size<<" "<<hf_ave_size<<" "<<total_sum_size<<" "<<hf_ave_size<<" "<<performance_result<<" ";
			ofs.close();
		}

		/**
		  * function used to output the analysis information
		  * hfpoint by computation
		  * written by wanjingyi
		  * */
		 void save_anaylysis_size(char* freqid_filename,char* write_filename,int hfRate)
		 {
			 //load hfPoint
			numOfHFpoint = 0;//first line is the number of HFpoints
			numOfHFpoint = static_cast<int> ( (double)numOfVertices*hfRate/(double)1000);
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

};

class DLabel : public Label {
	public:
	vector<index_t1> bindex_; // Backward labels.

	index_t1_p* bindex_p;
	
	
	two_index_t1_p* b_two_index_p;


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
		////const index_t1_p &idx_s = index_p[s];
		////const index_t1_p &idx_t = bindex_p[t];

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

		const index_t1_p &idx_s = index_p[s];
		const index_t1_p &idx_t = bindex_p[t];


		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];

			if (v1 == v2) {
				if (v1 == numOfVertices) break;  //two_index_t1_p Sentinel
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

		index_p = (index_t1_p*)memalign(64, numOfVertices * sizeof(index_t1_p));
		bindex_p = (index_t1_p*)memalign(64, numOfVertices * sizeof(index_t1_p));

		cout << numOfVertices << " vertices." << endl;

		for (NodeID v = 0; v < numOfVertices; ++v) {
			index_t1_p &idx = index_p[v];
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


			index_t1_p &bidx = bindex_p[v];
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

		/**
		 * @description: load each vertice's label size 
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */ 
		void load_label_size(char* load_filename){
			spt_v_num.resize(numOfVertices,0);
			spt_v_num_f.resize(numOfVertices,0);
			spt_v_num_r.resize(numOfVertices,0);
			ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
			NodeID isize_f=0,isize_r=0; //each line representing the label size of each vertice
			//read each line representing HFpoint to vector 
			char line[24];
			int i=0;
			while (in.getline(line,sizeof(line)))
			{
				stringstream ls(line);
				ls>>isize_f>>isize_r;
				spt_v_num_f[i]=isize_f;
				spt_v_num_r[i]=isize_r;
				spt_v_num[i++]=isize_f+isize_r;
			}
			in.close();
			if(i!=numOfVertices) cout<<"i!=numOfVertices"<<endl;
			cout<<"load_label_size: numOfVertices="<<spt_v_num.size()<<endl;
			return;
		}

	};



#endif