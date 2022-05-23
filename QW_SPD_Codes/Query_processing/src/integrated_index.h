/*
 * @Descripttion: index structure to store index for integrated solutions
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-04-13 17:06:36
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-02-26 10:30:38
 */
#ifndef _INTEGRATED_INDEX_H
#define _INTEGRATED_INDEX_H

#include <iostream>
#include <vector>
#include <bitset>
#include <limits>
#include <climits>
#include <malloc.h>
#include <sys/time.h>
#include <xmmintrin.h> 
#include "./paras.h"
#include "./time_util.h"
#include "./utils.h"
#include "./graph.h"
#include "./h2h_index.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define INF_WEIGHT SP_Constants::INF_WEIGHT

using namespace std;
using namespace time_util;

struct integrated_index_t{
	vector<NodeID> spt_v;
	vector<EdgeWeight> spt_d;
    vector<NodeID> pos;
    vector<EdgeWeight> dis;
    NodeID level=0;
    double siz=0;
};

struct integrated_index_t_p{
	NodeID* spt_v=nullptr;
	EdgeWeight* spt_d=nullptr;
    EdgeWeight* dis=nullptr;
	NodeID* pos=nullptr;
}__attribute__((aligned(64)));

//Add for directed
struct integrated_bindex_t{
	vector<NodeID> spt_v;
	vector<EdgeWeight> spt_d;
	vector<NodeID> r_spt_v;
	vector<EdgeWeight> r_spt_d;
    vector<NodeID> pos;
    vector<EdgeWeight> dis;
    vector<EdgeWeight> r_dis;
    NodeID level=0;
    double siz=0; //average label size as start
    double r_siz=0;//average label size as termination
};

//Add for directed
struct integrated_bindex_t_p{
	NodeID* spt_v=nullptr;
	EdgeWeight* spt_d=nullptr;
	NodeID* r_spt_v=nullptr;
	EdgeWeight* r_spt_d=nullptr;
    EdgeWeight* dis=nullptr;
    EdgeWeight* r_dis=nullptr;
	NodeID* pos=nullptr;
}__attribute__((aligned(64)));

/**
 * @description: base class for integrated index
 * @author: Wan Jingyi
 */
class Integrated_Index{
public:
    //***********variables*************
    vector<NodeID> inv;
    vector<NodeID> rank;
    vector<integrated_index_t> index_;
    vector<integrated_bindex_t> bindex_;//Add for directed

    vector<double> spt_v_num; //store the size of each vertice's labels
    vector<double> spt_v_num_f; //store the size of each vertice's forward labels
    vector<double> spt_v_num_r; //store the size of each vertice's reverse labels

    bool* isDeleted=nullptr;
    integrated_index_t_p* index_p=nullptr; 
    integrated_bindex_t_p* bindex_p=nullptr; //Add for directed
    LCA lca;//find the least common ancestor structure

    int _numOfOriginalVertices=0;
    int _numOfOverlayVertices=0;

    NodeID _max_overlay_hub{0};//the max overlay hub rank value
    NodeID _max_overlay_hub_r{0};//the max overlay hub rank value, add for directed
    //**********constructions and deconstructions*********
    Integrated_Index(){
        index_p=nullptr;
        bindex_p=nullptr;
    }

    ~Integrated_Index(){
        //clear index vector
        if(index_.size()!=0){
            for(int v=0;v<numOfVertices;++v) 
            index_.clear();
        }
    }

    //******************functions****************
	/**
  * @description: pll merge sort used for intergrated_query stage
  * @param {const integrated_index_t_p} &idx_s
  * @param {const integrated_index_t_p} &idx_t
  * @return {*}
  * @author: Wan Jingyi
  */ 
    EdgeWeight query_p(const integrated_index_t_p &idx_s, const integrated_index_t_p &idx_t) {
		EdgeWeight distance = INF_WEIGHT;
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

	/**
     * @description: pll merge sort used for integrated-index construction stage
     * @param {NodeID} s
     * @param {NodeID} t
     * @return {*}
     * @author: Wan Jingyi
     */ 
    EdgeWeight query_p(NodeID s, NodeID t) {
		EdgeWeight distance = INF_WEIGHT;
        const integrated_index_t &idx_s = index_[s];
		const integrated_index_t &idx_t = index_[t];
		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0;  i<idx_s.spt_v.size()-1&&j<idx_t.spt_v.size()-1;) {
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

    
    //Add fordirected
	/**
  * @description: pll merge sort used for intergrated_query stage
  * @param {const integrated_index_t_p} &idx_s
  * @param {const integrated_index_t_p} &idx_t
  * @return {*}
  * @author: Wan Jingyi
  */ 
    EdgeWeight query_p_directed(const integrated_bindex_t_p &idx_s, const integrated_bindex_t_p &idx_t) {
		EdgeWeight distance = INF_WEIGHT;
		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.r_spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.r_spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0; ; ) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.r_spt_v[j];

			if (v1 == numOfVertices) break;  // Sentinel

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.r_spt_d[j];
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
     * @description: pll merge sort used for integrated-index construction stage
     * @param {NodeID} s
     * @param {NodeID} t
     * @return {*}
     * @author: Wan Jingyi
     */ 
    EdgeWeight query_p_directed(NodeID s, NodeID t) {
		EdgeWeight distance = INF_WEIGHT;
        const integrated_bindex_t &idx_s = bindex_[s];
		const integrated_bindex_t &idx_t = bindex_[t];
		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.r_spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t.r_spt_d[0], _MM_HINT_T0);

		for (int i = 0, j = 0;  i<idx_s.spt_v.size()-1&&j<idx_t.r_spt_v.size()-1;) {
			NodeID v1 = idx_s.spt_v[i], v2 = idx_t.r_spt_v[j];

			if (v1 == v2) {
				EdgeWeight td = idx_s.spt_d[i] + idx_t.r_spt_d[j];
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
     * @description: used to l-l query for integrated-query
     * @param {const} vector
     * @param {const} vector
     * @return {*}
     * @author: Wan Jingyi
     */ 
    EdgeWeight query_p(const vector<EdgeWeight>& idx_s_vec, const vector<EdgeWeight>& idx_t_vec) {
		EdgeWeight distance = INF_WEIGHT;
		_mm_prefetch(&idx_s_vec[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t_vec[0], _MM_HINT_T0);
        size_t end_index=idx_s_vec.size();
		for (size_t i = 0; i<end_index;++i ) {
            if(idx_s_vec[i]+idx_t_vec[i]<distance) distance=idx_s_vec[i]+idx_t_vec[i];
        }
		return distance;
	}

	/**
     * @description: used to l-l query for integrated-query, add for directed
     * @param {const} vector
     * @param {const} vector
     * @return {*}
     * @author: Wan Jingyi
     */ 
    EdgeWeight query_p_directed(const vector<EdgeWeight>& idx_s_vec, const vector<EdgeWeight>& idx_t_vec) {
		EdgeWeight distance = INF_WEIGHT;
		_mm_prefetch(&idx_s_vec[0], _MM_HINT_T0);
		_mm_prefetch(&idx_t_vec[0], _MM_HINT_T0);
        size_t end_index=idx_s_vec.size();
		for (size_t i = 0; i<end_index;++i ) {
            if(idx_s_vec[i]+idx_t_vec[i]<distance) distance=idx_s_vec[i]+idx_t_vec[i];
        }
		return distance;
	}

	/**
     * @description: 
     * @param {size_t} end_index
     * @return {*}
     * @author: Wan Jingyi
     */ 
    EdgeWeight query_p(EdgeWeight* idx_s_vec, EdgeWeight* idx_t_vec,size_t end_index) {
		EdgeWeight distance = INF_WEIGHT;
		_mm_prefetch(idx_s_vec, _MM_HINT_T0);
		_mm_prefetch(idx_t_vec, _MM_HINT_T0);

		for (size_t i = 0; i<end_index;++i ) {
            if(idx_s_vec[i]+idx_t_vec[i]<distance) distance=idx_s_vec[i]+idx_t_vec[i];
        }
		return distance;
	}
    //*****************************unused********************************//

    EdgeWeight query_integrated(NodeID s,NodeID t){
        EdgeWeight distance = INF_WEIGHT,distance_t = INF_WEIGHT;

        const integrated_index_t_p &idx_s = index_p[s];
		const integrated_index_t_p &idx_t = index_p[t];

        if(!isDeleted[s]&&!isDeleted[t]){//h-h
            return query_p(idx_s,idx_t);
        }else{//l-h h-l l-l
            int lca_id=lca.lca_query(s,t);
            //cout<<"("<<s<<","<<t<<") "<<"lca="<<lca_id;//to be deleted
            if(lca_id==s){
                //cout<<" 0"<<endl;
                return idx_t.dis[idx_s.pos[1]];
            }else if(lca_id==t){
                //cout<<" 1"<<endl;
                return idx_s.dis[idx_t.pos[1]];
            }else{
                //cout<<" 2:";
                if(!isDeleted[s]&&isDeleted[t]){
                    //cout<<"sh-tl"<<endl;
                    const integrated_index_t_p &idx_lca = index_p[lca_id];
                    _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int end_len=idx_lca.pos[0];
                    for(int i=0;i<end_len;++i){
                         const integrated_index_t_p &idx_lca_t=index_p[idx_lca.pos[end_len+1+i]];
                        distance_t=query_p(idx_s,idx_lca_t)+idx_t.dis[idx_lca.pos[1+i]];
                        if(distance_t<distance) distance=distance_t;
                    }
                }else if(isDeleted[s]&&!isDeleted[t]){
                    //cout<<"sl-th"<<endl;
                    const integrated_index_t_p &idx_lca = index_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int end_len=idx_lca.pos[0];
                    for(int i=0;i<end_len;++i){
                        const integrated_index_t_p &idx_lca_t=index_p[idx_lca.pos[end_len+1+i]];
                        distance_t=idx_s.dis[idx_lca.pos[1+i]]+query_p(idx_lca_t,idx_t);
                        if(distance_t<distance) distance=distance_t;
                    }
                }else{
                    //cout<<"sl-tl"<<endl;
                    const integrated_index_t_p &idx_lca = index_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int end_len=idx_lca.pos[0];
                    for(int i=1;i<=end_len;++i){
                        int pos_i=idx_lca.pos[i];
                        distance_t=idx_s.dis[pos_i]+idx_t.dis[pos_i];
                        if(distance_t<distance) distance=distance_t;
                    }
                    return distance;
                }

            }
        }
        return distance;
    }

    EdgeWeight query_integrated_1(NodeID s,NodeID t){
        EdgeWeight distance = INF_WEIGHT,distance_t = INF_WEIGHT;
        //cout<<s<<"-"<<t<<":";//for debug
        if(lca._index[s]==-1&&lca._index[t]==-1){//h-h
            //cout<<"h-h"<<endl;//for debug
            const integrated_index_t_p &idx_s = index_p[s];
            const integrated_index_t_p &idx_t = index_p[t];
            return query_p(idx_s,idx_t);
        }else if(lca._index[s]==-1&&lca._index[t]!=-1){//h-l
            //cout<<"h-l ";//for debug
            const integrated_index_t_p &idx_t = index_p[t];
            const integrated_index_t_p &idx_root = index_p[idx_t.dis[0]];
            _mm_prefetch(&idx_root.spt_v[0], _MM_HINT_T0);
            //judge borders
            NodeID border_id;
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                if(border_id==s){
                    //cout<<"hash"<<endl;//for debug
                    return idx_t.dis[i+1];
                }           
            }
            //cout<<"pll"<<endl;//for debug
            //pll 2-hop
            const integrated_index_t_p &idx_s = index_p[s];
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                const integrated_index_t_p &idx_border = index_p[border_id];
                distance_t=query_p(idx_s,idx_border)+idx_t.dis[i+1];
                if(distance_t<distance) distance=distance_t;
            }
            return distance;
        }else if(lca._index[s]!=-1&&lca._index[t]==-1){//l-h
            //cout<<"l-h ";//for debug
            const integrated_index_t_p &idx_s = index_p[s];
            const integrated_index_t_p &idx_root = index_p[idx_s.dis[0]];
            _mm_prefetch(&idx_root.spt_v[0], _MM_HINT_T0);
            //judge borders
            NodeID border_id;
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                if(border_id==t){
                    //cout<<"hash"<<endl;//for debug
                    return idx_s.dis[i+1];
                }           
            }
            //cout<<"pll"<<endl;//for debug
            //pll 2-hop
            const integrated_index_t_p &idx_t = index_p[t];
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                const integrated_index_t_p &idx_border = index_p[border_id];
                distance_t=query_p(idx_border,idx_t)+idx_s.dis[i+1];
                if(distance_t<distance) distance=distance_t;
            }
            return distance;
        }else{//l-l
            //cout<<"l-l ";//for debug
            const integrated_index_t_p &idx_s = index_p[s];
            const integrated_index_t_p &idx_t = index_p[t];
            if(idx_s.dis[0]==idx_t.dis[0]){//in the same tree
               //cout<<"h2h ";//for debug
                int lca_id=lca.lca_query(s,t);
                //cout<<"lca="<<lca_id;//for debug
                if(lca_id==s){
                    //cout<<" 0"<<endl;//for debug
                   // _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    return idx_t.dis[idx_s.pos[0]+1];
                }else if(lca_id==t){
                    //cout<<" 1"<<endl;//for debug
                    //_mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    return idx_s.dis[idx_t.pos[0]+1];
                }else{
                    //cout<<" 2"<<endl;//for debug
                    const integrated_index_t_p &idx_lca = index_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int pos_i;
                    for(int i=0;(pos_i=idx_lca.pos[i])!=numOfVertices;++i){
                        distance_t=idx_s.dis[pos_i+1]+idx_t.dis[pos_i+1];
                        if(distance_t<distance) distance=distance_t;
                    }
                    return distance;
                }
            }else{//cover different trees
               //cout<<s<<"-"<<t<<" pll:"<<endl;//for debug
               double time_scan_add,time_scan=0;//for debug
               //variables 
               double time_variable=GetCurrentTimeSec();
                NodeID border_id;
                NodeID spt_v_t;
                EdgeWeight spt_d_t;
                int i,j;
                // EdgeWeight* idx_s_vec=new EdgeWeight[_max_overlay_hub+1];
                // EdgeWeight* idx_t_vec=new EdgeWeight[_max_overlay_hub+1];
                // for(size_t i=0;i<_max_overlay_hub+1;++i) idx_s_vec[i]=idx_t_vec[i]=INF_WEIGHT;
                
                //first scan the s node region
                const integrated_index_t_p &idx_root_s = index_p[idx_s.dis[0]];
                _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                _mm_prefetch(&idx_root_s.spt_v[0], _MM_HINT_T0);
                vector<EdgeWeight> idx_s_vec(_max_overlay_hub+1,INF_WEIGHT);
                for(i=0;(border_id=idx_root_s.spt_v[i])!=numOfVertices;++i){
                    const integrated_index_t_p &idx_border = index_p[border_id];
                    for(j=0;(spt_v_t=idx_border.spt_v[j])!=numOfVertices;++j){
                        spt_v_t=idx_border.spt_v[j];
                        spt_d_t=idx_border.spt_d[j]+idx_s.dis[i+1];
                        if(spt_d_t<idx_s_vec[spt_v_t]) idx_s_vec[spt_v_t]=spt_d_t;
                    }
                }
                //second scan the s node region
                const integrated_index_t_p &idx_root_t = index_p[idx_t.dis[0]]; 
                _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                _mm_prefetch(&idx_root_t.spt_v[0], _MM_HINT_T0);
                vector<EdgeWeight> idx_t_vec(_max_overlay_hub+1,INF_WEIGHT);
                for(i=0;(border_id=idx_root_t.spt_v[i])!=numOfVertices;++i){
                    const integrated_index_t_p &idx_border = index_p[border_id];
                    for(j=0;(spt_v_t=idx_border.spt_v[j])!=numOfVertices;++j){
                        spt_v_t=idx_border.spt_v[j];
                        spt_d_t=idx_border.spt_d[j]+idx_t.dis[i+1];
                        if(spt_d_t<idx_t_vec[spt_v_t]) idx_t_vec[spt_v_t]=spt_d_t;
                    }
                }
                //merge scane
                return query_p(idx_s_vec,idx_t_vec);
                //return query_p(idx_s_vec,idx_t_vec,_max_overlay_hub+1);
                // for(int i=0;(border_id_s=idx_root_s.spt_v[i])!=numOfVertices;++i){
                //     const integrated_index_t_p &idx_border_s = index_p[border_id_s];
                //     for(int j=0;(border_id_t=idx_root_t.spt_v[j])!=numOfVertices;++j){
                //         const integrated_index_t_p &idx_border_t = index_p[border_id_t];
                //         distance_t=query_p(idx_border_s,idx_border_t)+idx_s.dis[i+1]+idx_t.dis[j+1];
                //         if(distance_t<distance) distance=distance_t;
                //     }
                // }
                //return distance;
            }
        }
        return distance;
    }

    EdgeWeight query_integrated_1(NodeID s,NodeID t,int& type){
        EdgeWeight distance = INF_WEIGHT,distance_t = INF_WEIGHT;
        //cout<<s<<"-"<<t<<":";//for debug
        if(lca._index[s]==-1&&lca._index[t]==-1){//h-h
            //cout<<"h-h"<<endl;//for debug
            const integrated_index_t_p &idx_s = index_p[s];
            const integrated_index_t_p &idx_t = index_p[t];
            type=0;
            return query_p(idx_s,idx_t);
        }else if(lca._index[s]==-1&&lca._index[t]!=-1){//h-l
            //cout<<"h-l ";//for debug
            const integrated_index_t_p &idx_t = index_p[t];
            const integrated_index_t_p &idx_root = index_p[idx_t.dis[0]];
            _mm_prefetch(&idx_root.spt_v[0], _MM_HINT_T0);
            //judge borders
            NodeID border_id;
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                if(border_id==s){
                    //cout<<"hash"<<endl;//for debug
                    type=1;
                    return idx_t.dis[i+1];
                }           
            }
            //cout<<"pll"<<endl;//for debug
            //pll 2-hop
            const integrated_index_t_p &idx_s = index_p[s];
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                const integrated_index_t_p &idx_border = index_p[border_id];
                distance_t=query_p(idx_s,idx_border)+idx_t.dis[i+1];
                if(distance_t<distance) distance=distance_t;
            }
            type=2;
            return distance;
        }else if(lca._index[s]!=-1&&lca._index[t]==-1){//l-h
            //cout<<"l-h ";//for debug
            const integrated_index_t_p &idx_s = index_p[s];
            const integrated_index_t_p &idx_root = index_p[idx_s.dis[0]];
            _mm_prefetch(&idx_root.spt_v[0], _MM_HINT_T0);
            //judge borders
            NodeID border_id;
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                if(border_id==t){
                    //cout<<"hash"<<endl;//for debug
                    type=1;
                    return idx_s.dis[i+1];
                }           
            }
            //cout<<"pll"<<endl;//for debug
            //pll 2-hop
            const integrated_index_t_p &idx_t = index_p[t];
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                const integrated_index_t_p &idx_border = index_p[border_id];
                distance_t=query_p(idx_border,idx_t)+idx_s.dis[i+1];
                if(distance_t<distance) distance=distance_t;
            }
            type=2;
            return distance;
        }else{//l-l
            //cout<<"l-l ";//for debug
            const integrated_index_t_p &idx_s = index_p[s];
            const integrated_index_t_p &idx_t = index_p[t];
            if(idx_s.dis[0]==idx_t.dis[0]){//in the same tree
                type=3;
               //cout<<"h2h ";//for debug
                int lca_id=lca.lca_query(s,t);
                //cout<<"lca="<<lca_id;//for debug
                if(lca_id==s){
                    //cout<<" 0"<<endl;//for debug
                   // _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    return idx_t.dis[idx_s.pos[0]+1];
                }else if(lca_id==t){
                    //cout<<" 1"<<endl;//for debug
                    //_mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    return idx_s.dis[idx_t.pos[0]+1];
                }else{
                    //cout<<" 2"<<endl;//for debug
                    const integrated_index_t_p &idx_lca = index_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int pos_i;
                    for(int i=0;(pos_i=idx_lca.pos[i])!=numOfVertices;++i){
                        distance_t=idx_s.dis[pos_i+1]+idx_t.dis[pos_i+1];
                        if(distance_t<distance) distance=distance_t;
                    }
                    return distance;
                }
            }else{//cover different trees
                type=4;
               //cout<<s<<"-"<<t<<" pll:"<<endl;//for debug
               double time_scan_add,time_scan=0;//for debug
               //variables 
               double time_variable=GetCurrentTimeSec();
                NodeID border_id;
                NodeID spt_v_t;
                EdgeWeight spt_d_t;
                int i,j;
                //first scan the s node region
                const integrated_index_t_p &idx_root_s = index_p[idx_s.dis[0]];
                _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                _mm_prefetch(&idx_root_s.spt_v[0], _MM_HINT_T0);
                vector<EdgeWeight> idx_s_vec(_max_overlay_hub+1,INF_WEIGHT);
                for(i=0;(border_id=idx_root_s.spt_v[i])!=numOfVertices;++i){
                    const integrated_index_t_p &idx_border = index_p[border_id];
                    for(j=0;(spt_v_t=idx_border.spt_v[j])!=numOfVertices;++j){
                        spt_v_t=idx_border.spt_v[j];
                        spt_d_t=idx_border.spt_d[j]+idx_s.dis[i+1];
                        if(spt_d_t<idx_s_vec[spt_v_t]) idx_s_vec[spt_v_t]=spt_d_t;
                    }
                }
                //second scan the s node region
                const integrated_index_t_p &idx_root_t = index_p[idx_t.dis[0]]; 
                _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                _mm_prefetch(&idx_root_t.spt_v[0], _MM_HINT_T0);
                vector<EdgeWeight> idx_t_vec(_max_overlay_hub+1,INF_WEIGHT);
                for(i=0;(border_id=idx_root_t.spt_v[i])!=numOfVertices;++i){
                    const integrated_index_t_p &idx_border = index_p[border_id];
                    for(j=0;(spt_v_t=idx_border.spt_v[j])!=numOfVertices;++j){
                        spt_v_t=idx_border.spt_v[j];
                        spt_d_t=idx_border.spt_d[j]+idx_t.dis[i+1];
                        if(spt_d_t<idx_t_vec[spt_v_t]) idx_t_vec[spt_v_t]=spt_d_t;
                    }
                }
                //merge scane
                return query_p(idx_s_vec,idx_t_vec);
            }
        }
        return distance;
    }

    //Add for directed
    EdgeWeight query_integrated_1_directed(NodeID s,NodeID t){
        EdgeWeight distance = INF_WEIGHT, distance_t = INF_WEIGHT;
        //cout<<s<<"-"<<t<<":";//for debug
        if(lca._index[s]==-1&&lca._index[t]==-1){//h-h
            //cout<<"h-h"<<endl;//for debug
            const integrated_bindex_t_p &idx_s = bindex_p[s];
            const integrated_bindex_t_p &idx_t = bindex_p[t];
            return query_p_directed(idx_s,idx_t);
        }else if(lca._index[s]==-1&&lca._index[t]!=-1){//h-l
            //cout<<"h-l ";//for debug
            const integrated_bindex_t_p &idx_t = bindex_p[t];
            const integrated_bindex_t_p &idx_root = bindex_p[idx_t.pos[0]];
            _mm_prefetch(&idx_root.spt_v[0], _MM_HINT_T0);
            //judge borders
            NodeID border_id;
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                if(border_id==s){
                    //cout<<"hash"<<endl;//for debug
                    return idx_t.r_dis[i];
                }           
            }
            //cout<<"pll"<<endl;//for debug
            //pll 2-hop
            const integrated_bindex_t_p &idx_s = bindex_p[s];
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                const integrated_bindex_t_p &idx_border = bindex_p[border_id];
                distance_t=query_p_directed(idx_s,idx_border)+idx_t.r_dis[i];
                if(distance_t<distance) distance=distance_t;
            }
            return distance;
        }else if(lca._index[s]!=-1&&lca._index[t]==-1){//l-h
            //cout<<"l-h ";//for debug
            const integrated_bindex_t_p &idx_s = bindex_p[s];
            const integrated_bindex_t_p &idx_root = bindex_p[idx_s.pos[0]];
            _mm_prefetch(&idx_root.spt_v[0], _MM_HINT_T0);
            //judge borders
            NodeID border_id;
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                if(border_id==t){
                    //cout<<"hash"<<endl;//for debug
                    return idx_s.dis[i];
                }           
            }
            //cout<<"pll"<<endl;//for debug
            //pll 2-hop
            const integrated_bindex_t_p &idx_t = bindex_p[t];
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                const integrated_bindex_t_p &idx_border = bindex_p[border_id];
                distance_t=query_p_directed(idx_border,idx_t)+idx_s.dis[i];
                if(distance_t<distance) distance=distance_t;
            }
            return distance;
        }else{//l-l
            //cout<<"l-l ";//for debug
            const integrated_bindex_t_p &idx_s = bindex_p[s];
            const integrated_bindex_t_p &idx_t = bindex_p[t];
            if(idx_s.pos[0]==idx_t.pos[0]){//in the same tree
               //cout<<"h2h ";//for debug
                int lca_id=lca.lca_query(s,t);
                //cout<<"lca="<<lca_id;//for debug
                if(lca_id==s){
                    //cout<<" 0"<<endl;//for debug
                   // _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    return idx_t.r_dis[idx_s.pos[1]];
                }else if(lca_id==t){
                    //cout<<" 1"<<endl;//for debug
                    //_mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    return idx_s.dis[idx_t.pos[1]];
                }else{
                    //cout<<" 2"<<endl;//for debug
                    const integrated_bindex_t_p &idx_lca = bindex_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_t.r_dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int pos_i;
                    for(int i=1;(pos_i=idx_lca.pos[i])!=numOfVertices;++i){
                        distance_t=idx_s.dis[pos_i]+idx_t.r_dis[pos_i];
                        if(distance_t<distance) distance=distance_t;
                    }
                    return distance;
                }
            }else{//cover different trees
               //cout<<s<<"-"<<t<<" pll:"<<endl;//for debug
               double time_scan_add,time_scan=0;//for debug
               //variables 
               double time_variable=GetCurrentTimeSec();
                NodeID border_id;
                NodeID spt_v_t,r_spt_v_t;
                EdgeWeight spt_d_t,r_spt_d_t;
                int i,j;
                // EdgeWeight* idx_s_vec=new EdgeWeight[_max_overlay_hub+1];
                // EdgeWeight* idx_t_vec=new EdgeWeight[_max_overlay_hub+1];
                // for(size_t i=0;i<_max_overlay_hub+1;++i) idx_s_vec[i]=idx_t_vec[i]=INF_WEIGHT;
                
                //first scan the s node region
                const integrated_bindex_t_p &idx_root_s = bindex_p[idx_s.pos[0]];
                _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                _mm_prefetch(&idx_root_s.spt_v[0], _MM_HINT_T0);
                vector<EdgeWeight> idx_s_vec(_numOfOverlayVertices,INF_WEIGHT);
                for(i=0;(border_id=idx_root_s.spt_v[i])!=numOfVertices;++i){
                    if(idx_s.dis[i]==INF_WEIGHT) continue;
                    const integrated_bindex_t_p &idx_border = bindex_p[border_id];
                    for(j=0;(spt_v_t=idx_border.spt_v[j])!=numOfVertices;++j){
                        spt_v_t=idx_border.spt_v[j];
                        spt_d_t=idx_border.spt_d[j]+idx_s.dis[i];
                        if(spt_d_t<idx_s_vec[spt_v_t]) idx_s_vec[spt_v_t]=spt_d_t;
                    }
                }
                //second scan the t node region
                const integrated_bindex_t_p &idx_root_t = bindex_p[idx_t.pos[0]]; 
                _mm_prefetch(&idx_t.r_dis[0], _MM_HINT_T0);
                _mm_prefetch(&idx_root_t.r_spt_v[0], _MM_HINT_T0);
                vector<EdgeWeight> idx_t_vec(_numOfOverlayVertices,INF_WEIGHT);
                for(i=0;(border_id=idx_root_t.spt_v[i])!=numOfVertices;++i){
                    if(idx_t.r_dis[i]==INF_WEIGHT) continue;
                    const integrated_bindex_t_p &idx_border = bindex_p[border_id];
                    for(j=0;(spt_v_t=idx_border.r_spt_v[j])!=numOfVertices;++j){
                        r_spt_v_t=idx_border.r_spt_v[j];
                        r_spt_d_t=idx_border.r_spt_d[j]+idx_t.r_dis[i];
                        //std::cout<<"r_spt_v_t="<<r_spt_v_t<<" r_spt_d_t="<<r_spt_d_t<<std::endl;
                        if(r_spt_d_t<idx_t_vec[r_spt_v_t]) idx_t_vec[r_spt_v_t]=r_spt_d_t;
                    }
                }
                //merge scane
                return query_p_directed(idx_s_vec,idx_t_vec);
                //return query_p(idx_s_vec,idx_t_vec,_max_overlay_hub+1);
                // for(int i=0;(border_id_s=idx_root_s.spt_v[i])!=numOfVertices;++i){
                //     const integrated_index_t_p &idx_border_s = index_p[border_id_s];
                //     for(int j=0;(border_id_t=idx_root_t.spt_v[j])!=numOfVertices;++j){
                //         const integrated_index_t_p &idx_border_t = index_p[border_id_t];
                //         distance_t=query_p(idx_border_s,idx_border_t)+idx_s.dis[i+1]+idx_t.dis[j+1];
                //         if(distance_t<distance) distance=distance_t;
                //     }
                // }
                //return distance;
            }
        }
        return distance;
    }

    EdgeWeight query_integrated_1_directed(NodeID s,NodeID t,int& type){
        EdgeWeight distance = INF_WEIGHT, distance_t = INF_WEIGHT;
        //cout<<s<<"-"<<t<<":";//for debug
        if(lca._index[s]==-1&&lca._index[t]==-1){//h-h
            //cout<<"h-h"<<endl;//for debug
            type=0;
            const integrated_bindex_t_p &idx_s = bindex_p[s];
            const integrated_bindex_t_p &idx_t = bindex_p[t];
            return query_p_directed(idx_s,idx_t);
        }else if(lca._index[s]==-1&&lca._index[t]!=-1){//h-l
            //cout<<"h-l ";//for debug
            const integrated_bindex_t_p &idx_t = bindex_p[t];
            const integrated_bindex_t_p &idx_root = bindex_p[idx_t.pos[0]];
            _mm_prefetch(&idx_root.spt_v[0], _MM_HINT_T0);
            //judge borders
            NodeID border_id;
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                if(border_id==s){
                    //cout<<"hash"<<endl;//for debug
                    type=1;
                    return idx_t.r_dis[i];
                }           
            }
            //cout<<"pll"<<endl;//for debug
            //pll 2-hop
            const integrated_bindex_t_p &idx_s = bindex_p[s];
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                const integrated_bindex_t_p &idx_border = bindex_p[border_id];
                distance_t=query_p_directed(idx_s,idx_border)+idx_t.r_dis[i];
                if(distance_t<distance) distance=distance_t;
            }
            type=2;
            return distance;
        }else if(lca._index[s]!=-1&&lca._index[t]==-1){//l-h
            //cout<<"l-h ";//for debug
            const integrated_bindex_t_p &idx_s = bindex_p[s];
            const integrated_bindex_t_p &idx_root = bindex_p[idx_s.pos[0]];
            _mm_prefetch(&idx_root.spt_v[0], _MM_HINT_T0);
            //judge borders
            NodeID border_id;
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                if(border_id==t){
                    //cout<<"hash"<<endl;//for debug
                    type=1;
                    return idx_s.dis[i];
                }           
            }
            //cout<<"pll"<<endl;//for debug
            //pll 2-hop
            const integrated_bindex_t_p &idx_t = bindex_p[t];
            for(int i=0;(border_id=idx_root.spt_v[i])!=numOfVertices;++i){
                const integrated_bindex_t_p &idx_border = bindex_p[border_id];
                distance_t=query_p_directed(idx_border,idx_t)+idx_s.dis[i];
                if(distance_t<distance) distance=distance_t;
            }
            type=2;
            return distance;
        }else{//l-l
            //cout<<"l-l ";//for debug
            const integrated_bindex_t_p &idx_s = bindex_p[s];
            const integrated_bindex_t_p &idx_t = bindex_p[t];
            if(idx_s.pos[0]==idx_t.pos[0]){//in the same tree
                type=3;
               //cout<<"h2h ";//for debug
                int lca_id=lca.lca_query(s,t);
                //cout<<"lca="<<lca_id;//for debug
                if(lca_id==s){
                    //cout<<" 0"<<endl;//for debug
                   // _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    return idx_t.r_dis[idx_s.pos[1]];
                }else if(lca_id==t){
                    //cout<<" 1"<<endl;//for debug
                    //_mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    return idx_s.dis[idx_t.pos[1]];
                }else{
                    //cout<<" 2"<<endl;//for debug
                    const integrated_bindex_t_p &idx_lca = bindex_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_t.r_dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int pos_i;
                    for(int i=1;(pos_i=idx_lca.pos[i])!=numOfVertices;++i){
                        distance_t=idx_s.dis[pos_i]+idx_t.r_dis[pos_i];
                        if(distance_t<distance) distance=distance_t;
                    }
                    return distance;
                }
            }else{//cover different trees
                type=4;
               //cout<<s<<"-"<<t<<" pll:"<<endl;//for debug
               double time_scan_add,time_scan=0;//for debug
               //variables 
               double time_variable=GetCurrentTimeSec();
                NodeID border_id;
                NodeID spt_v_t,r_spt_v_t;
                EdgeWeight spt_d_t,r_spt_d_t;
                int i,j;
                //first scan the s node region
                const integrated_bindex_t_p &idx_root_s = bindex_p[idx_s.pos[0]];
                _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                _mm_prefetch(&idx_root_s.spt_v[0], _MM_HINT_T0);
                vector<EdgeWeight> idx_s_vec(_numOfOverlayVertices,INF_WEIGHT);
                for(i=0;(border_id=idx_root_s.spt_v[i])!=numOfVertices;++i){
                    if(idx_s.dis[i]==INF_WEIGHT) continue;
                    const integrated_bindex_t_p &idx_border = bindex_p[border_id];
                    for(j=0;(spt_v_t=idx_border.spt_v[j])!=numOfVertices;++j){
                        spt_v_t=idx_border.spt_v[j];
                        spt_d_t=idx_border.spt_d[j]+idx_s.dis[i];
                        if(spt_d_t<idx_s_vec[spt_v_t]) idx_s_vec[spt_v_t]=spt_d_t;
                    }
                }
                //second scan the t node region
                const integrated_bindex_t_p &idx_root_t = bindex_p[idx_t.pos[0]]; 
                _mm_prefetch(&idx_t.r_dis[0], _MM_HINT_T0);
                _mm_prefetch(&idx_root_t.r_spt_v[0], _MM_HINT_T0);
                vector<EdgeWeight> idx_t_vec(_numOfOverlayVertices,INF_WEIGHT);
                for(i=0;(border_id=idx_root_t.spt_v[i])!=numOfVertices;++i){
                    if(idx_t.r_dis[i]==INF_WEIGHT) continue;
                    const integrated_bindex_t_p &idx_border = bindex_p[border_id];
                    for(j=0;(spt_v_t=idx_border.r_spt_v[j])!=numOfVertices;++j){
                        r_spt_v_t=idx_border.r_spt_v[j];
                        r_spt_d_t=idx_border.r_spt_d[j]+idx_t.r_dis[i];
                        //std::cout<<"r_spt_v_t="<<r_spt_v_t<<" r_spt_d_t="<<r_spt_d_t<<std::endl;
                        if(r_spt_d_t<idx_t_vec[r_spt_v_t]) idx_t_vec[r_spt_v_t]=r_spt_d_t;
                    }
                }
                //merge scane
                return query_p_directed(idx_s_vec,idx_t_vec);
            }
        }
        return distance;
    }

    /**
     * @description: used to query processing for the third solution
     * @param {NodeID} s
     * @param {NodeID} t
     * @return {*}
     * @author: Wan Jingyi
     */    
    EdgeWeight query_integrated_2(NodeID s,NodeID t){
        EdgeWeight distance = INF_WEIGHT,distance_t = INF_WEIGHT;
        //cout<<s<<"-"<<t<<":";//for debug
        const integrated_index_t_p &idx_s = index_p[s];
        const integrated_index_t_p &idx_t = index_p[t];
        if(lca._index[s]==-1||lca._index[t]==-1){//one is h
            //for debug
            //if(lca._index[s]==-1&&lca._index[t]!=-1) cout<<"h-l ";
            //else if(lca._index[s]==-1&&lca._index[t]!=-1) cout<<"l-h ";
            //else cout<<"h-h"<<endl;
            return query_p(idx_s,idx_t);
        }else{//l-l
            //cout<<"l-l ";//for debug
            if(idx_s.dis[0]==idx_t.dis[0]){//in the same tree
                //cout<<"h2h ";//for debug
                int lca_id=lca.lca_query(s,t);
                //cout<<"lca="<<lca_id;//for debug
                if(lca_id==s){
                    //cout<<" 0"<<endl;//for debug
                    return idx_t.dis[idx_s.pos[1]+1];
                }else if(lca_id==t){
                    //cout<<" 1"<<endl;//for debug
                    return idx_s.dis[idx_t.pos[1]+1];
                }else{
                    //cout<<" 2"<<endl;//for debug
                    const integrated_index_t_p &idx_lca = index_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int pos_i;
                    if(idx_lca.pos[0]){
                        distance=query_p(idx_s,idx_t);
                    }
                    for(int i=1;(pos_i=idx_lca.pos[i])!=numOfVertices;++i){
                        distance_t=idx_s.dis[pos_i+1]+idx_t.dis[pos_i+1];
                        if(distance_t<distance) distance=distance_t;
                    }
                    return distance;
                }
            }else{
                //cout<<"pll"<<endl;//for debug
                return query_p(idx_s,idx_t);
            }
        }
        return distance;
    }

    EdgeWeight query_integrated_2(NodeID s,NodeID t,int& type){
        EdgeWeight distance = INF_WEIGHT,distance_t = INF_WEIGHT;
        //cout<<s<<"-"<<t<<":";//for debug
        const integrated_index_t_p &idx_s = index_p[s];
        const integrated_index_t_p &idx_t = index_p[t];
        if(lca._index[s]==-1||lca._index[t]==-1){//one is h
            //for debug
            //if(lca._index[s]==-1&&lca._index[t]!=-1) cout<<"h-l ";
            //else if(lca._index[s]==-1&&lca._index[t]!=-1) cout<<"l-h ";
            //else cout<<"h-h"<<endl;
            type=0;
            return query_p(idx_s,idx_t);
        }else{//l-l
            //cout<<"l-l ";//for debug
            if(idx_s.dis[0]==idx_t.dis[0]){//in the same tree
                type=1;
                //cout<<"h2h ";//for debug
                int lca_id=lca.lca_query(s,t);
                //cout<<"lca="<<lca_id;//for debug
                if(lca_id==s){
                    //cout<<" 0"<<endl;//for debug
                    return idx_t.dis[idx_s.pos[1]+1];
                }else if(lca_id==t){
                    //cout<<" 1"<<endl;//for debug
                    return idx_s.dis[idx_t.pos[1]+1];
                }else{
                    //cout<<" 2"<<endl;//for debug
                    const integrated_index_t_p &idx_lca = index_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int pos_i;
                    if(idx_lca.pos[0]){
                        distance=query_p(idx_s,idx_t);
                    }
                    for(int i=1;(pos_i=idx_lca.pos[i])!=numOfVertices;++i){
                        distance_t=idx_s.dis[pos_i+1]+idx_t.dis[pos_i+1];
                        if(distance_t<distance) distance=distance_t;
                    }
                    return distance;
                }
            }else{
                //cout<<"pll"<<endl;//for debug
                type=2;
                return query_p(idx_s,idx_t);
            }
        }
        return distance;
    }

    /**
     * @description: used to query processing for the third solution, add for directed
     * @param {NodeID} s
     * @param {NodeID} t
     * @return {*}
     * @author: Wan Jingyi
     */    
    EdgeWeight query_integrated_2_directed(NodeID s,NodeID t){
        EdgeWeight distance = INF_WEIGHT,distance_t = INF_WEIGHT;
        //cout<<s<<"-"<<t<<":";//for debug
        const integrated_bindex_t_p &idx_s = bindex_p[s];
        const integrated_bindex_t_p &idx_t = bindex_p[t];
        if(lca._index[s]==-1||lca._index[t]==-1){//one is h
            return query_p_directed(idx_s,idx_t);
        }else{//l-l
            //cout<<"l-l ";//for debug
            if(idx_s.pos[1]==idx_t.pos[1]){//in the same tree
                //cout<<"h2h ";//for debug
                int lca_id=lca.lca_query(s,t);
                //cout<<"lca="<<lca_id;//for debug
                if(lca_id==s){
                    //cout<<" 0"<<endl;//for debug
                    return idx_t.r_dis[idx_s.pos[2]];
                }else if(lca_id==t){
                    //cout<<" 1"<<endl;//for debug
                    return idx_s.dis[idx_t.pos[2]];
                }else{
                    //cout<<" 2"<<endl;//for debug
                    const integrated_bindex_t_p &idx_lca =bindex_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_t.r_dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int pos_i;
                    if(idx_lca.pos[0]){
                        distance=query_p_directed(idx_s,idx_t);
                    }
                    for(int i=2;(pos_i=idx_lca.pos[i])!=numOfVertices;++i){
                        distance_t=idx_s.dis[pos_i]+idx_t.r_dis[pos_i];
                        if(distance_t<distance) distance=distance_t;
                    }
                    return distance;
                }
            }else{
                //cout<<"pll"<<endl;//for debug
                return query_p_directed(idx_s,idx_t);
            }
        }
        return distance;
    }

    EdgeWeight query_integrated_2_directed(NodeID s,NodeID t,int& type){
        EdgeWeight distance = INF_WEIGHT,distance_t = INF_WEIGHT;
        //cout<<s<<"-"<<t<<":";//for debug
        const integrated_bindex_t_p &idx_s = bindex_p[s];
        const integrated_bindex_t_p &idx_t = bindex_p[t];
        if(lca._index[s]==-1||lca._index[t]==-1){//one is h
            type=0;
            return query_p_directed(idx_s,idx_t);
        }else{//l-l
            //cout<<"l-l ";//for debug
            if(idx_s.pos[1]==idx_t.pos[1]){//in the same tree
                type=1;
                //cout<<"h2h ";//for debug
                int lca_id=lca.lca_query(s,t);
                //cout<<"lca="<<lca_id;//for debug
                if(lca_id==s){
                    //cout<<" 0"<<endl;//for debug
                    return idx_t.r_dis[idx_s.pos[2]];
                }else if(lca_id==t){
                    //cout<<" 1"<<endl;//for debug
                    return idx_s.dis[idx_t.pos[2]];
                }else{
                    //cout<<" 2"<<endl;//for debug
                    const integrated_bindex_t_p &idx_lca =bindex_p[lca_id];
                    _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_t.r_dis[0], _MM_HINT_T0);
                    _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
                    int pos_i;
                    if(idx_lca.pos[0]){
                        distance=query_p_directed(idx_s,idx_t);
                    }
                    for(int i=2;(pos_i=idx_lca.pos[i])!=numOfVertices;++i){
                        distance_t=idx_s.dis[pos_i]+idx_t.r_dis[pos_i];
                        if(distance_t<distance) distance=distance_t;
                    }
                    return distance;
                }
            }else{
                type=2;
                //cout<<"pll"<<endl;//for debug
                return query_p_directed(idx_s,idx_t);
            }
        }
        return distance;
    }

    void load_is_deleted(char* isDeletedFIleName){
        ifstream in(isDeletedFIleName);
        if(!in.is_open()) cerr<<isDeletedFIleName<<" cannot be opened!"<<endl;
        in>>_numOfOriginalVertices>>_numOfOverlayVertices;
        if(isDeleted==nullptr) isDeleted=new bool[_numOfOriginalVertices];
        int i,u,flag;
        //build the map relationship 
        for (i = 0; in >> u >> flag;i++) {
            isDeleted[u]=flag;
        }
        if(i!=_numOfOriginalVertices) cerr<<"error:i!=_numOfOriginalVertices!"<<endl;	
        in.close();
        
    }

    void load_labels(char* load_filename){
        if(index_p){
            for(size_t i=0;i<numOfVertices;++i){
                if(index_p[i].spt_v!=nullptr) free(index_p[i].spt_v);
                if(index_p[i].spt_d!=nullptr) free(index_p[i].spt_d);
                if(index_p[i].pos!=nullptr) free(index_p[i].pos);
                if(index_p[i].dis!=nullptr) free(index_p[i].dis);
            }
            free(index_p);
        }
        index_p=nullptr;

        ifstream ifs(load_filename);
        if(!ifs.is_open()) std::cout<<"Cannot open "<<load_filename<<endl;
        NodeID isize;
        NodeID v_tmp;
        EdgeWeight w_tmp;
        //first line
        ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		cout<<"numOfVertices = "<<numOfVertices<<endl;
        
        index_p = (integrated_index_t_p*)memalign(64, numOfVertices * sizeof(integrated_index_t_p));

        for(NodeID v=0;v<numOfVertices;++v){
            integrated_index_t_p& idx=index_p[v];
            //pos
            ifs.read((char*)&isize, sizeof(isize));
			idx.pos = (NodeID*)memalign(64, isize * sizeof(NodeID));
			for (NodeID i = 0; i < isize; ++i) {
				ifs.read((char*)&v_tmp, sizeof(v_tmp));
				idx.pos[i] = v_tmp;
			}
            if(!isDeleted[v]){
                //spt-v,spt-d
                ifs.read((char*)&isize, sizeof(isize));
                idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
                idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.spt_v[i] = v_tmp;
                    idx.spt_d[i] = w_tmp;
                }
            }else{
                //dis
                ifs.read((char*)&isize, sizeof(isize));
                idx.dis = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.dis[i] = w_tmp;
                }
            }

        }

        ifs.close();
    }


    /**
     * @description: load labels for the second solution
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void load_labels_1(char* load_filename){
        if(index_p){
            // for(size_t i=0;i<numOfVertices;++i){
                // if(index_p[i].spt_v!=nullptr){
                //     free(index_p[i].spt_v);
                //     index_p[i].spt_v=nullptr;
                //     cout<<"spt_v ";//for debug
                // }
                // if(index_p[i].spt_d!=nullptr){
                //     free(index_p[i].spt_d);
                //     index_p[i].spt_d=nullptr;
                //     cout<<"spt_d ";//for debug
                // }
                // if(index_p[i].pos!=nullptr){
                //     free(index_p[i].pos);
                //     index_p[i].pos=nullptr;
                //     cout<<"pos ";//for debug
                // }
                // if(index_p[i].dis!=nullptr){
                //     free(index_p[i].dis);
                //     index_p[i].dis=nullptr;
                //     cout<<"dis ";//for debug
                // }
            // }
            free(index_p);
        }
        index_p=nullptr;
        ifstream ifs(load_filename);
        if(!ifs.is_open()) std::cout<<"Cannot open "<<load_filename<<endl;
        NodeID isize;
        NodeID v_tmp;
        EdgeWeight w_tmp;
        //first line
        ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		//cout<<"numOfVertices = "<<numOfVertices<<endl;
        
        index_p = (integrated_index_t_p*)memalign(64, numOfVertices * sizeof(integrated_index_t_p));

        for(NodeID v=0;v<numOfVertices;++v){
            integrated_index_t_p& idx=index_p[v];
            if(isDeleted[v]){
                //pos
                ifs.read((char*)&isize, sizeof(isize));
                idx.pos = (NodeID*)memalign(64, isize * sizeof(NodeID));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    idx.pos[i] = v_tmp;
                }
                //dis
                ifs.read((char*)&isize, sizeof(isize));
                idx.dis = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.dis[i] = w_tmp;
                }
                //spt-v
                ifs.read((char*)&isize, sizeof(isize));
                if(isize==0) continue;
                idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    idx.spt_v[i] = v_tmp;
                }
            }else{
                //spt-v,spt-d
                ifs.read((char*)&isize, sizeof(isize));
                idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
                idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.spt_v[i] = v_tmp;
                    idx.spt_d[i] = w_tmp;
                    if((i!=isize-1)&&(v_tmp>_max_overlay_hub)) _max_overlay_hub=v_tmp;
                }
            }

        }

        ifs.close();
    }

    /**
     * @description: load labels for the second solution,add for directed
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void load_labels_1_directed(char* load_filename){
        if(bindex_p){
            free(bindex_p);
        }
        bindex_p=nullptr;
        ifstream ifs(load_filename);
        if(!ifs.is_open()) std::cout<<"Cannot open "<<load_filename<<endl;
        NodeID isize;
        NodeID v_tmp;
        EdgeWeight w_tmp;
        //first line
        ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		//cout<<"numOfVertices = "<<numOfVertices<<endl;
        
        bindex_p = (integrated_bindex_t_p*)memalign(64, numOfVertices * sizeof(integrated_bindex_t_p));

        for(NodeID v=0;v<numOfVertices;++v){
            integrated_bindex_t_p& idx=bindex_p[v];
            if(isDeleted[v]){
                //pos
                ifs.read((char*)&isize, sizeof(isize));
                idx.pos = (NodeID*)memalign(64, isize * sizeof(NodeID));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    idx.pos[i] = v_tmp;
                }
                //dis and r_dis
                ifs.read((char*)&isize, sizeof(isize));
                idx.dis = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.dis[i] = w_tmp;
                }
                idx.r_dis = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.r_dis[i] = w_tmp;
                }
                //spt_v_siz
                ifs.read((char*)&isize, sizeof(isize));
                if(isize==0) continue;
                //spt_v
                idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    idx.spt_v[i] = v_tmp;
                }
            }else{
                //spt_v,spt_d
                ifs.read((char*)&isize, sizeof(isize));
                idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
                idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.spt_v[i] = v_tmp;
                    idx.spt_d[i] = w_tmp;
                    //if((i!=isize-1)&&(v_tmp>_max_overlay_hub)) _max_overlay_hub=v_tmp;
                }
                //r_spt_v,r_spt_d
                ifs.read((char*)&isize, sizeof(isize));
                idx.r_spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
                idx.r_spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.r_spt_v[i] = v_tmp;
                    idx.r_spt_d[i] = w_tmp;
                    //if((i!=isize-1)&&(v_tmp>_max_overlay_hub_r)) _max_overlay_hub_r=v_tmp;
                }
            }

        }

        ifs.close();
    }

    /**
     * @description: load labels for the third solution
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void load_labels_2(char* load_filename){
        if(index_p){
            free(index_p);
        }
        index_p=nullptr;

        ifstream ifs(load_filename);
        if(!ifs.is_open()) std::cout<<"Cannot open "<<load_filename<<endl;
        NodeID isize;
        NodeID v_tmp;
        EdgeWeight w_tmp;
        //first line
        ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		cout<<"numOfVertices = "<<numOfVertices<<endl;
        
        index_p = (integrated_index_t_p*)memalign(64, numOfVertices * sizeof(integrated_index_t_p));

        for(NodeID v=0;v<numOfVertices;++v){
            integrated_index_t_p& idx=index_p[v];
            if(isDeleted[v]){
                //pos
                ifs.read((char*)&isize, sizeof(isize));
                idx.pos = (NodeID*)memalign(64, isize * sizeof(NodeID));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    idx.pos[i] = v_tmp;
                }
                //dis
                ifs.read((char*)&isize, sizeof(isize));
                idx.dis = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.dis[i] = w_tmp;
                }
            }
            //spt-v,spt-d
            ifs.read((char*)&isize, sizeof(isize));
            idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
            idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
            for (NodeID i = 0; i < isize; ++i) {
                ifs.read((char*)&v_tmp, sizeof(v_tmp));
                ifs.read((char*)&w_tmp, sizeof(w_tmp));
                idx.spt_v[i] = v_tmp;
                idx.spt_d[i] = w_tmp;
            }
        }
        ifs.close();
    }

    /**
     * @description: load labels for the third solution, add for directed
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void load_labels_2_directed(char* load_filename){
        if(bindex_p){
            free(bindex_p);
        }
        bindex_p=nullptr;

        ifstream ifs(load_filename);
        if(!ifs.is_open()) std::cout<<"Cannot open "<<load_filename<<endl;
        NodeID isize;
        NodeID v_tmp,r_v_tmp;
        EdgeWeight w_tmp,r_w_tmp;
        //first line
        ifs.read((char*)&isize, sizeof(isize));
		numOfVertices = isize;
		cout<<"numOfVertices = "<<numOfVertices<<endl;
        
        bindex_p = (integrated_bindex_t_p*)memalign(64, numOfVertices * sizeof(integrated_bindex_t_p));

        for(NodeID v=0;v<numOfVertices;++v){
            integrated_bindex_t_p& idx=bindex_p[v];
            if(isDeleted[v]){
                //pos
                ifs.read((char*)&isize, sizeof(isize));
                idx.pos = (NodeID*)memalign(64, isize * sizeof(NodeID));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&v_tmp, sizeof(v_tmp));
                    idx.pos[i] = v_tmp;
                }
                //dis
                ifs.read((char*)&isize, sizeof(isize));
                idx.dis = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.dis[i] = w_tmp;
                }
                //rdis
                idx.r_dis = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
                for (NodeID i = 0; i < isize; ++i) {
                    ifs.read((char*)&w_tmp, sizeof(w_tmp));
                    idx.r_dis[i] = w_tmp;
                }
            }
            //spt-v,spt-d
            ifs.read((char*)&isize, sizeof(isize));
            idx.spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
            idx.spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
            for (NodeID i = 0; i < isize; ++i) {
                ifs.read((char*)&v_tmp, sizeof(v_tmp));
                ifs.read((char*)&w_tmp, sizeof(w_tmp));
                idx.spt_v[i] = v_tmp;
                idx.spt_d[i] = w_tmp;
            }
            //r-spt-v,r-spt-d
            ifs.read((char*)&isize, sizeof(isize));
            idx.r_spt_v = (NodeID*)memalign(64, isize * sizeof(NodeID));
            idx.r_spt_d = (EdgeWeight*)memalign(64, isize * sizeof(EdgeWeight));
            for (NodeID i = 0; i < isize; ++i) {
                ifs.read((char*)&v_tmp, sizeof(v_tmp));
                ifs.read((char*)&w_tmp, sizeof(w_tmp));
                idx.r_spt_v[i] = v_tmp;
                idx.r_spt_d[i] = w_tmp;
            }
        }
        ifs.close();
    }

    void load_lca_variables(char* load_filename){
        lca.initialize(load_filename);
    }

    void save_labels(char* labelFileName){
        ofstream ofs(labelFileName, ios::binary | ios::out);
        if(!ofs.is_open()) cout<<"Cannot open"<<labelFileName<<endl;
        ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
        NodeID isize;
        //save format:size-pos,size-(spt_v spt_d),size-dis
		for (NodeID v = 0; v < numOfVertices; ++v) {
			const integrated_index_t& idx=index_[v];
            if((isize=idx.pos.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.pos[i], sizeof(idx.pos[i]));
                }
            }
            if((isize=idx.spt_v.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                    ofs.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
                }
            }
            if((isize=idx.dis.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.dis[i], sizeof(idx.dis[i]));
                }
            }
		}
		ofs.close();
        std::cout<<"Save integrated_index labels to binary file successfully!"<<std::endl;
    }

    void save_labels_differ_analysis_1(char* outputDirName){
        string overlay_label_file_name(outputDirName);
        string lf_label_file_name(outputDirName);
        overlay_label_file_name.append("_overlay.label");
        lf_label_file_name.append("_lf.label");
        ofstream ofs_overlay(overlay_label_file_name, ios::binary | ios::out);
        ofstream ofs_lf(lf_label_file_name, ios::binary | ios::out);
        if(!ofs_overlay.is_open()) cout<<"Cannot open"<<overlay_label_file_name<<endl;
        if(!ofs_lf.is_open()) cout<<"Cannot open"<<lf_label_file_name<<endl;
        ofs_overlay.write((const char*)&numOfVertices, sizeof(numOfVertices));
        NodeID isize;
        for (NodeID v = 0; v < numOfVertices; ++v) {
			const integrated_index_t& idx=index_[v];
            if(idx.pos.size()!=0){//lf
                if((isize=idx.pos.size())!=0){
                    ofs_lf.write((const char*)&isize, sizeof(isize));
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs_lf.write((const char*)&idx.pos[i], sizeof(idx.pos[i]));
                    }
                }
                if((isize=idx.dis.size())!=0){
                    ofs_lf.write((const char*)&isize, sizeof(isize));
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs_lf.write((const char*)&idx.dis[i], sizeof(idx.dis[i]));
                    }
                }
                if((isize=idx.spt_v.size())==0){
                    ofs_lf.write((const char*)&isize, sizeof(isize));
                }else{
                    ofs_lf.write((const char*)&isize, sizeof(isize));
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs_lf.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                    }
                }
            }else{
                isize=idx.spt_v.size();
                ofs_overlay.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs_overlay.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                    ofs_overlay.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
                }
            }
        }
		ofs_lf.close();
        ofs_overlay.close();
        std::cout<<"Save differ labels analysis to binary file successfully!"<<std::endl;
    }

    void save_labels_differ_analysis_2(char* outputDirName){
        string overlay_label_file_name(outputDirName);
        string lf_pll_file_name(outputDirName);
        string lf_h2h_file_name(outputDirName);
        overlay_label_file_name.append("_overlay.label");
        lf_pll_file_name.append("_lf_pll.label");
        lf_h2h_file_name.append("_lf_h2h.label");
        ofstream ofs_overlay(overlay_label_file_name, ios::binary | ios::out);
        ofstream ofs_lf_pll(lf_pll_file_name, ios::binary | ios::out);
        ofstream ofs_lf_h2h(lf_h2h_file_name, ios::binary | ios::out);
        if(!ofs_overlay.is_open()) cout<<"Cannot open"<<overlay_label_file_name<<endl;
        if(!ofs_lf_pll.is_open()) cout<<"Cannot open"<<lf_pll_file_name<<endl;
        if(!ofs_lf_h2h.is_open()) cout<<"Cannot open"<<lf_h2h_file_name<<endl;
        ofs_overlay.write((const char*)&numOfVertices, sizeof(numOfVertices));
        NodeID isize;
        for (NodeID v = 0; v < numOfVertices; ++v) {
			const integrated_index_t& idx=index_[v];
            if(idx.pos.size()!=0){//lf
                if((isize=idx.pos.size())!=0){
                    ofs_lf_h2h.write((const char*)&isize, sizeof(isize));
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs_lf_h2h.write((const char*)&idx.pos[i], sizeof(idx.pos[i]));
                    }
                }
                if((isize=idx.dis.size())!=0){
                    ofs_lf_h2h.write((const char*)&isize, sizeof(isize));
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs_lf_h2h.write((const char*)&idx.dis[i], sizeof(idx.dis[i]));
                    }
                }
                if((isize=idx.spt_v.size())!=0){
                    ofs_lf_pll.write((const char*)&isize, sizeof(isize));
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs_lf_pll.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                        ofs_lf_pll.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
                    }
                }
            }else{
                isize=idx.spt_v.size();
                ofs_overlay.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs_overlay.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                    ofs_overlay.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
                }
            }
        }
		ofs_lf_h2h.close();
        ofs_lf_pll.close();
        ofs_overlay.close();
        std::cout<<"Save differ labels analysis to binary file successfully!"<<std::endl;
    }

    /**
     * @description: save labels for the second solution
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void save_labels_1(char* labelFileName){
        ofstream ofs(labelFileName, ios::binary | ios::out);
        if(!ofs.is_open()) cout<<"Cannot open"<<labelFileName<<endl;
        ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
        NodeID isize;
        //save format:
        //hf:(spt-v spt-d)
        //lf:pos dis spt-v
		for (NodeID v = 0; v < numOfVertices; ++v) {
			const integrated_index_t& idx=index_[v];
            if((isize=idx.pos.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.pos[i], sizeof(idx.pos[i]));
                }
            }
            if((isize=idx.dis.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.dis[i], sizeof(idx.dis[i]));
                }
            }
            if((isize=idx.spt_v.size())==0){
                ofs.write((const char*)&isize, sizeof(isize));
            }else{
                if(idx.spt_d.size()!=0){//hf
                    ofs.write((const char*)&isize, sizeof(isize));
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                        ofs.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
                    }
                }else{
                    ofs.write((const char*)&isize, sizeof(isize));
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                    }
                }  
            }
		}
		ofs.close();
        std::cout<<"Save integrated_index labels to binary file successfully!"<<std::endl;
    }
    
    /**
     * @description: save labels for the second solution, add for directed
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void save_labels_1_directed(char* labelFileName){
        ofstream ofs(labelFileName, ios::binary | ios::out);
        if(!ofs.is_open()) cout<<"Cannot open"<<labelFileName<<endl;
        ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
        NodeID isize;
        //save format:
        //hf:spt_v_siz spt_v spt_d r_spt_v_siz r_spt_v r_spt_d
        //lf:pos_siz pos dis_siz dis r_dis spt_v_siz + (spt_v)
		for (NodeID v = 0; v < numOfVertices; ++v) {
			const integrated_bindex_t& idx=bindex_[v];
            //pos
            if((isize=idx.pos.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.pos[i], sizeof(idx.pos[i]));
                }
            }
            if((isize=idx.dis.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                //dis
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.dis[i], sizeof(idx.dis[i]));
                }
                //r_dis
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.r_dis[i], sizeof(idx.r_dis[i]));
                }
            }
            //spt_v_siz
            if((isize=idx.spt_v.size())==0){
                ofs.write((const char*)&isize, sizeof(isize));
            }else{
                if(idx.spt_d.size()!=0){//hf
                    ofs.write((const char*)&isize, sizeof(isize));
                    //spt_v spt_d
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                        ofs.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
                    }
                    //r_spt_v_siz
                    isize=idx.r_spt_v.size();
                    ofs.write((const char*)&isize, sizeof(isize));
                    //r_spt_v r_spt_d
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs.write((const char*)&idx.r_spt_v[i], sizeof(idx.r_spt_v[i]));
                        ofs.write((const char*)&idx.r_spt_d[i], sizeof(idx.r_spt_d[i]));
                    }
                }else{
                    ofs.write((const char*)&isize, sizeof(isize));
                    //spt_v
                    for (NodeID i = 0; i < isize; ++i) {
                        ofs.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                    }
                }  
            }
		}
		ofs.close();
        std::cout<<"Save integrated_index directed labels to binary file successfully!"<<std::endl;
    }

    /**
     * @description: save labels for the third solution
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void save_labels_2(char* labelFileName){
        ofstream ofs(labelFileName, ios::binary | ios::out);
        if(!ofs.is_open()) cout<<"Cannot open"<<labelFileName<<endl;
        ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
        NodeID isize;
        //save format:pos dis (spt-v spt-d)
		for (NodeID v = 0; v < numOfVertices; ++v) {
			const integrated_index_t& idx=index_[v];
            if((isize=idx.pos.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.pos[i], sizeof(idx.pos[i]));
                }
            }
            if((isize=idx.dis.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.dis[i], sizeof(idx.dis[i]));
                }
            }
            if((isize=idx.spt_v.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                    ofs.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
                }
            } 
		}
		ofs.close();
        std::cout<<"Save integrated_index labels to binary file successfully!"<<std::endl;
    }

    /**
     * @description: save labels for the third solution,add for directed
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void save_labels_2_directed(char* labelFileName){
        ofstream ofs(labelFileName, ios::binary | ios::out);
        if(!ofs.is_open()) cout<<"Cannot open"<<labelFileName<<endl;
        ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
        NodeID isize;
        //save format:pos_size pos dis_size dis r_dis (spt-v-size spt-v spt-d r-spt-v-size r-spt-v r-spt-d)
		for (NodeID v = 0; v < numOfVertices; ++v) {
			const integrated_bindex_t& idx=bindex_[v];
            if((isize=idx.pos.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.pos[i], sizeof(idx.pos[i]));
                }
            }
            if((isize=idx.dis.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.dis[i], sizeof(idx.dis[i]));
                }
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.r_dis[i], sizeof(idx.r_dis[i]));
                }
            }
            if((isize=idx.spt_v.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.spt_v[i], sizeof(idx.spt_v[i]));
                    ofs.write((const char*)&idx.spt_d[i], sizeof(idx.spt_d[i]));
                }
            } 
            if((isize=idx.r_spt_v.size())!=0){
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < isize; ++i) {
                    ofs.write((const char*)&idx.r_spt_v[i], sizeof(idx.r_spt_v[i]));
                    ofs.write((const char*)&idx.r_spt_d[i], sizeof(idx.r_spt_d[i]));
                }
            } 
		}
		ofs.close();
        std::cout<<"Save integrated_index labels to binary file successfully!"<<std::endl;
    }

    /**
     * @description: used to save rank
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void save_rank(char* outputDirName){
        //get rank and inv
        rank.resize(numOfVertices,0);
        inv.resize(numOfVertices,numOfVertices);
        for(NodeID v=0;v<numOfVertices;++v){
            rank[v]=numOfVertices-1-index_[v].level;
            //rank[v]=index_[v].level;//revised
            inv[rank[v]]=v;
        }
        string order_file_name(outputDirName);
        order_file_name.append("_all.order");
        ofstream ofs(order_file_name);
        if(!ofs.is_open()) cout<<"Cannort open "<<order_file_name<<endl;
        for (int i = 0; i < numOfVertices; ++i) {
			ofs << inv[i] << endl;
		}
        ofs.close();
        std::cout<<"Save integrated_index order to file successfully!"<<endl;
    }

    /**
     * @description: used to save rank,add for directed
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void save_rank_directed(char* outputDirName){
        //get rank and inv
        rank.resize(numOfVertices,0);
        inv.resize(numOfVertices,numOfVertices);
        for(NodeID v=0;v<numOfVertices;++v){
            rank[v]=numOfVertices-1-bindex_[v].level;
            //rank[v]=bindex_[v].level;//revised
            inv[rank[v]]=v;
        }
        string order_file_name(outputDirName);
        order_file_name.append("_all.order");
        ofstream ofs(order_file_name);
        if(!ofs.is_open()) cout<<"Cannort open "<<order_file_name<<endl;
        for (int i = 0; i < numOfVertices; ++i) {
			ofs << inv[i] << endl;
		}
        ofs.close();
        std::cout<<"Save integrated_index order to file successfully!"<<endl;
    }

    /**
     * @description: save label size used to compute query cost
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void save_label_size(char* labelSizeFileName){
        ofstream ofs(labelSizeFileName);
		if(!ofs.is_open()) {cerr<<"Cannot open "<<labelSizeFileName<<endl;}
        ofs.setf(ios::fixed);
        ofs.precision(4);
        for (int i = 0; i < numOfVertices; ++i) {
			ofs <<i<<" "<< index_[i].siz<< endl;
		}
        ofs.close();
        std::cout<<"Save integrated_index label size to file successfully!"<<std::endl;
    }

    /**
     * @description: save label size used to compute query cost,add for directed
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void save_label_size_directed(char* labelSizeFileName){
        ofstream ofs(labelSizeFileName);
		if(!ofs.is_open()) {cerr<<"Cannot open "<<labelSizeFileName<<endl;}
        ofs.setf(ios::fixed);
        ofs.precision(4);
        for (int i = 0; i < numOfVertices; ++i) {
			ofs <<i<<" "<< bindex_[i].siz<<" "<< bindex_[i].r_siz<< endl;
		}
        ofs.close();
        std::cout<<"Save integrated_index label size to file successfully!"<<std::endl;
    }

    /**
     * @description: read label size from file
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void load_label_size(char* load_filename){
        spt_v_num.resize(numOfVertices,0);
        ifstream in(load_filename);//input HFPoint file to ifstream
        if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
        double label_size=0; //each line representing the label size of each vertice
        char line[24];
        int i=0;
        //read each line representing HFpoint to vector 
        while (in.getline(line,sizeof(line)))
        {
            stringstream ls(line);
            ls>>label_size;
            spt_v_num[i++]=label_size;
        }
        in.close();
        return;
    }

    /**
     * @description: read label size from file, add for directed
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void load_label_size_directed(char* load_filename){
        spt_v_num.resize(numOfVertices,0);
        spt_v_num_f.resize(numOfVertices,0);
		spt_v_num_r.resize(numOfVertices,0);
		ifstream in(load_filename);//input HFPoint file to ifstream
		if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
		double isize_f=0,isize_r=0; //each line representing the label size of each vertice
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
        return;
    }

    /**
     * @description: compute query cost used bu queryPairFile, add for directed
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    double get_query_cost_directed(char* load_filename){
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
			double isize_s = spt_v_num_f[v];
			double isize_t = spt_v_num_r[v];
			//compute the query performance function
			double ratio_s=(double)queryTime_s_tmp[v]/(double)DIVISION_FACTOR;
			double ratio_t=(double)queryTime_t_tmp[v]/(double)DIVISION_FACTOR;
			performance_result+=ratio_s*isize_s+ratio_t*isize_t;
		}
		//std::cout<<"get_query_cost successfully!"<<std::endl;
		return performance_result;
    }

    /**
     * @description: compute query cost used bu queryPairFile
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
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
            double isize = spt_v_num[v];
            double ratio=(double)queryTime_tmp[v]/(double)DIVISION_FACTOR;
            performance_result+=ratio*isize;
        }
        //std::cout<<"get_query_cost successfully!"<<std::endl;
        return performance_result;
    }

    void save_label_size_byOrder(char* outputDirName){
        string label_size_file(outputDirName);
        label_size_file.append("_byOrder.size");
        ofstream ofs(label_size_file);
		if(!ofs.is_open()) {cerr<<"Cannot open "<<label_size_file<<endl;}
		for (int i = 0; i < numOfVertices; ++i) {
			ofs<<i<<" "<<index_[inv[i]].siz<< endl;
		}		
        ofs.close();
        std::cout<<"Save integrated_index label size byOrder to file successfully!"<<std::endl;
    }

    void write_labels(char* outputDirName){
        string label_list_file_name(outputDirName);
        label_list_file_name.append(".list");
        ofstream ofs(label_list_file_name);
        if(!ofs.is_open()) cout<<"Cannot open "<<label_list_file_name<<endl;
        for(NodeID v=0;v<numOfVertices;++v){
            ofs<<v<<":"<<endl;
            if(index_[v].spt_v.size()!=0)
            {
                ofs<<"spt:";
                for(size_t i=0;i<index_[v].spt_v.size();++i) ofs<<" ("<<index_[v].spt_v[i]<<","<<index_[v].spt_d[i]<<")";
                ofs<<endl;
            }
            if(index_[v].pos.size()!=0)
            {
                ofs<<"pos:";
                for(size_t i=0;i<index_[v].pos.size();++i) ofs<<" "<<index_[v].pos[i];
                ofs<<endl;
            }
            if(index_[v].dis.size()!=0)
            {
                ofs<<"dis:";
                for(size_t i=0;i<index_[v].dis.size();++i) ofs<<" "<<index_[v].dis[i];
                ofs<<endl;
            }
        }
        ofs.close();

        // string label_list_file_name1(outputDirName);
        // label_list_file_name1.append("_original.list");
        // ofstream ofs1(label_list_file_name1);
        // if(!ofs1.is_open()) cout<<"Cannot open "<<label_list_file_name<<endl;
        // ofs1.close();
        std::cout<<"Write label list file successfully!"<<std::endl;
    }

    /**
     * @description: write labels for the second solution
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void write_labels_1(char* outputDirName){
        string label_list_file_name(outputDirName);
        label_list_file_name.append(".list");
        ofstream ofs(label_list_file_name);
        if(!ofs.is_open()) cout<<"Cannot open "<<label_list_file_name<<endl;
        //pos dis spt-v spt-d
        for(NodeID v=0;v<numOfVertices;++v){
            ofs<<v<<":"<<endl;
            if(index_[v].pos.size()!=0)
            {
                ofs<<"pos:";
                for(size_t i=0;i<index_[v].pos.size();++i) ofs<<" "<<index_[v].pos[i];
                ofs<<endl;
            }
            if(index_[v].dis.size()!=0)
            {
                ofs<<"dis:";
                for(size_t i=0;i<index_[v].dis.size();++i) ofs<<" "<<index_[v].dis[i];
                ofs<<endl;
            }
            if(index_[v].spt_v.size()!=0)
            {
                ofs<<"spt-v:";
                for(size_t i=0;i<index_[v].spt_v.size();++i) ofs<<" "<<index_[v].spt_v[i];
                ofs<<endl;
            }
            if(index_[v].spt_d.size()!=0)
            {
                ofs<<"spt-d:";
                for(size_t i=0;i<index_[v].spt_d.size();++i) ofs<<" "<<index_[v].spt_d[i];
                ofs<<endl;
            }
        }
        ofs.close();
        std::cout<<"Write label list file successfully!"<<std::endl;
    }

    /**
     * @description: write labels for the second solution, add for directed
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    void write_labels_1_directed(char* outputDirName){
        string label_list_file_name(outputDirName);
        label_list_file_name.append(".list");
        ofstream ofs(label_list_file_name);
        if(!ofs.is_open()) cout<<"Cannot open "<<label_list_file_name<<endl;
        //pos dis spt-v spt-d
        for(NodeID v=0;v<numOfVertices;++v){
            ofs<<v<<":"<<endl;
            if(bindex_[v].pos.size()!=0)
            {
                ofs<<"pos:";
                for(size_t i=0;i<bindex_[v].pos.size();++i) ofs<<" "<<bindex_[v].pos[i];
                ofs<<endl;
            }
            if(bindex_[v].dis.size()!=0)
            {
                ofs<<"dis:";
                for(size_t i=0;i<bindex_[v].dis.size();++i) ofs<<" "<<bindex_[v].dis[i];
                ofs<<endl;
                ofs<<"r_dis:";
                for(size_t i=0;i<bindex_[v].r_dis.size();++i) ofs<<" "<<bindex_[v].r_dis[i];
                ofs<<endl;
            }
            if(bindex_[v].spt_v.size()!=0)
            {
                ofs<<"spt-v:";
                for(size_t i=0;i<bindex_[v].spt_v.size();++i) ofs<<" "<<bindex_[v].spt_v[i];
                ofs<<endl;
            }
            if(bindex_[v].spt_d.size()!=0)
            {
                ofs<<"spt-d:";
                for(size_t i=0;i<bindex_[v].spt_d.size();++i) ofs<<" "<<bindex_[v].spt_d[i];
                ofs<<endl;
                ofs<<"r-spt-v:";
                for(size_t i=0;i<bindex_[v].r_spt_v.size();++i) ofs<<" "<<bindex_[v].r_spt_v[i];
                ofs<<endl;
                ofs<<"r-spt-d:";
                for(size_t i=0;i<bindex_[v].r_spt_d.size();++i) ofs<<" "<<bindex_[v].r_spt_d[i];
                ofs<<endl;
            }
        }
        ofs.close();
        std::cout<<"Write label list file successfully!"<<std::endl;
    }

};

#endif
