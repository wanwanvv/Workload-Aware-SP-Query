/*
 * @Descripttion: tree decomposition adn h2h-index construction
 * @version: 
 * @Author: wanjingyi
 * @Date: 2021-02-09 15:33:58
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-02-16 14:48:13
 */
#ifndef  H2H_CONSTURCTION_H
#define H2H_CONSTURCTION_H

#include <vector>
#include <fstream>
#include <algorithm>
#include <cstdint> 
#include <iostream>
#include <list>
#include <queue>
#include <iomanip>
#include <math.h>
#include <numeric>
#include <unordered_map>
#include <omp.h>
#include "./heap.h"
#include "./paras.h"
#include "./utils.h"
#include "./time_util.h"
#include "./integrated_index.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define INF_WEIGHT SP_Constants::INF_WEIGHT

using namespace std;
using namespace time_util;

//compare firs by id then by weight
bool cmp_edge(const pair<int,int> a, const pair<int,int> b) {
	if(a.first!=b.first){
		return a.first<b.first;
	}else{
		return a.second<b.second;
	}
}

bool cmp_edge_1(const pair<NodeID,EdgeWeight> a, const pair<NodeID,EdgeWeight> b) {
	if(a.first!=b.first){
		return a.first<b.first;
	}else{
		return a.second<b.second;
	}
}

/*
 *@description: node structure class
 *@author: gaoyongyong
 *@date: 2021-02-09
*/
class node
{
public:
	int* _pos{NULL};//pos
	int* _dis{NULL};//dis
	vector<node*> _ans;//store the ancestors
	int* _neighbour{NULL};//vertices in a node\curNum with descending order
	node** _neighbourPointer{NULL};//vertice pointer in a node\curNum with descending order
	int* _neighbourDis{NULL};//curNum distances to other vertices in a node with descending order
	vector<node*>_childs;
	node* _parent{ NULL };
	int _height{ 0 };//root height=0
	int _rank{ 0 };//order values
	int _curNum=0;
	int _parentNum;
	int _width=0;//node width
	int _repr=-1;//A hack to allow the LCA implementation to inject an offset into a node.
	double _siz=0;//average labels size

	//Add for directed
	int* _r_pos{NULL};//pos for directed
	int* _r_dis{NULL};//dis for directed
	int* _r_neighbourDis{NULL};
	double _r_siz=0;//average labels size as termination

	int _dis_height{ 0 };//dis_height

	int _tree_id{0};
	//for parallel
	size_t _nBorder_siz{0};//lf:store the num of neighbors are hf;hf:store the numOftrees
	int* _nBorder_index{NULL};//hf point: store the index poses in each forest used for multiple threads
	//for parallel compound pll labels
	vector<NodeID> pll_v;
	vector<EdgeWeight> pll_d;
	vector<NodeID> pll_v_r;//add for directed
	vector<EdgeWeight> pll_d_r;//add for directed
	//*****************constructions and deconstructions*************
	node(){}
	node(int currNum) :_curNum(currNum){};
    ~node(){
		if(_pos!=NULL) delete[] _pos;
		if(_dis!=NULL) delete[] _dis;
		if(_neighbourDis!=NULL) delete[] _neighbourDis;
		if(_neighbour!=NULL) delete[] _neighbour;
		if(_neighbourPointer!=NULL) delete[] _neighbourPointer;
		if(_r_pos!=NULL) delete[] _r_pos;
		if(_r_dis!=NULL) delete[] _r_dis;
		if(_r_neighbourDis!=NULL) delete[] _r_neighbourDis;
    }
	//*******************functions***********************
	/**
  * @description: set each node's parent to build the tree
  * @Author: wanjingyi
  * @Date: 2021-02-10 20:52:19
  * @param {const} vector
  * @param {const} vector
  * @return {*}
  */
	void setParent(vector<node>& nodes,vector<pair<int,int> >& neighbourInfo_v,const vector<int>& rank){
		_width=neighbourInfo_v.size();
		_rank=rank[_curNum];
		//sort and store
		vector<pair<int,int> > rank_tmp(_width);
		_neighbour=new int[_width];
		_neighbourDis=new int[_width];
		_neighbourPointer=new node*[_width];
		for(int i=0;i<_width;++i) rank_tmp[i]=make_pair(rank[neighbourInfo_v[i].first],i);
		sort(rank_tmp.rbegin(),rank_tmp.rend());//sort
		for(int i=0;i<_width;++i){
			int k=rank_tmp[i].second;
			_neighbour[i]=neighbourInfo_v[k].first;
			_neighbourDis[i]=neighbourInfo_v[k].second;
			_neighbourPointer[i]=&nodes[_neighbour[i]];
		}
		_parentNum=_neighbour[_width-1];//the smallest order value
		_parent=&nodes[_parentNum];
		//set the current node as childs to its parent node
		_parent->_childs.push_back(this);
	}

	/**
	 * @description: set parent node for directed graph
	 * @param {const} vector
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void setParent_directed(vector<node>& nodes,vector<pair<int,int> >& r_neighbourInfo_v,vector<pair<int,int> >& neighbourInfo_v,const vector<int>& rank){
		_rank=rank[_curNum];
		//remove the duplicated neighbors
		int nbr,nbr_w;
		unordered_map<int,pair<int,int> > nbr_map;
		for(size_t i=0;i<r_neighbourInfo_v.size();++i){
			nbr=r_neighbourInfo_v[i].first;
			nbr_w=r_neighbourInfo_v[i].second;
			nbr_map[nbr]=make_pair(nbr_w,INF_WEIGHT);
		}
		for(size_t i=0;i<neighbourInfo_v.size();++i){
			nbr=neighbourInfo_v[i].first;
			nbr_w=neighbourInfo_v[i].second;
			if(nbr_map.find(nbr)==nbr_map.end()) nbr_map[nbr]=make_pair(INF_WEIGHT,nbr_w);
			else nbr_map[nbr].second=nbr_w;
		}
		_width=nbr_map.size();
		//sort the rank
		vector<pair<int,int> > rank_tmp;
		for(unordered_map<int,pair<int,int> >::iterator it=nbr_map.begin();it!=nbr_map.end();++it){
			rank_tmp.push_back(make_pair(rank[it->first],it->first));
		}
		sort(rank_tmp.rbegin(),rank_tmp.rend());
		//store the neighbors and relevant distances
		_neighbour=new int[_width];
		_neighbourDis=new int[_width];
		_r_neighbourDis=new int[_width];
		_neighbourPointer=new node*[_width];
		for(size_t i=0;i<_width;++i){
			nbr=rank_tmp[i].second;
			_neighbour[i]=nbr;
			_neighbourPointer[i]=&nodes[nbr];
			_r_neighbourDis[i]=nbr_map[nbr].first;
			_neighbourDis[i]=nbr_map[nbr].second;
		}
		//the smallest order value
		_parentNum=_neighbour[_width-1];
		//set the current node as childs to its parent node
		_parent=&nodes[_parentNum];
		_parent->_childs.push_back(this);

		//release memory
		//nbr_map.clear();
		//rank_tmp.clear();
	}

	/**
  * @description: used for h2h-index construction
  * @param {const vector<node*>} ans
  * @param {int} h
  * @return {*}
  * @author: Wan Jingyi
  */ 
 	void compute_diis(const vector<node*> ans,int h){
		//initialize variables
		_ans=ans;
		_height=h;
		_pos=new int[_width+1];
		_dis=new int[_height+1];
		int i=0,j=0;
		for(i=0;i<=_height;++i) _dis[i]=INF_WEIGHT;
		//top down manner
		for(i=0;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_height;
			for(j=0;j<_ans.size();++j){
				if(j<=_pos[i]){
					int dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}else{
					int dis=_ans[j]->_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[i]=h;
		_dis[h]=0;
		//to be deleted
		// cout<<"*****************"<<_curNum<<"***************"<<endl;
		// cout<<"w="<<_width<<" h="<<_height<<endl;
		// cout<<"ans: ";
		// for (size_t i = 0; i <_height; i++) cout<<_ans[i]->_curNum<<" ";
		// cout<<endl;
		// cout<<"pos: ";
		// for(size_t i=0;i<=_width;++i) cout<<_pos[i]<<" ";
		// cout<<endl;
		// cout<<"dis: ";
		// for(size_t i=0;i<=_height;++i) cout<<_dis[i]<<" ";
		// cout<<endl;
	}

	/**
	 * @description: used for h2h-index construction, add for directed
	 * @param {const vector<node*>} ans
	 * @param {int} h
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
	void compute_diis_directed(const vector<node*> ans,int h){
		//initialize variables
		_ans=ans;
		_height=h;
		_pos=new int[_width+1];
		_dis=new int[_height+1];
		_r_dis=new int[_height+1];
		int i=0,j=0;
		for(i=0;i<=_height;++i) _dis[i]=_r_dis[i]=INF_WEIGHT;
		//top down manner
		for(i=0;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_height;
			for(j=0;j<_ans.size();++j){
				if(j<=_pos[i]){
					//forward
					int dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					int r_dis=neg_node->_r_dis[j]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}else{
					//forward
					int dis=_ans[j]->_r_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					int r_dis=_ans[j]->_dis[_pos[i]]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[i]=h;
		_dis[h]=0;
		_r_dis[h]=0;
	}

	/**
	 * @description: used for integrated-index second solution 
	 * @param {const} vector
	 * @param {int} node_height
	 * @param {int} border_size
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void compute_diis(const vector<node*>& ans,int node_height,int border_size){
		//initialize variables
		//need modification
		_ans=ans;
		_ans.erase(_ans.begin(),_ans.begin()+border_size);
		_height=node_height;
		_dis_height=node_height+border_size;
		_pos=new int[_width+1];
		_dis=new int[_dis_height+1];
		int i=0,j=0;
		for(i=0;i<=_dis_height;++i) _dis[i]=INF_WEIGHT;
		//top down manner
		for(i=0;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					int dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}else{
					int dis=ans[j]->_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[_width]=_dis_height;
		_dis[_dis_height]=0;
	}

	/**
	 * @description: used for integrated-index second solution, add fordirected
	 * @param {const} vector
	 * @param {int} node_height
	 * @param {int} border_size
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void compute_diis_directed(const vector<node*>& ans,int node_height,int border_size){
		//initialize variables
		//need modification
		_ans=ans;
		_ans.erase(_ans.begin(),_ans.begin()+border_size);
		_height=node_height;
		_dis_height=node_height+border_size;
		_pos=new int[_width+1];
		_dis=new int[_dis_height+1];
		_r_dis=new int[_dis_height+1];
		int i=0,j=0;
		for(i=0;i<=_dis_height;++i) _dis[i]=_r_dis[i]=INF_WEIGHT;
		//top down manner
		int dis,r_dis;
		for(i=0;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					//forward
					dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=neg_node->_r_dis[j]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}else{
					//forward
					dis=ans[j]->_r_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=ans[j]->_dis[_pos[i]]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[_width]=_dis_height;
		_dis[_dis_height]=0;
		_r_dis[_dis_height]=0;
	}

	/**
	 * @description: used to compute dis position for the second solution in multithreads
	 * @param {int} tree_index
	 * @param {const} vector
	 * @param {int} node_height
	 * @param {int} border_size
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void compute_diis_multiThreads(int tree_index,const vector<node*>& ans,int node_height,int border_size){
		//initialize variables
		//need modification
		_ans=ans;
		_ans.erase(_ans.begin(),_ans.begin()+border_size);
		_height=node_height;
		_dis_height=node_height+border_size;
		_pos=new int[_width+1];
		_dis=new int[_dis_height+1];
		int i=0,j=0;
		for(i=0;i<=_dis_height;++i) _dis[i]=INF_WEIGHT;
		//top down manner
		for(i=0;i<_nBorder_siz;++i){
			node* neg_node=ans[_neighbourPointer[i]->_nBorder_index[tree_index]];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					int dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}else{
					int dis=ans[j]->_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}
			}
		}
		for(;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					int dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}else{
					int dis=ans[j]->_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[_width]=_dis_height;
		_dis[_dis_height]=0;
	}

	/**
	 * @description: used to compute distance array for the second solution in mutiple threads, add for directed
	 * @param {int} tree_index
	 * @param {const} vector
	 * @param {int} node_height
	 * @param {int} border_size
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void compute_diis_multiThreads_directed(int tree_index,const vector<node*>& ans,int node_height,int border_size){
		//initialize variables
		//need modification
		_ans=ans;
		_ans.erase(_ans.begin(),_ans.begin()+border_size);
		_height=node_height;
		_dis_height=node_height+border_size;
		_pos=new int[_width+1];
		_dis=new int[_dis_height+1];
		_r_dis=new int[_dis_height+1];
		int i=0,j=0;
		for(i=0;i<=_dis_height;++i) _dis[i]=_r_dis[i]=INF_WEIGHT;
		//top down manner
		int dis,r_dis;
		for(i=0;i<_nBorder_siz;++i){
			node* neg_node=ans[_neighbourPointer[i]->_nBorder_index[tree_index]];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					//forward
					dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=neg_node->_r_dis[j]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}else{
					//forward
					dis=ans[j]->_r_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=ans[j]->_dis[_pos[i]]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}
			}
		}
		for(;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					//forward
					dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=neg_node->_r_dis[j]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}else{
					//forward
					dis=ans[j]->_r_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=ans[j]->_dis[_pos[i]]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[_width]=_dis_height;
		_dis[_dis_height]=0;
		_r_dis[_dis_height]=0;
	}

	/**
	 * @description: single thread for the third solution
	 * @param {const} vector
	 * @param {const} vector
	 * @param {int} node_height
	 * @param {int} border_size
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void compute_diis_and_pushPLL(Integrated_Index& integrated_index,const vector<int>& _tree_borders,const vector<node*>& ans,int node_height,int border_size){
		//initialize variables
		vector<integrated_index_t>& index_=integrated_index.index_;
		//need modification
		_ans=ans;
		_ans.erase(_ans.begin(),_ans.begin()+border_size);
		_height=node_height;
		_dis_height=node_height+border_size;
		_pos=new int[_width+1];
		_dis=new int[_dis_height+1];
		int i=0,j=0;
		for(i=0;i<=_dis_height;++i) _dis[i]=INF_WEIGHT;
		//top down manner
		for(i=0;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					int dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}else{
					int dis=ans[j]->_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[_width]=_dis_height;
		_dis[_dis_height]=0;
		//push down pll labels
		integrated_index_t& idx_curr=index_[_curNum];
		//***************BucketSort***************
		//double time_push=GetCurrentTimeSec();//for debug
		NodeID spt_v_t;
		EdgeWeight spt_d_t;
		size_t max_length=integrated_index._max_overlay_hub+1;
		vector<EdgeWeight> buckets(max_length,INF_WEIGHT);
		for(int i=0;i<border_size;++i){
			const integrated_index_t& idx_border=index_[_tree_borders[i]];
			for(int j=0;j<idx_border.spt_v.size()-1;++j){
				spt_v_t=idx_border.spt_v[j];
				spt_d_t=idx_border.spt_d[j]+_dis[i];
				if(spt_d_t<buckets[spt_v_t]) buckets[spt_v_t]=spt_d_t;
			}
		}
		idx_curr.spt_v.reserve(max_length+1);
		idx_curr.spt_d.reserve(max_length+1);
		for(NodeID i=0;i<max_length;++i){
			if((spt_d_t=buckets[i])!=INF_WEIGHT){
				idx_curr.spt_v.push_back(i);
				idx_curr.spt_d.push_back(spt_d_t);
			}
		}
		idx_curr.spt_v.push_back(numOfVertices);
		idx_curr.spt_d.push_back(INF_WEIGHT);
		//_siz=(double)(idx_curr.spt_v.size()-1);//temp store the spt_v size
		//time_push=GetCurrentTimeSec()-time_push;//for debug
		//***************direct sort*****************
		// vector<pair<NodeID,EdgeWeight> > tmp_idx_curr;
		// for(int i=0;i<border_size;++i){
		// 	const integrated_index_t& idx_border=index_[_tree_borders[i]];
		// 	for(int j=0;j<idx_border.spt_v.size()-1;++j){
		// 		tmp_idx_curr.push_back(make_pair(idx_border.spt_v[j],idx_border.spt_d[j]+_dis[i]));
		// 	}
		// }
		// size_t pll_cnt=tmp_idx_curr.size();
		// sort(tmp_idx_curr.begin(),tmp_idx_curr.end(),cmp_edge_1);
		// idx_curr.spt_v.reserve(pll_cnt);
		// idx_curr.spt_d.reserve(pll_cnt);
		// NodeID last_spt_v=numOfVertices,spt_v_t;
		// i=0;
		// while(i<pll_cnt){
		// 	for(j=i;j<pll_cnt&&(spt_v_t=tmp_idx_curr[j].first)==last_spt_v;++j);
		// 	if(j==pll_cnt) break;
		// 	idx_curr.spt_v.push_back(spt_v_t);
		// 	idx_curr.spt_d.push_back(tmp_idx_curr[j].second);
		// 	last_spt_v=spt_v_t;
		// 	i=++j;
		// }
	}


	/**
	 * @description: Single thread for the third solution, add for directed
	 * @param {const} vector
	 * @param {const} vector
	 * @param {int} node_height
	 * @param {int} border_size
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void compute_diis_and_pushPLL_directed(Integrated_Index& integrated_index,const vector<int>& _tree_borders,const vector<node*>& ans,int node_height,int border_size){
		//initialize variables
		vector<integrated_bindex_t>& bindex_=integrated_index.bindex_;
		//need modification
		_ans=ans;
		_ans.erase(_ans.begin(),_ans.begin()+border_size);
		_height=node_height;
		_dis_height=node_height+border_size;
		_pos=new int[_width+1];
		_dis=new int[_dis_height+1];
		_r_dis=new int[_dis_height+1];
		int i=0,j=0;
		for(i=0;i<=_dis_height;++i) _dis[i]=_r_dis[i]=INF_WEIGHT;
		//top down manner
		int dis,r_dis;
		for(i=0;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					//forward
					dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=neg_node->_r_dis[j]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}else{
					//forward
					dis=ans[j]->_r_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=ans[j]->_dis[_pos[i]]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[_width]=_dis_height;
		_dis[_dis_height]=0;
		_r_dis[_dis_height]=0;
		//push down pll labels
		integrated_bindex_t& idx_curr=bindex_[_curNum];
		//***************BucketSort***************//
		NodeID spt_v_t,r_spt_v_t;
		EdgeWeight spt_d_t,r_spt_d_t;
		size_t max_length=integrated_index._max_overlay_hub+1;
		size_t max_length_r=integrated_index._max_overlay_hub_r+1;
		vector<EdgeWeight> buckets(max_length,INF_WEIGHT);
		vector<EdgeWeight> buckets_r(max_length_r,INF_WEIGHT);
		for(int i=0;i<border_size;++i){
			const integrated_bindex_t& idx_border=bindex_[_tree_borders[i]];
			//forward
			if(_dis[i]!=INF_WEIGHT){
				for(int j=0;j<idx_border.spt_v.size()-1;++j){
					spt_v_t=idx_border.spt_v[j];
					spt_d_t=idx_border.spt_d[j]+_dis[i];
					if(spt_d_t<buckets[spt_v_t]) buckets[spt_v_t]=spt_d_t;
				}
			}
			//backward
			if(_r_dis[i]!=INF_WEIGHT){
				for(int j=0;j<idx_border.r_spt_v.size()-1;++j){
					r_spt_v_t=idx_border.r_spt_v[j];
					r_spt_d_t=idx_border.r_spt_d[j]+_r_dis[i];
					if(r_spt_d_t<buckets_r[r_spt_v_t]) buckets_r[r_spt_v_t]=r_spt_d_t;
				}
			}
		}
		//forward
		idx_curr.spt_v.reserve(max_length+1);
		idx_curr.spt_d.reserve(max_length+1);
		for(NodeID i=0;i<max_length;++i){
			if((spt_d_t=buckets[i])!=INF_WEIGHT){
				idx_curr.spt_v.push_back(i);
				idx_curr.spt_d.push_back(spt_d_t);
			}
		}
		//backward
		idx_curr.r_spt_v.reserve(max_length_r+1);
		idx_curr.r_spt_d.reserve(max_length_r+1);
		for(NodeID i=0;i<max_length_r;++i){
			if((r_spt_d_t=buckets_r[i])!=INF_WEIGHT){
				idx_curr.r_spt_v.push_back(i);
				idx_curr.r_spt_d.push_back(r_spt_d_t);
			}
		}
		//end sentinal forward and backward
		idx_curr.spt_v.push_back(numOfVertices);
		idx_curr.spt_d.push_back(INF_WEIGHT);
		idx_curr.r_spt_v.push_back(numOfVertices);
		idx_curr.r_spt_d.push_back(INF_WEIGHT);
		//forward and backward size
		//_r_siz=(double)(idx_curr.r_spt_v.size()-1);
		//_siz=(double)(idx_curr.spt_v.size()-1);
		return;
	}

	//Multiple threads for the third solution
	void compute_diis_and_pushPLL_multiThreads(int tree_index,const vector<integrated_index_t>& index_,const vector<int>& _tree_borders,const vector<node*>& ans,int node_height,int border_size,int max_overlay_hub){
		//initialize variables
		//need modification
		_ans=ans;
		_ans.erase(_ans.begin(),_ans.begin()+border_size);
		_height=node_height;
		_dis_height=node_height+border_size;
		_pos=new int[_width+1];
		_dis=new int[_dis_height+1];
		int i=0,j=0;
		for(i=0;i<=_dis_height;++i) _dis[i]=INF_WEIGHT;
		//top down manner
		for(i=0;i<_nBorder_siz;++i){
			node* neg_node=ans[_neighbourPointer[i]->_nBorder_index[tree_index]];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					int dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}else{
					int dis=ans[j]->_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}
			}
		}
		for(;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					int dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}else{
					int dis=ans[j]->_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[_width]=_dis_height;
		_dis[_dis_height]=0;
		NodeID spt_v_t;
		EdgeWeight spt_d_t;
		size_t max_length=max_overlay_hub+1;
		vector<EdgeWeight> buckets(max_length,INF_WEIGHT);
		for(int i=0;i<border_size;++i){
			const integrated_index_t& idx_border=index_[_tree_borders[i]];
			for(int j=0;j<idx_border.spt_v.size()-1;++j){
				spt_v_t=idx_border.spt_v[j];
				spt_d_t=idx_border.spt_d[j]+_dis[i];
				if(spt_d_t<buckets[spt_v_t]) buckets[spt_v_t]=spt_d_t;
			}
		}
		pll_v.reserve(max_length+1);
		pll_d.reserve(max_length+1);
		for(NodeID i=0;i<max_length;++i){
			if((spt_d_t=buckets[i])!=INF_WEIGHT){
				pll_v.push_back(i);
				pll_d.push_back(spt_d_t);
			}
		}
		pll_v.push_back(numOfVertices);
		pll_d.push_back(INF_WEIGHT);
		_siz=(double)(pll_v.size()-1);//temp store the spt_v size
	}

	//Multiple threads for the third solution, add for directed
	void compute_diis_and_pushPLL_multiThreads_directed(int tree_index,const vector<integrated_bindex_t>& bindex_,const vector<int>& _tree_borders,const vector<node*>& ans,int node_height,int border_size,int max_overlay_hub_r,int max_overlay_hub){
		//initialize variables
		//need modification
		_ans=ans;
		_ans.erase(_ans.begin(),_ans.begin()+border_size);
		_height=node_height;
		_dis_height=node_height+border_size;
		_pos=new int[_width+1];
		_dis=new int[_dis_height+1];
		_r_dis=new int[_dis_height+1];
		int i=0,j=0;
		for(i=0;i<=_dis_height;++i) _dis[i]=_r_dis[i]=INF_WEIGHT;
		//top down manner
		int dis,r_dis;
		for(i=0;i<_nBorder_siz;++i){
			node* neg_node=ans[_neighbourPointer[i]->_nBorder_index[tree_index]];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					//forward
					dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=neg_node->_r_dis[j]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}else{
					//forward
					dis=ans[j]->_r_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=ans[j]->_dis[_pos[i]]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}
			}
		}
		for(;i<_width;++i){
			node* neg_node=_neighbourPointer[i];
			//pos
			_pos[i]=neg_node->_dis_height;
			//dis
			for(j=0;j<_dis_height;++j){
				if(j<=_pos[i]){
					//forward
					dis=neg_node->_dis[j]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=neg_node->_r_dis[j]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}else{
					//forward
					dis=ans[j]->_r_dis[_pos[i]]+_neighbourDis[i];
					if(dis<_dis[j]) _dis[j]=dis;
					//backward
					r_dis=ans[j]->_dis[_pos[i]]+_r_neighbourDis[i];
					if(r_dis<_r_dis[j]) _r_dis[j]=r_dis;
				}
			}
		}
		//_curNum dis and pos
		_pos[_width]=_dis_height;
		_dis[_dis_height]=0;
		_r_dis[_dis_height]=0;
		NodeID spt_v_t,r_spt_v_t;
		EdgeWeight spt_d_t,r_spt_d_t;
		size_t max_length_r=max_overlay_hub_r+1;
		size_t max_length=max_overlay_hub+1;
		vector<EdgeWeight> buckets_r(max_length_r,INF_WEIGHT);
		vector<EdgeWeight> buckets(max_length,INF_WEIGHT);
		for(int i=0;i<border_size;++i){
			const integrated_bindex_t& idx_border=bindex_[_tree_borders[i]];
			//forward
			if(_dis[i]!=INF_WEIGHT){
				for(int j=0;j<idx_border.spt_v.size()-1;++j){
					spt_v_t=idx_border.spt_v[j];
					spt_d_t=idx_border.spt_d[j]+_dis[i];
					if(spt_d_t<buckets[spt_v_t]) buckets[spt_v_t]=spt_d_t;
				}
			}
			//backward
			if(_r_dis[i]!=INF_WEIGHT){
				for(int j=0;j<idx_border.r_spt_v.size()-1;++j){
					r_spt_v_t=idx_border.r_spt_v[j];
					r_spt_d_t=idx_border.r_spt_d[j]+_r_dis[i];
					if(r_spt_d_t<buckets_r[r_spt_v_t]) buckets_r[r_spt_v_t]=r_spt_d_t;
				}
			}
		}
		//forward
		pll_v.reserve(max_length+1);
		pll_d.reserve(max_length+1);
		for(NodeID i=0;i<max_length;++i){
			if((spt_d_t=buckets[i])!=INF_WEIGHT){
				pll_v.push_back(i);
				pll_d.push_back(spt_d_t);
			}
		}
		pll_v.push_back(numOfVertices);
		pll_d.push_back(INF_WEIGHT);
		//backward
		pll_v_r.reserve(max_length_r+1);
		pll_d_r.reserve(max_length_r+1);
		for(NodeID i=0;i<max_length_r;++i){
			if((r_spt_d_t=buckets_r[i])!=INF_WEIGHT){
				pll_v_r.push_back(i);
				pll_d_r.push_back(r_spt_d_t);
			}
		}
		pll_v_r.push_back(numOfVertices);
		pll_d_r.push_back(INF_WEIGHT);
		//_siz=(double)(pll_v.size()-1);//temp store the spt_v size
	}

};

/*
 *@description: tree structure class
 *@author: gaoyongyong
 *@date: 2021-02-09
*/
class tree
{
public:
	vector<node> _nodes;
	int _numOfNodes;
	node* _root;
	vector<int> _euler;
	vector<int> _level;
	vector<int> _index;//get the index by node num
	vector<int> _start;//store the dfs start pos
	vector<int> _end;//store the dfs end pos
	vector<int> _dfs;//tmp list to store the dfs id list for count
	vector<double> _labelSize;//store the 2-hop labels size
	int _tree_height;
	int _max_node_num;
	int _max_node_width;

	//*************constructions**************
	tree(){
		_numOfNodes=0;
		_root=NULL;
		_tree_height=-1;
		_max_node_num=-1;
		_max_node_width=-1;
	}
	
	tree(int n){
		_numOfNodes=n;
		_nodes.resize(_numOfNodes);
		_root=NULL;
		//initialize
		for(int i=0;i<_numOfNodes;++i){
			node _node(i);
			_nodes[i]=_node;
		}
	}
	//*************deconstructions****************
	~tree(){

	}
	void initialize(int n){
		_numOfNodes=n;
		_nodes.resize(_numOfNodes);
		_root=NULL;
		//initialize
		for(int i=0;i<_numOfNodes;++i){
			node _node(i);
			_nodes[i]=_node;
		}
	}

	void compute_rmq_variables()
	{
		int _n=numOfVertices*2-1;
		_euler.reserve(_n);
        _level.reserve(_n);
        _index.resize(numOfVertices,-1);
		dfs_preprocess(_root);
	}

	/**
	 * @Author: wanjingyi
	 * @description: used for euler and rmq
	 * @param {*}
	 * @return {*}
	 */ 
	void dfs_preprocess(node* curr_node){
        int id=curr_node->_curNum;
        _euler.push_back(id);
        _level.push_back(curr_node->_height);
        if(_index[id]==-1) _index[id]=_euler.size()-1;
        for(size_t i=0;i<curr_node->_childs.size();++i)
        {
            dfs_preprocess(curr_node->_childs[i]);
            _euler.push_back(id);
            _level.push_back(curr_node->_height);
        }
    }

	/**
	 * @Author: wanjingyi
	 * @description: used for label size and query cost computation
	 * @param {*}
	 * @return {*}
	 */ 
	void dfs_compute(node* curr_node){
		int id=curr_node->_curNum;
		_dfs.push_back(id);
		_start[id]=_dfs.size()-1;
		for(size_t i=0;i<curr_node->_childs.size();++i){
			dfs_compute(curr_node->_childs[i]);
		}
		_dfs.push_back(id);
		_end[id]=_dfs.size()-1;
	}

	void compute_node_label_size(){
		int ans_num,tree_num,left_num,ans_total_size;
		double label_size;
		//compute |T(vi)|
		_dfs.reserve(_numOfNodes);
		_start.resize(_numOfNodes);
		_end.resize(_numOfNodes);
		_labelSize.resize(_numOfNodes,0);
		dfs_compute(_root);
		vector<int> ().swap(_dfs);
		_dfs.clear();
		for(int v=0;v<_numOfNodes;++v){
			node* curr_node=&_nodes[v];
			ans_num=curr_node->_ans.size();
			if(ans_num==0) label_size=1;//root node
			else{
				label_size=0;
				ans_total_size=0;
				tree_num=(_end[v]-_start[v]-1)/2+1;
				left_num=_numOfNodes-ans_num-tree_num;
				for(int i=0;i<ans_num;++i){
					ans_total_size+=curr_node->_ans[i]->_width+1;
				}
				label_size=ans_num+tree_num+((double)ans_total_size/(double)ans_num)*left_num;
				label_size=label_size/(double)_numOfNodes;
			}
			curr_node->_siz=label_size;
			_labelSize[v]=label_size;
		}
		std::cout<<"Compute node label size successfully!"<<endl;
	}

	void analysis_and_save(char* save_fileDirname,vector<unsigned int>& query_time,vector<bool>& HFPointFlag,int numOfHfpoint,int max_query_time,long long total_query_time){
		//output each node query cost format: v query_cost
		string avgQueryFile(save_fileDirname);
		avgQueryFile.append("_h2h.cost");
		ofstream ofs_avg(avgQueryFile);
		if(!ofs_avg.is_open()) cout<<"Cannot open "<<avgQueryFile<<endl;
		ofs_avg.precision(6);
		//compute |T(vi)|
		_dfs.reserve(_numOfNodes);
		_start.resize(_numOfNodes);
		_end.resize(_numOfNodes);
		_labelSize.resize(_numOfNodes,0);
		dfs_compute(_root);
		vector<int> ().swap(_dfs);
		_dfs.clear();
		//compute query cost
		double total_ave_size=0.0,hf_ave_size=0.0,hf_sum_size=0.0,ratio_s,ratio;
		double total_performance=0,total_performance_s=0, ave_performance=0,ave_performance_s=0;//total performance function 
		double label_size;
		double total_sum_size=0;
		int ans_num,tree_num,left_num,ans_total_size;
		//bool isFirst=true;
		for(int v=0;v<_numOfNodes;++v){
			node* curr_node=&_nodes[v];
			ans_num=curr_node->_ans.size();
			if(ans_num==0){
				label_size=1;//root node
			}
			else{
				label_size=0;
				ans_total_size=0;
				tree_num=(_end[v]-_start[v]-1)/2+1;
				left_num=_numOfNodes-ans_num-tree_num;
				for(int i=0;i<ans_num;++i){
					ans_total_size+=curr_node->_ans[i]->_width+1;
				}
				label_size=ans_num+tree_num+((double)ans_total_size/(double)ans_num)*left_num;
				label_size=label_size/(double)_numOfNodes;
			}
			curr_node->_siz=label_size;
			total_sum_size+=label_size;
			_labelSize[v]=label_size;
			//cout<<v<<" ans_num"<<ans_num<<" tree_num="<<tree_num<<" left_num="<<left_num<<" label_size="<<label_size<<endl;//to be deleted
			if(HFPointFlag[v]) hf_sum_size+=label_size;
			ratio_s=(double)query_time[v]/(double)DIVISION_FACTOR;
			ratio=(double)query_time[v]/((double)(max_query_time*DIVISION_FACTOR));
			ave_performance=ratio*label_size;
			ave_performance_s=ratio_s*label_size;
			total_performance+=ave_performance;
			//cout<<v<<" ave_performance="<<ave_performance<<" total_performance="<<total_performance<<endl;//to be deleted
			//if(query_time[v]!=0) cout<<v<<" ratio="<<ratio<<" query_time="<<query_time[v]<<" total_query_time="<<total_query_time<<" ave_performance="<<ave_performance<<" total_performance="<<total_performance<<" ave_performance_s="<<ave_performance_s<<" total_performance_s="<<total_performance_s<<endl;
			ofs_avg<<v<<" "<<label_size<<" "<<ave_performance_s<<endl;
		}
		ofs_avg.close();
		total_ave_size= (double) total_sum_size/(double) _numOfNodes;
		hf_ave_size= hf_sum_size/(double) numOfHfpoint;
		total_performance_s=total_performance*((double)max_query_time);
		cout<<"_numOfNodes = "<<_numOfNodes<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<endl;
		cout<<"numOfHFpoint = "<<numOfHfpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
		cout<<"nomalization performance_result = "<<total_performance<<endl;
		cout<<"standard performance_result = "<<total_performance_s<<endl;
		//output all query cost
		string allQueryFile(save_fileDirname);
		allQueryFile.append("_analysis.cost");
		ofstream ofs_all(allQueryFile);
		if(!ofs_all.is_open()) cout<<"Cannot open "<<allQueryFile<<endl;
		//write performance function to file
		ofs_all<<"_numOfNodes = "<<_numOfNodes<<" total_ave_size = "<<total_ave_size<<endl;
		ofs_all<<"numOfHFpoint = "<<numOfHfpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
		ofs_all<<"nomalization performance_result = "<<total_performance<<endl;
		ofs_all<<"standard performance_result = "<<total_performance_s<<endl;
		if(_tree_height!=-1){
			ofs_all<<"tree_height = "<<_tree_height<<" max_node_num = "<<_max_node_num<<" max_node_width = "<<_max_node_width<<endl;
		}
		ofs_all.close();
		//output size for debug
		string treeLabelSizeFile(save_fileDirname);
		treeLabelSizeFile.append("_tree_ordered.size");
		ofstream ofs_tree_size(treeLabelSizeFile);
		if(!ofs_tree_size.is_open()) cout<<"Cannot open "<<treeLabelSizeFile<<endl;
		vector<pair<double,int> > tree_size(_numOfNodes);
		for(int v=0;v<_numOfNodes;++v){
			tree_size[v]=make_pair(_labelSize[v],v);
		}
		sort(tree_size.begin(),tree_size.end());
		for(int r=0;r<_numOfNodes;++r){
			ofs_tree_size<<tree_size[r].second<<" "<<tree_size[r].first<<" "<<_nodes[tree_size[r].second]._height<<endl;
		}
		ofs_tree_size.close();
	}

	void save_node_label_size(char* outputDirName){
		compute_node_label_size();
		string label_size_file(outputDirName);
		label_size_file.append(".size");
		ofstream ofs(label_size_file);
		if(!ofs.is_open()) {cerr<<"Cannot open "<<label_size_file<<endl;}
		ofs.precision(4);
		for(int v=0;v<_numOfNodes;++v){
			ofs << v<<" "<<_labelSize[v]<< endl;
		}
		ofs.close();
		cout<<"save_node_label_size successfully!"<<endl;
	}

	void merge_h2h_index(vector<integrated_index_t>& index_,vector<bool>& isDeleted){
		for(int v=0;v<numOfVertices;++v){
			int i,k;
			const node& curr_node=_nodes[v];
			integrated_index_t& idx=index_[v];
			if(!isDeleted[v]){
				//pos
				k=curr_node._width+1;//num of pos
				idx.pos.resize(k*2+1);
				idx.pos[0]=k;
				int s_pos=1,s_id=k+1;
				idx.pos[s_pos++]=curr_node._pos[k-1];
				idx.pos[s_id++]=curr_node._curNum;
				for(i=0;i<k-1;++i){
					int p=curr_node._pos[k-2-i];
					idx.pos[s_pos+i]=p;
					idx.pos[s_id+i]=curr_node._ans[p]->_curNum;
				}

			}else{
				//pos
				k=curr_node._width+1;
				idx.pos.resize(k+1);
				idx.pos[0]=k;
				for(i=k-1;i>=0;--i){
					idx.pos[k-i]=curr_node._pos[i];
				}
				//dis
				k=curr_node._height;
				idx.dis.resize(k+1);
				for(i=0;i<=k;++i){
					idx.dis[i]=curr_node._dis[i];
				}
				//level
				idx.level=curr_node._rank;
				//size
				idx.siz=curr_node._siz;
			}

		}
		std::cout<<"Merge h2h-index successfully!"<<std::endl;
	}

	void save_lca_variables(char* outputDirName){
		string lca_filename(outputDirName);
		lca_filename.append(".lca");
		ofstream out(lca_filename,ios::binary | ios::out);
		if(!out.is_open()) cout<<"Cannot open "<<lca_filename<<endl;
		//output _euler
		int esize=_euler.size();
		out.write((const char*)&esize,sizeof(esize));
		//cout<<"_euler:"<<endl;
		for(int i=0;i<esize;++i){
			out.write((const char*)&_euler[i],sizeof(_euler[i]));
			//cout<<_euler[i]<<" ";
		}
		//cout<<endl;
		//output _level
		int lsize=_level.size();
		out.write((const char*)&lsize,sizeof(lsize));
		//cout<<"_level:"<<endl;
		for(int i=0;i<lsize;++i){
			out.write((const char*)&_level[i],sizeof(_level[i]));
			//cout<<_level[i]<<" ";
		}
		//cout<<endl;
		//output _index
		int isize=_index.size();
		out.write((const char*)&isize,sizeof(isize));
		//cout<<"_index:"<<endl;
		for(int i=0;i<isize;++i){
			out.write((const char*)&_index[i],sizeof(_index[i]));
			//cout<<_index[i]<<" ";
		}
		//cout<<endl;
		out.close();
	}

	void save_index_and_lca_variables(char* save_fileDirname){
		string labelFile(save_fileDirname);
		labelFile.append(".h2hIndex");
		ofstream ofs(labelFile,ios::binary | ios::out);
		if(!ofs.is_open()) cout<<"Cannot open "<<labelFile<<endl;
		ofs.write((const char*)&_numOfNodes,sizeof(_numOfNodes));
		for(int v=0;v<_numOfNodes;++v){
			node& curr_node=_nodes[v];
			// cout<<v<<":"<<endl;//for debug
			// cout<<"pos:"<<endl;//for debug
			//pos
			ofs.write((const char*)&curr_node._width,sizeof(curr_node._width));
			for(int i=0;i<=curr_node._width;++i){
				ofs.write((const char*)&curr_node._pos[i],sizeof(curr_node._pos[i]));
				// cout<<curr_node._pos[i]<<" ";//for debug
			}
			//dis
			// cout<<endl;//for debug
			// cout<<"dis:"<<endl;//for debug
			ofs.write((const char*)&curr_node._height,sizeof(curr_node._height));
			for(int i=0;i<=curr_node._height;++i){
				ofs.write((const char*)&curr_node._dis[i],sizeof(curr_node._dis[i]));
				// cout<<curr_node._dis[i]<<" ";//for debug
			}
			// cout<<endl;//for debug
		}
		ofs.close();

		string lcaFile(save_fileDirname);
		lcaFile.append(".lca");
		ofstream out(lcaFile,ios::binary | ios::out);
		if(!out.is_open()) cout<<"Cannot open "<<lcaFile<<endl;
		//output _euler
		int esize=_euler.size();
		out.write((const char*)&esize,sizeof(esize));
		for(int i=0;i<esize;++i) out.write((const char*)&_euler[i],sizeof(_euler[i]));
		//output _level
		int lsize=_level.size();
		out.write((const char*)&lsize,sizeof(lsize));
		for(int i=0;i<lsize;++i) out.write((const char*)&_level[i],sizeof(_level[i]));
		//output _index
		int isize=_index.size();
		out.write((const char*)&isize,sizeof(isize));
		for(int i=0;i<isize;++i) out.write((const char*)&_index[i],sizeof(_index[i]));
		out.close();
		
	}

	void save_index_and_lca_variables_directed(char* save_fileDirname){
		//*********************save h2h-index label********************//
		string labelFile(save_fileDirname);
		labelFile.append(".h2hIndex");
		ofstream ofs(labelFile,ios::binary | ios::out);
		if(!ofs.is_open()) cout<<"Cannot open "<<labelFile<<endl;
		ofs.write((const char*)&_numOfNodes,sizeof(_numOfNodes));
		for(int v=0;v<_numOfNodes;++v){
			node& curr_node=_nodes[v];
			ofs.write((const char*)&curr_node._width,sizeof(curr_node._width));
			//pos
			for(int i=0;i<=curr_node._width;++i){
				ofs.write((const char*)&curr_node._pos[i],sizeof(curr_node._pos[i]));
				// cout<<curr_node._pos[i]<<" ";//for debug
			}
			//r_dis
			ofs.write((const char*)&curr_node._height,sizeof(curr_node._height));
			for(int i=0;i<=curr_node._height;++i){
				ofs.write((const char*)&curr_node._r_dis[i],sizeof(curr_node._r_dis[i]));
				// cout<<curr_node._dis[i]<<" ";//for debug
			}
			//dis
			//ofs.write((const char*)&curr_node._height,sizeof(curr_node._height));
			for(int i=0;i<=curr_node._height;++i){
				ofs.write((const char*)&curr_node._dis[i],sizeof(curr_node._dis[i]));
				// cout<<curr_node._dis[i]<<" ";//for debug
			}
		}
		ofs.close();

		//*****************************save lca*********************//
		string lcaFile(save_fileDirname);
		lcaFile.append(".lca");
		ofstream out(lcaFile,ios::binary | ios::out);
		if(!out.is_open()) cout<<"Cannot open "<<lcaFile<<endl;
		//output _euler
		int esize=_euler.size();
		out.write((const char*)&esize,sizeof(esize));
		for(int i=0;i<esize;++i) out.write((const char*)&_euler[i],sizeof(_euler[i]));
		//output _level
		int lsize=_level.size();
		out.write((const char*)&lsize,sizeof(lsize));
		for(int i=0;i<lsize;++i) out.write((const char*)&_level[i],sizeof(_level[i]));
		//output _index
		int isize=_index.size();
		out.write((const char*)&isize,sizeof(isize));
		for(int i=0;i<isize;++i) out.write((const char*)&_index[i],sizeof(_index[i]));
		out.close();

	}

	void compute_tree_top_down(){
		//variables
		vector<vector<node*> > ans_list;
		int height=1;
		node* curr_node;
		int k,new_k;
		//bfs search
		queue<pair<node*,int> > nqueue;
		ans_list.push_back(vector<node*>(1,_root));
		for(size_t i=0;i<_root->_childs.size();++i) nqueue.push(make_pair(_root->_childs[i],0));
		nqueue.push(make_pair(nullptr,-1)); 
		while(!nqueue.empty()){
			if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
			vector<vector<node*> > ans_tmp;
			while(nqueue.front().first!=nullptr){
				pair<node*,int> front_pair=nqueue.front();
				curr_node=front_pair.first;
				k=front_pair.second;
				nqueue.pop();
				vector<node*> curr_ans=ans_list[k];
				curr_node->compute_diis(curr_ans,height);
				curr_ans.push_back(curr_node);
				ans_tmp.push_back(curr_ans);
				new_k=ans_tmp.size()-1;
				for(size_t i=0;i<curr_node->_childs.size();++i){
					nqueue.push(make_pair(curr_node->_childs[i],new_k));
				}
			}
			++height;
			nqueue.pop();
			nqueue.push(make_pair(nullptr,-1));
			ans_list.swap(ans_tmp);
		}
	}

	void compute_tree_top_down_directed(){
		//variables
		vector<vector<node*> > ans_list;
		int height=1;
		node* curr_node;
		int k,new_k;
		//bfs search
		queue<pair<node*,int> > nqueue;
		ans_list.push_back(vector<node*>(1,_root));
		for(size_t i=0;i<_root->_childs.size();++i) nqueue.push(make_pair(_root->_childs[i],0));
		nqueue.push(make_pair(nullptr,-1)); 
		while(!nqueue.empty()){
			if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
			vector<vector<node*> > ans_tmp;
			while(nqueue.front().first!=nullptr){
				pair<node*,int> front_pair=nqueue.front();
				curr_node=front_pair.first;
				k=front_pair.second;
				nqueue.pop();
				vector<node*> curr_ans=ans_list[k];
				curr_node->compute_diis_directed(curr_ans,height);
				curr_ans.push_back(curr_node);
				ans_tmp.push_back(curr_ans);
				new_k=ans_tmp.size()-1;
				for(size_t i=0;i<curr_node->_childs.size();++i){
					nqueue.push(make_pair(curr_node->_childs[i],new_k));
				}
			}
			++height;
			nqueue.pop();
			nqueue.push(make_pair(nullptr,-1));
			ans_list.swap(ans_tmp);
		}
	}

	void compute_tree_performance(){
		int max_node_width=0;
		int max_node_num=0;
		int node_num;
		int tree_height=0;
		int num_nodes=0;
		queue<node*> nqueue;
		nqueue.push(_root);
		nqueue.push(NULL);
		while(!nqueue.empty()){
			if(nqueue.size()==1&&nqueue.front()==NULL) break;
			node_num=0;
			while(nqueue.front()!=NULL){
				node_num++;
				node* curr_node=nqueue.front();
				nqueue.pop();
				if(curr_node->_width>max_node_width) max_node_width=curr_node->_width;
				//push all childs
				for(size_t i=0;i<curr_node->_childs.size();++i) nqueue.push(curr_node->_childs[i]);
			}
			tree_height++;
			num_nodes+=node_num;
			if(node_num>max_node_num) max_node_num=node_num;
			nqueue.pop();
			nqueue.push(NULL);
		}
		_tree_height=tree_height;
		_max_node_num=max_node_num;
		_max_node_width=max_node_width;
		if(num_nodes!=_numOfNodes){
			std::cout<<"num_nodes="<<std::endl;
			_numOfNodes=num_nodes;
		}
		std::cout<<"tree_height="<<tree_height<<endl;
		std::cout<<"max_node_num="<<max_node_num<<" max_node_width="<<max_node_width<<endl;
		std::cout<<"Compute tree performance successfully!!"<<std::endl;
	}

	void output_tree(char* outputDirName){
		string tree_output_filename(outputDirName);
		tree_output_filename.append("_tree.structure");
		ofstream tree_ofs(tree_output_filename);
		if(!tree_ofs.is_open()) cout<<"Cannot open "<<tree_output_filename<<endl;
		int max_node_width=0;
		int max_node_num=0;
		int node_num;
		int tree_height=0;
		queue<node*> nqueue;
		nqueue.push(_root);
		nqueue.push(NULL);
		while(!nqueue.empty()){
			if(nqueue.size()==1&&nqueue.front()==NULL) break;
			node_num=0;
			while(nqueue.front()!=NULL){
				node_num++;
				node* curr_node=nqueue.front();
				nqueue.pop();
				tree_ofs<<"("<<curr_node->_width+1<<":";
				if(curr_node->_width>max_node_width) max_node_width=curr_node->_width;
				for(size_t i=0;i<curr_node->_width;++i) tree_ofs<<curr_node->_neighbour[i]<<"-"<<curr_node->_neighbourDis[i]<<" ";
				tree_ofs<<curr_node->_curNum;
				tree_ofs<<") ";
				//push all childs
				for(size_t i=0;i<curr_node->_childs.size();++i) nqueue.push(curr_node->_childs[i]);
			}
			tree_height++;
			if(node_num>max_node_num) max_node_num=node_num;
			tree_ofs<<"["<<node_num<<"]";
			tree_ofs<<endl;
			nqueue.pop();
			nqueue.push(NULL);
		}
		tree_ofs.close();
		_tree_height=tree_height;
		_max_node_num=max_node_num;
		_max_node_width=max_node_width;
		std::cout<<"tree_height="<<tree_height<<endl;
		std::cout<<"max_node_num="<<max_node_num<<" max_node_width="<<max_node_width<<endl;
	}

	void build_dp_tree(vector<vector<pair<int,int> > >& neighborsInfo,vector<int>& rank){
		for(int v=0;v<_numOfNodes;++v) 
		{
			if(neighborsInfo[v].size()==0){
				cout<<"root="<<v<<endl;//to be deleted
				_root=&_nodes[v];
				_root->_parentNum=-1;
				_root->_rank=rank[v];
				_root->_height=0;
				_root->_width=0;
				_root->_pos=new int[1];
				_root->_dis=new int[1];
				_root->_pos[0]=0;
				_root->_dis[0]=0;
				continue;
			}
			_nodes[v].setParent(_nodes,neighborsInfo[v],rank);
		}
	}

	void build_dp_tree_directed(vector<vector<pair<int,int> > >& r_neighborsInfo,vector<vector<pair<int,int> > >& neighborsInfo,vector<int>& rank){
		for(int v=0;v<_numOfNodes;++v) 
		{
			if(neighborsInfo[v].size()==0&&r_neighborsInfo[v].size()==0){
				cout<<"root="<<v<<endl;
				_root=&_nodes[v];
				_root->_parentNum=-1;
				_root->_rank=rank[v];
				_root->_height=0;
				_root->_width=0;
				_root->_pos=new int[1];
				_root->_dis=new int[1];
				_root->_pos[0]=0;
				_root->_dis[0]=0;
				//Add for directed
				_root->_r_dis=new int[1];
				_root->_r_dis[0]=0;
				continue;
			}
			_nodes[v].setParent_directed(_nodes,r_neighborsInfo[v],neighborsInfo[v],rank);
		}
	}

//**********************utils functions*****************
	EdgeWeight query_pll(const integrated_index_t &idx_s, const integrated_index_t &idx_t) {
		EdgeWeight distance = INF_WEIGHT;
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

    EdgeWeight query_pll_directed(const integrated_bindex_t &idx_s, const integrated_bindex_t &idx_t) {
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

};

/**
 * @description: the forest class
 * @return {*}
 * @author: Wan Jingyi
 */
class forest:public tree{
	public:
		vector<node*> _trees;//store the root node  of each tree
		vector<int> _trees_num;//store the num of nodes of each tree
		vector<int> _trees_height;//store the height of each tree
		vector<int> _trees_width;//store the width of each tree
		vector<int> _trees_nodesNum;//store the maximal node num of each tree
		vector<vector<int> > _trees_borders;//store the borders id of each tree

		int _numOfTrees{0};//num of trees
		int _numOfForestNodes;//=numOfDeletedVertices
		int _max_allTree_num=-1;
		int _min_allTree_num=INF_WEIGHT;
		double _avg_allTree_num=0;
		int _max_allTree_height=-1;
		int _min_allTree_height=INF_WEIGHT;
		double _avg_allTree_height=0;
		int _max_allTree_width=-1;
		int _min_allTree_width=INF_WEIGHT;
		double _avg_allTree_width=0;
		int _max_allTree_nodesNum=-1;

		//for the second solution
		double _avg_borders_num=0;
		int _max_borders_num=-1;
		int _min_borders_num=INF_WEIGHT;
		//*******************constructions and deconstructions**********
		forest(){
			_numOfNodes=0;
			_numOfTrees=0;
			_numOfForestNodes=0;

		}

		forest(int numOfNodes,int numOfForestNodes){
			_numOfTrees=0;
			_numOfNodes=numOfNodes;
			_numOfForestNodes=numOfForestNodes;
			//reserve
			_trees.reserve(_numOfForestNodes);
			_trees_borders.reserve(_numOfForestNodes);
			_nodes.resize(_numOfNodes);
			//initialize allNodes
			for(int i=0;i<_numOfNodes;++i){
				node _node(i);
				_nodes[i]=_node;
			}
		}

		~forest(){}
		//***********************functions******************
		void initialize(int numOfNodes,int numOfForestNodes){
			_numOfTrees=0;
			_numOfNodes=numOfNodes;
			_numOfForestNodes=numOfForestNodes;
			//reserve
			_trees.reserve(_numOfForestNodes);
			_trees_borders.reserve(_numOfForestNodes);
			_nodes.resize(_numOfNodes);
			//initialize allNodes
			for(int i=0;i<_numOfNodes;++i){
				node _node(i);
				_nodes[i]=_node;
			}
		}

		/**
		 * @description: find and set the parent to create the structure for the second solution
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void build_dp_forest(vector<vector<pair<int,int> > >& neighborsInfo,vector<int>& rank,const vector<unsigned int>& queryTime,bool isSorted){
			//set parent for all deleted vertices
			for(int v=0;v<_numOfNodes;++v){
				_nodes[v]._rank=rank[v];
				if(rank[v]==-1) continue;
				if(neighborsInfo[v].size()==0){
					std::cout<<v<<"neighbor size=0!"<<std::endl;
					continue;
				}
				setParent_forest(v,neighborsInfo[v],rank,queryTime,isSorted);
			}
			//resize num vector
			_trees_num.resize(_numOfTrees,0);
		}

		//Add for directed
		/**
		 * @description: find and set the parent to create the structure for the second solution
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void build_dp_forest_directed(vector<vector<pair<int,int> > >& r_neighborsInfo,vector<vector<pair<int,int> > >& neighborsInfo,vector<int>& rank,const vector<unsigned int>& queryTime,bool isSorted){
			//set parent for all deleted vertices
			for(int v=0;v<_numOfNodes;++v){
				_nodes[v]._rank=rank[v];
				if(rank[v]==-1) continue;
				if(r_neighborsInfo[v].size()==0&&neighborsInfo[v].size()==0){
					std::cout<<v<<"neighbor size=0!"<<std::endl;
					continue;
				}
				setParent_forest_directed(v,r_neighborsInfo[v],neighborsInfo[v],rank,queryTime,isSorted);
			}
			//resize num vector
			_trees_num.resize(_numOfTrees,0);
		}

		/**
		 * @description: find and set the parent to create the structure for the third solution
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void build_dp_forest(vector<vector<pair<int,int> > >& neighborsInfo,vector<int>& rank){
			//set parent for all deleted vertices
			for(int v=0;v<_numOfNodes;++v){
				_nodes[v]._rank=rank[v];
				if(rank[v]==-1) continue;
				if(neighborsInfo[v].size()==0){
					std::cout<<v<<"neighbor size=0!"<<std::endl;
					continue;
				}
				setParent_forest(v,neighborsInfo[v],rank);
			}
			//resize num vector
			_trees_num.resize(_numOfTrees,0);
		}

		/**
		 * @description: find and set the parent to create the structure for the third solution,add for directed
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void build_dp_forest_directed(vector<vector<pair<int,int> > >& r_neighborsInfo,vector<vector<pair<int,int> > >& neighborsInfo,vector<int>& rank){
			//set parent for all deleted vertices
			for(int v=0;v<_numOfNodes;++v){
				_nodes[v]._rank=rank[v];
				if(rank[v]==-1) continue;
				if(r_neighborsInfo[v].size()==0&&neighborsInfo[v].size()==0){
					std::cout<<v<<"neighbor size=0!"<<std::endl;
					continue;
				}
				setParent_forest_directed(v,r_neighborsInfo[v],neighborsInfo[v],rank);
			}
			//resize num vector
			_trees_num.resize(_numOfTrees,0);
		}

		void compute_borders_poses(){
			for(int i=0;i<_numOfTrees;++i){
				for(int j=0;j<_trees_borders[i].size();++j){
					node* curr_node=&_nodes[_trees_borders[i][j]];
					if(curr_node->_nBorder_siz==0){
						curr_node->_nBorder_siz=_numOfTrees;
						curr_node->_nBorder_index=new int[_numOfTrees];
					}
					curr_node->_nBorder_index[i]=j;
				}
			}
		}

		//Single thread for the second solution (undirected)
		void compute_and_combinePLL_forest_top_down(Integrated_Index& integrated_index){
			for(int i=0;i<_numOfTrees;++i){
				compute_and_combinePLL_tree_top_down(i,integrated_index);
			}
		}

		//Single thread, Add for directed
		void compute_and_combinePLL_forest_top_down_directed(Integrated_Index& integrated_index){
			for(int i=0;i<_numOfTrees;++i){
				compute_and_combinePLL_tree_top_down_directed(i,integrated_index);
			}
		}

		//multiThreads version
		void compute_and_combinePLL_forest_top_down(const vector<integrated_index_t>& index_,int num_threads){
			omp_set_num_threads(num_threads);//设置线程数
			cout<<"start multithread compute..."<<std::endl;//to be deleted
			#pragma omp parallel for schedule(dynamic)
			//#pragma omp parallel for
				for(int i=0;i<_numOfTrees;++i){
					compute_and_combinePLL_tree_top_down(i,index_);
				}
			//#pragma omp parallel end
		}

		//multiThreads version, add for directed
		void compute_and_combinePLL_forest_top_down_directed(const vector<integrated_bindex_t>& bindex_,int num_threads){
			omp_set_num_threads(num_threads);//设置线程数
			cout<<"start multithread compute..."<<std::endl;//to be deleted
			#pragma omp parallel for schedule(dynamic)
			//#pragma omp parallel for
				for(int i=0;i<_numOfTrees;++i){
					compute_and_combinePLL_tree_top_down_directed(i,bindex_);
				}
			//#pragma omp parallel end
		}

		/**
		 * @description: single thread for the third solution
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void compute_and_pushPLL_forest_top_down(Integrated_Index& integrated_index){
			for(int i=0;i<_numOfTrees;++i){
				compute_and_pushPLL_tree_top_down(i,integrated_index);
			}
		}

		/**
		 * @description: single thread for the third solution,add for directed
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void compute_and_pushPLL_forest_top_down_directed(Integrated_Index& integrated_index){
			for(int i=0;i<_numOfTrees;++i){
				compute_and_pushPLL_tree_top_down_directed(i,integrated_index);
			}
		}

		//multiple threads version
		void compute_and_pushPLL_forest_top_down(const vector<integrated_index_t>& index_,int max_overlay_hub,int num_threads){
			omp_set_num_threads(num_threads);//设置线程数
			#pragma omp parallel for schedule(dynamic)
			//#pragma omp parallel for
				for(int i=0;i<_numOfTrees;++i){
					compute_and_pushPLL_tree_top_down(i,index_,max_overlay_hub);
				}
			cout<<"start multithread pushPll and compute finished!"<<std::endl;//to be deleted
			//#pragma omp parallel end
		}

		//multiple threads version for the third solution, add for directed
		void compute_and_pushPLL_forest_top_down_directed(const vector<integrated_bindex_t>& bindex_,int max_overlay_hub_r,int max_overlay_hub,int num_threads){
			omp_set_num_threads(num_threads);//设置线程数
			#pragma omp parallel for schedule(dynamic)
			//#pragma omp parallel for
				for(int i=0;i<_numOfTrees;++i){
					compute_and_pushPLL_tree_top_down_directed(i,bindex_,max_overlay_hub_r,max_overlay_hub);
				}
			cout<<"start multithread pushPll and compute finished!"<<std::endl;//to be deleted
			//#pragma omp parallel end
		}

		void compute_rmq_variables_forest(){
			int _n=_numOfForestNodes*2-1;
			_euler.reserve(_n);
			_level.reserve(_n);
			_index.resize(_numOfNodes,-1);
			for(int i=0;i<_numOfTrees;++i){
				dfs_preprocess_forest(_trees[i]);
			}
			// //to be deleted
			// cout<<"_euler:"<<_euler.size()<<endl;
			// for(size_t i=0;i<_euler.size();++i) cout<<_euler[i]<<" ";
			// cout<<endl;
			// cout<<"_level:"<<_level.size()<<endl;
			// for(size_t i=0;i<_level.size();++i) cout<<_level[i]<<" ";
			// cout<<endl;
			// cout<<"_index:"<<_index.size()<<endl;
			// for(size_t i=0;i<_index.size();++i) cout<<_index[i]<<" ";
			// cout<<endl;
		}

		/**
		 * @description: merge h2h_index for the third solution
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void merge_h2h_index_pushPLL_forest(vector<integrated_index_t>& index_){
			int i,j,k,v;
			for(v=0;v<_numOfNodes;++v){
				const node& curr_node=_nodes[v];
				if(curr_node._rank==-1) continue;
				integrated_index_t& idx=index_[v];
				//pos
				bool isHF=false;
				for(i=0;i<curr_node._width;++i){
					if(curr_node._neighbourPointer[i]->_rank==-1){
						isHF=true;
						continue;
					}else{
						break;
					}
				}
				k=curr_node._width-i+1;
				idx.pos.resize(k+2);
				idx.pos[0]=isHF;
				i=curr_node._dis_height-curr_node._height;//num of bordrs
				for(j=0;j<k;++j) idx.pos[j+1]=curr_node._pos[curr_node._width-j]-i;
				idx.pos[k+1]=_numOfNodes;
				//dis
				k=curr_node._height+1;
				j=curr_node._dis_height-curr_node._height;//num of bordrs
				idx.dis.resize(k+1);
				idx.dis[0]=_trees[curr_node._tree_id]->_curNum;
				for(i=0;i<k;++i){
					idx.dis[i+1]=curr_node._dis[j+i];
				}
				//level
				idx.level=numOfVertices-1-curr_node._rank;
				//size
				idx.siz=curr_node._siz;
			}
		}

		/**
		 * @description: merge h2h_index for the third solution,add for directed
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void merge_h2h_index_pushPLL_forest_directed(vector<integrated_bindex_t>& bindex_){
			int i,j,k,v;
			for(v=0;v<_numOfNodes;++v){
				const node& curr_node=_nodes[v];
				if(curr_node._rank==-1) continue;
				integrated_bindex_t& idx=bindex_[v];
				//pos
				bool isHF=false;
				for(i=0;i<curr_node._width;++i){
					if(curr_node._neighbourPointer[i]->_rank==-1){
						isHF=true;
						continue;
					}else{
						break;
					}
				}
				//pos
				k=curr_node._width-i+1;//pos size
				j=curr_node._dis_height-curr_node._height;//num of bordrs
				idx.pos.resize(k+3);
				idx.pos[0]=isHF;
				idx.pos[1]=_trees[curr_node._tree_id]->_curNum;
				for(i=0;i<k;++i) idx.pos[i+2]=curr_node._pos[curr_node._width-i]-j;
				idx.pos[k+2]=_numOfNodes;
				//dis and r_dis
				k=curr_node._height+1;
				idx.dis.resize(k);
				idx.r_dis.resize(k);
				for(i=0;i<k;++i){
					idx.dis[i]=curr_node._dis[j+i];
					idx.r_dis[i]=curr_node._r_dis[j+i];
				}
				//level
				idx.level=numOfVertices-1-curr_node._rank;
				//size
				idx.siz=curr_node._siz;
				idx.r_siz=curr_node._r_siz;
			}
		}

		//multiple threads version
		void merge_h2h_index_pushPLL_forest_multiThreads(vector<integrated_index_t>& index_){
			int i,j,k,v;
			for(v=0;v<_numOfNodes;++v){
				const node& curr_node=_nodes[v];
				if(curr_node._rank==-1) continue;
				integrated_index_t& idx=index_[v];
				k=curr_node._width-curr_node._nBorder_siz+1;
				idx.pos.resize(k+2);
				if(curr_node._nBorder_siz!=0) idx.pos[0]=1;
				else idx.pos[0]=0;
				i=curr_node._dis_height-curr_node._height;//num of bordrs
				for(j=0;j<k;++j) idx.pos[j+1]=curr_node._pos[curr_node._width-j]-i;
				idx.pos[k+1]=_numOfNodes;
				//dis
				k=curr_node._height+1;
				j=curr_node._dis_height-curr_node._height;//num of bordrs
				idx.dis.resize(k+1);
				idx.dis[0]=_trees[curr_node._tree_id]->_curNum;
				for(i=0;i<k;++i){
					idx.dis[i+1]=curr_node._dis[j+i];
				}
				//spt_v spt_d
				size_t spt_siz=curr_node.pll_v.size();
				idx.spt_v.resize(spt_siz);
				idx.spt_d.resize(spt_siz);
				for(size_t i=0;i<spt_siz;++i){
					idx.spt_v[i]=curr_node.pll_v[i];
					idx.spt_d[i]=curr_node.pll_d[i];
				}
				//level
				idx.level=numOfVertices-1-curr_node._rank;
				//size
				idx.siz=curr_node._siz;
			}
		}

		//multiple threads version for the third solution, add for directed
		void merge_h2h_index_pushPLL_forest_multiThreads_directed(vector<integrated_bindex_t>& bindex_){
			int i,j,k,v;
			for(v=0;v<_numOfNodes;++v){
				const node& curr_node=_nodes[v];
				if(curr_node._rank==-1) continue;
				integrated_bindex_t& idx=bindex_[v];
				k=curr_node._width-curr_node._nBorder_siz+1;
				idx.pos.resize(k+3);
				if(curr_node._nBorder_siz!=0) idx.pos[0]=1;
				else idx.pos[0]=0;
				idx.pos[1]=_trees[curr_node._tree_id]->_curNum;
				i=curr_node._dis_height-curr_node._height;//num of bordrs
				for(j=0;j<k;++j) idx.pos[j+2]=curr_node._pos[curr_node._width-j]-i;
				idx.pos[k+2]=_numOfNodes;
				//dis and r_dis
				k=curr_node._height+1;
				j=curr_node._dis_height-curr_node._height;//num of bordrs
				idx.dis.resize(k);
				idx.r_dis.resize(k);
				for(i=0;i<k;++i){
					idx.dis[i]=curr_node._dis[j+i];
					idx.r_dis[i]=curr_node._r_dis[j+i];
				}
				//spt_v spt_d
				size_t spt_siz=curr_node.pll_v.size();
				idx.spt_v.resize(spt_siz);
				idx.spt_d.resize(spt_siz);
				for(size_t i=0;i<spt_siz;++i){
					idx.spt_v[i]=curr_node.pll_v[i];
					idx.spt_d[i]=curr_node.pll_d[i];
				}
				//r_spt_v r_spt_d
				size_t r_spt_siz=curr_node.pll_v_r.size();
				idx.r_spt_v.resize(r_spt_siz);
				idx.r_spt_d.resize(r_spt_siz);
				for(size_t i=0;i<r_spt_siz;++i){
					idx.r_spt_v[i]=curr_node.pll_v_r[i];
					idx.r_spt_d[i]=curr_node.pll_d_r[i];
				}
				//level
				idx.level=numOfVertices-1-curr_node._rank;
				//size r_size
				idx.siz=curr_node._siz;
				idx.r_siz=curr_node._r_siz;
			}
		}

		/**
		 * @description: single thread merge to integrated_index
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void merge_h2h_index_forest(vector<integrated_index_t>& index_){
			int i,j,k,v;
			//root node
			for(i=0;i<_numOfTrees;++i){
				node* root=_trees[i];
				v=root->_curNum;
				integrated_index_t& idx=index_[v];
				//spt-v
				k=_trees_borders[i].size();
				idx.spt_v.resize(k+1);
				for(j=0;j<k;++j) idx.spt_v[j]=_trees_borders[i][j];
				idx.spt_v[k]=_numOfNodes;
				//pos
				k=root->_width+1;
				idx.pos.resize(k+1);
				for(j=0;j<k;++j) idx.pos[j]=root->_pos[k-1-j];
				idx.pos[k]=_numOfNodes;
				//dis
				k=root->_dis_height+1;
				idx.dis.resize(k+1);
				idx.dis[0]=v;
				for(j=0;j<k;++j) idx.dis[j+1]=root->_dis[j];
				//level
				idx.level=numOfVertices-1-root->_rank;
				//size
				idx.siz=root->_siz;
			}
			//non-root node
			for(v=0;v<_numOfNodes;++v){
				const node& curr_node=_nodes[v];
				if(curr_node._rank==-1||curr_node._parent==NULL) continue;
				integrated_index_t& idx=index_[v];
				//pos
				k=curr_node._width+1;
				idx.pos.resize(k+1);
				for(i=0;i<k;++i) idx.pos[i]=curr_node._pos[k-1-i];
				idx.pos[k]=_numOfNodes;
				//dis
				k=curr_node._dis_height+1;
				idx.dis.resize(k+1);
				idx.dis[0]=_trees[curr_node._tree_id]->_curNum;
				for(i=0;i<k;++i) idx.dis[i+1]=curr_node._dis[i];
				//level
				idx.level=numOfVertices-1-curr_node._rank;
				//size
				idx.siz=curr_node._siz;
			}
		}

		/**
		 * @description: single thread merge to integrated_index, add for directed
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void merge_h2h_index_forest_directed(vector<integrated_bindex_t>& bindex_){
			int i,j,k,v;
			//root node
			for(i=0;i<_numOfTrees;++i){
				node* root=_trees[i];
				v=root->_curNum;
				integrated_bindex_t& idx=bindex_[v];
				//spt-v
				k=_trees_borders[i].size();
				idx.spt_v.resize(k+1);
				for(j=0;j<k;++j) idx.spt_v[j]=_trees_borders[i][j];
				idx.spt_v[k]=_numOfNodes;
				//pos
				k=root->_width+1;
				idx.pos.resize(k+2);
				idx.pos[0]=v;
				for(j=0;j<k;++j) idx.pos[j+1]=root->_pos[k-1-j];
				idx.pos[k+1]=_numOfNodes;
				//dis and r_dis
				k=root->_dis_height+1;
				idx.dis.resize(k);
				idx.r_dis.resize(k);
				for(j=0;j<k;++j){
					idx.dis[j]=root->_dis[j];
					idx.r_dis[j]=root->_r_dis[j];
				}
				//level
				idx.level=numOfVertices-1-root->_rank;
				//size
				idx.siz=root->_siz;
				//r_size
				idx.r_siz=root->_r_siz;
			}
			//non-root node
			for(v=0;v<_numOfNodes;++v){
				const node& curr_node=_nodes[v];
				if(curr_node._rank==-1||curr_node._parent==NULL) continue;
				integrated_bindex_t& idx=bindex_[v];
				//pos
				k=curr_node._width+1;
				idx.pos.resize(k+2);
				idx.pos[0]=_trees[curr_node._tree_id]->_curNum;
				for(i=0;i<k;++i) idx.pos[i+1]=curr_node._pos[k-1-i];
				idx.pos[k+1]=_numOfNodes;
				//dis
				k=curr_node._dis_height+1;
				idx.dis.resize(k);
				idx.r_dis.resize(k);
				for(i=0;i<k;++i){
					idx.dis[i]=curr_node._dis[i];
					idx.r_dis[i]=curr_node._r_dis[i];
				}
				//level
				idx.level=numOfVertices-1-curr_node._rank;
				//size
				idx.siz=curr_node._siz;
				//r_size
				idx.siz=curr_node._r_siz;
			}
		}

		/**
		 * @description: save lca variables to binary file
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
		void save_lca_variables_forest(char* outputDirName){
			string lca_filename(outputDirName);
			lca_filename.append(".lca");
			ofstream out(lca_filename,ios::binary | ios::out);
			if(!out.is_open()) cout<<"Cannot open "<<lca_filename<<endl;
			int esize=_euler.size();
			int lsize=_level.size();
			int isize=_index.size();
			int _n=_numOfForestNodes*2-1;
			std::cout<<"_numOfForestNodes*2-1="<<_n<<" total_euler_size="<<esize<<" total_level_size="<<lsize<<" total_index_size="<<lsize<<std::endl;
			//_euler, _level,_index
			out.write((const char*)&esize,sizeof(esize));
			for(int i=0;i<esize;++i) out.write((const char*)&_euler[i],sizeof(_euler[i]));
			out.write((const char*)&lsize,sizeof(lsize));
			for(int i=0;i<lsize;++i) out.write((const char*)&_level[i],sizeof(_level[i]));
			out.write((const char*)&isize,sizeof(isize)); 
			for(int i=0;i<isize;++i) out.write((const char*)&_index[i],sizeof(_index[i]));
			out.close();
			//std::cout << "Save forest lca variables successfully!" << std::endl;
		}

		//for the second integrated solution
		void compute_node_label_size_forest_2(vector<integrated_index_t>& index_,const vector<int>& rank){
			//initialize
			_start.resize(_numOfNodes);
			_end.resize(_numOfNodes);
			for(int i=0;i<_numOfTrees;++i){
				node* root=_trees[i];
				_dfs.reserve(_trees_num[i]);
				dfs_compute(root);
				vector<int> ().swap(_dfs);
				_dfs.clear();
			}
			int cnt=0;//for debug
			for(int v=0;v<_numOfNodes;++v){
				if(rank[v]==-1) continue;
				++cnt;
				int ans_num,tree_num,left_num,border_num,ans_total_size;
				double label_size_1=1,label_size_2=0,label_size_3=0,label_size_4=0,label_size=0;
				node* curr_node=&_nodes[v];
				int tree_index=curr_node->_tree_id;
				ans_num=curr_node->_ans.size();
				tree_num=(_end[v]-_start[v]-1)/2+1;
				border_num=_trees_borders[tree_index].size();
				left_num=_trees_num[tree_index]-ans_num-tree_num;
				ans_total_size=0;
				for(int i=0;i<ans_num;++i){
					ans_total_size+=curr_node->_ans[i]->_width+1;
				}
				if(ans_num!=0) label_size_2=(double)ans_total_size/(double)ans_num;
				label_size_3=(double)((1+border_num)*border_num)/2.0;
				for(int i=0;i<border_num;++i){
					label_size_4+=index_[_trees_borders[tree_index][i]].spt_v.size();
				}
				label_size=label_size_1*(ans_num+tree_num)+label_size_2*left_num+label_size_3*border_num+label_size_4*(_numOfNodes-_trees_num[tree_index]-border_num);
				label_size=label_size/(double)_numOfNodes;
				curr_node->_siz=label_size;
			}
			if(cnt!=_numOfForestNodes) cout<<"Error:cnt!=_numOfForestNodes!"<<" cnt="<<cnt<<" _numOfForestNodes="<<_numOfForestNodes<<endl;
			std::cout<<"Compute node label size forest successfully!"<<endl;
		}

		//for the second integrated solution, add for directed
		void compute_node_label_size_forest_2_directed(vector<integrated_bindex_t>& bindex_,const vector<int>& rank){
			//initialize
			_start.resize(_numOfNodes);
			_end.resize(_numOfNodes);
			for(int i=0;i<_numOfTrees;++i){
				node* root=_trees[i];
				_dfs.reserve(_trees_num[i]);
				dfs_compute(root);
				vector<int> ().swap(_dfs);
				_dfs.clear();
			}
			int cnt=0;//for debug
			for(int v=0;v<_numOfNodes;++v){
				if(rank[v]==-1) continue;
				++cnt;
				int ans_num,tree_num,left_num,border_num,ans_total_size;
				double label_size_1=1,label_size_2=0,label_size_3=0,label_size_4=0,label_size=0;
				node* curr_node=&_nodes[v];
				int tree_index=curr_node->_tree_id;
				ans_num=curr_node->_ans.size();
				tree_num=(_end[v]-_start[v]-1)/2+1;
				border_num=_trees_borders[tree_index].size();
				left_num=_trees_num[tree_index]-ans_num-tree_num;
				ans_total_size=0;
				for(int i=0;i<ans_num;++i){
					ans_total_size+=curr_node->_ans[i]->_width+1;
				}
				if(ans_num!=0) label_size_2=(double)ans_total_size/(double)ans_num;
				label_size_3=(double)border_num;
				//***********s and t is different only for label_size_4************
				//as start
				for(int i=0;i<border_num;++i){
					if(curr_node->_dis[i]==INF_WEIGHT) continue;
					label_size_4+=bindex_[_trees_borders[tree_index][i]].spt_v.size()-1;
				}
				label_size=label_size_1*(ans_num+tree_num)+label_size_2*left_num+label_size_3*border_num+label_size_4*(_numOfNodes-_trees_num[tree_index]-border_num);
				label_size=label_size/(double)_numOfNodes;
				curr_node->_siz=label_size;
				//as termination
				label_size_4=0;
				label_size=0;
				for(int i=0;i<border_num;++i){
					if(curr_node->_r_dis[i]==INF_WEIGHT) continue;
					label_size_4+=bindex_[_trees_borders[tree_index][i]].r_spt_v.size()-1;
				}
				label_size=label_size_1*(ans_num+tree_num)+label_size_2*left_num+label_size_3*border_num+label_size_4*(_numOfNodes-_trees_num[tree_index]-border_num);
				label_size=label_size/(double)_numOfNodes;
				curr_node->_r_siz=label_size;
			}
			if(cnt!=_numOfForestNodes) cout<<"Error:cnt!=_numOfForestNodes!"<<" cnt="<<cnt<<" _numOfForestNodes="<<_numOfForestNodes<<endl;
			
		}

		//for the third integrated solution
		void compute_node_label_size_forest_3(vector<integrated_index_t>& index_,const vector<int>& rank,bool isMultiple=false){
			//initialize
			_start.resize(_numOfNodes);
			_end.resize(_numOfNodes);
			for(int i=0;i<_numOfTrees;++i){
				node* root=_trees[i];
				_dfs.reserve(_trees_num[i]);
				dfs_compute(root);
				vector<int> ().swap(_dfs);
				_dfs.clear();
			}
			int cnt=0;//for debug
			for(int v=0;v<_numOfNodes;++v){
				if(rank[v]==-1) continue;
				++cnt;
				int ans_num,tree_num,left_num,border_num,ans_total_size;
				double label_size_1=1,label_size_2=0,label_size_3=0,label_size=0;
				node* curr_node=&_nodes[v];
				int tree_index=curr_node->_tree_id;
				ans_num=curr_node->_ans.size();
				tree_num=(_end[v]-_start[v]-1)/2+1;
				border_num=_trees_borders[tree_index].size();
				left_num=_trees_num[tree_index]-ans_num-tree_num;
				if(isMultiple==true) label_size_3=curr_node->pll_v.size()-1;
				else label_size_3=index_[v].spt_v.size()-1;
				ans_total_size=0;
				for(int i=0;i<ans_num;++i){
					ans_total_size+=curr_node->_ans[i]->_width+1;
					if(curr_node->_ans[i]->_nBorder_siz!=0) ans_total_size+=label_size_3;
				}
				if(ans_num!=0) label_size_2=(double)ans_total_size/(double)ans_num;
				label_size=label_size_1*(ans_num+tree_num)+label_size_2*left_num+label_size_3*(_numOfNodes-_trees_num[tree_index]);
				label_size=label_size/(double)_numOfNodes;
				curr_node->_siz=label_size;
			}
			if(cnt!=_numOfForestNodes) cout<<"Error:cnt!=_numOfForestNode!"<<endl;
			std::cout<<"Compute node label size forest successfully!"<<endl;
		}

		//for the third integrated solution,add for directed
		void compute_node_label_size_forest_3_directed(vector<integrated_bindex_t>& bindex_,const vector<int>& rank,bool isMultiple=false){
			//initialize
			_start.resize(_numOfNodes);
			_end.resize(_numOfNodes);
			for(int i=0;i<_numOfTrees;++i){
				node* root=_trees[i];
				_dfs.reserve(_trees_num[i]);
				dfs_compute(root);
				vector<int> ().swap(_dfs);
				_dfs.clear();
			}
			int cnt=0;//for debug
			for(int v=0;v<_numOfNodes;++v){
				if(rank[v]==-1) continue;
				++cnt;
				int ans_num,tree_num,left_num,border_num,ans_total_size;
				double label_size_1=1,label_size_2,label_size_3,label_size=0,r_label_size=0;
				node* curr_node=&_nodes[v];
				int tree_index=curr_node->_tree_id;
				ans_num=curr_node->_ans.size();
				tree_num=(_end[v]-_start[v]-1)/2+1;
				border_num=_trees_borders[tree_index].size();
				left_num=_trees_num[tree_index]-ans_num-tree_num;
				//forward
				label_size_2=0;
				if(isMultiple) label_size_3=curr_node->pll_v.size()-1;
				else label_size_3=bindex_[v].spt_v.size()-1;
				ans_total_size=0;
				for(int i=0;i<ans_num;++i){
					ans_total_size+=curr_node->_ans[i]->_width+1;
					if(curr_node->_ans[i]->_nBorder_siz!=0) ans_total_size+=label_size_3;
				}
				if(ans_num!=0) label_size_2=(double)ans_total_size/(double)ans_num;
				label_size=label_size_1*(ans_num+tree_num)+label_size_2*left_num+label_size_3*(_numOfNodes-_trees_num[tree_index]);
				label_size=label_size/(double)_numOfNodes;
				curr_node->_siz=label_size;
				//backward
				label_size_2=0;
				if(isMultiple) label_size_3=curr_node->pll_v_r.size()-1;
				else label_size_3=bindex_[v].r_spt_v.size()-1;
				ans_total_size=0;
				for(int i=0;i<ans_num;++i){
					ans_total_size+=curr_node->_ans[i]->_width+1;
					if(curr_node->_ans[i]->_nBorder_siz!=0) ans_total_size+=label_size_3;
				}
				if(ans_num!=0) label_size_2=(double)ans_total_size/(double)ans_num;
				r_label_size=label_size_1*(ans_num+tree_num)+label_size_2*left_num+label_size_3*(_numOfNodes-_trees_num[tree_index]);
				r_label_size=r_label_size/(double)_numOfNodes;
				curr_node->_r_siz=r_label_size;
			}
			if(cnt!=_numOfForestNodes) cout<<"Error:cnt!=_numOfForestNode!"<<endl;
			std::cout<<"Compute node label size forest successfully!"<<endl;
		}

		//*****************utils functions****************
		/**
		 * @description: output each tree structure to txt file
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void output_forest(char* write_filename){
			string forest_output_filename(write_filename);
			forest_output_filename.append("_forest.structure");
			ofstream forest_ofs(forest_output_filename);
			if(!forest_ofs.is_open()) std::cout<<"Cannot open "<<forest_output_filename<<std::endl;
			for(int i=0;i<_numOfTrees;++i){
				output_tree_forest(i,forest_ofs);
			}
			forest_ofs.close();
		}

		/**
		 * @description: 
		 * @param {*}
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void count_and_analyze_forest(){
			//resize
			_trees_height.resize(_numOfTrees,0);
			_trees_width.resize(_numOfTrees,0);
			_trees_nodesNum.resize(_numOfTrees,0);
			long long sum_allTree_num=0,sum_allTree_height=0,sum_allTree_width=0,nodes_cnt=0,sum_borders_num=0;
			for(int i=0;i<_numOfTrees;++i) 
			{
				count_and_analyze_tree(i);
				sum_borders_num+=_trees_borders[i].size()*_trees_num[i];
				_min_borders_num=min(_min_borders_num,int(_trees_borders[i].size()));
				_max_borders_num=max(_max_borders_num,int(_trees_borders[i].size()));
				nodes_cnt+=_trees_num[i];
				_min_allTree_num=min(_min_allTree_num,_trees_num[i]);
				_max_allTree_num=max(_max_allTree_num,_trees_num[i]);
				sum_allTree_num+=_trees_num[i];
				_min_allTree_height=min(_min_allTree_height,_trees_height[i]);
				_max_allTree_height=max(_max_allTree_height,_trees_height[i]);
				sum_allTree_height+=_trees_height[i];
				_min_allTree_width=min(_min_allTree_num,_trees_width[i]);
				_max_allTree_width=max(_max_allTree_num,_trees_width[i]);
				sum_allTree_width+=_trees_width[i];
				_max_allTree_nodesNum=max(_max_allTree_nodesNum,_trees_nodesNum[i]);
			}
			if(nodes_cnt!=_numOfForestNodes) std::cout<<"nodes_cnt="<<nodes_cnt<<" _numOfForestNodes="<<_numOfForestNodes<<std::endl;
			//count average
			_avg_borders_num=(double)sum_borders_num/(double)_numOfForestNodes;
			_avg_allTree_num=(double)sum_allTree_num/(double)_numOfTrees;
			_avg_allTree_height=(double)sum_allTree_height/(double)_numOfTrees;
			_avg_allTree_width=(double)sum_allTree_width/(double)_numOfTrees;
			// _min_allTree_num=*min_element(_trees_num.begin(),_trees_num.end());
			// _max_allTree_num=*max_element(_trees_num.begin(),_trees_num.end());
			// _max_allTree_height=*max_element(_trees_height.begin(),_trees_height.end());
			// _max_allTree_width=*max_element(_trees_width.begin(),_trees_width.end());
			// _max_allTree_nodesNum=*max_element(_trees_nodesNum.begin(),_trees_nodesNum.end());
			std::cout<<"_min_borders_num = "<<_min_borders_num<<" _max_borders_num = "<<_max_borders_num<<" _avg_borders_num = "<<_avg_borders_num<<std::endl;
			std::cout<<"_min_allTree_num = "<<_min_allTree_num<<" _max_allTree_num = "<<_max_allTree_num<<" _avg_allTree_num = "<<_avg_allTree_num<<std::endl;
			std::cout<<"_min_allTree_height = "<<_min_allTree_height<<" _max_allTree_height = "<<_max_allTree_height<<" _avg_allTree_height = "<<_avg_allTree_height<<std::endl;
			std::cout<<"_min_allTree_width = "<<_min_allTree_width<<" _max_allTree_width = "<<_max_allTree_width<<" _avg_allTree_width = "<<_avg_allTree_width<<std::endl;
			std::cout<<"_max_allTree_node_num = "<<_max_allTree_nodesNum<<std::endl;
			std::cout<<"count_and_analyze_forest successfully!!"<<std::endl;
		}

	protected:
		void output_tree_forest(int tree_index,ofstream& tree_ofs){
			int node_num=0;
			node* root=_trees[tree_index];
			queue<node*> nqueue;
			nqueue.push(root);
			nqueue.push(NULL);
			tree_ofs<<_trees_num[tree_index]<<":"<<endl;
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front()==NULL) break;
				node_num=0;
				while(nqueue.front()!=NULL){
					node_num++;
					node* curr_node=nqueue.front();
					nqueue.pop();
					tree_ofs<<"("<<curr_node->_width+1<<":";
					if(curr_node!=root) for(size_t i=0;i<curr_node->_width;++i) tree_ofs<<curr_node->_neighbour[i]<<"-"<<curr_node->_neighbourDis[i]<<" ";
					tree_ofs<<curr_node->_curNum;
					tree_ofs<<") ";
					//push all childs
					for(size_t i=0;i<curr_node->_childs.size();++i) nqueue.push(curr_node->_childs[i]);
				}
				tree_ofs<<"["<<node_num<<"]";
				tree_ofs<<endl;
				nqueue.pop();
				nqueue.push(NULL);
			}
			tree_ofs<<"*****************************************************************"<<endl;
		}

		void count_and_analyze_tree(int tree_index){
			node* root=_trees[tree_index];
			int max_node_width=0;
			int max_node_num=0;
			int node_num;
			int tree_height=0;
			int num_nodes=0;
			queue<node*> nqueue;
			nqueue.push(root);
			nqueue.push(NULL);
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front()==NULL) break;
				node_num=0;
				while(nqueue.front()!=NULL){
					node_num++;
					node* curr_node=nqueue.front();
					nqueue.pop();
					if(curr_node!=root&&curr_node->_width>max_node_width) max_node_width=curr_node->_width;
					//push all childs
					for(size_t i=0;i<curr_node->_childs.size();++i) nqueue.push(curr_node->_childs[i]);
				}
				tree_height++;
				num_nodes+=node_num;
				if(node_num>max_node_num) max_node_num=node_num;
				nqueue.pop();
				nqueue.push(NULL);
			}
			_trees_height[tree_index]=tree_height;
			_trees_width[tree_index]=max_node_width;
			_trees_nodesNum[tree_index]=max_node_num;
			_trees_num[tree_index]=num_nodes;
		}

		/**
		 * @description: set parent for each tree in the forest for the second solution
		 * @param {const} vector
		 * @param {const} vector
		 * @return {*}
		 * @author: Wan Jingyi
		 */ 
 		void setParent_forest(int v,const vector<pair<int,int> >& neighbourInfo_v,const vector<int>& rank,const vector<unsigned int>& queryTime,bool isSorted){
			node* curr_node=&_nodes[v];
			int width=neighbourInfo_v.size();
			//sort and store
			vector<pair<int,int> > rank_in_tmp;
			vector<int > rank_border_tmp;

			for(int i=0;i<width;++i){
				int rank_v=rank[neighbourInfo_v[i].first];
				if(rank_v!=-1) rank_in_tmp.push_back(make_pair(rank_v,i));
				else rank_border_tmp.push_back(i);
			}
			//set root
			if(rank_in_tmp.size()==0){
				curr_node->_parentNum=-1;
				curr_node->_width=width;
				curr_node->_neighbour=new int[width];
				curr_node->_neighbourDis=new int[width];
				curr_node->_neighbourPointer=new node*[width];
				curr_node->_nBorder_siz=width;//added
				//reorder borders
				if(isSorted){
					vector<pair<unsigned int,int> > borders_rank(width);
					for(int i=0;i<width;++i) borders_rank[i]=make_pair(queryTime[neighbourInfo_v[i].first],i);
					sort(borders_rank.rbegin(),borders_rank.rend());
					for(int i=0;i<width;++i) rank_border_tmp[i]=borders_rank[i].second;
				}
				vector<int> _borders(width);
				for(int i=0;i<width;++i){
					int k=rank_border_tmp[i];
					_borders[i]=neighbourInfo_v[k].first;
					curr_node->_neighbour[i]=_borders[i];
					curr_node->_neighbourDis[i]=neighbourInfo_v[k].second;
					curr_node->_neighbourPointer[i]=&_nodes[_borders[i]];
				}
				_trees.push_back(curr_node);
				_trees_borders.push_back(_borders);
				++_numOfTrees;
			}else{
				curr_node->_width=width;
				curr_node->_neighbour=new int[width];
				curr_node->_neighbourDis=new int[width];
				curr_node->_neighbourPointer=new node*[width];
				curr_node->_nBorder_siz=rank_border_tmp.size();//added
				sort(rank_in_tmp.rbegin(),rank_in_tmp.rend());//sort
				int i;
				for(i=0;i<rank_border_tmp.size();++i){
					int k=rank_border_tmp[i];
					curr_node->_neighbour[i]=neighbourInfo_v[k].first;
					curr_node->_neighbourDis[i]=neighbourInfo_v[k].second;
					curr_node->_neighbourPointer[i]=&_nodes[curr_node->_neighbour[i]];
				}
				for(;i<width;++i){
					int k=rank_in_tmp[i-rank_border_tmp.size()].second;
					curr_node->_neighbour[i]=neighbourInfo_v[k].first;
					curr_node->_neighbourDis[i]=neighbourInfo_v[k].second;
					curr_node->_neighbourPointer[i]=&_nodes[curr_node->_neighbour[i]];
				}
				curr_node->_parentNum=curr_node->_neighbour[width-1];//the smallest order value
				curr_node->_parent=&_nodes[curr_node->_parentNum];
				//set the current node as childs to its parent node
				curr_node->_parent->_childs.push_back(curr_node);
			}
		}

		//Add for directed
		/**
		 * @description: set parent for each tree in the forest for the second solution
		 * @param {const} vector
		 * @param {const} vector
		 * @return {*}
		 * @author: Wan Jingyi
		 */ 
 		void setParent_forest_directed(int v,const vector<pair<int,int> >& r_neighbourInfo_v,const vector<pair<int,int> >& neighbourInfo_v,const vector<int>& rank,const vector<unsigned int>& queryTime,bool isSorted){
			node* curr_node=&_nodes[v];
			//remove the duplicated neighbors
			int nbr,nbr_w;
			unordered_map<int,pair<int,int> > nbr_map;
			for(size_t i=0;i<r_neighbourInfo_v.size();++i){
				nbr=r_neighbourInfo_v[i].first;
				nbr_w=r_neighbourInfo_v[i].second;
				nbr_map[nbr]=make_pair(nbr_w,INF_WEIGHT);
			}
			for(size_t i=0;i<neighbourInfo_v.size();++i){
				nbr=neighbourInfo_v[i].first;
				nbr_w=neighbourInfo_v[i].second;
				if(nbr_map.find(nbr)==nbr_map.end()) nbr_map[nbr]=make_pair(INF_WEIGHT,nbr_w);
				else nbr_map[nbr].second=nbr_w;
			}
			int width=nbr_map.size();
			//sort by rank
			vector<pair<int,int> > rank_in_tmp;
			vector<int > rank_border_tmp;
			for(unordered_map<int,pair<int,int> >::iterator it=nbr_map.begin();it!=nbr_map.end();++it){
				int rank_v=rank[it->first];
				if(rank_v!=-1) rank_in_tmp.push_back(make_pair(rank_v,it->first));
				else rank_border_tmp.push_back(it->first);
			}
			//set root
			if(rank_in_tmp.size()==0){
				curr_node->_parentNum=-1;
				curr_node->_width=width;
				curr_node->_height=0;
				curr_node->_neighbour=new int[width];
				curr_node->_neighbourDis=new int[width];
				curr_node->_r_neighbourDis=new int[width];
				curr_node->_neighbourPointer=new node*[width];
				curr_node->_nBorder_siz=width;//added
				//curr_node->_dis_height=width;
				//curr_node->_pos=new int[width+1];
				//curr_node->_dis=new int[width+1];
				//curr_node->_r_dis=new int[width+1];
				//reorder borders
				if(isSorted){
					vector<pair<unsigned int,int> > borders_rank(width);
					for(int i=0;i<width;++i) borders_rank[i]=make_pair(queryTime[rank_border_tmp[i]],rank_border_tmp[i]); //pair(rank_v,v)
					sort(borders_rank.rbegin(),borders_rank.rend());
					for(int i=0;i<width;++i) rank_border_tmp[i]=borders_rank[i].second;
				}
				vector<int> _borders(width);
				for(int i=0;i<width;++i){
					nbr=rank_border_tmp[i];
					_borders[i]=nbr;
					curr_node->_neighbour[i]=nbr;
					curr_node->_r_neighbourDis[i]=nbr_map[rank_border_tmp[i]].first;
					curr_node->_neighbourDis[i]=nbr_map[rank_border_tmp[i]].second;
					curr_node->_neighbourPointer[i]=&_nodes[nbr];
					//curr_node->_pos[i]=i;
					//curr_node->_r_dis[i]=nbr_map[rank_border_tmp[i]].first;
					//curr_node->_dis[i]=nbr_map[rank_border_tmp[i]].second;
				}
				//curr_node->_dis[width]=0;
				//curr_node->_r_dis[width]=0;
				//curr_node->_pos[width]=width;
				_trees.push_back(curr_node);
				_trees_borders.push_back(_borders);
				++_numOfTrees;
			}else{
				curr_node->_width=width;
				curr_node->_neighbour=new int[width];
				curr_node->_neighbourDis=new int[width];
				curr_node->_r_neighbourDis=new int[width];
				curr_node->_neighbourPointer=new node*[width];
				curr_node->_nBorder_siz=rank_border_tmp.size();//added
				sort(rank_in_tmp.rbegin(),rank_in_tmp.rend());//sort
				int i;
				for(i=0;i<rank_border_tmp.size();++i){
					nbr=rank_border_tmp[i];
					curr_node->_neighbour[i]=nbr;
					curr_node->_r_neighbourDis[i]=nbr_map[nbr].first;
					curr_node->_neighbourDis[i]=nbr_map[nbr].second;
					curr_node->_neighbourPointer[i]=&_nodes[nbr];
				}
				for(;i<width;++i){
					nbr=rank_in_tmp[i-rank_border_tmp.size()].second;
					curr_node->_neighbour[i]=nbr;
					curr_node->_r_neighbourDis[i]=nbr_map[nbr].first;
					curr_node->_neighbourDis[i]=nbr_map[nbr].second;
					curr_node->_neighbourPointer[i]=&_nodes[nbr];
				}
				curr_node->_parentNum=curr_node->_neighbour[width-1];//the smallest order value
				curr_node->_parent=&_nodes[curr_node->_parentNum];
				//set the current node as childs to its parent node
				curr_node->_parent->_childs.push_back(curr_node);
			}
		}

		/**
		 * @description: set parent for each tree in the forest for the third solution
		 * @param {const} vector
		 * @param {const} vector
		 * @return {*}
		 * @author: Wan Jingyi
		 */ 
 		void setParent_forest(int v,const vector<pair<int,int> >& neighbourInfo_v,const vector<int>& rank){
			node* curr_node=&_nodes[v];
			int width=neighbourInfo_v.size();
			//sort and store
			vector<pair<int,int> > rank_in_tmp;
			vector<int > rank_border_tmp;

			for(int i=0;i<width;++i){
				int rank_v=rank[neighbourInfo_v[i].first];
				if(rank_v!=-1) rank_in_tmp.push_back(make_pair(rank_v,i));
				else rank_border_tmp.push_back(i);
			}
			curr_node->_width=width;
			curr_node->_neighbour=new int[width];
			curr_node->_neighbourDis=new int[width];
			curr_node->_neighbourPointer=new node*[width];
			curr_node->_nBorder_siz=rank_border_tmp.size();//added
			//set root
			if(rank_in_tmp.size()==0){
				curr_node->_parentNum=-1;
				curr_node->_width=width;
				curr_node->_height=0;
				vector<int> _borders(width);
				for(int i=0;i<width;++i){
					int k=rank_border_tmp[i];
					_borders[i]=neighbourInfo_v[k].first;
					curr_node->_neighbour[i]=neighbourInfo_v[k].first;
					curr_node->_neighbourDis[i]=neighbourInfo_v[k].second;
					curr_node->_neighbourPointer[i]=&_nodes[curr_node->_neighbour[i]];
				}
				_trees.push_back(curr_node);
				_trees_borders.push_back(_borders);
				++_numOfTrees;
			}else{
				sort(rank_in_tmp.rbegin(),rank_in_tmp.rend());//sort
				int i;
				for(i=0;i<rank_border_tmp.size();++i){
					int k=rank_border_tmp[i];
					curr_node->_neighbour[i]=neighbourInfo_v[k].first;
					curr_node->_neighbourDis[i]=neighbourInfo_v[k].second;
					curr_node->_neighbourPointer[i]=&_nodes[curr_node->_neighbour[i]];
				}
				for(;i<width;++i){
					int k=rank_in_tmp[i-rank_border_tmp.size()].second;
					curr_node->_neighbour[i]=neighbourInfo_v[k].first;
					curr_node->_neighbourDis[i]=neighbourInfo_v[k].second;
					curr_node->_neighbourPointer[i]=&_nodes[curr_node->_neighbour[i]];
				}
				curr_node->_parentNum=curr_node->_neighbour[width-1];//the smallest order value
				curr_node->_parent=&_nodes[curr_node->_parentNum];
				//set the current node as childs to its parent node
				curr_node->_parent->_childs.push_back(curr_node);
			}
		}

		/**
		 * @description: set parent for each tree in the forest for the third solution, add for directed
		 * @param {const} vector
		 * @param {const} vector
		 * @return {*}
		 * @author: Wan Jingyi
		 */ 
 		void setParent_forest_directed(int v,const vector<pair<int,int> >& r_neighbourInfo_v,const vector<pair<int,int> >& neighbourInfo_v,const vector<int>& rank){
			node* curr_node=&_nodes[v];
			//remove the duplicated neighbors
			int nbr,nbr_w;
			unordered_map<int,pair<int,int> > nbr_map;
			for(size_t i=0;i<r_neighbourInfo_v.size();++i){
				nbr=r_neighbourInfo_v[i].first;
				nbr_w=r_neighbourInfo_v[i].second;
				nbr_map[nbr]=make_pair(nbr_w,INF_WEIGHT);
			}
			for(size_t i=0;i<neighbourInfo_v.size();++i){
				nbr=neighbourInfo_v[i].first;
				nbr_w=neighbourInfo_v[i].second;
				if(nbr_map.find(nbr)==nbr_map.end()) nbr_map[nbr]=make_pair(INF_WEIGHT,nbr_w);
				else nbr_map[nbr].second=nbr_w;
			}
			int width=nbr_map.size();
			//sort by rank
			vector<pair<int,int> > rank_in_tmp;
			vector<int > rank_border_tmp;
			for(unordered_map<int,pair<int,int> >::iterator it=nbr_map.begin();it!=nbr_map.end();++it){
				int rank_v=rank[it->first];
				if(rank_v!=-1) rank_in_tmp.push_back(make_pair(rank_v,it->first));
				else rank_border_tmp.push_back(it->first);
			}
			curr_node->_width=width;
			curr_node->_neighbour=new int[width];
			curr_node->_neighbourDis=new int[width];
			curr_node->_r_neighbourDis=new int[width];
			curr_node->_neighbourPointer=new node*[width];
			curr_node->_nBorder_siz=rank_border_tmp.size();//added
			//set root
			if(rank_in_tmp.size()==0){
				curr_node->_parentNum=-1;
				curr_node->_width=width;
				curr_node->_height=0;
				vector<int> _borders(width);
				for(int i=0;i<width;++i){
					nbr=rank_border_tmp[i];
					_borders[i]=nbr;
					curr_node->_neighbour[i]=nbr;
					curr_node->_r_neighbourDis[i]=nbr_map[nbr].first;
					curr_node->_neighbourDis[i]=nbr_map[nbr].second;
					curr_node->_neighbourPointer[i]=&_nodes[nbr];
				}
				_trees.push_back(curr_node);
				_trees_borders.push_back(_borders);
				++_numOfTrees;
			}else{
				sort(rank_in_tmp.rbegin(),rank_in_tmp.rend());//sort
				int i;
				for(i=0;i<rank_border_tmp.size();++i){
					nbr=rank_border_tmp[i];
					curr_node->_neighbour[i]=nbr;
					curr_node->_r_neighbourDis[i]=nbr_map[nbr].first;
					curr_node->_neighbourDis[i]=nbr_map[nbr].second;
					curr_node->_neighbourPointer[i]=&_nodes[nbr];
				}
				for(;i<width;++i){
					nbr=rank_in_tmp[i-rank_border_tmp.size()].second;
					curr_node->_neighbour[i]=nbr;
					curr_node->_r_neighbourDis[i]=nbr_map[nbr].first;
					curr_node->_neighbourDis[i]=nbr_map[nbr].second;
					curr_node->_neighbourPointer[i]=&_nodes[nbr];
				}
				curr_node->_parentNum=curr_node->_neighbour[width-1];//the smallest order value
				curr_node->_parent=&_nodes[curr_node->_parentNum];
				//set the current node as childs to its parent node
				curr_node->_parent->_childs.push_back(curr_node);
			}
		}

	//Single thread
		void compute_and_combinePLL_tree_top_down(int tree_index,Integrated_Index& integrated_index){
			node* root=_trees[tree_index];
			root->_tree_id=tree_index;
			//variables
			int tree_num_of_nodes=0;
			vector<vector<node*> > ans_list;
			int height=0;
			node* curr_node;
			int k,new_k;
			size_t borders_siz=_trees_borders[tree_index].size();
			vector<node*> ans_initial(borders_siz);
			//compute borders dis and pos
			for(int i=0;i<borders_siz;++i){
				int border_id=_trees_borders[tree_index][i];
				curr_node=&_nodes[border_id];
				ans_initial[i]=curr_node;
				curr_node->_dis_height=i;
				//release pointer space
				if(curr_node->_pos!=NULL) delete []curr_node->_pos;
				if(curr_node->_dis!=NULL) delete []curr_node->_dis;
				curr_node->_pos=new int[i+1];
				curr_node->_dis=new int[i+1];
				//int  dpg_distance=INF_WEIGHT;
				for(int j=0;j<i;j++){
					curr_node->_pos[j]=j;
					int border_ans_id=_trees_borders[tree_index][j];
					// for(list<edge>::const_iterator iter=es[border_id].begin();iter!=es[border_id].end();++iter){
					// 	if(iter->v==border_ans_id){
					// 		dpg_distance=iter->weight;
					// 		break;
					// 	}
					// }
					int pll_distance=integrated_index.query_p(border_ans_id,border_id);
					//if(dpg_distance!=pll_distance) std::cout<<"dpg_distance="<<dpg_distance<<" pll_distance="<<pll_distance<<std::endl;
					curr_node->_dis[j]=pll_distance;
				}
				curr_node->_pos[i]=i;
				curr_node->_dis[i]=0;
			}
			ans_list.push_back(ans_initial);
			//bfs search
			queue<pair<node*,int> > nqueue;
			nqueue.push(make_pair(root,0));
			nqueue.push(make_pair(nullptr,-1)); 
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
				vector<vector<node*> > ans_tmp;
				while(nqueue.front().first!=nullptr){
					pair<node*,int> front_pair=nqueue.front();
					curr_node=front_pair.first;
					//add tree_id to node
					curr_node->_tree_id=tree_index;
					k=front_pair.second;
					nqueue.pop();
					++tree_num_of_nodes;
					vector<node*> curr_ans=ans_list[k];
					curr_node->compute_diis(curr_ans,height,borders_siz);//add the borders size
					curr_ans.push_back(curr_node);
					ans_tmp.push_back(curr_ans);
					new_k=ans_tmp.size()-1;
					for(size_t i=0;i<curr_node->_childs.size();++i){
						nqueue.push(make_pair(curr_node->_childs[i],new_k));
					}
				}
				++height;
				nqueue.pop();
				nqueue.push(make_pair(nullptr,-1));
				ans_list.swap(ans_tmp);
			}
			if(_trees_num[tree_index]==-1||_trees_num[tree_index]!=tree_num_of_nodes){
				_trees_num[tree_index]=tree_num_of_nodes;
				//std::cout<<"error:_trees_num[tree_index]="<<_trees_num[tree_index]<<" tree_num_of_nodes="<<tree_num_of_nodes<<std::endl;//to be deleted
			}
		}

		//Single thread, add for directed
		void compute_and_combinePLL_tree_top_down_directed(int tree_index,Integrated_Index& integrated_index){
			node* root=_trees[tree_index];
			root->_tree_id=tree_index;
			//variables
			int tree_num_of_nodes=0;
			vector<vector<node*> > ans_list;
			int height=0;
			node* curr_node;
			int k,new_k;
			size_t borders_siz=_trees_borders[tree_index].size();
			vector<node*> ans_initial(borders_siz);
			//compute borders dis and pos
			for(int i=0;i<borders_siz;++i){
				int border_id=_trees_borders[tree_index][i];
				curr_node=&_nodes[border_id];
				ans_initial[i]=curr_node;
				curr_node->_dis_height=i;
				//release pointer space
				if(curr_node->_pos!=NULL) delete []curr_node->_pos;
				if(curr_node->_dis!=NULL) delete []curr_node->_dis;
				if(curr_node->_r_dis!=NULL) delete []curr_node->_r_dis;
				curr_node->_pos=new int[i+1];
				curr_node->_dis=new int[i+1];
				curr_node->_r_dis=new int[i+1];
				//int  dpg_distance=INF_WEIGHT;
				for(int j=0;j<i;j++){
					curr_node->_pos[j]=j;
					int border_ans_id=_trees_borders[tree_index][j];
					curr_node->_dis[j]=integrated_index.query_p_directed(border_id,border_ans_id);
					curr_node->_r_dis[j]=integrated_index.query_p_directed(border_ans_id,border_id);
				}
				curr_node->_pos[i]=i;
				curr_node->_dis[i]=0;
				curr_node->_r_dis[i]=0;
			}
			ans_list.push_back(ans_initial);
			//bfs search
			queue<pair<node*,int> > nqueue;
			nqueue.push(make_pair(root,0));
			nqueue.push(make_pair(nullptr,-1)); 
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
				vector<vector<node*> > ans_tmp;
				while(nqueue.front().first!=nullptr){
					pair<node*,int> front_pair=nqueue.front();
					curr_node=front_pair.first;
					//add tree_id to node
					curr_node->_tree_id=tree_index;
					k=front_pair.second;
					nqueue.pop();
					++tree_num_of_nodes;
					vector<node*> curr_ans=ans_list[k];
					curr_node->compute_diis_directed(curr_ans,height,borders_siz);//add the borders size
					curr_ans.push_back(curr_node);
					ans_tmp.push_back(curr_ans);
					new_k=ans_tmp.size()-1;
					for(size_t i=0;i<curr_node->_childs.size();++i){
						nqueue.push(make_pair(curr_node->_childs[i],new_k));
					}
				}
				++height;
				nqueue.pop();
				nqueue.push(make_pair(nullptr,-1));
				ans_list.swap(ans_tmp);
			}
			if(_trees_num[tree_index]==-1||_trees_num[tree_index]!=tree_num_of_nodes){
				_trees_num[tree_index]=tree_num_of_nodes;
				//std::cout<<"error:_trees_num[tree_index]="<<_trees_num[tree_index]<<" tree_num_of_nodes="<<tree_num_of_nodes<<std::endl;//to be deleted
			}
		}

		//multiple threads version
		void compute_and_combinePLL_tree_top_down(int tree_index,const vector<integrated_index_t>& index_){
			node* root=_trees[tree_index];
			root->_tree_id=tree_index;
			//variables
			int tree_num_of_nodes=0;
			vector<vector<node*> > ans_list;
			int height=0;
			node* curr_node;
			int k,new_k;
			size_t borders_siz=_trees_borders[tree_index].size();
			//new borders node vector
			vector<node> borders_nodes_vec(borders_siz);
			vector<node*> ans_initial(borders_siz);
			//compute borders dis and pos
			for(int i=0;i<borders_siz;++i){
				int border_id=_trees_borders[tree_index][i];
				borders_nodes_vec[i]._curNum=border_id;
				curr_node=&borders_nodes_vec[i];
				ans_initial[i]=curr_node;
				curr_node->_dis_height=i;
				curr_node->_pos=new int[i+1];
				curr_node->_dis=new int[i+1];
				//int  dpg_distance=INF_WEIGHT;
				for(int j=0;j<i;j++){
					curr_node->_pos[j]=j;
					int border_ans_id=_trees_borders[tree_index][j];
					const integrated_index_t &idx_s=index_[border_ans_id];
					const integrated_index_t &idx_t=index_[border_id];
					//int pll_distance=integrated_index.query_p(border_ans_id,border_id);
					int pll_distance=query_pll(idx_s,idx_t);
					curr_node->_dis[j]=pll_distance;
				}
				curr_node->_pos[i]=i;
				curr_node->_dis[i]=0;
			}
			ans_list.push_back(ans_initial);
			//bfs search
			queue<pair<node*,int> > nqueue;
			nqueue.push(make_pair(root,0));
			nqueue.push(make_pair(nullptr,-1)); 
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
				vector<vector<node*> > ans_tmp;
				while(nqueue.front().first!=nullptr){
					pair<node*,int> front_pair=nqueue.front();
					curr_node=front_pair.first;
					//add tree_id to node
					curr_node->_tree_id=tree_index;
					k=front_pair.second;
					nqueue.pop();
					++tree_num_of_nodes;
					vector<node*> curr_ans=ans_list[k];
					curr_node->compute_diis_multiThreads(tree_index,curr_ans,height,borders_siz);//add the borders size
					curr_ans.push_back(curr_node);
					ans_tmp.push_back(curr_ans);
					new_k=ans_tmp.size()-1;
					for(size_t i=0;i<curr_node->_childs.size();++i){
						nqueue.push(make_pair(curr_node->_childs[i],new_k));
					}
				}
				++height;
				nqueue.pop();
				nqueue.push(make_pair(nullptr,-1));
				ans_list.swap(ans_tmp);
			}
			if(_trees_num[tree_index]==-1||_trees_num[tree_index]!=tree_num_of_nodes){
				_trees_num[tree_index]=tree_num_of_nodes;
				//std::cout<<"error:_trees_num[tree_index]="<<_trees_num[tree_index]<<" tree_num_of_nodes="<<tree_num_of_nodes<<std::endl;//to be deleted
			}
		}

		//multiple threads version, add for directed
		void compute_and_combinePLL_tree_top_down_directed(int tree_index,const vector<integrated_bindex_t>& bindex_){
			node* root=_trees[tree_index];
			root->_tree_id=tree_index;
			//variables
			int tree_num_of_nodes=0;
			vector<vector<node*> > ans_list;
			int height=0;
			node* curr_node;
			int k,new_k;
			size_t borders_siz=_trees_borders[tree_index].size();
			//new borders node vector
			vector<node> borders_nodes_vec(borders_siz);
			vector<node*> ans_initial(borders_siz);
			//compute borders dis and pos
			for(int i=0;i<borders_siz;++i){
				int border_id=_trees_borders[tree_index][i];
				borders_nodes_vec[i]._curNum=border_id;
				curr_node=&borders_nodes_vec[i];
				ans_initial[i]=curr_node;
				curr_node->_dis_height=i;
				curr_node->_pos=new int[i+1];
				curr_node->_dis=new int[i+1];
				curr_node->_r_dis=new int[i+1];
				//int  dpg_distance=INF_WEIGHT;
				for(int j=0;j<i;j++){
					curr_node->_pos[j]=j;
					int border_ans_id=_trees_borders[tree_index][j];
					const integrated_bindex_t &idx_s=bindex_[border_ans_id];
					const integrated_bindex_t &idx_t=bindex_[border_id];
					//int pll_distance=integrated_index.query_p(border_ans_id,border_id);
					curr_node->_dis[j]=query_pll_directed(idx_t,idx_s);
					curr_node->_r_dis[j]=query_pll_directed(idx_s,idx_t);
				}
				curr_node->_pos[i]=i;
				curr_node->_dis[i]=0;
				curr_node->_r_dis[i]=0;
			}
			ans_list.push_back(ans_initial);
			//bfs search
			queue<pair<node*,int> > nqueue;
			nqueue.push(make_pair(root,0));
			nqueue.push(make_pair(nullptr,-1)); 
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
				vector<vector<node*> > ans_tmp;
				while(nqueue.front().first!=nullptr){
					pair<node*,int> front_pair=nqueue.front();
					curr_node=front_pair.first;
					//add tree_id to node
					curr_node->_tree_id=tree_index;
					k=front_pair.second;
					nqueue.pop();
					++tree_num_of_nodes;
					vector<node*> curr_ans=ans_list[k];
					curr_node->compute_diis_multiThreads_directed(tree_index,curr_ans,height,borders_siz);//add the borders size
					curr_ans.push_back(curr_node);
					ans_tmp.push_back(curr_ans);
					new_k=ans_tmp.size()-1;
					for(size_t i=0;i<curr_node->_childs.size();++i){
						nqueue.push(make_pair(curr_node->_childs[i],new_k));
					}
				}
				++height;
				nqueue.pop();
				nqueue.push(make_pair(nullptr,-1));
				ans_list.swap(ans_tmp);
			}
			if(_trees_num[tree_index]==-1||_trees_num[tree_index]!=tree_num_of_nodes){
				_trees_num[tree_index]=tree_num_of_nodes;
				//std::cout<<"error:_trees_num[tree_index]="<<_trees_num[tree_index]<<" tree_num_of_nodes="<<tree_num_of_nodes<<std::endl;//to be deleted
			}
		}

		/**
		 * @description: single thread for the third solution
		 * @param {int} tree_index
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void compute_and_pushPLL_tree_top_down(int tree_index,Integrated_Index& integrated_index){
			node* root=_trees[tree_index];
			//variables
			int tree_num_of_nodes=0;
			vector<vector<node*> > ans_list;
			int height=0;
			node* curr_node;
			int k,new_k;
			size_t borders_siz=_trees_borders[tree_index].size();
			vector<node*> ans_initial(borders_siz);
			//compute borders dis and pos
			for(int i=0;i<borders_siz;++i){
				int border_id=_trees_borders[tree_index][i];
				curr_node=&_nodes[border_id];
				ans_initial[i]=curr_node;
				curr_node->_dis_height=i;
				//release pointer space
				if(curr_node->_pos!=NULL) delete []curr_node->_pos;
				if(curr_node->_dis!=NULL) delete []curr_node->_dis;
				curr_node->_pos=new int[i+1];
				curr_node->_dis=new int[i+1];
				//int  dpg_distance=INF_WEIGHT;
				for(int j=0;j<i;j++){
					curr_node->_pos[j]=j;
					int border_ans_id=_trees_borders[tree_index][j];
					int pll_distance=integrated_index.query_p(border_ans_id,border_id);
					curr_node->_dis[j]=pll_distance;
				}
				curr_node->_pos[i]=i;
				curr_node->_dis[i]=0;
			}
			ans_list.push_back(ans_initial);
			//bfs search
			queue<pair<node*,int> > nqueue;
			nqueue.push(make_pair(root,0));
			nqueue.push(make_pair(nullptr,-1)); 
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
				vector<vector<node*> > ans_tmp;
				while(nqueue.front().first!=nullptr){
					pair<node*,int> front_pair=nqueue.front();
					curr_node=front_pair.first;
					//add tree_id to node
					curr_node->_tree_id=tree_index;
					k=front_pair.second;
					nqueue.pop();
					++tree_num_of_nodes;
					vector<node*> curr_ans=ans_list[k];
					curr_node->compute_diis_and_pushPLL(integrated_index,_trees_borders[tree_index],curr_ans,height,borders_siz);
					curr_ans.push_back(curr_node);
					ans_tmp.push_back(curr_ans);
					new_k=ans_tmp.size()-1;
					for(size_t i=0;i<curr_node->_childs.size();++i){
						nqueue.push(make_pair(curr_node->_childs[i],new_k));
					}
				}
				++height;
				nqueue.pop();
				nqueue.push(make_pair(nullptr,-1));
				ans_list.swap(ans_tmp);
			}
			if(_trees_num[tree_index]==-1||_trees_num[tree_index]!=tree_num_of_nodes){
				_trees_num[tree_index]=tree_num_of_nodes;
				//std::cout<<"error:_trees_num[tree_index]="<<_trees_num[tree_index]<<" tree_num_of_nodes="<<tree_num_of_nodes<<std::endl;//to be deleted
			}
		}

		/**
		 * @description: single thread for the third solution, add for directed
		 * @param {int} tree_index
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void compute_and_pushPLL_tree_top_down_directed(int tree_index,Integrated_Index& integrated_index){
			node* root=_trees[tree_index];
			//variables
			int tree_num_of_nodes=0;
			vector<vector<node*> > ans_list;
			int height=0;
			node* curr_node;
			int k,new_k;
			size_t borders_siz=_trees_borders[tree_index].size();
			vector<node*> ans_initial(borders_siz);
			//compute borders dis and pos
			for(int i=0;i<borders_siz;++i){
				int border_id=_trees_borders[tree_index][i];
				curr_node=&_nodes[border_id];
				ans_initial[i]=curr_node;
				curr_node->_dis_height=i;
				//release pointer space
				if(curr_node->_pos!=NULL) delete []curr_node->_pos;
				if(curr_node->_dis!=NULL) delete []curr_node->_dis;
				if(curr_node->_r_dis!=NULL) delete []curr_node->_r_dis;
				curr_node->_pos=new int[i+1];
				curr_node->_dis=new int[i+1];
				curr_node->_r_dis=new int[i+1];
				//int  dpg_distance=INF_WEIGHT;
				for(int j=0;j<i;j++){
					curr_node->_pos[j]=j;
					int border_ans_id=_trees_borders[tree_index][j];
					curr_node->_dis[j]=integrated_index.query_p_directed(border_id,border_ans_id);
					curr_node->_r_dis[j]=integrated_index.query_p_directed(border_ans_id,border_id);
				}
				curr_node->_pos[i]=i;
				curr_node->_dis[i]=0;
				curr_node->_r_dis[i]=0;
			}
			ans_list.push_back(ans_initial);
			//bfs search
			queue<pair<node*,int> > nqueue;
			nqueue.push(make_pair(root,0));
			nqueue.push(make_pair(nullptr,-1)); 
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
				vector<vector<node*> > ans_tmp;
				while(nqueue.front().first!=nullptr){
					pair<node*,int> front_pair=nqueue.front();
					curr_node=front_pair.first;
					//add tree_id to node
					curr_node->_tree_id=tree_index;
					k=front_pair.second;
					nqueue.pop();
					++tree_num_of_nodes;
					vector<node*> curr_ans=ans_list[k];
					curr_node->compute_diis_and_pushPLL_directed(integrated_index,_trees_borders[tree_index],curr_ans,height,borders_siz);
					curr_ans.push_back(curr_node);
					ans_tmp.push_back(curr_ans);
					new_k=ans_tmp.size()-1;
					for(size_t i=0;i<curr_node->_childs.size();++i){
						nqueue.push(make_pair(curr_node->_childs[i],new_k));
					}
				}
				++height;
				nqueue.pop();
				nqueue.push(make_pair(nullptr,-1));
				ans_list.swap(ans_tmp);
			}
			if(_trees_num[tree_index]==-1||_trees_num[tree_index]!=tree_num_of_nodes){
				_trees_num[tree_index]=tree_num_of_nodes;
				//std::cout<<"error:_trees_num[tree_index]="<<_trees_num[tree_index]<<" tree_num_of_nodes="<<tree_num_of_nodes<<std::endl;//to be deleted
			}
		}

		//multiple threads version
		void compute_and_pushPLL_tree_top_down(int tree_index,const vector<integrated_index_t>& index_,int max_overlay_hub){
			node* root=_trees[tree_index];
			//variables
			int tree_num_of_nodes=0;
			vector<vector<node*> > ans_list;
			int height=0;
			node* curr_node;
			int k,new_k;
			size_t borders_siz=_trees_borders[tree_index].size();
			//new borders node vector
			vector<node> borders_nodes_vec(borders_siz);
			vector<node*> ans_initial(borders_siz);
			//compute borders dis and pos
			for(int i=0;i<borders_siz;++i){
				int border_id=_trees_borders[tree_index][i];
				borders_nodes_vec[i]._curNum=border_id;
				curr_node=&borders_nodes_vec[i];
				ans_initial[i]=curr_node;
				curr_node->_dis_height=i;
				curr_node->_pos=new int[i+1];
				curr_node->_dis=new int[i+1];
				//int  dpg_distance=INF_WEIGHT;
				for(int j=0;j<i;j++){
					curr_node->_pos[j]=j;
					int border_ans_id=_trees_borders[tree_index][j];
					const integrated_index_t &idx_s=index_[border_ans_id];
					const integrated_index_t &idx_t=index_[border_id];
					int pll_distance=query_pll(idx_s,idx_t);
					curr_node->_dis[j]=pll_distance;
				}
				curr_node->_pos[i]=i;
				curr_node->_dis[i]=0;
			}
			ans_list.push_back(ans_initial);
			//bfs search
			queue<pair<node*,int> > nqueue;
			nqueue.push(make_pair(root,0));
			nqueue.push(make_pair(nullptr,-1)); 
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
				vector<vector<node*> > ans_tmp;
				while(nqueue.front().first!=nullptr){
					pair<node*,int> front_pair=nqueue.front();
					curr_node=front_pair.first;
					//add tree_id to node
					curr_node->_tree_id=tree_index;
					k=front_pair.second;
					nqueue.pop();
					++tree_num_of_nodes;
					vector<node*> curr_ans=ans_list[k];
					curr_node->compute_diis_and_pushPLL_multiThreads(tree_index,index_,_trees_borders[tree_index],curr_ans,height,borders_siz,max_overlay_hub);
					curr_ans.push_back(curr_node);
					ans_tmp.push_back(curr_ans);
					new_k=ans_tmp.size()-1;
					for(size_t i=0;i<curr_node->_childs.size();++i){
						nqueue.push(make_pair(curr_node->_childs[i],new_k));
					}
				}
				++height;
				nqueue.pop();
				nqueue.push(make_pair(nullptr,-1));
				ans_list.swap(ans_tmp);
			}
			if(_trees_num[tree_index]==-1||_trees_num[tree_index]!=tree_num_of_nodes){
				_trees_num[tree_index]=tree_num_of_nodes;
				//std::cout<<"error:_trees_num[tree_index]="<<_trees_num[tree_index]<<" tree_num_of_nodes="<<tree_num_of_nodes<<std::endl;//to be deleted
			}
		}

		//multiple threads version for the third solution, add for directed
		void compute_and_pushPLL_tree_top_down_directed(int tree_index,const vector<integrated_bindex_t>& bindex_,int max_overlay_hub_r,int max_overlay_hub){
			node* root=_trees[tree_index];
			//variables
			int tree_num_of_nodes=0;
			vector<vector<node*> > ans_list;
			int height=0;
			node* curr_node;
			int k,new_k;
			size_t borders_siz=_trees_borders[tree_index].size();
			//new borders node vector
			vector<node> borders_nodes_vec(borders_siz);
			vector<node*> ans_initial(borders_siz);
			//compute borders dis and pos
			for(int i=0;i<borders_siz;++i){
				int border_id=_trees_borders[tree_index][i];
				borders_nodes_vec[i]._curNum=border_id;
				curr_node=&borders_nodes_vec[i];
				ans_initial[i]=curr_node;
				curr_node->_dis_height=i;
				curr_node->_pos=new int[i+1];
				curr_node->_dis=new int[i+1];
				curr_node->_r_dis=new int[i+1];
				//int  dpg_distance=INF_WEIGHT;
				for(int j=0;j<i;j++){
					curr_node->_pos[j]=j;
					int border_ans_id=_trees_borders[tree_index][j];
					const integrated_bindex_t &idx_s=bindex_[border_ans_id];
					const integrated_bindex_t &idx_t=bindex_[border_id];
					curr_node->_dis[j]=query_pll_directed(idx_t,idx_s);
					curr_node->_r_dis[j]=query_pll_directed(idx_s,idx_t);
				}
				curr_node->_pos[i]=i;
				curr_node->_dis[i]=0;
				curr_node->_r_dis[i]=0;
			}
			ans_list.push_back(ans_initial);
			//bfs search
			queue<pair<node*,int> > nqueue;
			nqueue.push(make_pair(root,0));
			nqueue.push(make_pair(nullptr,-1)); 
			while(!nqueue.empty()){
				if(nqueue.size()==1&&nqueue.front().first==nullptr) break;
				vector<vector<node*> > ans_tmp;
				while(nqueue.front().first!=nullptr){
					pair<node*,int> front_pair=nqueue.front();
					curr_node=front_pair.first;
					//add tree_id to node
					curr_node->_tree_id=tree_index;
					k=front_pair.second;
					nqueue.pop();
					++tree_num_of_nodes;
					vector<node*> curr_ans=ans_list[k];
					curr_node->compute_diis_and_pushPLL_multiThreads_directed(tree_index,bindex_,_trees_borders[tree_index],curr_ans,height,borders_siz,max_overlay_hub_r,max_overlay_hub);
					curr_ans.push_back(curr_node);
					ans_tmp.push_back(curr_ans);
					new_k=ans_tmp.size()-1;
					for(size_t i=0;i<curr_node->_childs.size();++i){
						nqueue.push(make_pair(curr_node->_childs[i],new_k));
					}
				}
				++height;
				nqueue.pop();
				nqueue.push(make_pair(nullptr,-1));
				ans_list.swap(ans_tmp);
			}
			if(_trees_num[tree_index]==-1||_trees_num[tree_index]!=tree_num_of_nodes){
				_trees_num[tree_index]=tree_num_of_nodes;
				//std::cout<<"error:_trees_num[tree_index]="<<_trees_num[tree_index]<<" tree_num_of_nodes="<<tree_num_of_nodes<<std::endl;//to be deleted
			}
		}

		void dfs_preprocess_forest(node* curr_node){
			int id=curr_node->_curNum;
			_euler.push_back(id);
			_level.push_back(curr_node->_height);
			if(_index[id]==-1) _index[id]=_euler.size()-1;
			for(size_t i=0;i<curr_node->_childs.size();++i)
			{
				dfs_preprocess_forest(curr_node->_childs[i]);
				_euler.push_back(id);
				_level.push_back(curr_node->_height);
			}
		}

};

/*
 *@description: construct h2h-Indexing class
 *@author: gaoyongyong
 *@date: 2021-02-09
*/
class H2H_construction{
public:
        //********************variables********************
        vector<list<edge> > es;//update overlay_grph
		vector<list<edge> > es_original;//original graph

        vector<list<edge> > r_es;//update overlay_grph for directed
		vector<list<edge> > r_es_original;//original graph for directed
	
		vector<NodeID> HFPoint;//store the HFpoint
		vector<bool> HFPointFlag;//hfpoint flag by node index
		
        vector<vector<int> > adj;//store the adjacent vertices
		vector< vector<int> > adj_weight;//store the adjacent weight

        vector<vector<int> > r_adj;//store the adjacent vertices for directed
		vector< vector<int> > r_adj_weight;//store the adjacent weight for directed
        
        int numOfHfpoint;//original num of HFpoints
		int numOfOriginalVertices{0};
		int numOfOriginalEdges{0};

		long long total_query_time;//total query time for overlay left vertices
		int max_query_time;//max querytime for overlay left vertices
		vector<unsigned int> query_time;//get query time for overlay left vertices

		vector<int> dp_rank;//fetch the delete order calue by index
		vector<int> dp_inv;//fetach the vertex index by order value
		vector<vector<pair<int,int> > > neighborsInfo;//store the node neighbors
		vector<vector<pair<int,int> > > r_neighborsInfo;//store the node neighbors for directed

		tree _dp_tree;
		forest _dp_forest;

		double time_delete,avg_time_delete;
		double time_build_tree,avg_time_build_tree;
		double time_compute,avg_time_compute;
		double time_rmq,avg_time_rmq;
		double time_index,avg_time_index;

        //****************constructions and deconstructions**************
		 /**
		* @description: used for the third integration solution
		* @param {*}
		* @return {*}
		* @author: Wan Jingyi
		*/   
		H2H_construction(WGraph& overlay_wgraph,char* graphFileName,const vector<bool>& isHFPoint,
												vector<bool>& isDeleted,vector<int>& newToOriginal,vector<int>& originalToNew,
												int& numOfOverlayVertices,int& numOfDeletedVertices,int& numOfOverlayEdges,
												int thresholdDegree,int bestDegree)
		{
			//load graph
			load_original_graph(graphFileName);
			//delete node
			create_node_hf_overlay(overlay_wgraph,newToOriginal,originalToNew,isDeleted,isHFPoint,numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,thresholdDegree,bestDegree);
			std::cout<<"******************create_node_hf_overlay successfully!*****************"<<std::endl;
			//build tree
			create_forest(numOfOriginalVertices,numOfDeletedVertices);
			std::cout<<"********************Create_forest finished!**********************"<<std::endl;
		}

		 /**
		* @description: used for the third integration solution,add for directed
		* @param {*}
		* @return {*}
		* @author: Wan Jingyi
		*/   
		H2H_construction(WGraph& overlay_wgraph,char* graphFileName,const vector<bool>& isHFPoint,
												vector<bool>& isDeleted,vector<int>& newToOriginal,vector<int>& originalToNew,
												int& numOfOverlayVertices,int& numOfDeletedVertices,int& numOfOverlayEdges,
												int graphType,int thresholdDegree,int bestDegree)
		{
			//load graph
			load_original_graph(graphFileName,graphType);
			//delete node
			create_node_hf_overlay_directed(overlay_wgraph,newToOriginal,originalToNew,isDeleted,isHFPoint,numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,thresholdDegree,bestDegree);
			std::cout<<"******************create_node_hf_overlay successfully!*****************"<<std::endl;
			//build tree
			create_forest_directed(numOfOriginalVertices,numOfDeletedVertices);
			std::cout<<"********************Create_forest_directed finished!**********************"<<std::endl;
		}

		 /**
		* @description: used for the second integration solution
		* @param {*}
		* @return {*}
		* @author: Wan Jingyi
		*/   
		H2H_construction(char* graphFileName,WGraph& overlay_wgraph,const vector<bool>& isHFPoint,
												vector<bool>& isDeleted,vector<int>& newToOriginal,vector<int>& originalToNew,
												int& numOfOverlayVertices,int& numOfDeletedVertices,int& numOfOverlayEdges,
												int thresholdDegree,const vector<unsigned int>& queryTime,int bestDegree,
												bool isSorted=false)
		{
			//load graph
			load_original_graph(graphFileName);
			//delete node
			create_node_hf_overlay(overlay_wgraph,newToOriginal,originalToNew,isDeleted,isHFPoint,numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,thresholdDegree,bestDegree);
			std::cout<<"******************create_node_hf_overlay successfully!*****************"<<std::endl;
			//build tree
			create_forest(numOfOriginalVertices,numOfDeletedVertices,queryTime,isSorted);
			std::cout<<"********************Create_forest finished!**********************"<<std::endl;
		}

		/**
	 * @description: used for the second integration solution, Add for directed
	 * @param {*}
	 * @return {*}
	 * @author: Wan Jingyi
	 */  
  		H2H_construction(char* graphFileName,WGraph& overlay_wgraph,const vector<bool>& isHFPoint,
												vector<bool>& isDeleted,vector<int>& newToOriginal,vector<int>& originalToNew,
												int& numOfOverlayVertices,int& numOfDeletedVertices,int& numOfOverlayEdges,
												int thresholdDegree,const vector<unsigned int>& queryTime,int bestDegree,
												int graphType,bool isSorted=false)//delete outputDirName
		{
			//load graph
			load_original_graph(graphFileName,graphType);
			//delete node
			create_node_hf_overlay_directed(overlay_wgraph,newToOriginal,originalToNew,isDeleted,isHFPoint,numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,thresholdDegree,bestDegree);
			std::cout<<"******************create_node_hf_overlay_directed successfully!*****************"<<std::endl;
			//build tree
			create_forest_directed(numOfOriginalVertices,numOfDeletedVertices,queryTime,isSorted);
			std::cout<<"********************Create_forest_directed finished!**********************"<<std::endl;
		}

		 /**
		* @description: used for the first integration solution
		* @param {*}
		* @return {*}
		* @author: Wan Jingyi
		*/   
   		H2H_construction(vector<integrated_index_t>& index_,WGraph& overlay_wgraph,char* graphFileName,
		   										const vector<bool>& isHFPoint,double& time_tree_index,double& time_compute_rmq,
												vector<bool>& isDeleted,vector<int>& newToOriginal,vector<int>& originalToNew,
												int& numOfOverlayVertices,int& numOfDeletedVertices,int& numOfOverlayEdges,
												int thresholdDegree)
		{
			//load graph
			load_original_graph(graphFileName);
            //delete node
			create_node_hf_all(overlay_wgraph,newToOriginal,originalToNew,isDeleted,isHFPoint,numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,thresholdDegree);
			//build tree
			time_build_tree=GetCurrentTimeSec();
			create_tree();
			time_build_tree=GetCurrentTimeSec()-time_build_tree;
			//release memory
			vector<vector<pair<int,int> > >().swap(neighborsInfo);
			neighborsInfo.clear();
			//compute tree
			time_compute=GetCurrentTimeSec();
			_dp_tree.compute_tree_top_down();
			time_compute=GetCurrentTimeSec()-time_compute;
			//compute rmq variables needed
			time_rmq=GetCurrentTimeSec();
			_dp_tree.compute_rmq_variables();
			time_rmq=GetCurrentTimeSec()-time_rmq;
			time_index=time_delete+time_build_tree+time_compute;
			//merge index(".h2hIndex")
			_dp_tree.compute_node_label_size();
			_dp_tree.merge_h2h_index(index_,isDeleted);
			//output time
			avg_time_index=time_index/(double) numOfVertices;
			avg_time_delete=time_delete/(double) numOfVertices;
			avg_time_build_tree=time_build_tree/(double) numOfVertices;
			avg_time_compute=time_compute/(double) numOfVertices;
			avg_time_rmq=time_rmq/(double) numOfVertices;
			cout << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
			cout << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
			cout << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
			cout << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
			cout << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
			time_compute_rmq=time_rmq;
			time_tree_index=time_index;
		}

        ~H2H_construction(){}

        //****************functions**************
protected:
	void create_node_by_input_order(char *inputOrderFile){
		int u,v,w;//temp variables u-start v-end
		int cnt=0;
		int curr_id,curr_deg;
		edge temp;//edge temp variable
		dp_rank.resize(numOfVertices);
		dp_inv.resize(numOfVertices);
		neighborsInfo.resize(numOfVertices);
		ifstream in(inputOrderFile);
		if(!in.is_open()) cout<<"Cannot open "<<inputOrderFile<<endl;
		dp_rank.resize(numOfVertices);
		for(int i=0;i<numOfVertices;++i){
			in>>u;
			dp_rank[u]=i;
			dp_inv[i]=u;
		}
		in.close();
		//for each vertex to construct node
		for(int r=0;r<numOfVertices;++r){
			curr_id=dp_inv[r];
			curr_deg=es[curr_id].size();
			++cnt;
			if(curr_deg==0){
				continue;
			}
			else if(curr_deg==1){
					u=es[curr_id].front().v;
					w=es[curr_id].front().weight;
					neighborsInfo[curr_id].push_back(make_pair(u,w));//store the neighbors information
					es[curr_id].clear();
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==v){
							es[u].erase(it);
							break;
						}
					}
			}else{//degree>=2
				//deal with each neighbor pair
				vector<pair<int,int> > neighbors;//store the neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
				neighborsInfo[curr_id]=neighbors;
				//u-v menas neighbor pair
				for(int i=0;i<neighbors.size()-1; i++){
					u=neighbors[i].first;
					for(int j=i+1;j<neighbors.size();j++){
						v=neighbors[j].first;
						isLink=false;
						w=neighbors[i].second+neighbors[j].second;
						for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
							//if exist edge(u,v),update the weight(u,v) to minimun
							if(iter_u->v==v){
								if(w<iter_u->weight){
									//update two directions
									iter_u->weight=w;
									for (list<edge>::iterator iter_v = es[v].begin(); iter_v != es[v].end(); iter_v++){
										if(iter_v->v==u){
											iter_v->weight=w;
											break;
										}
									}
								}
								isLink=true;
								break;
							}
						}
						//if no edge(u,v),then insert
						if(!isLink){
							temp.u=u;
							temp.v=v;
							temp.weight=w;
							es[u].push_back(temp);
							temp.u=v;
							temp.v=u;
							es[v].push_back(temp);
						}
					}
				}
				//process neighbor's edges and update dqueue
				for(int i=0;i<neighbors.size();++i){
					u=neighbors[i].first;
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
						if(it->v== v){
							es[u].erase(it);
							break;
						}
					}
				}
				//remove curr_id and its adjacent edges from H
				es[curr_id].clear();
			}
		}//while ends
	}

	/**
     * @description: create each node by ascending degree order
	 * to solve line problem, add some little revisions
     * @Author: wanjingyi
     * @Date: 2021-02-10 20:53:42
     * @param {*}
     * @return {*}
     */
    void create_node_deg_revise_line(){
		std::cout<<"****************************createNode by deg  revise line start!***************************"<<std::endl;
		//used variables
		benchmark::heap<2, int, int> dqueue(numOfVertices);
		int u,v,w;//temp variables u-start v-end
		int cnt=0;
		int curr_id,curr_deg;
		int line_end=-1;//the node adjacent to last deg1 node deleted
		edge temp;//edge temp variable
		dp_rank.resize(numOfVertices);
		dp_inv.resize(numOfVertices);
		neighborsInfo.resize(numOfVertices);
		//push degree
		for(int i=0;i<numOfVertices;++i){
				int deg=es[i].size();
				if(deg==0){
					cout<<"error: isolated vertex "<<i<<endl;
					continue;
				}
				dqueue.update(i,deg);	
		}
		//for each vertex to construct node
		while(!dqueue.empty()){
			dqueue.extract_min(curr_id,curr_deg);
			dp_rank[curr_id]=cnt;
			dp_inv[cnt]=curr_id;
			++cnt;
			if(curr_deg==0){
				continue;
			}
			else if(curr_deg==1){
					if(curr_deg>=dqueue.size()-1){//process line
						if(curr_id==line_end&&dqueue.top()==curr_deg){
							cout<<"curr_id="<<curr_id<<" curr_deg="<<curr_deg<<" line_end="<<line_end<<" next:("<<dqueue.top_value()<<","<<dqueue.top()<<")"<<endl;
							dqueue.update(curr_id,curr_deg);
							continue;
						}
						line_end=curr_id;
					}
					u=es[curr_id].front().v;
					w=es[curr_id].front().weight;
					neighborsInfo[curr_id].push_back(make_pair(u,w));//store the neighbors information
					es[curr_id].clear();
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==v){
							es[u].erase(it);
							dqueue.update(u,es[u].size());
							break;
						}
					}
			}else{//degree>=2
				line_end=-1;
				//deal with each neighbor pair
				vector<pair<int,int> > neighbors;//store the neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
				neighborsInfo[curr_id]=neighbors;
				//u-v menas neighbor pair
				for(int i=0;i<neighbors.size()-1; i++){
					u=neighbors[i].first;
					for(int j=i+1;j<neighbors.size();j++){
						v=neighbors[j].first;
						isLink=false;
						w=neighbors[i].second+neighbors[j].second;
						for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
							//if exist edge(u,v),update the weight(u,v) to minimun
							if(iter_u->v==v){
								if(w<iter_u->weight){
									//update two directions
									iter_u->weight=w;
									for (list<edge>::iterator iter_v = es[v].begin(); iter_v != es[v].end(); iter_v++){
										if(iter_v->v==u){
											iter_v->weight=w;
											break;
										}
									}
								}
								isLink=true;
								break;
							}
						}
						//if no edge(u,v),then insert
						if(!isLink){
							temp.u=u;
							temp.v=v;
							temp.weight=w;
							es[u].push_back(temp);
							temp.u=v;
							temp.v=u;
							es[v].push_back(temp);
						}
					}
				}
				//process neighbor's edges and update dqueue
				for(int i=0;i<neighbors.size();++i){
					u=neighbors[i].first;
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
						if(it->v== v){
							es[u].erase(it);
							dqueue.update(u,es[u].size());//update degree
							break;
						}
					}
				}
				//remove curr_id and its adjacent edges from H
				es[curr_id].clear();
			}
		}//while ends
		//if(cnt!=numOfVertices) std::cout<<"Error:cnt!=numOfVertices!"<<std::endl;
		//std::cout<<"****************************createNode by deg finished!***************************"<<endl;
    }

	/**
     * @description: create each node by ascending degree order
     * @Author: wanjingyi
     * @Date: 2021-02-10 20:53:42
     * @param {*}
     * @return {*}
     */
    void create_node_deg(){
		//std::cout<<"****************************createNode by deg start!***************************"<<std::endl;
		//used variables
		benchmark::heap<2, int, int> dqueue(numOfVertices);
		int u,v,w;//temp variables u-start v-end
		int cnt=0;
		int curr_id,curr_deg;
		edge temp;//edge temp variable
		dp_rank.resize(numOfVertices);
		dp_inv.resize(numOfVertices);
		neighborsInfo.resize(numOfVertices);
		//push degree
		for(int i=0;i<numOfVertices;++i){
				int deg=es[i].size();
				if(deg==0){
					cout<<"error: isolated vertex "<<i<<endl;
					continue;
				}
				dqueue.update(i,deg);	
		}
		//for each vertex to construct node
		while(!dqueue.empty()){
			dqueue.extract_min(curr_id,curr_deg);
			dp_rank[curr_id]=cnt;
			dp_inv[cnt]=curr_id;
			++cnt;
			if(curr_deg==0){
				continue;
			}
			else if(curr_deg==1){
					u=es[curr_id].front().v;
					w=es[curr_id].front().weight;
					neighborsInfo[curr_id].push_back(make_pair(u,w));//store the neighbors information
					es[curr_id].clear();
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==v){
							es[u].erase(it);
							dqueue.update(u,es[u].size());
							break;
						}
					}
			}else{//degree>=2
				//deal with each neighbor pair
				vector<pair<int,int> > neighbors;//store the neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
				neighborsInfo[curr_id]=neighbors;
				//u-v menas neighbor pair
				for(int i=0;i<neighbors.size()-1; i++){
					u=neighbors[i].first;
					for(int j=i+1;j<neighbors.size();j++){
						v=neighbors[j].first;
						isLink=false;
						w=neighbors[i].second+neighbors[j].second;
						for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
							//if exist edge(u,v),update the weight(u,v) to minimun
							if(iter_u->v==v){
								if(w<iter_u->weight){
									//update two directions
									iter_u->weight=w;
									for (list<edge>::iterator iter_v = es[v].begin(); iter_v != es[v].end(); iter_v++){
										if(iter_v->v==u){
											iter_v->weight=w;
											break;
										}
									}
								}
								isLink=true;
								break;
							}
						}
						//if no edge(u,v),then insert
						if(!isLink){
							temp.u=u;
							temp.v=v;
							temp.weight=w;
							es[u].push_back(temp);
							temp.u=v;
							temp.v=u;
							es[v].push_back(temp);
						}
					}
				}
				//process neighbor's edges and update dqueue
				for(int i=0;i<neighbors.size();++i){
					u=neighbors[i].first;
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
						if(it->v== v){
							es[u].erase(it);
							dqueue.update(u,es[u].size());//update degree
							break;
						}
					}
				}
				//remove curr_id and its adjacent edges from H
				es[curr_id].clear();
			}
		}//while ends
		//if(cnt!=numOfVertices) std::cout<<"Error:cnt!=numOfVertices!"<<std::endl;
		//std::cout<<"****************************createNode by deg finished!***************************"<<endl;
    }

	//Add for directed
	void create_node_deg_directed(){
		//used variables
		benchmark::heap<2, int, int> dqueue(numOfVertices);
		int curr_id,curr_deg,deg;
		int u,v,u_w,v_w,w;//temp variables u-start v-end
		int cnt=0;//current num of deleted nodes
		edge temp;//edge temp variable
		dp_rank.resize(numOfVertices);
		dp_inv.resize(numOfVertices);
		neighborsInfo.resize(numOfVertices);
		r_neighborsInfo.resize(numOfVertices);
		//get the initial degree
		for(int i=0;i<numOfVertices;++i){
			if(es[i].size()==0&&r_es[i].size()==0){
				cout<<"erroe: isolated vertex "<<i<<std::endl;
				continue;
			}
			deg=es[i].size()*r_es[i].size();
			if(deg==1&&es[i].front().v==r_es[i].front().u) deg=0;
			//to be deleted
			//std::cout<<i<<":es_size="<<es[i].size()<<" r_es_size="<<r_es[i].size()<<" deg="<<deg<<std::endl;
			dqueue.update(i,deg);
		}
		//for each vertex to construct node
		while(!dqueue.empty()){
			dqueue.extract_min(curr_id,curr_deg);
			dp_rank[curr_id]=cnt;
			dp_inv[cnt]=curr_id;
			++cnt;
			if(curr_deg==0){//line
				if(es[curr_id].size()==0&&r_es[curr_id].size()==0){//only left one vertice
					std::cout<<"Isolated vertex!"<<std::endl;
				}else if(es[curr_id].size()==0&&r_es[curr_id].size()!=0){
					for(list<edge>::iterator it=r_es[curr_id].begin();it!=r_es[curr_id].end();++it){
						u=r_es[curr_id].front().u;
						u_w=r_es[curr_id].front().weight;
						r_neighborsInfo[curr_id].push_back(make_pair(u,u_w));
						for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
							if(it->v==curr_id){
								es[u].erase(it);
								deg=es[u].size()*r_es[u].size();
								if(deg==1&&es[u].front().v==r_es[u].front().u) deg=0;
								dqueue.update(u,deg);
								break;
							}
						}
					}
					r_es[curr_id].clear();
				}else if(r_es[curr_id].size()==0&&es[curr_id].size()!=0){//all edges are in frontward
					for(list<edge>::iterator it=es[curr_id].begin();it!=es[curr_id].end();++it){
						v=es[curr_id].front().v;
						v_w=es[curr_id].front().weight;
						neighborsInfo[curr_id].push_back(make_pair(v,v_w));
						for (list<edge>::iterator it = r_es[v].begin(); it != r_es[v].end(); it++) {
							if(it->u==curr_id){
								r_es[v].erase(it);
								deg=es[v].size()*r_es[v].size();
								if(deg==1&&es[v].front().v==r_es[v].front().u) deg=0;
								dqueue.update(v,deg);
								break;
							}
						}
					}
					es[curr_id].clear();
				}else{
					v=es[curr_id].front().v;
					v_w=es[curr_id].front().weight;
					u=r_es[curr_id].front().u;
					u_w=r_es[curr_id].front().weight;
					r_neighborsInfo[curr_id].push_back(make_pair(u,u_w));
					neighborsInfo[curr_id].push_back(make_pair(v,v_w));
					//to be deleted
					//if(u!=v) std::cout<<"Error: deg=0 one vertice!"<<std::endl;
					for (list<edge>::iterator it = r_es[v].begin(); it != r_es[v].end(); it++) {
						if(it->u==curr_id){
							r_es[v].erase(it);
							break;
						}
					}
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==curr_id){
							es[u].erase(it);
							break;
						}
					}
					deg=es[u].size()*r_es[u].size();
					if(deg==1&&es[u].front().v==r_es[u].front().u) deg=0;
					dqueue.update(u,deg);
					es[curr_id].clear();
					r_es[curr_id].clear();
				}
			}else{
				vector<pair<int,int> > neighbors;//store the forward neighbor and its distance
				vector<pair<int,int> > r_neighbors;//store the backward neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++){
					neighbors.push_back(make_pair(it->v,it->weight));
				}
				for (list<edge>::iterator it = r_es[curr_id].begin(); it != r_es[curr_id].end(); it++){
					r_neighbors.push_back(make_pair(it->u,it->weight));
				}
				r_neighborsInfo[curr_id]=r_neighbors;
				neighborsInfo[curr_id]=neighbors;
				//u-v menas neighbor pair
				for(int i=0;i<r_neighbors.size();++i){
					u=r_neighbors[i].first;
					u_w=r_neighbors[i].second;
					for(int j=0;j<neighbors.size();++j){
						v=neighbors[j].first;
						v_w=neighbors[j].second;
						if(u==v) continue;//not need to deal
						isLink=false;
						w=u_w+v_w;
						for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
							//if exist edge(u,v),update the weight(u,v) to minimun
							if(iter_u->v==v){
								if(w<iter_u->weight){
									iter_u->weight=w;//update forward direction
									for (list<edge>::iterator iter_v = r_es[v].begin(); iter_v != r_es[v].end(); iter_v++){
										if(iter_v->u==u){
											iter_v->weight=w;//update backward direction
											break;
										}
									}
									break;
								}
								isLink=true;
								break;
							}
						}
						//if no edge(u,v),then insert
						if(!isLink&&u!=v){
							temp.u=u;
							temp.v=v;
							temp.weight=w;
							es[u].push_back(temp);
							r_es[v].push_back(temp);
						}
					}
				}
				//process neighbor's edges and update dqueue
				for(int i=0;i<r_neighbors.size();++i){
					u=r_neighbors[i].first;
					//u_w=r_neighbors[i].second;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
						if(it->v==curr_id){
							es[u].erase(it);
							deg=es[u].size()*r_es[u].size();
							if(deg==1&&es[u].front().v==r_es[u].front().u) deg=0;
							dqueue.update(u,deg);
							break;
						}
					}
				}
				for(int j=0;j<neighbors.size();++j){
					v=neighbors[j].first;
					//v_w=neighbors[j].second;
					for (list<edge>::iterator it = r_es[v].begin(); it != r_es[v].end(); it++){
						if(it->u==curr_id){
							r_es[v].erase(it);
							deg=es[v].size()*r_es[v].size();
							if(deg==1&&es[v].front().v==r_es[v].front().u) deg=0;
							dqueue.update(v,deg);
							break;
						}
					}
				}
				//remove curr_id and its adjacent edges from H
				es[curr_id].clear();
				r_es[curr_id].clear();
			}
		}

	}

	void load_hf_query_setting(int hfRate,char* HFPoint_file){
		if(hfRate==0) numOfHfpoint=numOfVertices;
		else numOfHfpoint = static_cast<int> ( (double)numOfVertices*hfRate/(double)HF_DIVIDION);
		query_time.resize(numOfVertices,0);
		HFPointFlag.resize(numOfVertices,false);
        load_hfpoint_and_qt(HFPoint_file,numOfHfpoint,HFPoint,HFPointFlag,query_time,max_query_time,total_query_time);
		if(numOfHfpoint<=0) cout<<"error:numOfHFpoint<=0"<<endl;
	}

	void create_node_hf_overlay(WGraph& overlay_wgraph,vector<int>& newToOriginal,vector<int>& originalToNew,vector<bool>& isDeleted,const vector<bool>& isHFPoint,int& numOfOverlayVertices,int& numOfDeletedVertices,int& numOfOverlayEdges,int thresholdDegree,int bestDegree){
		//initialize variables
		double time_delete_add;
		time_delete=0;
		numOfHfpoint=0;
		time_delete_add=GetCurrentTimeSec();
		dp_rank.resize(numOfVertices,-1);
		dp_inv.reserve(numOfVertices);
		neighborsInfo.resize(numOfVertices);
		benchmark::heap<2, int, int> lqueue(numOfVertices);
		int u,v,w;//temp variables u-start v-end
		edge temp;//edge temp variable
		int d_cnt=0;//current num of deleted nodes
		int curr_id,curr_deg;
		//get the initial degree
		for(int i=0;i<numOfVertices;++i){
			//judge isolated nodes
			int deg=es[i].size();
			if(deg==0){
				cout<<"error: isolated vertex "<<i<<std::endl;
				continue;
			}
			if(isHFPoint[i]){
				++numOfHfpoint;
				continue;
			}
			lqueue.update(i,deg);
		}
		numOfDeletedVertices=numOfVertices-numOfHfpoint;
		//**********************the first stage**********************
		while(!lqueue.empty()){
			//cout<<"lqueue size="<<lqueue.size()<<std::endl;//to be deleted
			if(d_cnt>=numOfDeletedVertices) break;
			lqueue.extract_min(curr_id,curr_deg);
			if(curr_deg>bestDegree){
				//std::cout<<"break  "<<d_cnt<<":"<<curr_id<<"-"<<curr_deg<<" "<<std::endl;
				break;
			}
			if(curr_deg==0){
				continue;
			}
			isDeleted[curr_id]=true;
			dp_rank[curr_id]=d_cnt++;
			dp_inv.push_back(curr_id);
			if(curr_deg==1){//degree 1
				u=es[curr_id].front().v;
				w=es[curr_id].front().weight;
				neighborsInfo[curr_id].push_back(make_pair(u,w));//store the neighbors information
				es[curr_id].clear();
				v=curr_id;
				for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
					if(it->v==v){
						es[u].erase(it);
						if(!isHFPoint[u]) lqueue.update(u,es[u].size());
						//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
						break;
					}
				}
			}else{//degree >=2
				//deal with each neighbor pair
				vector<pair<int,int> > neighbors;//store the neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
				neighborsInfo[curr_id]=neighbors;//store the neighborinfo
				//u-v menas neighbor pair
				for(int i=0;i<neighbors.size()-1; i++){
					u=neighbors[i].first;
					for(int j=i+1;j<neighbors.size();j++){
						v=neighbors[j].first;
						isLink=false;
						w=neighbors[i].second+neighbors[j].second;
						for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
							//if exist edge(u,v),update the weight(u,v) to minimun
							if(iter_u->v==v){
								if(w<iter_u->weight){
									//update two directions
									iter_u->weight=w;
									for (list<edge>::iterator iter_v = es[v].begin(); iter_v != es[v].end(); iter_v++){
										if(iter_v->v==u){
											iter_v->weight=w;
											break;
										}
									}
								}
								isLink=true;
								break;
							}
						}
						//if no edge(u,v),then insert
						if(!isLink){
							temp.u=u;
							temp.v=v;
							temp.weight=w;
							es[u].push_back(temp);
							temp.u=v;
							temp.v=u;
							es[v].push_back(temp);
						}
					}
				}
				//process and delete each neighbor's edge to curr_id
				for(int i=0;i<neighbors.size();++i){
					u=neighbors[i].first;
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
						if(it->v== v){
							es[u].erase(it);
							if(!isHFPoint[u]) lqueue.update(u,es[u].size());
							//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
							break;
						}
					}
				}
				//delete curr_id and clear all edges
				es[curr_id].clear();
			}
		}//while ends
		time_delete_add=GetCurrentTimeSec()-time_delete_add;
		time_delete+=time_delete_add;
		//collect information and store overlay_graph
		numOfDeletedVertices=0;
		numOfOverlayEdges=0;
		int degree_cnt=0,degree_cnt_new=0,max_degree=0,max_degree_new=0;
		for(NodeID i=0;i<numOfVertices;++i){
			int isize=es_original[i].size();
			int isize_new=es[i].size();
			degree_cnt_new+=isize_new;
			degree_cnt+=isize;
			if(isize>max_degree) max_degree=isize;
			if(isize_new>max_degree_new) max_degree_new=isize_new;
			numOfOverlayEdges+=isize_new;
			if (isize_new == 0) {
				numOfDeletedVertices++;
			}
			//sort the edges by neighbor id
			es[i].sort([](const edge e1, const edge e2) {return e1.v < e2.v; });
		}
		numOfOverlayVertices=numOfOriginalVertices-numOfDeletedVertices;
		double avg_degree=(double)degree_cnt/(double)numOfOriginalVertices;
		double avg_degree_new=(double)degree_cnt_new/(double)numOfOverlayVertices;
		originalToNew.resize(numOfOriginalVertices, -1);
		newToOriginal.reserve(numOfOverlayVertices);
		for(NodeID i=0;i<numOfOriginalVertices;++i){
			if (0 != es[i].size() && !isDeleted[i]) {
				newToOriginal.push_back(i);
				originalToNew[i] = newToOriginal.size() - 1;
			}
		}
		//load graph
		numOfVertices=0;
		numOfEdges=0;
		overlay_wgraph.load_graph(es,originalToNew,numOfOriginalVertices);
		//output for debug
		std::cout<<"numOfDeletedVertices="<<numOfDeletedVertices<<" numOfOverlayVertices="<<numOfOverlayVertices<<" numOfOverlayEdges="<<numOfOverlayEdges<<std::endl;
		std::cout<<"degree_cnt="<<degree_cnt<<" avg_degree="<<avg_degree<<std::endl;
		std::cout<<"degree_cnt_new="<<degree_cnt_new<<" avg_degree_new="<<avg_degree_new<<std::endl;
		//*******************************second stage**************************
		if(thresholdDegree>bestDegree) 
		{
			numOfDeletedVertices=numOfOriginalVertices-numOfHfpoint;
			time_delete_add=GetCurrentTimeSec();
			while(!lqueue.empty()){
				//cout<<"lqueue size="<<lqueue.size()<<std::endl;//to be deleted
				if(d_cnt>=numOfDeletedVertices) break;
				lqueue.extract_min(curr_id,curr_deg);
				if(curr_deg>thresholdDegree){
					//std::cout<<"break  "<<d_cnt<<":"<<curr_id<<"-"<<curr_deg<<" "<<std::endl;
					break;
				}
				if(curr_deg==0){
					continue;
				}
				isDeleted[curr_id]=true;
				dp_rank[curr_id]=d_cnt++;
				dp_inv.push_back(curr_id);
				if(curr_deg==1){//degree 1
					u=es[curr_id].front().v;
					w=es[curr_id].front().weight;
					neighborsInfo[curr_id].push_back(make_pair(u,w));//store the neighbors information
					es[curr_id].clear();
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==v){
							es[u].erase(it);
							if(!isHFPoint[u]) lqueue.update(u,es[u].size());
							//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
							break;
						}
					}
				}else{//degree >=2
					//deal with each neighbor pair
					vector<pair<int,int> > neighbors;//store the neighbor and its distance
					bool isLink=false;//indicate whether exist edge(u,v)
					for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
					neighborsInfo[curr_id]=neighbors;//store the neighborinfo
					//u-v menas neighbor pair
					for(int i=0;i<neighbors.size()-1; i++){
						u=neighbors[i].first;
						for(int j=i+1;j<neighbors.size();j++){
							v=neighbors[j].first;
							isLink=false;
							w=neighbors[i].second+neighbors[j].second;
							for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
								//if exist edge(u,v),update the weight(u,v) to minimun
								if(iter_u->v==v){
									if(w<iter_u->weight){
										//update two directions
										iter_u->weight=w;
										for (list<edge>::iterator iter_v = es[v].begin(); iter_v != es[v].end(); iter_v++){
											if(iter_v->v==u){
												iter_v->weight=w;
												break;
											}
										}
									}
									isLink=true;
									break;
								}
							}
							//if no edge(u,v),then insert
							if(!isLink){
								temp.u=u;
								temp.v=v;
								temp.weight=w;
								es[u].push_back(temp);
								temp.u=v;
								temp.v=u;
								es[v].push_back(temp);
							}
						}
					}
					//process and delete each neighbor's edge to curr_id
					for(int i=0;i<neighbors.size();++i){
						u=neighbors[i].first;
						v=curr_id;
						for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
							if(it->v== v){
								es[u].erase(it);
								if(!isHFPoint[u]) lqueue.update(u,es[u].size());
								//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
								break;
							}
						}
					}
					//delete curr_id and clear all edges
					es[curr_id].clear();
				}
			}//while ends
			time_delete_add=GetCurrentTimeSec()-time_delete_add;
			time_delete+=time_delete_add;
			//collect information and store overlay_graph
			numOfDeletedVertices=0;
			numOfOverlayEdges=0;
			degree_cnt_new=0,max_degree_new=0;
			for(NodeID i=0;i<numOfOriginalVertices;++i){
				int isize_new=es[i].size();
				degree_cnt_new+=isize_new;
				if(isize_new>max_degree_new) max_degree_new=isize_new;
				numOfOverlayEdges+=isize_new;
				if (isize_new == 0) {
					numOfDeletedVertices++;
				}
				//sort the edges by neighbor id
				//es[i].sort([](const edge e1, const edge e2) {return e1.v < e2.v; });
			}
			numOfOverlayVertices=numOfOriginalVertices-numOfDeletedVertices;
			avg_degree_new=(double)degree_cnt_new/(double)numOfOverlayVertices;
			//output for debug
			std::cout<<"numOfDeletedVertices="<<numOfDeletedVertices<<" numOfOverlayVertices="<<numOfOverlayVertices<<" numOfOverlayEdges="<<numOfOverlayEdges<<std::endl;
			std::cout<<"degree_cnt_new="<<degree_cnt_new<<" avg_degree_new="<<avg_degree_new<<std::endl;
		}

	}

	void create_node_hf_overlay_directed(WGraph& overlay_wgraph,vector<int>& newToOriginal,vector<int>& originalToNew,vector<bool>& isDeleted,const vector<bool>& isHFPoint,int& numOfOverlayVertices,int& numOfDeletedVertices,int& numOfOverlayEdges,int thresholdDegree,int bestDegree){
		std::cout<<"******************create_node_hf_overlay_directed start!*****************"<<std::endl;
		//initialize variables
		double time_delete_add;
		time_delete=0;
		numOfHfpoint=0;
		time_delete_add=GetCurrentTimeSec();
		dp_rank.resize(numOfOriginalVertices,-1);
		dp_inv.reserve(numOfOriginalVertices);
		neighborsInfo.resize(numOfOriginalVertices);
		r_neighborsInfo.resize(numOfOriginalVertices);
		benchmark::heap<2, int, int> lqueue(numOfOriginalVertices);
		int u,v,u_w,v_w,w;//temp variables u-start v-end
		int curr_id,curr_deg,deg;
		edge temp;//edge temp variable
		int d_cnt=0;//current num of deleted nodes
		//get the initial degree
		for(int i=0;i<numOfOriginalVertices;++i){
			if(es[i].size()==0&&r_es[i].size()==0){
				cout<<"erroe: isolated vertex "<<i<<std::endl;
				continue;
			}
			if(isHFPoint[i]){
				++numOfHfpoint;
				continue;
			}
			deg=es[i].size()*r_es[i].size();
			if(deg==1&&es[i].front().v==r_es[i].front().u) deg=0;
			//to be deleted
			//std::cout<<i<<":es_size="<<es[i].size()<<" r_es_size="<<r_es[i].size()<<" deg="<<deg<<std::endl;
			lqueue.update(i,deg);
		}
		numOfDeletedVertices=numOfOriginalVertices-numOfHfpoint;
		//**********************the first stage**********************
		while(!lqueue.empty()){
			//cout<<"lqueue size="<<lqueue.size()<<std::endl;//to be deleted
			if(d_cnt>=numOfDeletedVertices) break;
			lqueue.extract_min(curr_id,curr_deg);
			if(curr_deg>bestDegree){
				std::cout<<"bestDegree break  "<<d_cnt<<":"<<curr_id<<"-"<<curr_deg<<" "<<std::endl;
				break;
			}
			isDeleted[curr_id]=true;
			dp_rank[curr_id]=d_cnt++;
			dp_inv.push_back(curr_id);
			if(curr_deg==0){//line
				if(es[curr_id].size()==0&&r_es[curr_id].size()==0){//only left one vertice
					std::cout<<"Isolated vertex!"<<std::endl;
				}else if(es[curr_id].size()==0&&r_es[curr_id].size()!=0){
					for(list<edge>::iterator it=r_es[curr_id].begin();it!=r_es[curr_id].end();++it){
						u=r_es[curr_id].front().u;
						u_w=r_es[curr_id].front().weight;
						r_neighborsInfo[curr_id].push_back(make_pair(u,u_w));
						for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
							if(it->v==curr_id){
								es[u].erase(it);
								if(!isHFPoint[u]){
									deg=es[u].size()*r_es[u].size();
									if(deg==1&&es[u].front().v==r_es[u].front().u) deg=0;
									lqueue.update(u,deg);
								}
								break;
							}
						}
					}
					r_es[curr_id].clear();
				}else if(r_es[curr_id].size()==0&&es[curr_id].size()!=0){//all edges are in frontward
					for(list<edge>::iterator it=es[curr_id].begin();it!=es[curr_id].end();++it){
						v=es[curr_id].front().v;
						v_w=es[curr_id].front().weight;
						neighborsInfo[curr_id].push_back(make_pair(v,v_w));
						for (list<edge>::iterator it = r_es[v].begin(); it != r_es[v].end(); it++) {
							if(it->u==curr_id){
								r_es[v].erase(it);
								if(!isHFPoint[v]){
									deg=es[v].size()*r_es[v].size();
									if(deg==1&&es[v].front().v==r_es[v].front().u) deg=0;
									lqueue.update(v,deg);
								}
								break;
							}
						}
					}
					es[curr_id].clear();
				}else{
					v=es[curr_id].front().v;
					v_w=es[curr_id].front().weight;
					u=r_es[curr_id].front().u;
					u_w=r_es[curr_id].front().weight;
					r_neighborsInfo[curr_id].push_back(make_pair(u,u_w));
					neighborsInfo[curr_id].push_back(make_pair(v,v_w));
					//to be deleted
					//if(u!=v) std::cout<<"Error: deg=0 one vertice!"<<std::endl;
					for (list<edge>::iterator it = r_es[v].begin(); it != r_es[v].end(); it++) {
						if(it->u==curr_id){
							r_es[v].erase(it);
							break;
						}
					}
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==curr_id){
							es[u].erase(it);
							break;
						}
					}
					if(!isHFPoint[u]){
						deg=es[u].size()*r_es[u].size();
						if(deg==1&&es[u].front().v==r_es[u].front().u) deg=0;
						lqueue.update(u,deg);
					}
					es[curr_id].clear();
					r_es[curr_id].clear();
				}
			}else{
				vector<pair<int,int> > neighbors;//store the forward neighbor and its distance
				vector<pair<int,int> > r_neighbors;//store the backward neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++){
					neighbors.push_back(make_pair(it->v,it->weight));
				}
				for (list<edge>::iterator it = r_es[curr_id].begin(); it != r_es[curr_id].end(); it++){
					r_neighbors.push_back(make_pair(it->u,it->weight));
				}
				//store the information
				r_neighborsInfo[curr_id]=r_neighbors;
				neighborsInfo[curr_id]=neighbors;
				//u-v menas neighbor pair
				for(int i=0;i<r_neighbors.size();++i){
					u=r_neighbors[i].first;
					u_w=r_neighbors[i].second;
					for(int j=0;j<neighbors.size();++j){
						v=neighbors[j].first;
						v_w=neighbors[j].second;
						if(u==v) continue;//not need to deal
						isLink=false;
						w=u_w+v_w;
						for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
							//if exist edge(u,v),update the weight(u,v) to minimun
							if(iter_u->v==v){
								if(w<iter_u->weight){
									iter_u->weight=w;//update forward direction
									for (list<edge>::iterator iter_v = r_es[v].begin(); iter_v != r_es[v].end(); iter_v++){
										if(iter_v->u==u){
											iter_v->weight=w;//update backward direction
											break;
										}
									}
									break;
								}
								isLink=true;
								break;
							}
						}
						//if no edge(u,v),then insert
						if(!isLink&&u!=v){
							temp.u=u;
							temp.v=v;
							temp.weight=w;
							es[u].push_back(temp);
							r_es[v].push_back(temp);
						}
					}
				}
				//process neighbor's edges and update dqueue
				for(int i=0;i<r_neighbors.size();++i){
					u=r_neighbors[i].first;
					//u_w=r_neighbors[i].second;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
						if(it->v==curr_id){
							es[u].erase(it);
							if(!isHFPoint[u]){
								deg=es[u].size()*r_es[u].size();
								if(deg==1&&es[u].front().v==r_es[u].front().u) deg=0;
								lqueue.update(u,deg);						
							}
							break;
						}
					}
				}
				for(int j=0;j<neighbors.size();++j){
					v=neighbors[j].first;
					//v_w=neighbors[j].second;
					for (list<edge>::iterator it = r_es[v].begin(); it != r_es[v].end(); it++){
						if(it->u==curr_id){
							r_es[v].erase(it);
							if(!isHFPoint[v]){
								deg=es[v].size()*r_es[v].size();
								if(deg==1&&es[v].front().v==r_es[v].front().u) deg=0;
								lqueue.update(v,deg);
							}
							break;
						}
					}
				}
				//remove curr_id and its adjacent edges from H
				es[curr_id].clear();
				r_es[curr_id].clear();
			}
		}
		time_delete_add=GetCurrentTimeSec()-time_delete_add;
		time_delete+=time_delete_add;
		//collect and output information
		numOfDeletedVertices=0;
		numOfOverlayEdges=0;
		int degree_cnt=0,degree_cnt_new=0,max_degree=0,max_degree_new=0;
		originalToNew.resize(numOfOriginalVertices, -1);
		newToOriginal.reserve(numOfOverlayVertices);
		for (int i = 0; i < numOfOriginalVertices; i++) {
			int isize=es_original[i].size()*r_es_original[i].size();
			degree_cnt+=isize;
			if(isize>max_degree) max_degree=isize;
			if ((0 != es[i].size()||0!=r_es[i].size())&& !isDeleted[i]) {
				newToOriginal.push_back(i);
				originalToNew[i] = newToOriginal.size() - 1;
				int isize_new=es[i].size()*r_es[i].size();
				degree_cnt_new+=isize_new;
				if(isize_new>max_degree_new) max_degree_new=isize_new;
			}else{
				numOfDeletedVertices++;
			}
		}
		numOfOverlayVertices=numOfOriginalVertices-numOfDeletedVertices;
		double avg_degree=(double)degree_cnt/(double)numOfOriginalVertices;
		double avg_degree_new=(double)degree_cnt_new/(double)numOfOverlayVertices;
		std::cout<<"numOfOriginalVertices="<<numOfOriginalVertices<<" numOfDeletedVertices="<<numOfDeletedVertices<<" numOfOverlayVertices="<<numOfOverlayVertices<<std::endl;
		std::cout<<"degree_cnt="<<degree_cnt<<" avg_degree="<<avg_degree<<std::endl;
		std::cout<<"degree_cnt_new="<<degree_cnt_new<<" avg_degree_new="<<avg_degree_new<<std::endl;
		//load overlay_graph
		overlay_wgraph.load_graph(es,originalToNew,numOfOriginalVertices);
		numOfOriginalEdges=numOfEdges;
		std::cout<<"numOfOverlayEdges="<<numOfOverlayEdges<<std::endl;
		//**********************the second stage**********************
		if(thresholdDegree<bestDegree){
			std::cout<<"thresholdDegree<bestDegree!"<<std::endl;
		}
		return;
	}

	/**
	 * @Author: wanjingyi
	 * @description: create node by two stages modified by wanjingyi
	 * @param {*}
	 * @return {*}
	 */ 
	void create_node_hf_all(WGraph& overlay_wgraph,vector<int>& newToOriginal,vector<int>& originalToNew,vector<bool>& isDeleted,const vector<bool>& isHFPoint,int& numOfOverlayVertices,int& numOfDeletedVertices,int& numOfOverlayEdges,int thresholdDegree){
		//initialize variables
		numOfHfpoint=0;
		double time_delete_1=GetCurrentTimeSec();
		//isDeleted.resize(numOfVertices,0);
		dp_rank.resize(numOfVertices,0);
		dp_inv.reserve(numOfVertices);
		neighborsInfo.resize(numOfVertices);
		benchmark::heap<2, int, int> lqueue(numOfVertices);
		int u,v,w;//temp variables u-start v-end
		edge temp;//edge temp variable
		int d_cnt=0;//current num of deleted nodes
		int curr_id,curr_deg;
		bool isOne=false;
		//get the initial degree
		for(int i=0;i<numOfVertices;++i){
			//judge isolated nodes
			int deg=es[i].size();
			if(deg==0){
				cout<<"error: isolated vertex "<<i<<std::endl;
				continue;
			}
			if(isHFPoint[i]){
				++numOfHfpoint;
				continue;
			}
			lqueue.update(i,deg);
		}
		numOfDeletedVertices=numOfVertices-numOfHfpoint;
		//for each vertex to construct node
		//**********************the first stage**********************
		while(!lqueue.empty()){
			//cout<<"lqueue size="<<lqueue.size()<<std::endl;//to be deleted
			if(d_cnt>=numOfDeletedVertices) break;
			lqueue.extract_min(curr_id,curr_deg);
			if(curr_deg>thresholdDegree){
				std::cout<<"break  "<<d_cnt<<":"<<curr_id<<"-"<<curr_deg<<" "<<std::endl;
				break;
			}
			if(curr_deg==0){
				continue;
			}
			isDeleted[curr_id]=true;
			dp_rank[curr_id]=d_cnt++;
			dp_inv.push_back(curr_id);
			if(curr_deg==1){//degree 1
				u=es[curr_id].front().v;
				w=es[curr_id].front().weight;
				neighborsInfo[curr_id].push_back(make_pair(u,w));//store the neighbors information
				es[curr_id].clear();
				v=curr_id;
				for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
					if(it->v==v){
						es[u].erase(it);
						if(!isHFPoint[u]) lqueue.update(u,es[u].size());
						//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
						break;
					}
				}
			}else{//degree >=2
				//deal with each neighbor pair
				vector<pair<int,int> > neighbors;//store the neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
				neighborsInfo[curr_id]=neighbors;//store the neighborinfo
				//u-v menas neighbor pair
				for(int i=0;i<neighbors.size()-1; i++){
					u=neighbors[i].first;
					for(int j=i+1;j<neighbors.size();j++){
						v=neighbors[j].first;
						isLink=false;
						w=neighbors[i].second+neighbors[j].second;
						for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
							//if exist edge(u,v),update the weight(u,v) to minimun
							if(iter_u->v==v){
								if(w<iter_u->weight){
									//update two directions
									iter_u->weight=w;
									for (list<edge>::iterator iter_v = es[v].begin(); iter_v != es[v].end(); iter_v++){
										if(iter_v->v==u){
											iter_v->weight=w;
											break;
										}
									}
								}
								isLink=true;
								break;
							}
						}
						//if no edge(u,v),then insert
						if(!isLink){
							temp.u=u;
							temp.v=v;
							temp.weight=w;
							es[u].push_back(temp);
							temp.u=v;
							temp.v=u;
							es[v].push_back(temp);
						}
					}
				}
				//process and delete each neighbor's edge to curr_id
				for(int i=0;i<neighbors.size();++i){
					u=neighbors[i].first;
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
						if(it->v== v){
							es[u].erase(it);
							if(!isHFPoint[u]) lqueue.update(u,es[u].size());
							//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
							break;
						}
					}
				}
				//delete curr_id and clear all edges
				es[curr_id].clear();
			}
		}//while ends
		time_delete_1=GetCurrentTimeSec()-time_delete_1;
		//collect information and store overlay_graph
		numOfDeletedVertices=0;
		numOfOverlayEdges=0;
		int degree_cnt=0,degree_cnt_new=0,max_degree=0,max_degree_new=0;
		for(NodeID i=0;i<numOfVertices;++i){
			int isize=es_original[i].size();
			int isize_new=es[i].size();
			degree_cnt_new+=isize_new;
			degree_cnt+=isize;
			if(isize>max_degree) max_degree=isize;
			if(isize_new>max_degree_new) max_degree_new=isize_new;
			numOfOverlayEdges+=isize_new;
			if (isize_new == 0) {
				numOfDeletedVertices++;
			}
			//sort the edges by neighbor id
			es[i].sort([](const edge e1, const edge e2) {return e1.v < e2.v; });
		}
		numOfOverlayVertices=numOfVertices-numOfDeletedVertices;
		double avg_degree=(double)degree_cnt/(double)numOfVertices;
		double avg_degree_new=(double)degree_cnt_new/(double)numOfOverlayVertices;
		originalToNew.resize(numOfVertices, -1);
		newToOriginal.reserve(numOfOverlayVertices);
		for(NodeID i=0;i<numOfVertices;++i){
			if (0 != es[i].size() && !isDeleted[i]) {
				newToOriginal.push_back(i);
				originalToNew[i] = newToOriginal.size() - 1;
			}
		}
		int numOfOriginalVertices=numOfVertices;
		int numOfOriginalEdges=numOfEdges;
		numOfVertices=0;
		overlay_wgraph.load_graph(es,originalToNew,numOfOriginalVertices);
		numOfVertices=numOfOriginalVertices;
		numOfEdges=numOfOriginalEdges;
		std::cout<<"numOfHfpoint="<<numOfHfpoint<<std::endl;
		std::cout<<"numOfOverlayVertices="<<numOfOverlayVertices<<std::endl;
		std::cout<<"numOfDeletedVertices="<<numOfDeletedVertices<<std::endl;
		std::cout<<"numOfOverlayEdges="<<numOfOverlayEdges<<std::endl;
		std::cout<<"degree_cnt="<<degree_cnt<<" avg_degree="<<avg_degree<<std::endl;
		std::cout<<"degree_cnt_new="<<degree_cnt_new<<" avg_degree_new="<<avg_degree_new<<std::endl;
		//********************second stage********************
		double time_delete_2=GetCurrentTimeSec();
		lqueue.clear();
		//push degree
		for(int i=0;i<numOfVertices;++i){
			if(isDeleted[i]) continue;
			int deg=es[i].size();
			if(deg==0){
				cout<<"error: isolated vertex "<<i<<endl;
				continue;
			}
			lqueue.update(i,deg);	
		}
		//for each vertex to construct node
		while(!lqueue.empty()){
			lqueue.extract_min(curr_id,curr_deg);
			dp_rank[curr_id]=d_cnt++;
			dp_inv.push_back(curr_id);
			if(curr_deg==0){
				continue;
			}
			else if(curr_deg==1){
					u=es[curr_id].front().v;
					w=es[curr_id].front().weight;
					neighborsInfo[curr_id].push_back(make_pair(u,w));//store the neighbors information
					es[curr_id].clear();
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==v){
							es[u].erase(it);
							lqueue.update(u,es[u].size());
							break;
						}
					}
			}else{//degree>=2
				//deal with each neighbor pair
				vector<pair<int,int> > neighbors;//store the neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
				neighborsInfo[curr_id]=neighbors;
				//u-v menas neighbor pair
				for(int i=0;i<neighbors.size()-1; i++){
					u=neighbors[i].first;
					for(int j=i+1;j<neighbors.size();j++){
						v=neighbors[j].first;
						isLink=false;
						w=neighbors[i].second+neighbors[j].second;
						for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
							//if exist edge(u,v),update the weight(u,v) to minimun
							if(iter_u->v==v){
								if(w<iter_u->weight){
									//update two directions
									iter_u->weight=w;
									for (list<edge>::iterator iter_v = es[v].begin(); iter_v != es[v].end(); iter_v++){
										if(iter_v->v==u){
											iter_v->weight=w;
											break;
										}
									}
								}
								isLink=true;
								break;
							}
						}
						//if no edge(u,v),then insert
						if(!isLink){
							temp.u=u;
							temp.v=v;
							temp.weight=w;
							es[u].push_back(temp);
							temp.u=v;
							temp.v=u;
							es[v].push_back(temp);
						}
					}
				}
				//process neighbor's edges and update dqueue
				for(int i=0;i<neighbors.size();++i){
					u=neighbors[i].first;
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
						if(it->v== v){
							es[u].erase(it);
							lqueue.update(u,es[u].size());//update degree
							break;
						}
					}
				}
				//remove curr_id and its adjacent edges from H
				es[curr_id].clear();
			}
		}//while ends
		time_delete_2=GetCurrentTimeSec()-time_delete_2;
		time_delete=time_delete_1+time_delete_2;
	}

	/**
     * @Author: wanjingyi
     * @description: create node by frequency
     * @param {int} hfRate
     * @return {*}
     */
    void create_node_hf(){
		//std::cout<<"****************************createNode by hf start!***************************"<<std::endl;
		benchmark::my_2heap<2, int,int,int> pqueue(numOfVertices);
		int u,v,w;//temp variables u-start v-end
		int cnt=0;
		int curr_id,curr_deg,curr_freq;
		edge temp;//edge temp variable
		dp_rank.resize(numOfVertices);
		dp_inv.resize(numOfVertices);
		neighborsInfo.resize(numOfVertices);
		//push degree and frequency
		for(int i=0;i<numOfVertices;++i){
			int deg=es[i].size();
			if(deg==0){
				cout<<"error: isolated vertex "<<i<<endl;
				continue;
			}
			pqueue.update(i,deg,query_time[i]);	
		}
		//to be deleted
		// string debug_file_name("../experiment/NewYork/h2h_index/hf/delete.order");
		// ofstream ofs_debug(debug_file_name);
		//if(!ofs_debug.is_open()) cout<<"Cannot open "<<debug_file_name<<endl; 
		//for each vertex to construct node
		while(!pqueue.empty()){
			pqueue.extract_min(curr_id,curr_deg,curr_freq);
			//ofs_debug<<curr_id<<" "<<curr_deg<<" "<<curr_freq<<endl;//to be deleted
			dp_rank[curr_id]=cnt;
			dp_inv[cnt]=curr_id;
			++cnt;
			if(curr_deg==0){
				continue;
			}
			else if(curr_deg==1){
					u=es[curr_id].front().v;
					w=es[curr_id].front().weight;
					neighborsInfo[curr_id].push_back(make_pair(u,w));//store the neighbors information
					es[curr_id].clear();
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==v){
							es[u].erase(it);
							pqueue.update(u,es[u].size(),curr_freq);
							break;
						}
					}
			}else{//degree>=2
				//deal with each neighbor pair
				vector<pair<int,int> > neighbors;//store the neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
				neighborsInfo[curr_id]=neighbors;
				//u-v menas neighbor pair
				for(int i=0;i<neighbors.size()-1; i++){
					u=neighbors[i].first;
					for(int j=i+1;j<neighbors.size();j++){
						v=neighbors[j].first;
						isLink=false;
						w=neighbors[i].second+neighbors[j].second;
						for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
							//if exist edge(u,v),update the weight(u,v) to minimun
							if(iter_u->v==v){
								if(w<iter_u->weight){
									//update two directions
									iter_u->weight=w;
									for (list<edge>::iterator iter_v = es[v].begin(); iter_v != es[v].end(); iter_v++){
										if(iter_v->v==u){
											iter_v->weight=w;
											break;
										}
									}
								}
								isLink=true;
								break;
							}
						}
						//if no edge(u,v),then insert
						if(!isLink){
							temp.u=u;
							temp.v=v;
							temp.weight=w;
							es[u].push_back(temp);
							temp.u=v;
							temp.v=u;
							es[v].push_back(temp);
						}
					}
				}
				//process neighbor's edges and update dqueue
				for(int i=0;i<neighbors.size();++i){
					u=neighbors[i].first;
					v=curr_id;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
						if(it->v== v){
							es[u].erase(it);
							pqueue.update(u,es[u].size(),curr_freq);//update degree
							break;
						}
					}
				}
				//remove curr_id and its adjacent edges from H
				es[curr_id].clear();
			}
		}//while ends
		//std::cout<<"****************************createNode by hf finished!***************************"<<std::endl;
    }

	/**
	 * @description: create tree for all vertices
	 * @param {*}
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
	void create_tree(){
		_dp_tree.initialize(numOfVertices);
		_dp_tree.build_dp_tree(neighborsInfo,dp_rank);
	}

	/**
	 * @description: create tree for all vertices
	 * @param {*}
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
	void create_tree_directed(){
		_dp_tree.initialize(numOfVertices);
		_dp_tree.build_dp_tree_directed(r_neighborsInfo,neighborsInfo,dp_rank);
	}

	/**
	 * @description: create forest for the second solution
	 * @param {int} numOfOverlayVertices
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void create_forest(int numOfNodes,int numOfForestNodes,const vector<unsigned int>& queryTime,bool isSorted){
		//build forest
		time_build_tree=GetCurrentTimeSec();
		_dp_forest.initialize(numOfNodes,numOfForestNodes);
		_dp_forest.build_dp_forest(neighborsInfo,dp_rank,queryTime,isSorted);
		time_build_tree=GetCurrentTimeSec()-time_build_tree;
		std::cout<<"_numOfTrees="<<_dp_forest._numOfTrees<<std::endl;
		//release memory
		vector<vector<pair<int,int> > >().swap(neighborsInfo);
		neighborsInfo.clear();
		_dp_forest.count_and_analyze_forest();
	}

	/**
	 * @description: create forest for the second solution, Add for directed
	 * @param {int} numOfNodes
	 * @param {int} numOfForestNodes
	 * @param {const} vector
	 * @param {bool} isSorted
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void create_forest_directed(int numOfNodes,int numOfForestNodes,const vector<unsigned int>& queryTime,bool isSorted){
		time_build_tree=GetCurrentTimeSec();
		_dp_forest.initialize(numOfNodes,numOfForestNodes);
		_dp_forest.build_dp_forest_directed(r_neighborsInfo,neighborsInfo,dp_rank,queryTime,isSorted);
		time_build_tree=GetCurrentTimeSec()-time_build_tree;
		std::cout<<"_numOfTrees="<<_dp_forest._numOfTrees<<std::endl;
		//release memory
		vector<vector<pair<int,int> > >().swap(neighborsInfo);
		neighborsInfo.clear();
		vector<vector<pair<int,int> > >().swap(r_neighborsInfo);
		r_neighborsInfo.clear();
		//count and analyze the forest
		_dp_forest.count_and_analyze_forest();
		return;
	}

	/**
	 * @description: create forest for the third solution
	 * @param {int} numOfOverlayVertices
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void create_forest(int numOfNodes,int numOfForestNodes){
		//build forest
		time_build_tree=GetCurrentTimeSec();
		_dp_forest.initialize(numOfNodes,numOfForestNodes);
		_dp_forest.build_dp_forest(neighborsInfo,dp_rank);
		time_build_tree=GetCurrentTimeSec()-time_build_tree;
		std::cout<<"_numOfTrees="<<_dp_forest._numOfTrees<<std::endl;
		//release memory
		vector<vector<pair<int,int> > >().swap(neighborsInfo);
		neighborsInfo.clear();
		_dp_forest.count_and_analyze_forest();
	}

	/**
	 * @description: create forest for the third solution,add for directed
	 * @param {int} numOfOverlayVertices
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void create_forest_directed(int numOfNodes,int numOfForestNodes){
		//build forest
		time_build_tree=GetCurrentTimeSec();
		_dp_forest.initialize(numOfNodes,numOfForestNodes);
		_dp_forest.build_dp_forest_directed(r_neighborsInfo,neighborsInfo,dp_rank);
		time_build_tree=GetCurrentTimeSec()-time_build_tree;
		std::cout<<"_numOfTrees="<<_dp_forest._numOfTrees<<std::endl;
		//release memory
		vector<vector<pair<int,int> > >().swap(r_neighborsInfo);
		r_neighborsInfo.clear();
		vector<vector<pair<int,int> > >().swap(neighborsInfo);
		neighborsInfo.clear();
		_dp_forest.count_and_analyze_forest();
	}
	
		/**
		 * @description: 
		 * @param {const} char
		 * @param {int} graph_type
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
		//Add for directed
  		void load_original_graph(const char* graph_file,int graph_type=0){
			if(numOfVertices!=0){
				numOfVertices=0;
				numOfEdges=0;
			}
			vector<edge> Edges;
			ifstream ifs(graph_file);
			//ifs >> numOfVertices >> numOfEdges;
			edge temp;
			int u=0,v=0,w,i=0;
			int self_loop_cnt=0,duplicated_cnt_forward=0,duplicated_cnt_backward=0,insert_edged_cnt_forward=0,insert_edged_cnt_backward=0;
			//read the original graph
			std::cout<<"************************load original graph start!*****************"<<std::endl;
			for (i = 0; ifs >> u >> v >> w;i++) {
				numOfVertices = max(numOfVertices, max(u, v) + 1);
				if(u==v){
					self_loop_cnt++;
					continue;
				}
				temp.u = u;
				temp.v = v;
				temp.weight = w;
				Edges.push_back(temp);
				numOfEdges++;
			}

			if (ifs.bad()) return;
			ifs.close();

			es_original.resize(numOfVertices);
			if(DIRECTED_FLAG) r_es_original.resize(numOfVertices);

			std::cout<<"read lines = "<<i<<std::endl;
			std::cout<<"self_loop_cnt = "<<self_loop_cnt<<std::endl;
			std::cout<<"numOfVertices = "<<numOfVertices<<std::endl;
			std::cout<<"original numOfEdges = "<<numOfEdges<<std::endl;

			adj.resize(numOfVertices);
			adj_weight.resize(numOfVertices);
			assert(adj.size() == numOfVertices && adj_weight.size() == numOfVertices && "No enough memory for adj lists!");
			if (DIRECTED_FLAG == true) {
				r_adj.resize(numOfVertices);
				r_adj_weight.resize(numOfVertices);
				assert(r_adj.size() == numOfVertices && r_adj_weight.size() == numOfVertices  && "No enough memory for adj lists!");
			}

			if(graph_type==1){
				//adjacent list
				for(i=0;i<Edges.size();++i){
					u=Edges[i].u;
					v=Edges[i].v;
					w=Edges[i].weight;
					adj[u].push_back(v);
					adj_weight[u].push_back(w);
					if(DIRECTED_FLAG==true){
						r_adj[v].push_back(u);
						r_adj_weight[v].push_back(w);
					}
				}
			}else{
				//adjacent list
				for(i=0;i<Edges.size();++i){
					u=Edges[i].u;
					v=Edges[i].v;
					w=Edges[i].weight;
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

			//sort the adj and adj_weight
			for(v=0;v<numOfVertices;++v){
				vector<int>& adj_v=adj[v];
				vector<int>& adj_weight_v=adj_weight[v];
				vector<pair<int,int> > adj_vw_v(adj_v.size());
				for (i = 0; i < adj_v.size(); ++i) {
					u = adj_v[i];
					w = adj_weight_v[i];
					adj_vw_v[i] = make_pair(u, w);
				}
				sort(adj_vw_v.begin(), adj_vw_v.end(),cmp_edge_1);//cmp1
				for (i = 0; i < adj_v.size(); ++i) {
					adj_v[i] = adj_vw_v[i].first;
					adj_weight_v[i] = adj_vw_v[i].second;
				}
				adj_vw_v.clear();

				if(DIRECTED_FLAG==true){
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

			//remove the duplicated edges
			int last_v;
			for(u=0;u<numOfVertices;++u){
				last_v=numOfVertices;
				for(i=0;i<adj[u].size();++i){
					if(adj[u][i]==last_v){
						duplicated_cnt_forward++;
						continue;
					}
					//add edges
					temp.u=u;
					temp.v=adj[u][i];
					temp.weight=adj_weight[u][i];
					//es[u].push_back(temp);
					es_original[u].push_back(temp);
					last_v=adj[u][i];
					insert_edged_cnt_forward++;
				}
			}
			std::cout<<"duplicated_cnt_forward = "<<duplicated_cnt_forward<<std::endl;
			std::cout<<"insert_edged_cnt_forward = "<<insert_edged_cnt_forward<<std::endl;

			if(DIRECTED_FLAG==true){
				int last_u;
				for(v=0;v<numOfVertices;++v){
					last_u=numOfVertices;
					for(i=0;i<r_adj[v].size();++i){
						if(r_adj[v][i]==last_u){
							duplicated_cnt_backward++;
							continue;
						}
						//add edges
						temp.u=r_adj[v][i];
						temp.v=v;
						temp.weight=r_adj_weight[v][i];
						//es[u].push_back(temp);
						r_es_original[v].push_back(temp);
						last_u=r_adj[v][i];
						insert_edged_cnt_backward++;
					}
				}
				std::cout<<"duplicated_cnt_backward = "<<duplicated_cnt_backward<<std::endl;
				std::cout<<"insert_edged_cnt_backward = "<<insert_edged_cnt_backward<<std::endl;
			}

			numOfEdges=insert_edged_cnt_forward+insert_edged_cnt_backward;
			std::cout<<"numOfEdges = "<<numOfEdges<<std::endl;
			std::cout<<"numOfVertices = "<<numOfVertices<<std::endl;

			adj.clear();
			r_adj.clear();
			adj.shrink_to_fit();
			r_adj.shrink_to_fit();

			adj_weight.clear();
			r_adj_weight.clear();
			adj_weight.shrink_to_fit();
			r_adj_weight.shrink_to_fit();

			//copy
            es.reserve(numOfVertices);
			es=es_original;
			if(DIRECTED_FLAG==true){
				r_es.reserve(numOfVertices);
				r_es=r_es_original;
			}
			numOfOriginalVertices=numOfVertices;
			numOfOriginalEdges=numOfEdges;
			std::cout<<"************************load original graph finished!*****************"<<std::endl;
			return;
		}

public:
	//***************************tree*****************************
	void outputOrder(char* writeFileName){
		string order_file_name(writeFileName);
		order_file_name.append("_delete.order");
		ofstream ofs(order_file_name);
		if(!ofs.is_open()) cout<<"Cannot open "<<order_file_name<<endl;
		for(size_t i=0;i<numOfVertices;++i) ofs<<dp_inv[i]<<endl;
		ofs.close();
		std::cout<<"Output tree delete order to file finished!"<<std::endl;
	}

	void outputTreeStructure(char* writeFileName){
		_dp_tree.output_tree(writeFileName);
		std::cout<<"Output tree structure to file finished!"<<std::endl;
	}

	void saveLcaVariables(char* writeFileName){
		_dp_tree.save_lca_variables(writeFileName);
		std::cout<<"***********************Save lca variables to binary file succeffully!***********************"<<std::endl;
	}

	void outputAnalysis(char* wirteFileName){
		string tree_analysis_file(wirteFileName);
		tree_analysis_file.append("_tree.analysis");
		ofstream ofs(tree_analysis_file);
		if(!ofs.is_open()) cout<<"Cannot open "<<tree_analysis_file<<endl;
		ofs.setf(ios::fixed);
		ofs.precision(4);
		ofs<<"tree_height="<<_dp_tree._tree_height<<endl;
		ofs<<"max_node_num="<<_dp_tree._max_node_num<<endl;
		ofs<<" max_node_width="<<_dp_tree._max_node_width<<endl;
		ofs << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		ofs << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		ofs << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		ofs << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		ofs << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
		ofs.close();
		std::cout<<"Save analysis result to file succeffully!"<<std::endl;
	}

	//********************************forest*******************************
	void compute_forest(Integrated_Index& integrated_index,double& time_tree_index,double& time_compute_rmq){
		//compute tree
		time_compute=GetCurrentTimeSec();
		_dp_forest.compute_and_combinePLL_forest_top_down(integrated_index);
		time_compute=GetCurrentTimeSec()-time_compute;
		std::cout<<"**********************Compute_forest_top_down finished!**********************"<<std::endl;
		time_rmq=GetCurrentTimeSec();
		_dp_forest.compute_rmq_variables_forest();
		time_rmq=GetCurrentTimeSec()-time_rmq;
		std::cout<<"**********************Compute_rmq_variables_forest finished!**********************"<<std::endl;
		time_index=time_delete+time_build_tree+time_compute;
		time_tree_index=time_index;
		time_compute_rmq=time_rmq;
		//compute label size
		_dp_forest.compute_node_label_size_forest_2(integrated_index.index_,dp_rank);
		//merge to integrated index
		_dp_forest.merge_h2h_index_forest(integrated_index.index_);
		std::cout<<"**********************Merge_h2h_index_forest finished!**********************"<<std::endl;
		//output time
		avg_time_index=time_index/(double) numOfVertices;
		avg_time_delete=time_delete/(double) numOfVertices;
		avg_time_build_tree=time_build_tree/(double) numOfVertices;
		avg_time_compute=time_compute/(double) numOfVertices;
		avg_time_rmq=time_rmq/(double) numOfVertices;
		cout << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		cout << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		cout << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		cout << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		cout << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
	}

	//Add for directed (single thread)
	void compute_forest_directed(Integrated_Index& integrated_index,double& time_tree_index,double& time_compute_rmq){
		//compute tree
		time_compute=GetCurrentTimeSec();
		_dp_forest.compute_and_combinePLL_forest_top_down_directed(integrated_index);
		time_compute=GetCurrentTimeSec()-time_compute;
		std::cout<<"**********************Compute_forest_top_down finished!**********************"<<std::endl;
		time_rmq=GetCurrentTimeSec();
		_dp_forest.compute_rmq_variables_forest();
		time_rmq=GetCurrentTimeSec()-time_rmq;
		std::cout<<"**********************Compute_rmq_variables_forest finished!**********************"<<std::endl;
		time_index=time_delete+time_build_tree+time_compute;
		time_tree_index=time_index;
		time_compute_rmq=time_rmq;
		//compute label size
		_dp_forest.compute_node_label_size_forest_2_directed(integrated_index.bindex_,dp_rank);
		std::cout<<"**********************Compute node label size forest successfully!**********************"<<endl;
		//merge to integrated index
		_dp_forest.merge_h2h_index_forest_directed(integrated_index.bindex_);
		std::cout<<"**********************Merge_h2h_index_forest finished!**********************"<<std::endl;
		//output time
		avg_time_index=time_index/(double) numOfVertices;
		avg_time_delete=time_delete/(double) numOfVertices;
		avg_time_build_tree=time_build_tree/(double) numOfVertices;
		avg_time_compute=time_compute/(double) numOfVertices;
		avg_time_rmq=time_rmq/(double) numOfVertices;
		cout << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		cout << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		cout << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		cout << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		cout << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
	}

	//multiple threads version
	void compute_forest(Integrated_Index& integrated_index,double& time_tree_index,double& time_compute_rmq,int num_threads){
		//compute tree
		time_compute=GetCurrentTimeSec();
		_dp_forest.compute_borders_poses();
		_dp_forest.compute_and_combinePLL_forest_top_down(integrated_index.index_,num_threads);
		time_compute=GetCurrentTimeSec()-time_compute;
		std::cout<<"**********************Compute_forest_top_down finished!**********************"<<std::endl;
		time_rmq=GetCurrentTimeSec();
		_dp_forest.compute_rmq_variables_forest();
		time_rmq=GetCurrentTimeSec()-time_rmq;
		std::cout<<"**********************Compute_rmq_variables_forest finished!**********************"<<std::endl;
		time_index=time_delete+time_build_tree+time_compute;
		time_tree_index=time_index;
		time_compute_rmq=time_rmq;
		//compute label size
		_dp_forest.compute_node_label_size_forest_2(integrated_index.index_,dp_rank);
		//merge to integrated index
		_dp_forest.merge_h2h_index_forest(integrated_index.index_);
		std::cout<<"**********************Merge_h2h_index_forest finished!**********************"<<std::endl;
		//output time
		avg_time_index=time_index/(double) numOfVertices;
		avg_time_delete=time_delete/(double) numOfVertices;
		avg_time_build_tree=time_build_tree/(double) numOfVertices;
		avg_time_compute=time_compute/(double) numOfVertices;
		avg_time_rmq=time_rmq/(double) numOfVertices;
		cout << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		cout << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		cout << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		cout << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		cout << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
	}

	//multiple threads version, add for directed
	void compute_forest_directed(Integrated_Index& integrated_index,double& time_tree_index,double& time_compute_rmq,int num_threads){
		//compute tree
		time_compute=GetCurrentTimeSec();
		_dp_forest.compute_borders_poses();
		_dp_forest.compute_and_combinePLL_forest_top_down_directed(integrated_index.bindex_,num_threads);
		time_compute=GetCurrentTimeSec()-time_compute;
		std::cout<<"**********************Compute_forest_top_down finished!**********************"<<std::endl;
		time_rmq=GetCurrentTimeSec();
		_dp_forest.compute_rmq_variables_forest();
		time_rmq=GetCurrentTimeSec()-time_rmq;
		std::cout<<"**********************Compute_rmq_variables_forest finished!**********************"<<std::endl;
		time_index=time_delete+time_build_tree+time_compute;
		time_tree_index=time_index;
		time_compute_rmq=time_rmq;
		//compute label size
		_dp_forest.compute_node_label_size_forest_2_directed(integrated_index.bindex_,dp_rank);
		//merge to integrated index
		_dp_forest.merge_h2h_index_forest_directed(integrated_index.bindex_);
		std::cout<<"**********************Merge_h2h_index_forest finished!**********************"<<std::endl;
		//output time
		avg_time_index=time_index/(double) numOfVertices;
		avg_time_delete=time_delete/(double) numOfVertices;
		avg_time_build_tree=time_build_tree/(double) numOfVertices;
		avg_time_compute=time_compute/(double) numOfVertices;
		avg_time_rmq=time_rmq/(double) numOfVertices;
		cout << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		cout << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		cout << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		cout << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		cout << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
	}

	/**
	 * @description: single thread push forest for the third solution
	 * @param {*}
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void push_forest(Integrated_Index& integrated_index,double& time_tree_index,double& time_compute_rmq){
		//compute tree
		time_compute=GetCurrentTimeSec();
		_dp_forest.compute_and_pushPLL_forest_top_down(integrated_index);
		time_compute=GetCurrentTimeSec()-time_compute;
		std::cout<<"**********************Compute_forest_top_down finished!**********************"<<std::endl;
		time_rmq=GetCurrentTimeSec();
		_dp_forest.compute_rmq_variables_forest();
		time_rmq=GetCurrentTimeSec()-time_rmq;
		std::cout<<"**********************Compute_rmq_variables_forest finished!**********************"<<std::endl;
		time_index=time_delete+time_build_tree+time_compute;
		time_tree_index=time_index;
		time_compute_rmq=time_rmq;
		//compute label size
		_dp_forest.compute_node_label_size_forest_3(integrated_index.index_,dp_rank);
		//merge to integrated index
		_dp_forest.merge_h2h_index_pushPLL_forest(integrated_index.index_);
		std::cout<<"**********************Merge_h2h_index_forest finished!**********************"<<std::endl;
		//output time
		avg_time_index=time_index/(double) numOfVertices;
		avg_time_delete=time_delete/(double) numOfVertices;
		avg_time_build_tree=time_build_tree/(double) numOfVertices;
		avg_time_compute=time_compute/(double) numOfVertices;
		avg_time_rmq=time_rmq/(double) numOfVertices;
		cout << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		cout << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		cout << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		cout << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		cout << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
	}

	/**
	 * @description: single thread push forest for the third solution, add for directed
	 * @param {*}
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void push_forest_directed(Integrated_Index& integrated_index,double& time_tree_index,double& time_compute_rmq){
		//compute tree
		time_compute=GetCurrentTimeSec();
		_dp_forest.compute_and_pushPLL_forest_top_down_directed(integrated_index);
		time_compute=GetCurrentTimeSec()-time_compute;
		std::cout<<"**********************Compute_forest_top_down finished!**********************"<<std::endl;
		time_rmq=GetCurrentTimeSec();
		_dp_forest.compute_rmq_variables_forest();
		time_rmq=GetCurrentTimeSec()-time_rmq;
		std::cout<<"**********************Compute_rmq_variables_forest finished!**********************"<<std::endl;
		time_index=time_delete+time_build_tree+time_compute;
		time_tree_index=time_index;
		time_compute_rmq=time_rmq;
		//compute label size
		_dp_forest.compute_node_label_size_forest_3_directed(integrated_index.bindex_,dp_rank);
		//merge to integrated index
		_dp_forest.merge_h2h_index_pushPLL_forest_directed(integrated_index.bindex_);
		std::cout<<"**********************Merge_h2h_index_forest finished!**********************"<<std::endl;
		//output time
		avg_time_index=time_index/(double) numOfVertices;
		avg_time_delete=time_delete/(double) numOfVertices;
		avg_time_build_tree=time_build_tree/(double) numOfVertices;
		avg_time_compute=time_compute/(double) numOfVertices;
		avg_time_rmq=time_rmq/(double) numOfVertices;
		cout << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		cout << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		cout << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		cout << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		cout << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
	}

	//multiple threads version
	void push_forest(Integrated_Index& integrated_index,double& time_tree_index,double& time_compute_rmq,int num_threads){
		double time_push=0,time_dis=0;
		//compute tree
		time_compute=GetCurrentTimeSec();
		_dp_forest.compute_borders_poses();
		_dp_forest.compute_and_pushPLL_forest_top_down(integrated_index.index_,integrated_index._max_overlay_hub,num_threads);
		time_compute=GetCurrentTimeSec()-time_compute;
		std::cout<<"**********************Compute_forest_top_down finished!**********************"<<std::endl;
		time_rmq=GetCurrentTimeSec();
		_dp_forest.compute_rmq_variables_forest();
		time_rmq=GetCurrentTimeSec()-time_rmq;
		std::cout<<"**********************Compute_rmq_variables_forest finished!**********************"<<std::endl;
		time_index=time_delete+time_build_tree+time_compute;
		time_tree_index=time_index;
		time_compute_rmq=time_rmq;
		//compute label size
		_dp_forest.compute_node_label_size_forest_3(integrated_index.index_,dp_rank,true);
		//merge to integrated index
		_dp_forest.merge_h2h_index_pushPLL_forest_multiThreads(integrated_index.index_);
		std::cout<<"**********************Merge_h2h_index_forest finished!**********************"<<std::endl;
		//output time
		avg_time_index=time_index/(double) numOfVertices;
		avg_time_delete=time_delete/(double) numOfVertices;
		avg_time_build_tree=time_build_tree/(double) numOfVertices;
		avg_time_compute=time_compute/(double) numOfVertices;
		avg_time_rmq=time_rmq/(double) numOfVertices;
		cout << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		cout << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		cout << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		cout << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		cout << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
	}

	//multiple threads version for the third solution, add for directed
	void push_forest_directed(Integrated_Index& integrated_index,double& time_tree_index,double& time_compute_rmq,int num_threads){
		double time_push=0,time_dis=0;
		//compute tree
		time_compute=GetCurrentTimeSec();
		_dp_forest.compute_borders_poses();
		_dp_forest.compute_and_pushPLL_forest_top_down_directed(integrated_index.bindex_,integrated_index._max_overlay_hub_r,integrated_index._max_overlay_hub,num_threads);
		time_compute=GetCurrentTimeSec()-time_compute;
		std::cout<<"**********************Compute_forest_top_down finished!**********************"<<std::endl;
		time_rmq=GetCurrentTimeSec();
		_dp_forest.compute_rmq_variables_forest();
		time_rmq=GetCurrentTimeSec()-time_rmq;
		std::cout<<"**********************Compute_rmq_variables_forest finished!**********************"<<std::endl;
		time_index=time_delete+time_build_tree+time_compute;
		time_tree_index=time_index;
		time_compute_rmq=time_rmq;
		//compute label size
		_dp_forest.compute_node_label_size_forest_3_directed(integrated_index.bindex_,dp_rank,true);
		//merge to integrated index
		_dp_forest.merge_h2h_index_pushPLL_forest_multiThreads_directed(integrated_index.bindex_);
		std::cout<<"**********************Merge_h2h_index_forest finished!**********************"<<std::endl;
		//output time
		avg_time_index=time_index/(double) numOfVertices;
		avg_time_delete=time_delete/(double) numOfVertices;
		avg_time_build_tree=time_build_tree/(double) numOfVertices;
		avg_time_compute=time_compute/(double) numOfVertices;
		avg_time_rmq=time_rmq/(double) numOfVertices;
		cout << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		cout << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		cout << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		cout << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		cout << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
	}

	void outputForestStructure(char* writeFileName){
		_dp_forest.output_forest(writeFileName);
		std::cout<<"Output forest structure to file finished!"<<std::endl;
	}

	void saveLcaVariables_forest(char* writeFileName){
		_dp_forest.save_lca_variables_forest(writeFileName);
		std::cout<<"Save lca variables to binary file succeffully!"<<std::endl;
	}

	void outputAnalysis_forest(char* wirteFileName){
		string tree_analysis_file(wirteFileName);
		tree_analysis_file.append("_forest.analysis");
		ofstream ofs(tree_analysis_file);
		if(!ofs.is_open()) cout<<"Cannot open "<<tree_analysis_file<<endl;
		ofs.setf(ios::fixed);
		ofs.precision(4);
		ofs<<"_max_allTree_numOfNodes="<<_dp_forest._max_allTree_num<<" _min_allTree_numOfNodes="<<_dp_forest._min_allTree_num<<" _avg_allTree_numOfNodes="<<_dp_forest._avg_allTree_num<<std::endl;
		ofs<<"_max_allTree_height="<<_dp_forest._max_allTree_height<<" _min_allTree_height="<<_dp_forest._min_allTree_height<<" _avg_allTree_height="<<_dp_forest._avg_allTree_height<<std::endl;
		ofs<<"_max_allTree_width="<<_dp_forest._max_allTree_width<<" _min_allTree_width="<<_dp_forest._min_allTree_width<<" _avg_allTree_width="<<_dp_forest._avg_allTree_width<<std::endl;
		ofs<<"_max_allTree_node_num="<<_dp_forest._max_allTree_nodesNum<<std::endl;
		ofs << "indexing_time=" << time_index *1e6 << " avg_indexing_time=" << avg_time_index *1e6 << endl;
		ofs << "deleting_time=" << time_delete *1e6 << " avg_deleting_time=" << avg_time_delete *1e6  << endl;
		ofs << "build_tree_time=" << time_build_tree *1e6 << " avg_build_tree_time=" << avg_time_build_tree *1e6 << endl;
		ofs << "compute_time=" << time_compute *1e6  << " avg_compute_time=" << avg_time_compute *1e6 << endl;
		ofs << "rmq_dfs_time=" << time_rmq *1e6 << " avg_rmq_dfs_time=" << avg_time_rmq *1e6  << endl;
		ofs.close();
		std::cout<<"Save forest analysis result to file succeffully!"<<std::endl;
	}

	void appendAnalysis_forest(char* wirteFileName){
		ofstream ofs(wirteFileName,ios::app|ios::out);
		if(!ofs.is_open()) cout<<"Cannot open "<<wirteFileName<<endl;
		ofs.setf(ios::fixed);
		ofs.precision(4);
		ofs<<_dp_forest._trees.size()<<" ";
		ofs<<_dp_forest._max_allTree_num<<" "<<_dp_forest._avg_allTree_num<<" ";
		ofs<<_dp_forest._max_allTree_width<<" "<<_dp_forest._avg_allTree_width<<" ";
		ofs<<_dp_forest._max_allTree_height<<" "<<_dp_forest._avg_allTree_height<<" ";
		ofs<<_dp_forest._max_borders_num<<" "<<_dp_forest._avg_borders_num<<" ";
		//ofs<<time_index *1e6<<" "<<time_rmq *1e6<<" ";
		ofs.close();
		std::cout<<"Append forest analysis result to file succeffully!"<<std::endl;
	}

};

#endif // ! H2H_CONSTURCTION_H

