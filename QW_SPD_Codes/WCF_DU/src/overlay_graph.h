/*
 * @Description:  overlay graph and its functions
 * @Author: wanjingyi
 * @Date: 2021-01-26 17:33:32
 * @LastEditTime: 2021-10-16 09:46:14
 */
#ifndef OVERLAY_GRAPH_H
#define OVERLAY_GRAPH_H

#include <vector>
#include <fstream>
#include <algorithm>
#include <cstdint> 
#include <iostream>
#include <list>
#include <map>
#include <omp.h>
#include <queue>
#include <unordered_map>
#include "./heap.h"
#include "./paras.h"
#include "./graph.h"
#include "./utils.h"
#include "./labels.h"
#include "./time_util.h"
#include "./ordering.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
using namespace std;
using namespace time_util;

//compare firs by id then by weight
bool cmp_edge_id(const pair<int,int> a, const pair<int,int> b) {
	if(a.first!=b.first){
		return a.first<b.first;
	}else{
		return a.second<b.second;
	}
}

bool cmp_neg(const pair<EdgeWeight,NodeID> a, const pair<EdgeWeight,NodeID> b){
	if(a.first!=b.first){
		return a.first<b.first;
	}else{
		return a.second<b.second;
	}
}

struct overlay_index_t1 {
	vector<NodeID> spt_v;
	vector<EdgeWeight> spt_d;

	NodeID size() {
		return spt_v.size();
	}
};

struct index_t11_p {
	NodeID* spt_v;
	EdgeWeight* spt_d;
}__attribute__((aligned(64)));  // Aligned for cache lines;

/*
 *@description: containing left constraction graph and deleted vertices
 *@author: wanjingyi
 *@date: 2021-01-25
*/
class Overlay_graph{
public:
		struct edge {
			int u;
			int v;
			int weight;
		};
		//********************variables********************
		vector<list<edge> > es;//update overlay_grph
		vector<list<edge> > es_original;//original graph

		vector<list<edge> > r_es;//update overlay_graph for directed
		vector<list<edge> > r_es_original;//original graph	for directed
		
		vector<NodeID> HFPoint;//store the HFpoint
		vector<bool> HFPointFlag;//hfpoint flag by node index
		vector<int> newToOriginal;// get original ids by new ids
		vector<int> originalToNew;//get original ids by new ids
		
		vector<vector<int> > adj;//store the adjacent vertices
		vector< vector<int> > adj_weight;//store the adjacent weight
		
		vector<vector<int> > r_adj;//store the adjacent vertices for directed
		vector< vector<int> > r_adj_weight;//store the adjacent weight for directed
		
		vector<bool> isDeleted;//indicate whether the node has been deleted
		vector<NodeID> deleteToOriginal;//fetch the original idx by sequential index

		long long total_query_time;//total query time for overlay left vertices
		long long total_query_time_hf;//total query time for overlay left vertices
		int max_query_time;//max querytime for overlay left vertices
		int min_query_time;//min querytime for overlay left vertices
		vector<unsigned int> queryTime;//get query time for overlay left vertices

		int numOfHfpoint;//original num of HFpoints
		int numOfOverlayVertices;//num of vertices in overlay graph
		int numOfOriginalVertices;//num of vertices in original graph
		int _numOfEdges;//num of edges in overlay graph
		int _numOfEdges_r;//num of backward edges in overlay graph
		int _numOfEdges_f;//num of foward edges in overlay graph
		int numOfDeletedVertices;//num of deleted vertices

		//********************constructions and deconstructions********************
		Overlay_graph(){
			numOfVertices = 0;
			numOfEdges = 0;
			numOfHfpoint=0;
			numOfOverlayVertices=0;
			numOfOriginalVertices=0;
			_numOfEdges=0;
			HFPointFlag.clear();
			HFPoint.clear();
			adj.clear();
			adj_weight.clear();
			es.clear();
			es_original.clear();
			isDeleted.clear();
			deleteToOriginal.clear();
		}
		~Overlay_graph(){}

		//********************functions****************
		/**
		 * @Author: wanjingyi
		 * @description: read from file indicates whether the node is deleted (true-deleted false-left)
		 * @param {*}
		 * @return {*}
		 */
		void load_is_deleted_and_map(char* isDeletedFIleName){
			numOfDeletedVertices=0;
			ifstream in(isDeletedFIleName);
			if(!in.is_open()) cerr<<isDeletedFIleName<<" cannot be opened!"<<std::endl;
			in>>numOfOriginalVertices>>numOfOverlayVertices;
			numOfDeletedVertices=numOfOriginalVertices-numOfOverlayVertices;
			cout<<"numOfDeletedVertices="<<numOfDeletedVertices<<" numOfOriginalVertices="<<numOfOriginalVertices<<" numOfOverlayVerticese="<<numOfOverlayVertices<<std::endl;
			isDeleted.resize(numOfOriginalVertices,0);
			deleteToOriginal.resize(numOfDeletedVertices);
			int i,u,flag,d_i=0;
			//build the map relationship 
			for (i = 0; in >> u >> flag;i++) {
				isDeleted[u]=flag;
				if(flag){
					deleteToOriginal[d_i++]=u;
				}
			}
			if(i!=numOfOriginalVertices) cerr<<"error:i!=numOfOriginalVertices!"<<std::endl;	
			in.close();
		}

		/**
		 * @Author: wanjingyi
		 * @description: read from file indicates whether the node is deleted (true-deleted false-left)
		 * @param {*}
		 * @return {*}
		 */
		void load_is_deleted(char* isDeletedFIleName){
			numOfDeletedVertices=0;
			ifstream in(isDeletedFIleName);
			if(!in.is_open()) cerr<<isDeletedFIleName<<" cannot be opened!"<<std::endl;
			in>>numOfOriginalVertices>>numOfOverlayVertices;
			numOfVertices=numOfOriginalVertices;
			isDeleted.resize(numOfOriginalVertices,0);
			int i,u,flag;
			//build the map relationship 
			for (i = 0; in >> u >> flag;i++) {
				isDeleted[u]=flag;
				if(flag) ++numOfDeletedVertices;
			}
			if(i!=numOfOriginalVertices) cerr<<"error:i!=numOfOriginalVertices!"<<std::endl;	
			cout<<"numOfVertices="<<numOfVertices<<" numOfDeletedVertices="<<numOfDeletedVertices<<" numOfOriginalVertices="<<numOfOriginalVertices<<" numOfOverlayVerticese="<<numOfOverlayVertices<<std::endl;
			in.close();
		}

		/**
		 * @Author: wanjingyi
		 * @description: load new node ids and original ids
		 * @param {*}
		 * @return {*}
		 */
		void load_id_map(char* newToOriginalFileName){
			ifstream in(newToOriginalFileName);
			if(!in.is_open()) cerr<<newToOriginalFileName<<" cannot be opened!"<<std::endl;
			in>>numOfOriginalVertices>>numOfOverlayVertices;
			numOfVertices=numOfOriginalVertices;
			cout<<"numOfVertices="<<numOfVertices<<" numOfOriginalVertices="<<numOfOriginalVertices<<" numOfOverlayVertices="<<numOfOverlayVertices<<std::endl;
			originalToNew.resize(numOfOriginalVertices, -1);
			newToOriginal.resize(numOfOverlayVertices);
			int i,u,v;
			//build the map relationship 
			for (i = 0; in >> u >> v;i++) {
				newToOriginal[u]=v;
				originalToNew[v]=u;
			}
			if(i!=numOfOverlayVertices) cerr<<"error:i!=numOfOverlayVertices!"<<std::endl;
			in.close();
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
				sort(adj_vw_v.begin(), adj_vw_v.end(),cmp_edge_id);//cmp1
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
			std::cout<<"************************load original graph finished!*****************"<<std::endl;
			return;
		}

		void bfs_sp_search(int source,vector<int>& dist_table,int hop_limit,benchmark::heap<2, int, int>& pqueue,vector<int>& hops,vector<bool>& vis){
			//clear tmp list
			pqueue.clear();
			for(int i=0;i<numOfVertices;++i){
				hops[i]=0;
				vis[i]=false;
			}
			int u,v,v_w,u_w,v_d;
			pqueue.update(source,0);
			dist_table[source]=0;
			while(!pqueue.empty()){
				pqueue.extract_min(u,u_w);
				vis[u]=true;
				if(hops[u]>hop_limit) break;
				for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
					v=it->v;
					if(vis[v]) continue;
					v_w=it->weight;
					v_d=u_w+v_w;
					if(v_d<dist_table[v]){
						dist_table[v]=v_d;
						hops[v]=hops[u]+1;
						pqueue.update(v,v_d);
					}
				}
			}
		}
};

/*
 *@description: contruct labelings on top of the node deletion
 *@author: wanjingyi
 *@date: 2021-01-25
*/
class Construction_Overlay:public Overlay_graph{
	public:
		//********************variables********************
            vector<NodeID> original_rank;
            vector<NodeID> original_inv;
            vector<NodeID> overlay_rank;
            vector<NodeID> overlay_inv;

			vector<unsigned int> betweenessValues;//get bet values for overlay left vertices
			vector<unsigned int> coverageValues;//get cov values for overlay left vertices
			vector<unsigned int> depthValues;//get dep values for overlay left vertices
		//********************constructions and deconstructions********************

		//********************functions********************
		/**
		 * @description: load left hfpoint setted before node deletion
		 * @Author: wanjingyi 
		 * @Date: 2021-02-01 11:49:02
		 * @param {*}
		 * @return {*}
		 */
		void load_hfpoint_setting(char* HFPoint_file,int hfRate){
			cout<<"**************load_hfpoint_setting start!****************"<<std::endl;
			if(hfRate==0){
				numOfHfpoint=numOfOriginalVertices;
			}else{
				numOfHfpoint= 0;//first line is the number of HFpoints
				numOfHfpoint= static_cast<int> ( (double)numOfOriginalVertices*hfRate/(double)HF_DIVIDION);
			}
			if(numOfHfpoint<=0) cout<<"error:numOfHfpoint<=0"<<std::endl;
			HFPointFlag.resize(numOfOverlayVertices,0);
			queryTime.resize(numOfOverlayVertices,0);
			min_query_time=INF_WEIGHT;
			max_query_time=0;
			total_query_time=0;
			total_query_time_hf=0;
			cout<<"hfpoint resize successfully!"<<std::endl;
			load_hfpoint_overlay(HFPoint_file);
			cout<<"**************load_hfpoint_setting finished!****************"<<std::endl;
		}

		/*
		*@description: load hfpoint
		*@author: wanjingyi
		*@date: 2021-01-16
		*/
		void load_hfpoint_overlay(char* load_filename){
			std::cout<<"initial numOfHfpoint = "<<numOfHfpoint<<std::endl;
			std::ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int id;int t;
			int cnt_hf=0;
			char line[24];
			//read each line representing HFpoint to vector 
			while (in.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>id>>t;		
				total_query_time+=t;
				
				if(cnt_hf>=numOfHfpoint) continue;
				total_query_time_hf+=t;
				if(t>max_query_time) max_query_time=t;
				if(t<min_query_time) min_query_time=t;
				NodeID v=originalToNew[id];
				queryTime[v]=t;
				HFPoint.push_back(v);
				HFPointFlag[v]=true;
				++cnt_hf;
			}
			if(cnt_hf<numOfHfpoint){
				numOfHfpoint=cnt_hf;
				std::cout<<"real numOfHfpoint = "<<numOfHfpoint<<std::endl;
			}
			std::cout<<"max_query_time = "<<max_query_time<<std::endl;
			std::cout<<"min_query_time = "<<min_query_time<<std::endl;
			std::cout<<"toal_query_time = "<<total_query_time<<std::endl;
			std::cout<<"toal_query_time_hf = "<<total_query_time_hf<<std::endl;
		}	

		/**
	 * @description: extract order values for left vertices
	 * @Author: wanjingyi
	 * @Date: 2021-02-02 16:08:52
	 * @param {*}
	 * @return {*}
	 */
		void extract_order_values(char* values_filename,vector<unsigned int>& order_values){
			order_values.resize(numOfOverlayVertices,0);
			ifstream in(values_filename);
			if(!in.is_open()) cout<<"Cannot open "<<values_filename<<std::endl;
			NodeID i,value,u,v,j=0;
			for(i=0;in>>v>>value;i++){
				u=originalToNew[v];
				if(u!=-1){
					order_values[u]=value;
					++j;
				}
			}
			if(i!=numOfOriginalVertices) cout<<"i!=numOfOriginalVertices"<<std::endl;
			if(j!=numOfOverlayVertices) cout<<"j!=numOfOverlayVertices"<<std::endl;
			in.close();
		}

		/**
	 * @description: extract order from already-computed order
	 * @Author: wanjingyi
	 * @Date: 2021-02-01 09:59:59
	 * @param {*}
	 * @return {*}
	 */
		void extract_order_from_original_order(char* order_filename){
			original_inv.resize(numOfOriginalVertices);
			original_rank.resize(numOfOriginalVertices);
			overlay_inv.resize(numOfOverlayVertices);
			overlay_rank.resize(numOfOverlayVertices);
			//read original order
			get_order_from_file(order_filename,original_rank,original_inv,numOfOriginalVertices);
			//extract order
			NodeID rank=0;
			for(int r=0;r<numOfOriginalVertices;++r){
				int v=originalToNew[original_inv[r]];
				if(v==-1) continue;
				overlay_inv[rank]=v;
				overlay_rank[v]=rank;
				rank++;
			}
			if(rank!=numOfOverlayVertices) cout<<"error:rank!=numOfOverlayVertices!"<<std::endl;
			cout<<"rank = "<<rank<<std::endl;
			
		}

		/*
		 *@description: 
		 *@author: wanjingyi
		 *@date: 2021-01-25
		*/
		void compareOrder( char* originalOrderFileName, char* output_dir_name){
			ifstream ifs2(originalOrderFileName);
			if(!ifs2.is_open()) cerr<<originalOrderFileName<<" cannot be opend!"<<std::endl;
			string output_order_filename(output_dir_name);
			output_order_filename.append("original_overlay.order");
			ofstream out(output_order_filename);
			original_rank.resize(numOfOriginalVertices);
			original_inv.resize(numOfOriginalVertices);
			get_order_from_file(originalOrderFileName, original_rank,original_inv,numOfOriginalVertices);
			//extract the left order from original order
			int r,i=0;
			for(r=0;r<numOfOriginalVertices;++r){
				int v=originalToNew[original_inv[r]];
				if(v!=-1){
					out<<v<<std::endl;
					++i;
				}
			}
			out.close();
			if(i!=numOfOverlayVertices) cerr<<"error:i!=numOfOverlayVertices!"<<std::endl;
		}

		void get_betweenness(){

		}

};

/*
 *@description: contruct labelings for the deleted vertices
 *@author: wanjingyi
 *@date: 2021-02-03
*/
class Construction_Deleted:public Overlay_graph{
public:
	//*********variables**********
	//HFLabel labels;
	vector<overlay_index_t1> index_overlay;
	vector<NodeID> overlay_rank;//fetch rank by original index
	vector<NodeID> overlay_inv;//fetch original index by rank
	double _labeling_time,_ave_labeling_time;

	//*************constructions and deconstructions***************
	Construction_Deleted(WGraph& original_wgraph,char* isDeletedFileName,char* newToOriginalFileName,char*overlayLabelsFileName,char* queryFreqFileName,char* outputDirName,char* deletedLabelsFileName,char* labelSizeFileName,int hfRate,int num_threads){
		load_is_deleted(isDeletedFileName);
		std::cout<<"load is deleted successfully!"<<std::endl;
		index_overlay.resize(numOfVertices);//numOfvertices=numOfOriginalVertices
		construct_labels(original_wgraph,isDeletedFileName,num_threads);
		load_id_map(newToOriginalFileName);
		std::cout<<"load overlay map successfully!"<<std::endl;
		load_overlay_labels(overlayLabelsFileName);
		std::cout<<"load overlay labels successfully!"<<std::endl;
		save_and_analyze_index(outputDirName,deletedLabelsFileName,labelSizeFileName,queryFreqFileName,hfRate);
		std::cout<<"save and analysis index successfully!"<<std::endl;
	}

	//experiment model
	Construction_Deleted(WGraph& original_wgraph,char* isDeletedFileName,char* newToOriginalFileName,char*overlayLabelsFileName,char* queryFreqFileName,char* outputDirName,char* deletedLabelsFileName,char* labelSizeFileName,int hfRate,int num_threads,bool isExperiment){
		load_is_deleted(isDeletedFileName);
		std::cout<<"load is deleted successfully!"<<std::endl;
		index_overlay.resize(numOfVertices);
		construct_labels(original_wgraph,isDeletedFileName,num_threads);
		load_id_map(newToOriginalFileName);
		std::cout<<"load overlay map successfully!"<<std::endl;
		load_overlay_labels(overlayLabelsFileName);
		std::cout<<"load overlay labels successfully!"<<std::endl;
		save_and_analyze_index(outputDirName,deletedLabelsFileName,labelSizeFileName,queryFreqFileName,hfRate,isExperiment);
		std::cout<<"save and analysis index successfully!"<<std::endl;
	}

	//Construction_Deleted(){}
	~Construction_Deleted(){

	}


protected:
	//*************functions****************
	void construct_labels(WGraph& original_wgraph,char* isDeletedFileName,int num_threads){
		//construction
		_labeling_time=GetCurrentTimeSec();
		bfs_local_search(original_wgraph,num_threads);
		_labeling_time = GetCurrentTimeSec() - _labeling_time;
		_ave_labeling_time=_labeling_time/(double) numOfOverlayVertices;
		cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << std::endl;
		cout<<"average labeling time:"<<_ave_labeling_time*1e6 <<  " microseconds" << std::endl;
	}

	void append_expriment_result(char* write_filename){
		double performance_result=0,standard_performance_result=0;//total performance function 
		long long total_sum_size=0,real_total_sum_size=0,delete_sum_size=0,real_delete_sum_size=0,hf_sum_size=0,overlay_sum_size=0;
		double total_ave_size=0,real_total_ave_size=0,delete_ave_size=0,real_delete_ave_size=0,hf_ave_size=0,overlay_ave_size=0;
		NodeID isize;
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			isize = index_overlay[v].size()-1;
			total_sum_size+=isize;
			if(isDeleted[v]) delete_sum_size+=isize;
			else overlay_sum_size+=isize;
			if(HFPointFlag[v]) hf_sum_size+=isize;
			//compute the query performance function
			double ratio=(double)queryTime[v]/((double)(max_query_time*DIVISION_FACTOR));
			performance_result+=ratio*(double)isize;	
		}
		standard_performance_result=performance_result*((double)max_query_time);
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			if(isDeleted[v]){
				isize=0;
				int i;
				for(i=1;i<index_overlay[v].spt_v[0];++i) isize+=index_overlay[index_overlay[v].spt_v[i]].size()-1;
				isize+=index_overlay[v].spt_d[0]-index_overlay[v].spt_v[0];
				real_delete_sum_size+=isize;
			}else{
				isize = index_overlay[v].size()-1;
			}
			real_total_sum_size+=isize;
		}
		//print
		total_ave_size= (double) total_sum_size/(double) numOfVertices;
		real_total_ave_size= (double) real_total_sum_size/(double) numOfVertices;
		hf_ave_size= (double) hf_sum_size/(double) numOfHfpoint;
		delete_ave_size=(double) delete_sum_size/(double) numOfDeletedVertices;
		real_delete_ave_size=(double) real_delete_sum_size/(double) numOfDeletedVertices;
		overlay_ave_size=(double) overlay_sum_size/(double) numOfOverlayVertices; 
		std::cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<" real_total_sum_size = "<<real_total_sum_size<<" real_total_ave_size = "<<real_total_ave_size<<std::endl;
		std::cout<<"numOfHFpoint = "<<numOfHfpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<std::endl;
		std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" overlay_sum_size = "<<overlay_sum_size<<" overlay_ave_size = "<<overlay_ave_size<<std::endl;
		std::cout<<"numOfDeletedVertices = "<<numOfDeletedVertices<<" delete_sum_size = "<<delete_sum_size<<" delete_ave_size = "<<delete_ave_size<<" real_delete_sum_size = "<<real_delete_sum_size<<" real_delete_ave_size = "<<real_delete_ave_size<<std::endl;
		std::cout<<"nomalization performance_result = "<<performance_result<<std::endl;
		std::cout<<"standard performance_result = "<<performance_result<<std::endl;
		//append to file
		ofstream ofs(write_filename,ios::app|ios::out);//append way
		if(!ofs.is_open()) cout<<"Cannot open "<<write_filename<<endl;
		ofs.setf(ios::fixed);
		ofs.precision(4);
		//ofs out
		ofs<<_labeling_time* 1e6<<" ";
		//ofs<<standard_performance_result<<" ";
		//ofs<<hf_sum_size<<" "<<hf_ave_size<<" ";
		//ofs<<overlay_sum_size<<" "<<overlay_ave_size<<" ";
		ofs<<delete_sum_size<<" "<<delete_ave_size<<" "<<real_delete_sum_size<<" "<<real_delete_ave_size<<" ";
		ofs<<total_sum_size<<" "<<total_ave_size<<" "<<real_total_sum_size<<" "<<real_total_ave_size<<" ";
		ofs.close();
	}

	void save_and_analyze_index(char* outputDirName,char* deletedLabelsFileName,char* labelSizeFileName,char* queryFreqFileName,int hfRate){
		//save all labels to binaty file
		save_all_labels(deletedLabelsFileName);
		//save label size
		save_all_label_size(labelSizeFileName);
		//write labels_lias and size for debug
		write_labels(outputDirName);
		//load hfpoint for analysis
		load_hfpoint_setting(queryFreqFileName,hfRate);
		//save_analysis
		save_analysis(outputDirName);
		//save_analysis_overlay
		save_analysis_hierarchy_construction(outputDirName);
	}

	void save_and_analyze_index(char* outputDirName,char* deletedLabelsFileName,char* labelSizeFileName,char* queryFreqFileName,int hfRate,bool isExperiment){
		//save all labels to binaty file
		save_all_labels(deletedLabelsFileName);
		//save label size
		save_all_label_size(labelSizeFileName,isExperiment);
		//load hfpoint for analysis
		load_hfpoint_setting(queryFreqFileName,hfRate);
		//save_analysis
		append_expriment_result(outputDirName);
	}

	//*************functions****************
	void save_all_label_size(char* label_size_file){
		ofstream ofs(label_size_file); 
		if(!ofs.is_open()) {cerr<<"Cannot open "<<label_size_file<<endl;}
		for (int v = 0; v < numOfVertices; ++v) {
			ofs << index_overlay[v].size()-1 << endl;
		}
		ofs.close();
		string label_size_file_name(label_size_file);
		label_size_file_name.append("_real");
		ofstream ofs_real(label_size_file_name); 
		if(!ofs_real.is_open()) {cerr<<"Cannot open "<<label_size_file_name<<endl;}
		NodeID isize;
		for (int v = 0; v < numOfVertices; ++v) {
			if(isDeleted[v]){
				isize=0;
				int i;
				for(i=1;i<index_overlay[v].spt_v[0];++i) isize+=index_overlay[index_overlay[v].spt_v[i]].size()-1;
				isize+=index_overlay[v].spt_d[0]-index_overlay[v].spt_v[0];
			}else{
				isize = index_overlay[v].size()-1;
			}
			ofs_real<<isize<<endl;
		}
		ofs_real.close();
	}

	void save_all_label_size(char* label_size_file,bool isExperiment){
		ofstream ofs(label_size_file); 
		if(!ofs.is_open()) {cerr<<"Cannot open "<<label_size_file<<endl;}
		NodeID isize;
		for (int v = 0; v < numOfVertices; ++v) {
			if(isDeleted[v]){
				isize=0;
				int i;
				for(i=1;i<index_overlay[v].spt_v[0];++i) isize+=index_overlay[index_overlay[v].spt_v[i]].size()-1;
				isize+=index_overlay[v].spt_d[0]-index_overlay[v].spt_v[0];
			}else{
				isize = index_overlay[v].size()-1;
			}
			ofs<<isize<<endl;
		}
		ofs.close();
	}

	void save_analysis(char* write_filename){
		//*************output  analysis size*******************
		cout<<"**************analysis****************"<<std::endl;
		long long total_sum_size=0,hf_sum_size=0,delete_sum_size=0,overlay_sum_size=0;
		double total_ave_size=0,hf_ave_size=0,delete_ave_size=0,overlay_ave_size=0;
		double performance_result=0;//total performance function 
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			NodeID isize = index_overlay[v].size()-1;
			total_sum_size+=isize;
			if(isDeleted[v]) delete_sum_size+=isize;
			else overlay_sum_size+=isize;
			if(HFPointFlag[v]) hf_sum_size+=isize;
			//compute the query performance function
			double ratio=(double)queryTime[v]/((double)(max_query_time*DIVISION_FACTOR));
			performance_result+=ratio*(double)isize;
		}
		total_ave_size= (double) total_sum_size/(double) numOfVertices;
		hf_ave_size= (double) hf_sum_size/(double) numOfHfpoint;
		delete_ave_size=(double) delete_sum_size/(double) numOfDeletedVertices;
		overlay_ave_size=(double) overlay_sum_size/(double) numOfOverlayVertices; 
		cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<std::endl;
		cout<<"numOfHFpoint = "<<numOfHfpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<std::endl;
		cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" overlay_sum_size = "<<overlay_sum_size<<" overlay_ave_size = "<<overlay_ave_size<<std::endl;
		cout<<"numOfDeletedVertices = "<<numOfDeletedVertices<<" delete_sum_size = "<<delete_sum_size<<" delete_ave_size = "<<delete_ave_size<<std::endl;
		string output_asizeFileName(write_filename);
		output_asizeFileName.append(".analysis");
		ofstream ofs(output_asizeFileName);
		if(!ofs.is_open()) {cerr<<"Cannot open "<<output_asizeFileName<<std::endl;}
		ofs<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<std::endl;
		ofs<<"numOfHFpoint = "<<numOfHfpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<std::endl;
		ofs<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" overlay_sum_size = "<<overlay_sum_size<<" overlay_ave_size = "<<overlay_ave_size<<std::endl;
		ofs<<"numOfDeletedVertices = "<<numOfDeletedVertices<<" delete_sum_size = "<<delete_sum_size<<" delete_ave_size = "<<delete_ave_size<<std::endl;
		//write performance function to file
		std::cout<<"nomalization performance_result = "<<performance_result<<std::endl;
		ofs<<"nomalization performance_result = "<<performance_result<<std::endl;
		performance_result=performance_result*((double)max_query_time);
		std::cout<<"standard performance_result = "<<performance_result<<std::endl;
		ofs<<"standard performance_result = "<<performance_result<<std::endl;
		ofs.close();
		cout<<"**************analysis****************"<<std::endl;
	}

	void save_analysis_hierarchy_construction(char* write_filename){
		//*************output  overlay analysis size*******************
		cout<<"**************analysis overlay****************"<<std::endl;
		long long total_sum_size=0,hf_sum_size=0,delete_sum_size=0,overlay_sum_size=0;
		double total_ave_size=0,hf_ave_size=0,delete_ave_size=0,overlay_ave_size=0;
		double performance_result=0;//total performance function 
		NodeID isize;
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			if(isDeleted[v]){
				isize=0;
				int i;
				for(i=1;i<index_overlay[v].spt_v[0];++i) isize+=index_overlay[index_overlay[v].spt_v[i]].size()-1;
				isize+=index_overlay[v].spt_d[0]-index_overlay[v].spt_v[0];
				delete_sum_size+=isize;
			}else{
				isize = index_overlay[v].size()-1;
				overlay_sum_size+=isize;
			}
			total_sum_size+=isize;
			if(HFPointFlag[v]) hf_sum_size+=isize;
			//compute the query performance function
			double ratio=(double)queryTime[v]/((double)(max_query_time*DIVISION_FACTOR));
			performance_result+=ratio*(double)isize;
		}
		total_ave_size= (double) total_sum_size/(double) numOfVertices;
		hf_ave_size= (double) hf_sum_size/(double) numOfHfpoint;
		delete_ave_size=(double) delete_sum_size/(double) numOfDeletedVertices;
		overlay_ave_size=(double) overlay_sum_size/(double) numOfOverlayVertices; 
		cout<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<std::endl;
		cout<<"numOfHFpoint = "<<numOfHfpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<std::endl;
		cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" overlay_sum_size = "<<overlay_sum_size<<" overlay_ave_size = "<<overlay_ave_size<<std::endl;
		cout<<"numOfDeletedVertices = "<<numOfDeletedVertices<<" delete_sum_size = "<<delete_sum_size<<" delete_ave_size = "<<delete_ave_size<<std::endl;
		string output_asizeFileName_overlay(write_filename);
		output_asizeFileName_overlay.append("_overlay.analysis");
		ofstream ofs_overlay(output_asizeFileName_overlay.c_str());
		if(!ofs_overlay.is_open()) {cerr<<"Cannot open "<<output_asizeFileName_overlay<<std::endl;}
		ofs_overlay<<"numOfVertices = "<<numOfVertices<<" total_sum_size = "<<total_sum_size<<" total_ave_size = "<<total_ave_size<<std::endl;
		ofs_overlay<<"numOfHFpoint = "<<numOfHfpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<std::endl;
		ofs_overlay<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" overlay_sum_size = "<<overlay_sum_size<<" overlay_ave_size = "<<overlay_ave_size<<std::endl;
		ofs_overlay<<"numOfDeletedVertices = "<<numOfDeletedVertices<<" delete_sum_size = "<<delete_sum_size<<" delete_ave_size = "<<delete_ave_size<<std::endl;
		//write performance function to file
		std::cout<<"nomalization performance_result = "<<performance_result<<std::endl;
		ofs_overlay<<"nomalization performance_result = "<<performance_result<<std::endl;
		performance_result=performance_result*((double)max_query_time);
		std::cout<<"standard performance_result = "<<performance_result<<std::endl;
		ofs_overlay<<"standard performance_result = "<<performance_result<<std::endl;
		ofs_overlay.close();
		cout<<"**************analysis overlay****************"<<std::endl;
	}

	void save_all_labels(char* save_filename){
		ofstream ofs(save_filename, ios::binary | ios::out);
		if(!ofs.is_open()) cout<<"Cnnot open"<<save_filename<<std::endl;
		ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
		for (NodeID v = 0; v < numOfVertices; ++v) {
			NodeID isize = index_overlay[v].size();
			ofs.write((const char*)&isize, sizeof(isize));
			for (NodeID i = 0; i < index_overlay[v].size(); ++i) {
				ofs.write((const char*)&index_overlay[v].spt_v[i], sizeof(index_overlay[v].spt_v[i]));
				ofs.write((const char*)&index_overlay[v].spt_d[i], sizeof(index_overlay[v].spt_d[i]));
			}
		}
		ofs.close();
	}

	void write_labels(char* write_filename){
		string write_filename_prefix(write_filename);
		string write_filename1=write_filename_prefix.append(".list");
		ofstream ofs(write_filename1.c_str());		
		string labelSizefile(write_filename);
		labelSizefile.append(".size");
		ofstream ofs_size(labelSizefile);
		for (NodeID v = 0; v < numOfVertices; ++v) 
		{
			NodeID isize = index_overlay[v].size();
			ofs_size << isize-1 << std::endl;
			ofs <<isize<<" "<<v;
			for (NodeID i = 0; i < index_overlay[v].size(); ++i) {
				ofs<<" "<<'('<<index_overlay[v].spt_v[i]<<","<<index_overlay[v].spt_d[i]<<")";
			}
			ofs<<std::endl;
		}
		ofs.close();
		ofs_size.close();
	}

	/**
	 * @description: load left hfpoint setted before node deletion
	 * @Author: wanjingyi 
	 * @Date: 2021-02-01 11:49:02
	 * @param {*}
	 * @return {*}
	 */
	void load_hfpoint_setting(char* HFPoint_file,int hfRate){
		cout<<"**************load_hfpoint_setting start!****************"<<std::endl;
		if(hfRate==0){
			numOfHfpoint=numOfOriginalVertices;
		}else{
			numOfHfpoint= 0;//first line is the number of HFpoints
			numOfHfpoint= static_cast<int> ( (double)numOfOriginalVertices*hfRate/(double)HF_DIVIDION);
		}
		if(numOfHfpoint<=0) cout<<"error:numOfHfpoint<=0"<<std::endl;
		HFPointFlag.resize(numOfOriginalVertices,0);
		queryTime.resize(numOfOriginalVertices,0);
		max_query_time=0;
		total_query_time=0;
		cout<<"hfpoint resize successfully!"<<std::endl;
		load_hfpoint_overlay(HFPoint_file);
		cout<<"**************load_hfpoint_setting finished!****************"<<std::endl;
	}

	/*
	*@description: load hfpoint
	*@author: wanjingyi
	*@date: 2021-01-16
	*/
	void load_hfpoint_overlay(char* load_filename){
		std::cout<<"initial numOfHfpoint = "<<numOfHfpoint<<std::endl;
		std::ifstream in(load_filename);//input HFPoint file to ifstream
		if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
		int id;int t;
		int cnt_hf=0;
		char line[24];
		//read each line representing HFpoint to vector 
		while (in.getline(line,sizeof(line))){
			if(cnt_hf>=numOfHfpoint) break;
			std::stringstream hp(line);
			hp>>id>>t;		
			//cnt
			queryTime[id]=t;
			if(t>max_query_time) max_query_time=t;
			total_query_time+=t;
			HFPointFlag[id]=true;
			++cnt_hf;
		}
		if(cnt_hf<numOfHfpoint){
			numOfHfpoint=cnt_hf;
			std::cout<<"real numOfHfpoint = "<<numOfHfpoint<<std::endl;
		}
		std::cout<<"max_query_time = "<<max_query_time<<std::endl;
		std::cout<<"toal_query_time = "<<total_query_time<<std::endl;
	}	

	/**
	 * @description: read overlay labels from file
	 * @Author: wanjingyi
	 * @Date: 2021-02-03 19:01:34
	 * @param {const} char
	 * @return {*}
	 */
	void load_overlay_labels(const char* load_filename){
		ifstream ifs(load_filename);
		NodeID isize = 0;
		ifs.read((char*)&isize, sizeof(isize));
		if(isize!=numOfOverlayVertices) cout<<"error: isize="<<isize<<std::endl;
		for (NodeID i= 0; i < numOfOverlayVertices; ++i) {
			NodeID v=newToOriginal[i];
			overlay_index_t1& idx = index_overlay[v];
			ifs.read((char*)&isize, sizeof(isize));
			//cout<<"i="<<i<<" v="<<v<<" isize="<<isize<<" size()="<<idx.spt_d.size()<<std::endl;
			//if(idx.spt_v.size()!=0) cout<<" has size="<<idx.spt_v.size()<<std::endl;
			//else cout<<std::endl;
			idx.spt_v.resize(isize);
			idx.spt_d.resize(isize);
			//read hubs for each vertex
			for (NodeID j = 0; j < isize; ++j) {
				NodeID hub;
				EdgeWeight hub_weight;
				ifs.read((char*)&hub, sizeof(hub));
				ifs.read((char*)&hub_weight, sizeof(hub_weight));
				idx.spt_v[j] = hub;
				idx.spt_d[j] = hub_weight;
			}
			//adjust the last one
			idx.spt_v[isize-1]=numOfOriginalVertices;
		}
		ifs.close();
	}

	/**
	 * @description: load overlay vertices indexing order
	 * @Author: wanjingyi
	 * @Date: 2021-02-03 20:28:32
	 * @param {const} char
	 * @return {*}
	 */
	void load_overlay_rank(const char*load_filename){
		//resize
		overlay_rank.resize(numOfOriginalVertices);
		overlay_inv.resize(numOfOriginalVertices);
		ifstream in(load_filename);
		if(!in.is_open()) cout<<"Cannot open "<<load_filename<<std::endl;
		NodeID v,i;
		for(i=0;in>>v;++i){
			overlay_rank[newToOriginal[v]]=i;
			overlay_inv[i]=newToOriginal[v];
		}
		if(i!=numOfOverlayVertices) cout<<"error:i!=numOfOverlayVertices i="<<i<<std::endl;
		in.close();
	}

	/**
	 * @description: bfs for each vertex
	 * @Author: wanjingyi
	 * @Date: 2021-02-04 17:17:27
	 * @param {NodeID} source
	 * @return {*}
	 */
	void source_bfs(NodeID source,WGraph& wgraph,vector<bool>& vis,benchmark::heap<2, EdgeWeight, NodeID>& pqueue,vector<EdgeWeight>& distances,vector<NodeID>& precedants,overlay_index_t1& idx_source){
		vector<pair<NodeID,EdgeWeight> > hf_neighbors;
		vector<pair<NodeID,EdgeWeight> > lf_neighbors;
		lf_neighbors.push_back(make_pair(source,0));
		pqueue.update(source,0);
		distances[source]=0;
		NodeID u,v; EdgeWeight u_w,v_w,v_d;
		while(!pqueue.empty()){
			bool isPruned=false;
			pqueue.extract_min(u,u_w);
			vis[u]=true;
			if(u!=source){
				if(!isDeleted[u]){
					hf_neighbors.push_back(make_pair(u,u_w));
					isPruned=true;
				}
				else{
					if(isDeleted[precedants[u]]){
						lf_neighbors.push_back(make_pair(u,u_w));
					}else{
						isPruned=true;
					}
				} 
			}
			//update distances and precedants
			if(isPruned){
				for (EdgeID eid = wgraph.vertices[u]; eid < wgraph.vertices[u + 1]; ++eid) {
					v = wgraph.edges[eid].first;
					v_w=wgraph.edges[eid].second;
					v_d=distances[u]+v_w;
					if(!vis[v]&&v_d<=distances[v]){//<= further pruned
						distances[v]=v_d;
						precedants[v]=u;
						if(pqueue.contains(v)) pqueue.update(v,v_d);
					}
				}
			}else{
				for (EdgeID eid = wgraph.vertices[u]; eid < wgraph.vertices[u + 1]; ++eid) {
					v = wgraph.edges[eid].first;
					v_w=wgraph.edges[eid].second;
					v_d=distances[u]+v_w;
					if(!vis[v]){
						if(v_d<distances[v]){
							distances[v]=v_d;
							precedants[v]=u;
							pqueue.update(v,v_d);
						}
					}
				}
			}
		}
		EdgeWeight last_w=INF_WEIGHT;
		int i=0,j,hf_len=hf_neighbors.size()+1,len=lf_neighbors.size()+hf_neighbors.size()+1;
		idx_source.spt_v.resize(len);
		idx_source.spt_d.resize(len);
		idx_source.spt_v[i]=hf_len;
		idx_source.spt_d[i++]=len;
		//sort
		sort(hf_neighbors.begin(),hf_neighbors.end());
		for(j=0;i<hf_len;++i){
			idx_source.spt_v[i]=hf_neighbors[j].first;
			idx_source.spt_d[i]=hf_neighbors[j++].second;
		}
		for(j=0;i<len;++i){
			idx_source.spt_v[i]=lf_neighbors[j].first;
			idx_source.spt_d[i]=lf_neighbors[j++].second;
		}
		//clear tmp list
		for(i=0;i<numOfOriginalVertices;++i){
			vis[i]=false;
			distances[i]=INF_WEIGHT;
			precedants[i]=numOfOriginalVertices;
		}
		hf_neighbors.clear();
		lf_neighbors.clear();
		pqueue.clear();
	}

	// void source_bfs(NodeID source,WGraph& wgraph,vector<bool>& vis,benchmark::heap<2, EdgeWeight, NodeID>& pqueue,vector<EdgeWeight>& distances,overlay_index_t1& idx_source){
	// 	vector<pair<EdgeWeight,NodeID> > hf_neighbors;
	// 	vector<pair<EdgeWeight,NodeID> > lf_neighbors;
	// 	pqueue.update(source,0);
	// 	distances[source]=0;
	// 	NodeID u,v; EdgeWeight u_w,v_w,v_d;
	// 	while(!pqueue.empty()){
	// 		pqueue.extract_min(u,u_w);
	// 		vis[u]=true;
	// 		if(u!=source){
	// 			if(!isDeleted[u]) hf_neighbors.push_back(make_pair(u_w,u));
	// 			else lf_neighbors.push_back(make_pair(u_w,u));
	// 		}
	// 		if(!isDeleted[u]) goto pruned;
	// 		for (EdgeID eid = wgraph.vertices[u]; eid < wgraph.vertices[u + 1]; ++eid) {
	// 			v = wgraph.edges[eid].first;
    //             v_w=wgraph.edges[eid].second;
    //             v_d=distances[u]+v_w;
	// 			if(!vis[v]){
	// 				if(v_d<distances[v]){
	// 					distances[v]=v_d;
	// 					pqueue.update(v,v_d);
	// 				}
	// 			}
	// 		}
	// 		pruned:{}
	// 	}
	// 	//****************store index*****************
	// 	//sort(neighbors.begin(),neighbors.end(),cmp_neg);
	// 	EdgeWeight last_w=INF_WEIGHT;
	// 	// for(size_t i=0;i<neighbors.size();++i){
	// 	// 	v=neighbors[i].second;
	// 	// 	v_d=neighbors[i].first;
	// 	// 	if(v_d==last_w) continue;
	// 	// 	idx_source.spt_v.push_back(overlay_inv[v]);
	// 	// 	idx_source.spt_d.push_back(v_d);
	// 	// 	last_w=v_d;
	// 	// }
	// 	//store the (lf_staart_position,lf_size)
	// 	idx_source.spt_v.push_back(hf_neighbors.size()+1);
	// 	idx_source.spt_d.push_back(lf_neighbors.size()+hf_neighbors.size()+1);
	// 	for(size_t i=0;i<hf_neighbors.size();++i){
	// 		idx_source.spt_v.push_back(hf_neighbors[i].second);
	// 		idx_source.spt_d.push_back(hf_neighbors[i].first);
	// 	}
	// 	for(size_t i=0;i<lf_neighbors.size();++i){
	// 		idx_source.spt_v.push_back(lf_neighbors[i].second);
	// 		idx_source.spt_d.push_back(lf_neighbors[i].first);
	// 	}
	// 	//clear tmp list
	// 	for(NodeID i=0;i<numOfOriginalVertices;++i){
	// 		vis[i]=false;
	// 		distances[i]=INF_WEIGHT;
	// 	}
	// 	hf_neighbors.clear();
	// 	lf_neighbors.clear();
	// 	pqueue.clear();
	// }

	/**
	 * @description: bfs search for deleted vertices to find the left neigbors 
	 * @Author: wanjingyi
	 * @Date: 2021-02-03 19:33:03
	 * @param {*}
	 * @return {*}
  	*/ 
 	void bfs_local_search(WGraph& original_wgrah,int _num_threads=5){
		int num_threads = _num_threads;
		//variables
		vector<vector<bool> > vis(num_threads, vector<bool>(numOfOriginalVertices,false)); //vis for omp
		vector<vector<EdgeWeight> > distances(num_threads, vector<EdgeWeight>(numOfOriginalVertices, INF_WEIGHT));
		vector<vector<NodeID> > precedants(num_threads,vector<NodeID>(numOfOriginalVertices,numOfOriginalVertices));
        vector<benchmark::heap<2, EdgeWeight, NodeID> > pqueue(num_threads, benchmark::heap<2, EdgeWeight, NodeID>(numOfOriginalVertices) );
		omp_set_num_threads(num_threads);//设置线程数
		cout<<"start multithread compute..."<<std::endl;//to be deleted
		#pragma omp parallel for schedule(dynamic)
		//#pragma omp parallel for
			for (NodeID v = 0; v<numOfOriginalVertices; ++v) {
				if(!isDeleted[v]) continue;
				source_bfs(v,original_wgrah,vis[omp_get_thread_num()],pqueue[omp_get_thread_num()],distances[omp_get_thread_num()],precedants[omp_get_thread_num()],index_overlay[v]);
			}
		cout<<"start multithread check..."<<std::endl;//to be deleted
		double check_time=GetCurrentTimeSec();
		omp_set_num_threads(num_threads);//设置线程数
		#pragma omp parallel for schedule(dynamic)
			for (NodeID v = 0; v<numOfOriginalVertices; ++v) {
				if(!isDeleted[v]) continue;
				check_region_sp_distances(index_overlay[v]);
			}
		check_time=GetCurrentTimeSec()-check_time;
		cout<<"check_time="<<check_time*1e-6<<endl;
	}

	/**
  * @Descripttion: 
  * @version: 
  * @Author: Wan Jingyi
  * @Date: 2021-04-11 21:04:46
  * @LastEditors: Wan Jingyi
  * @LastEditTime: Do not Edit
  * @param {NodeID} s
  */ 
 void check_region_sp_distances(overlay_index_t1& idx_s){
	 	EdgeWeight distance,distance_tmp,td;
		int s_end=idx_s.spt_v[0];
		_mm_prefetch(&idx_s.spt_v[0], _MM_HINT_T0);
		_mm_prefetch(&idx_s.spt_d[0], _MM_HINT_T0);
		for(int k=s_end;k<idx_s.spt_d[0];++k){
			NodeID t=idx_s.spt_v[k];
			distance=idx_s.spt_d[k];
			const overlay_index_t1& idx_t=index_overlay[t];
			int t_end=idx_t.spt_v[0];
			distance_tmp=INF_WEIGHT;
			_mm_prefetch(&idx_t.spt_v[0], _MM_HINT_T0);
			_mm_prefetch(&idx_t.spt_d[0], _MM_HINT_T0);
			for (int i = 1, j = 1; ; ) {
				if (i == s_end||j==t_end) break;  // Sentinel
				NodeID v1 = idx_s.spt_v[i], v2 = idx_t.spt_v[j];
				if (v1 == v2) {
					td = idx_s.spt_d[i] + idx_t.spt_d[j];
					if (td < distance_tmp) distance_tmp = td;
					++i;
					++j;
				} 
				else {
					i += v1 < v2 ? 1 : 0;
					j += v1 > v2 ? 1 : 0;
				}
			}
			//update
			if(distance_tmp<distance) idx_s.spt_d[k]=distance_tmp;
		}
	}

};

class Overlay_labels: public Overlay_graph{
	public:
		//****************variables****************
		//HFLabel labels;
		vector<overlay_index_t1> index_overlay;//all index
		vector<NodeID> overlay_rank;//fetch rank by original index
		vector<NodeID> overlay_inv;//fetch original index by rank

		//**************constructions****************
		Overlay_labels(char* newToOriginalFileName,char* isDeletedFileName,char* deletedLabelsFileName,char* overlayOrderFileName,char* overlayLabelsFileName,char* outputDirName,char* queryFreqFileName,int hfRate){
			load_is_deleted_and_map(isDeletedFileName);
			cout<<"load deleted map successfully!"<<std::endl;
			load_id_map(newToOriginalFileName);
			cout<<"load overlay map successfully!"<<std::endl;
			//load overlay rank
			load_overlay_rank(overlayOrderFileName);
			cout<<"load overlay rank successfully!"<<std::endl;
			index_overlay.resize(numOfVertices);
			cout<<"index_overlay resize "<<numOfVertices<<" successfully!"<<std::endl;
			//load overlay labels
			load_overlay_labels(overlayLabelsFileName);
			cout<<"load overlay labels successfully!"<<std::endl;
			//load deleted load_overlay_label
			load_deleted_labels(deletedLabelsFileName);
			cout<<"load deleted labels successfully!"<<std::endl;
			//save_and_analyze_index(outputDirName,queryFreqFileName,hfRate);
		}

		~Overlay_labels(){
		}

	protected:
		void merge_labels(char* newToOriginalFileName,char* isDeletedFileName,char* deletedLabelsFileName,char* overlayOrderFileName,char* overlayLabelsFileName){
			//load map relationship
			//numOfOriginalVertices, numOfOverlayVertices and numOfDeletedVertices read
			load_is_deleted_and_map(isDeletedFileName);
			cout<<"load deleted map successfully!"<<std::endl;
			load_id_map(newToOriginalFileName);
			cout<<"load overlay map successfully!"<<std::endl;
			//load overlay rank
			load_overlay_rank(overlayOrderFileName);
			cout<<"load overlay rank successfully!"<<std::endl;
			index_overlay.resize(numOfVertices);
			cout<<"index_overlay resize successfully!"<<std::endl;
			//load overlay labels
			load_overlay_labels(overlayLabelsFileName);
			cout<<"load overlay labels successfully!"<<std::endl;
			//load deleted load_overlay_label
			load_deleted_labels(deletedLabelsFileName);
			cout<<"load deleted labels successfully!"<<std::endl;
		}

		// void save_and_analyze_index(char* outputDirName,char* queryFreqFileName,int hfRate){
		// 	//save all labels to binaty file
		// 	save_all_labels(outputDirName);
		// 	//write labels_lias and size for debug
		// 	write_labels(outputDirName);
		// 	//load hfpoint for analysis
		// 	load_hfpoint_setting(queryFreqFileName,hfRate);
		// 	//save_analysis
		// 	save_analysis(outputDirName);
		// 	//save_analysis_overlay
		// 	save_analysis_overlay(outputDirName);
		// }

		/**
		 * @description: load left hfpoint setted before node deletion
		 * @Author: wanjingyi 
		 * @Date: 2021-02-01 11:49:02
		 * @param {*}
		 * @return {*}
		 */
		void load_hfpoint_setting(char* HFPoint_file,int hfRate){
			cout<<"**************load_hfpoint_setting start!****************"<<std::endl;
			if(hfRate==0){
				numOfHfpoint=numOfOriginalVertices;
			}else{
				numOfHfpoint= 0;//first line is the number of HFpoints
				numOfHfpoint= static_cast<int> ( (double)numOfOriginalVertices*hfRate/(double)HF_DIVIDION);
			}
			if(numOfHfpoint<=0) cout<<"error:numOfHfpoint<=0"<<std::endl;
			HFPointFlag.resize(numOfOriginalVertices,0);
			queryTime.resize(numOfOriginalVertices,0);
			max_query_time=0;
			total_query_time=0;
			cout<<"hfpoint resize successfully!"<<std::endl;
			load_hfpoint_overlay(HFPoint_file);
			cout<<"**************load_hfpoint_setting finished!****************"<<std::endl;
		}

		/*
		*@description: load hfpoint
		*@author: wanjingyi
		*@date: 2021-01-16
		*/
		void load_hfpoint_overlay(char* load_filename){
			std::cout<<"initial numOfHfpoint = "<<numOfHfpoint<<std::endl;
			std::ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int id;int t;
			int cnt=0,cnt_hf=0;
			char line[24];
			//read each line representing HFpoint to vector 
			while (in.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>id>>t;		
				//cnt
				queryTime[id]=t;
				if(cnt==0) max_query_time=t;
				total_query_time+=t;
				if(cnt_hf<numOfHfpoint)
				{
					//HFPoint.push_back(id);
					HFPointFlag[id]=true;
					++cnt_hf;
				}
				++cnt;
			}
			if(cnt_hf<numOfHfpoint){
				numOfHfpoint=cnt_hf;
				std::cout<<"real numOfHfpoint = "<<numOfHfpoint<<std::endl;
			}
			std::cout<<"max_query_time = "<<max_query_time<<std::endl;
			std::cout<<"toal_query_time = "<<total_query_time<<std::endl;
		}	

		/**
		 * @description: read overlay labels from file
		 * @Author: wanjingyi
		 * @Date: 2021-02-03 19:01:34
		 * @param {const} char
		 * @return {*}
		 */
		void load_overlay_labels(const char* load_filename){
			ifstream ifs(load_filename);
			NodeID isize = 0;
			ifs.read((char*)&isize, sizeof(isize));
			if(isize!=numOfOverlayVertices) cout<<"error: isize="<<isize<<std::endl;
			for (NodeID i= 0; i < numOfOverlayVertices; ++i) {
				NodeID v=newToOriginal[i];
				overlay_index_t1& idx = index_overlay[v];
				ifs.read((char*)&isize, sizeof(isize));
				cout<<"i="<<i<<" v="<<v<<" isize="<<isize<<" size()="<<idx.spt_d.size()<<std::endl;
				//if(idx.spt_v.size()!=0) cout<<" has size="<<idx.spt_v.size()<<std::endl;
				//else cout<<std::endl;
				idx.spt_v.resize(isize);
				idx.spt_d.resize(isize);
				//read hubs for each vertex
				for (NodeID j = 0; j < isize; ++j) {
					NodeID hub;
					EdgeWeight hub_weight;
					ifs.read((char*)&hub, sizeof(hub));
					ifs.read((char*)&hub_weight, sizeof(hub_weight));
					idx.spt_v[j] = hub;
					idx.spt_d[j] = hub_weight;
				}
				//adjust the last one
				idx.spt_v[isize-1]=numOfOriginalVertices;
			}
			ifs.close();
		}

		/**
	 * @Author: wanjingyi
	 * @description: load deleted labels from binary file
	 * @param {const} char
	 * @return {*}
	 */  
    void load_deleted_labels(const char* load_filename){
		   //vector<overlay_index_t1> & index_overlay=labels.index_;
			ifstream ifs(load_filename);
			NodeID isize = 0,v;
			ifs.read((char*)&isize, sizeof(isize));
			if(isize!=numOfDeletedVertices) cout<<"error: isize="<<isize<<std::endl;
			for (NodeID i= 0; i < numOfDeletedVertices;++i) {
				ifs.read((char*)&isize, sizeof(isize));
				v=deleteToOriginal[i];
				index_overlay[v].spt_v.resize(isize);
				index_overlay[v].spt_d.resize(isize);
				//read hubs for each vertex
				for (NodeID j = 0; j < isize; ++j) {
					NodeID hub;
					EdgeWeight hub_weight;
					ifs.read((char*)&hub, sizeof(hub));
					ifs.read((char*)&hub_weight, sizeof(hub_weight));
					index_overlay[v].spt_v[j] = hub;
					index_overlay[v].spt_d[j] = hub_weight;
				}
			}
			ifs.close();
		}

		/**
		 * @description: load overlay vertices indexing order
		 * @Author: wanjingyi
		 * @Date: 2021-02-03 20:28:32
		 * @param {const} char
		 * @return {*}
		 */
		void load_overlay_rank(const char*load_filename){
			//resize
			overlay_rank.resize(numOfOriginalVertices);
			overlay_inv.resize(numOfOriginalVertices);
			ifstream in(load_filename);
			if(!in.is_open()) cout<<"Cannot open "<<load_filename<<std::endl;
			NodeID v,i;
			for(i=0;in>>v;++i){
				overlay_rank[newToOriginal[v]]=i;
				overlay_inv[i]=newToOriginal[v];
			}
			if(i!=numOfOverlayVertices) cout<<"error:i!=numOfOverlayVertices i="<<i<<std::endl;
			in.close();
		}
};

/*
 *@description: Distance Preserved graph for deleteNode
 *@author: wanjingyi
 *@date: 2021-01-20
*/
class Node_deletion:public Overlay_graph
{
	public:
		//variables
		vector<list<edge>> removedEdge;//deleted node degree in original graph
		vector<bool> isDeletedLFNode; //indicate whether the node need to delete in LFPointx
		vector<int> level;//delete order by vertex index
		vector<NodeID> _freq_inv;//fetch the index by rank
        vector<NodeID> _freq_rank;//fetch the rank  by index
		int _numOfDeletedVertices;//num of deleted vertices
		//int _numOfDeletedEdges;//num of deleted edges
		bool isDeleteLHPoint;//whether to delete hfPoint

		//variables for hfpoint delete way
		vector<int> _delete_rank;//fetch the delete rank by vertex index
		vector<int> _delete_inv;//fetch the vertex index by delete rank
		int max_degree;//original max_degree
		int _max_degree;//overlay_graph max_degree
		int thresholdDegree;//degree threshold for delete Lfpoint

		double dtime;//delete_time

		//consturctions abd deconstructions
		Node_deletion(){
			_numOfDeletedVertices=0;
			isDeleteLHPoint = false;
			removedEdge.clear();
			isDeletedLFNode.clear();
			level.clear();
			_delete_rank.clear();
			_delete_inv.clear();
			thresholdDegree=1;
		}
		~Node_deletion(){}

		/*
		 *@description: load original graph and bankup to es for deletion
		 *@author: wanjingyi
		 *@date: 2021-01-25
		*/
		// void load_and_copy_original_graph(const char* graph_file){
		// 	load_original_graph(graph_file,es_original);
		// 	es.reserve(numOfVertices);
		// 	//es.assign(es_original.begin(),es_original.end());//clear and deep copy
		// 	es=es_original;
		// }

		/*
		 *@description: load original graph and bankup to es for deletion
		 *@author: wanjingyi
		 *@date: 2021-01-25
		*/
		void load_and_copy_original_graph(const char* graph_file,int graph_type){
			load_original_graph(graph_file,graph_type);
			es.reserve(numOfVertices);
			//es.assign(es_original.begin(),es_original.end());//clear and deep copy
			es=es_original;
			if(DIRECTED_FLAG==true){
				r_es.reserve(numOfVertices);
				r_es=r_es_original;
			}
		}

		//Add for directed
		/**
		 * @description: 
		 * @param {int} deleteModel
		 * @param {int} MAXDEGREE
		 * @param {int} deleteNodeRate
		 * @param {int} hfRate
		 * @return {*}
		 * @author: Wan Jingyi
		 */  
  		void deleteNode_directed(int deleteModel,int MAXDEGREE,int deleteNodeRate,int hfRate,char* queryPair_file){
			dtime=GetCurrentTimeSec();
			if(MAXDEGREE==0) thresholdDegree=INF_WEIGHT;//means delete all non-hfpoint
			else thresholdDegree=MAXDEGREE;
			if(hfRate==0){
					numOfHfpoint=numOfVertices;
			}else{
				numOfHfpoint= 0;//first line is the number of HFpoints
				numOfHfpoint= static_cast<int> ( (double)numOfVertices*hfRate/(double)HF_DIVIDION);
			}
			if(numOfHfpoint<=0) cout<<"error:numOfHfpoint<=0"<<std::endl;
			if(deleteModel==0){//delete node by degree
				std::cout<<"deletemodel==0 minDegree..."<<std::endl;
				deleteNodeByDegree_directed();
			}else if(deleteModel==1){ //delete node by degree and freq (no local search)
				std::cout<<"deletemodel==1 minDegreeFreq..."<<std::endl;
			}
			return;
		}

		/*
		 *@description: deleteNode by degreeOrder and frequency
		 *@author: gaoyongyong
		 *@date: 2021-01-21
		*/
		void deleteNode(int deleteModel,int MAXDEGREE,int deleteNodeRate,int hfRate,char* HFPoint_file,int bfsHop_num,int searchThreshold){
			dtime=GetCurrentTimeSec();
			if(MAXDEGREE==0) thresholdDegree=INF_WEIGHT;//means delete all non-hfpoint
			else thresholdDegree=MAXDEGREE;
			if(hfRate==0){
					numOfHfpoint=numOfVertices;
			}else{
				numOfHfpoint= 0;//first line is the number of HFpoints
				numOfHfpoint= static_cast<int> ( (double)numOfVertices*hfRate/(double)HF_DIVIDION);
			}
			if(numOfHfpoint<=0) cout<<"error:numOfHfpoint<=0"<<std::endl;
			if(deleteModel==0){//delete node by degree
				//std::cout<<"Delete node by degree..."<<std::endl;
				initIterms();//resize variables
				//first delete the vertices with degree 1
				for (int i = 0; i < numOfVertices; i++) {
					if (1 == es[i].size()) {
						degreeOneDelete(es, removedEdge, isDeleted, i);
					}
				}
				//std::cout<<"Delete all vertices with degree 1 finished!"<<std::endl;
				//deleteNode by ascending degree order
				for (int degree = 2; degree <= MAXDEGREE; degree++) {
					for (int i = 0; i < numOfVertices; i++) {
						if (degree == es[i].size()) {
							degreeDelete(es, removedEdge, isDeleted, i, degree);
							//cout << "delete cur node " << i << std::endl;
						}
					}
				}
				//std::cout<<"deleteNode by ascending degree order(from 2) finished!"<<std::endl;		
			}else if (deleteModel==1){//delete node by freq degree
				//std::cout<<"Delete Hfpoint by hfpoint..."<<std::endl;
				//init variables
				HFPointFlag.resize(numOfVertices,0);//resize
				isDeleted.resize(numOfVertices,0);
				_delete_inv.reserve(numOfVertices);
				//load hfpoint
				load_hfpoint(HFPoint_file,numOfHfpoint,HFPoint,HFPointFlag);
				if(deleteNodeRate==0){//means left HFPoint
					_numOfDeletedVertices=numOfVertices-numOfHfpoint;
				}else{
					double _deleteNodeRate=(double)deleteNodeRate/(double)1000;
					_numOfDeletedVertices=static_cast<int>(numOfVertices* _deleteNodeRate)+1;
				}
				//cout<<"initial numOfDeletedVertices = "<<_numOfDeletedVertices<<std::endl;
				//rankVerticesByFreq();
				int result=deleteNodeByFreqDegree();
			}else if(deleteModel==2){//delete node by freq and local search
				//std::cout<<"Delete Hfpoint by hfpoint and local search..."<<" bfsHop_num="<<bfsHop_num<<" searchThreshold="<<searchThreshold<<std::endl;
				//init variables
				HFPointFlag.resize(numOfVertices,0);//resize
				isDeleted.resize(numOfVertices,0);
				_delete_inv.reserve(numOfVertices);
				//load hfpoint
				load_hfpoint(HFPoint_file,numOfHfpoint,HFPoint,HFPointFlag);
				if(deleteNodeRate==0){//means left HFPoint
					_numOfDeletedVertices=numOfVertices-numOfHfpoint;
				}else{
					double _deleteNodeRate=(double)deleteNodeRate/(double)1000;
					_numOfDeletedVertices=static_cast<int>(numOfVertices* _deleteNodeRate)+1;
				}
				int result=deleteNodeByFreqDegree(bfsHop_num);
			}else if(deleteModel==3){//delete node by limited local serach
				//std::cout<<"Delete Hfpoint by hfpoint and local search..."<<" bfsHop_num="<<bfsHop_num<<" searchThreshold="<<searchThreshold<<std::endl;
				//init variables
				HFPointFlag.resize(numOfVertices,0);//resize
				isDeleted.resize(numOfVertices,0);
				_delete_inv.reserve(numOfVertices);
				//load hfpoint
				load_hfpoint(HFPoint_file,numOfHfpoint,HFPoint,HFPointFlag);
				if(deleteNodeRate==0){//means left HFPoint
					_numOfDeletedVertices=numOfVertices-numOfHfpoint;
				}else{
					double _deleteNodeRate=(double)deleteNodeRate/(double)1000;
					_numOfDeletedVertices=static_cast<int>(numOfVertices* _deleteNodeRate)+1;
				}
				int result=deleteNodeByFreqDegree(bfsHop_num,searchThreshold);
			}
			dtime=GetCurrentTimeSec()-dtime;

			//analyze and count the overlay_graph
			int numOfDeletedVertices=0;//num of deleted nodes
			_numOfEdges=0;
			int degreeOne=0;//num of nodes with left degree 1
			//vector<int> sumOfBlanks(MAXDEGREE, 0);
			for (int i = 0; i < numOfVertices; i++) {
				int isize=es[i].size();
				_numOfEdges+=isize;
				if (isize == 0) {
					numOfDeletedVertices++;
					//sumOfBlanks[es_original[i].size()]++;
				}
				if (isize == 1) {	
					degreeOne++;
				}
			}

			numOfOverlayVertices=numOfVertices-numOfDeletedVertices;
			std::cout<<"numOfHfpoint = "<<numOfHfpoint<<std::endl;
			std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<std::endl;
			std::cout<<"numOfDeletedVertices = "<<numOfDeletedVertices<<std::endl;
			std::cout<<"_numOfEdges = "<<_numOfEdges<<std::endl;
			double ave_dtime=dtime/(double)numOfDeletedVertices;
			std::cout<<"deleteNode time = "<<dtime*1e6<<" microseconds"<<std::endl;
			std::cout<<"average deleteNode time = "<<ave_dtime*1e6<<" microseconds"<<std::endl;
			//std::cout<<"numOfDeletedEdges = "<<numOfEdges-_numOfEdges<<std::endl;
		}

		/*
		 *@description: deleteNode by degreeOrder and frequency
		 *@author: gaoyongyong
		 *@date: 2021-01-21
		*/
		void deleteNode(int deleteModel,int MAXDEGREE,int deleteNodeRate,int hfRate,char* HFPoint_file,int bfsHop_num,int searchThreshold,char* analysisFileName){
			dtime=GetCurrentTimeSec();
			if(MAXDEGREE==0) thresholdDegree=INF_WEIGHT;//means delete all non-hfpoint
			else thresholdDegree=MAXDEGREE;
			if(hfRate==0){
					numOfHfpoint=numOfVertices;
			}else{
				numOfHfpoint= 0;//first line is the number of HFpoints
				numOfHfpoint= static_cast<int> ( (double)numOfVertices*hfRate/(double)HF_DIVIDION);
			}
			if(numOfHfpoint<=0) cout<<"error:numOfHfpoint<=0"<<std::endl;
			if(deleteModel==0){//delete node by degree
				//std::cout<<"Delete node by degree..."<<std::endl;
				initIterms();//resize variables
				//first delete the vertices with degree 1
				for (int i = 0; i < numOfVertices; i++) {
					if (1 == es[i].size()) {
						degreeOneDelete(es, removedEdge, isDeleted, i);
					}
				}
				//std::cout<<"Delete all vertices with degree 1 finished!"<<std::endl;
				//deleteNode by ascending degree order
				for (int degree = 2; degree <= MAXDEGREE; degree++) {
					for (int i = 0; i < numOfVertices; i++) {
						if (degree == es[i].size()) {
							degreeDelete(es, removedEdge, isDeleted, i, degree);
							//cout << "delete cur node " << i << std::endl;
						}
					}
				}
				//std::cout<<"deleteNode by ascending degree order(from 2) finished!"<<std::endl;		
			}else if (deleteModel==1){//delete node by freq degree
				//std::cout<<"Delete Hfpoint by hfpoint..."<<std::endl;
				//init variables
				HFPointFlag.resize(numOfVertices,0);//resize
				isDeleted.resize(numOfVertices,0);
				_delete_inv.reserve(numOfVertices);
				//load hfpoint
				load_hfpoint(HFPoint_file,numOfHfpoint,HFPoint,HFPointFlag);
				if(deleteNodeRate==0){//means left HFPoint
					_numOfDeletedVertices=numOfVertices-numOfHfpoint;
				}else{
					double _deleteNodeRate=(double)deleteNodeRate/(double)1000;
					_numOfDeletedVertices=static_cast<int>(numOfVertices* _deleteNodeRate)+1;
				}
				//cout<<"initial numOfDeletedVertices = "<<_numOfDeletedVertices<<std::endl;
				//rankVerticesByFreq();
				int result=deleteNodeByFreqDegree();
			}else if(deleteModel==2){//delete node by freq and local search
				//std::cout<<"Delete Hfpoint by hfpoint and local search..."<<" bfsHop_num="<<bfsHop_num<<" searchThreshold="<<searchThreshold<<std::endl;
				//init variables
				HFPointFlag.resize(numOfVertices,0);//resize
				isDeleted.resize(numOfVertices,0);
				_delete_inv.reserve(numOfVertices);
				//load hfpoint
				load_hfpoint(HFPoint_file,numOfHfpoint,HFPoint,HFPointFlag);
				if(deleteNodeRate==0){//means left HFPoint
					_numOfDeletedVertices=numOfVertices-numOfHfpoint;
				}else{
					double _deleteNodeRate=(double)deleteNodeRate/(double)1000;
					_numOfDeletedVertices=static_cast<int>(numOfVertices* _deleteNodeRate)+1;
				}
				int result=deleteNodeByFreqDegree(bfsHop_num);
			}else if(deleteModel==3){//delete node by limited local serach
				//std::cout<<"Delete Hfpoint by hfpoint and local search..."<<" bfsHop_num="<<bfsHop_num<<" searchThreshold="<<searchThreshold<<std::endl;
				//init variables
				HFPointFlag.resize(numOfVertices,0);//resize
				isDeleted.resize(numOfVertices,0);
				_delete_inv.reserve(numOfVertices);
				//load hfpoint
				load_hfpoint(HFPoint_file,numOfHfpoint,HFPoint,HFPointFlag);
				if(deleteNodeRate==0){//means left HFPoint
					_numOfDeletedVertices=numOfVertices-numOfHfpoint;
				}else{
					double _deleteNodeRate=(double)deleteNodeRate/(double)1000;
					_numOfDeletedVertices=static_cast<int>(numOfVertices* _deleteNodeRate)+1;
				}
				int result=deleteNodeByFreqDegree(bfsHop_num,searchThreshold);
			}
			dtime=GetCurrentTimeSec()-dtime;

			//analyze and count the overlay_graph
			int numOfDeletedVertices=0;//num of deleted nodes
			_numOfEdges=0;
			int degreeOne=0;//num of nodes with left degree 1
			//vector<int> sumOfBlanks(MAXDEGREE, 0);
			for (int i = 0; i < numOfVertices; i++) {
				int isize=es[i].size();
				_numOfEdges+=isize;
				if (isize == 0) {
					numOfDeletedVertices++;
					//sumOfBlanks[es_original[i].size()]++;
				}
				if (isize == 1) {	
					degreeOne++;
				}
			}
			//output experiment result  to file
			numOfOverlayVertices=numOfVertices-numOfDeletedVertices;
			ofstream ofs(analysisFileName,ios::app|ios::out);
			if(!ofs.is_open()) cout<<"Cannot open "<<analysisFileName<<std::endl;
			ofs.setf(ios::fixed);
			ofs.precision(4);
			double ave_dtime=dtime/(double)numOfDeletedVertices;
			std::cout<<"deleteNode time = "<<dtime*1e6<<" microseconds"<<std::endl;
			std::cout<<"average deleteNode time = "<<ave_dtime*1e6<<" microseconds"<<std::endl;
			ofs<<dtime*1e6<<" "<<ave_dtime*1e6<<" "<<numOfOverlayVertices<<" "<<_numOfEdges<<" ";
			ofs.close();
		}

		/*
		 *@description: generate and ouptu overlay graph
		 *@author: wanjingyi
		 *@date: 2021-01-21
		*/
		void generate_overlay_graph(char* output_dir_name,char* analysisFIleName){
			std::cout<<"******************generate and output overlay graph start!************"<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				es[i].sort([](const edge e1, const edge e2) {return e1.v < e2.v; });//sort the edges by neighbor id
			}
			originalToNew.resize(numOfVertices, -1);
			newToOriginal.reserve(numOfOverlayVertices);
			//relabel the undeleted vertices form 0
			//build the map relationship 

			//************output the new to original vertice map************
			for (int i = 0; i < numOfVertices; i++) {
				if (0 != es[i].size() && !isDeleted[i]) {
					newToOriginal.push_back(i);
					originalToNew[i] = newToOriginal.size() - 1;
				}
			}			
			string output_map_file(output_dir_name);
			output_map_file.append(".newToOriginal");
			ofstream ofsmap(output_map_file);
			if(!ofsmap.is_open()) std::cerr<<output_map_file<<" cannot be opened!"<<std::endl;
			ofsmap<<numOfVertices<<" "<<numOfOverlayVertices<<std::endl;
			for (int i = 0; i < newToOriginal.size(); i++) {
				ofsmap <<i<<" "<< newToOriginal[i] << std::endl;
			}
			ofsmap.close();
			std::cout<<"write "<<output_map_file<<" successfully!"<<std::endl;

			//************output whether the original vertice is deleeted 1-deleted 0-no deleted************
			string outputIsDeleted_file(output_dir_name);
			outputIsDeleted_file.append(".isDeleted");
			ofstream ofsIsDeleted(outputIsDeleted_file);
			if(!ofsIsDeleted.is_open()) std::cerr<<outputIsDeleted_file<<" cannot be opened!"<<std::endl;
			ofsIsDeleted<<numOfVertices<<" "<<numOfOverlayVertices<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				ofsIsDeleted <<i<<" "<< isDeleted[i] << std::endl;
			}
			ofsIsDeleted.close();
			std::cout<<"write "<<outputIsDeleted_file<<" successfully!"<<std::endl;

			//************output new graph file using new ids************
			string output_file(output_dir_name);
			output_file.append("_overlay_graph.new");
			ofstream ofs(output_file);
			if(!ofs.is_open()) std::cerr<<output_file<<" cannot be opened!"<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				if (es[i].size() != 0) {
					for (list<edge>::iterator it = es[i].begin(); it != es[i].end(); it++) {
						if (originalToNew[it->u] != -1 && originalToNew[it->v] != -1)
							ofs << originalToNew[it->u] << " " << originalToNew[it->v] << " " << it->weight << std::endl;
					}
				}
			}
			ofs.close();
			std::cout<<"write "<<output_file<<" successfully!"<<std::endl;

			//************output new graph file using original ids************
			string output_file_original(output_dir_name);
			output_file_original.append("_overlay_graph.original");
			ofstream ofs1(output_file_original);
			if(!ofs1.is_open()) std::cerr<<output_file_original<<" cannot be opened!"<<std::endl;
			ofs1<<numOfVertices<<" "<<numOfOverlayVertices<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				if (es[i].size() != 0) {
					for (list<edge>::iterator it = es[i].begin(); it != es[i].end(); it++) {
						ofs1 << it->u << " " << it->v << " " << it->weight << std::endl;
					}
				}
			}
			ofs1.close();
			std::cout<<"write "<<output_file_original<<" successfully!"<<std::endl;

			//************output degree analysis************
			string output_degree_file(output_dir_name);
			int degree_cnt=0,degree_cnt_new=0;
			max_degree=_max_degree=0;
			for (int i = 0; i < numOfVertices; i++) {
				degree_cnt+=es_original[i].size();
				degree_cnt_new+=es[i].size();
				if(es_original[i].size()>max_degree) max_degree=es_original[i].size();
				if(es[i].size()>_max_degree) _max_degree=es[i].size();
			}		
			double avg_degree=(double)degree_cnt/(double)numOfVertices;
			double avg_degree_new=(double)degree_cnt_new/(double)numOfOverlayVertices;

			ofstream ofs_result(analysisFIleName,ios::app|ios::out);
			if(!ofs_result.is_open()) cout<<"Cannot open "<<analysisFIleName<<std::endl;
			ofs_result.setf(ios::fixed);
			ofs_result.precision(4);
			ofs_result<<avg_degree_new<<" "<<_max_degree<<endl;
			ofs_result.close();
			return;
		}

		//Add for directed (debug)
		void generate_overlay_graph_directed(char* output_dir_name,int graph_type=0){
			std::cout<<"******************generate and output directed overlay graph start!************"<<std::endl;
			
			//****************sort and build the map relationship************//
			// for (int i = 0; i < numOfVertices; i++) {
			// 	//sort the edges by neighbor id
			// 	es[i].sort([](const edge e1, const edge e2) {return e1.v < e2.v; });
			// 	r_es[i].sort([](const edge e1, const edge e2) {return e1.u < e2.u; });
			// }

			//************deal with the left edges************//
			// string output_test_file(output_dir_name);
			// output_test_file.append(".test_graph");
			// ofstream ofs_test(output_test_file);
			// if(!ofs_test.is_open()) std::cout<<"Cannot open the "<<output_test_file<<std::endl;

			originalToNew.resize(numOfVertices, -1);
			newToOriginal.reserve(numOfOverlayVertices);
			//build the map relationship
			for (int i = 0; i < numOfVertices; i++) {
				if ((0 != es[i].size()||0!=r_es[i].size())&& !isDeleted[i]) {
					newToOriginal.push_back(i);
					originalToNew[i] = newToOriginal.size() - 1;
				}
			}			
			adj.resize(numOfOverlayVertices);
			adj_weight.resize(numOfOverlayVertices);
			for (int i = 0; i < numOfVertices; i++) {//collect all left directed edges
				if (r_es[i].size() != 0) {
					for (list<edge>::iterator it = r_es[i].begin(); it != r_es[i].end(); it++) {
						if (originalToNew[it->u] != -1 && originalToNew[it->v] != -1){
							adj[originalToNew[it->u]].push_back(originalToNew[it->v]);
							adj_weight[originalToNew[it->u]].push_back(it->weight);
							++_numOfEdges_r;
							//to be deleted
							//ofs_test<<it->u<<" "<<it->v<<" "<<it->weight<<endl;
						}
					}
				}
				if (es[i].size() != 0) {
					for (list<edge>::iterator it = es[i].begin(); it != es[i].end(); it++) {
						if (originalToNew[it->u] != -1 && originalToNew[it->v] != -1){
							adj[originalToNew[it->u]].push_back(originalToNew[it->v]);
							adj_weight[originalToNew[it->u]].push_back(it->weight);
							++_numOfEdges_f;
							//to be deleted
							//ofs_test<<it->u<<" "<<it->v<<" "<<it->weight<<endl;
						}
					}
				}
			}
			//remove the duplicated save edges
			for(int u=0;u<numOfOverlayVertices;++u){
				vector<int>& adj_u=adj[u];
				vector<int>& adj_weight_u=adj_weight[u];
				vector<pair<int,int> > adj_uw_u(adj_u.size());
				for (int i = 0; i < adj_u.size(); ++i) {
					int v = adj_u[i];
					int w = adj_weight_u[i];
					adj_uw_u[i] = make_pair(v, w);
				}
				sort(adj_uw_u.begin(), adj_uw_u.end(),cmp_edge_id);//cmp1
				for (int i = 0; i < adj_u.size(); ++i) {
					adj_u[i] = adj_uw_u[i].first;
					adj_weight_u[i] = adj_uw_u[i].second;
				}
				adj_uw_u.clear();
			}

			//************output the new to original vertice map************//
			string output_map_file(output_dir_name);
			output_map_file.append(".newToOriginal");
			ofstream ofsmap(output_map_file);
			if(!ofsmap.is_open()) std::cerr<<output_map_file<<" cannot be opened!"<<std::endl;
			ofsmap<<numOfOriginalVertices<<" "<<numOfOverlayVertices<<std::endl;
			for (int i = 0; i < newToOriginal.size(); i++) {
				ofsmap <<i<<" "<< newToOriginal[i] << std::endl;
			}
			ofsmap.close();
			std::cout<<"write "<<output_map_file<<" successfully!"<<std::endl;

			//************output whether the original vertice is deleeted 1-deleted 0-no deleted************//
			string outputIsDeleted_file(output_dir_name);
			outputIsDeleted_file.append(".isDeleted");
			ofstream ofsIsDeleted(outputIsDeleted_file);
			if(!ofsIsDeleted.is_open()) std::cerr<<outputIsDeleted_file<<" cannot be opened!"<<std::endl;
			ofsIsDeleted<<numOfVertices<<" "<<numOfOverlayVertices<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				ofsIsDeleted <<i<<" "<< isDeleted[i] << std::endl;
			}
			ofsIsDeleted.close();
			std::cout<<"write "<<outputIsDeleted_file<<" successfully!"<<std::endl;

			//************output new graph file using new ids************//
			string output_file(output_dir_name);
			output_file.append("_overlay_graph.new");
			ofstream ofs(output_file);
			if(!ofs.is_open()) std::cerr<<output_file<<" cannot be opened!"<<std::endl;
			if(graph_type==0){
				for(int u=0;u<numOfOverlayVertices;++u){
					int last_v=numOfOverlayVertices;
					for(int i=0;i<adj[u].size();++i){
						if(adj[u][i]==last_v||adj[u][i]<=u){
							continue;
						}
						ofs<<u<<" "<<adj[u][i]<<" "<<adj_weight[u][i]<<endl;
						++_numOfEdges;
						last_v=adj[u][i];
					}
				}
				_numOfEdges=_numOfEdges*2;
			}else if(graph_type==1){
				for(int u=0;u<numOfOverlayVertices;++u){
					int last_v=numOfOverlayVertices;
					for(int i=0;i<adj[u].size();++i){
						if(adj[u][i]==last_v){
							continue;
						}
						ofs<<u<<" "<<adj[u][i]<<" "<<adj_weight[u][i]<<endl;
						++_numOfEdges;
						last_v=adj[u][i];
					}
				}
			}

			ofs.close();
			std::cout<<"_numOfEdges_f = "<<_numOfEdges_f<<std::endl;
			std::cout<<"_numOfEdges_r = "<<_numOfEdges_r<<std::endl;
			std::cout<<"_numOfEdges = "<<_numOfEdges<<std::endl;
			std::cout<<"write "<<output_file<<" successfully!"<<std::endl;

			//************output new graph file using original ids************
			string output_file_original(output_dir_name);
			output_file_original.append("_overlay_graph.original");
			ofstream ofs1(output_file_original);
			if(!ofs1.is_open()) std::cerr<<output_file_original<<" cannot be opened!"<<std::endl;
			ofs1<<numOfVertices<<" "<<numOfOverlayVertices<<std::endl;
			if(graph_type==0){
				for(int u=0;u<numOfOverlayVertices;++u){
					int last_v=numOfOverlayVertices;
					for(int i=0;i<adj[u].size();++i){
						if(adj[u][i]==last_v||adj[u][i]<=u){
							continue;
						}
						ofs1<<newToOriginal[u]<<" "<<newToOriginal[adj[u][i]]<<" "<<adj_weight[u][i]<<endl;
						last_v=adj[u][i];
					}
				}
			}else if(graph_type==1){
				for(int u=0;u<numOfOverlayVertices;++u){
					int last_v=numOfOverlayVertices;
					for(int i=0;i<adj[u].size();++i){
						if(adj[u][i]==last_v){
							continue;
						}
						ofs1<<newToOriginal[u]<<" "<<newToOriginal[adj[u][i]]<<" "<<adj_weight[u][i]<<endl;
						last_v=adj[u][i];
					}
				}
			}
			ofs1.close();
			std::cout<<"write "<<output_file_original<<" successfully!"<<std::endl;

			//************output delete order************//
			string output_delete_order_file(output_dir_name);
			output_delete_order_file.append(".delete_order");
			ofstream ofsDorder(output_delete_order_file);				
			if(!ofsDorder.is_open()) std::cerr<<output_delete_order_file<<" cannot be opened!"<<std::endl;
			for(int i=0;i<_delete_inv.size();++i){
				ofsDorder<<_delete_inv[i]<<std::endl;
			}
			ofsDorder.close();
			std::cout<<"write "<<output_delete_order_file<<" successfully!"<<std::endl;

			std::cout<<"******************generate and output directed overlay graph finished!************"<<std::endl;
		
			adj.clear();
			adj_weight.clear();
			adj.shrink_to_fit();
			adj_weight.shrink_to_fit();
			return;
		}

		/*
		 *@description: generate and ouptu overlay graph
		 *@author: wanjingyi
		 *@date: 2021-01-21
		*/
		void generate_overlay_graph(char* output_dir_name){
			std::cout<<"******************generate and output overlay graph start!************"<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				es[i].sort([](const edge e1, const edge e2) {return e1.v < e2.v; });//sort the edges by neighbor id
			}
			originalToNew.resize(numOfVertices, -1);
			newToOriginal.reserve(numOfOverlayVertices);
			//relabel the undeleted vertices form 0
			//build the map relationship 

			//************output the new to original vertice map************
			for (int i = 0; i < numOfVertices; i++) {
				if (0 != es[i].size() && !isDeleted[i]) {
					newToOriginal.push_back(i);
					originalToNew[i] = newToOriginal.size() - 1;
				}
			}			
			string output_map_file(output_dir_name);
			output_map_file.append(".newToOriginal");
			ofstream ofsmap(output_map_file);
			if(!ofsmap.is_open()) std::cerr<<output_map_file<<" cannot be opened!"<<std::endl;
			ofsmap<<numOfVertices<<" "<<numOfOverlayVertices<<std::endl;
			for (int i = 0; i < newToOriginal.size(); i++) {
				ofsmap <<i<<" "<< newToOriginal[i] << std::endl;
			}
			ofsmap.close();
			std::cout<<"write "<<output_map_file<<" successfully!"<<std::endl;

			//************output degree analysis************
			string output_degree_file(output_dir_name);
			int degree_cnt=0,degree_cnt_new=0;
			max_degree=_max_degree=0;
			output_degree_file.append(".degree");
			ofstream ofsdeg(output_degree_file);
			if(!ofsdeg.is_open()) std::cerr<<output_degree_file<<" cannot be opened!"<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				degree_cnt+=es_original[i].size();
				degree_cnt_new+=es[i].size();
				if(es_original[i].size()>max_degree) max_degree=es_original[i].size();
				if(es[i].size()>_max_degree) _max_degree=es[i].size();
			}		
			double avg_degree=(double)degree_cnt/(double)numOfVertices;
			double avg_degree_new=(double)degree_cnt_new/(double)numOfOverlayVertices;
			//count the degree distribution
			vector<int> deg_distribute(max_degree/5+1,0);
			vector<int> _deg_distribute(_max_degree/5+1,0);
			for (int i = 0; i < numOfVertices; i++) {
				if(es_original[i].size()!= 0) deg_distribute[es_original[i].size()/5]++;
				if(es[i].size()!= 0) _deg_distribute[es[i].size()/5]++;
			}
			std::cout<<"Original Average degree = "<<avg_degree<<std::endl;
			std::cout<<"Original Max degree = "<<max_degree<<std::endl;
			std::cout<<"Overlay_graph Average degree = "<<avg_degree_new<<std::endl;
			std::cout<<"Overlay_graph Max degree = "<<_max_degree<<std::endl;
			ofsdeg<<"Original Average degree = "<<avg_degree<<std::endl;
			ofsdeg<<"Original Max degree = "<<max_degree<<std::endl;
			ofsdeg<<"Overlay_graph Average degree = "<<avg_degree_new<<std::endl;
			ofsdeg<<"Overlay_graph Max degree = "<<_max_degree<<std::endl;
			//output degree distribution
			ofsdeg<<"Original degree distribution:"<<std::endl;
			for(int i=0;i<deg_distribute.size();i++) ofsdeg<<"["<<5*i<<","<<5*(i+1)<<") "<<deg_distribute[i]<<std::endl;
			ofsdeg<<"Overlay_graph degree distribution:"<<std::endl;
			for(int i=0;i<_deg_distribute.size();i++) ofsdeg<<"["<<5*i<<","<<5*(i+1)<<") "<<_deg_distribute[i]<<std::endl;
			ofsdeg.close();
			std::cout<<"write "<<output_degree_file<<" successfully!"<<std::endl;

			//************output whether the original vertice is deleeted 1-deleted 0-no deleted************
			string outputIsDeleted_file(output_dir_name);
			outputIsDeleted_file.append(".isDeleted");
			ofstream ofsIsDeleted(outputIsDeleted_file);
			if(!ofsIsDeleted.is_open()) std::cerr<<outputIsDeleted_file<<" cannot be opened!"<<std::endl;
			ofsIsDeleted<<numOfVertices<<" "<<numOfOverlayVertices<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				ofsIsDeleted <<i<<" "<< isDeleted[i] << std::endl;
			}
			ofsIsDeleted.close();
			std::cout<<"write "<<outputIsDeleted_file<<" successfully!"<<std::endl;

			//************output new graph file using new ids************
			string output_file(output_dir_name);
			output_file.append("_overlay_graph.new");
			ofstream ofs(output_file);
			if(!ofs.is_open()) std::cerr<<output_file<<" cannot be opened!"<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				if (es[i].size() != 0) {
					for (list<edge>::iterator it = es[i].begin(); it != es[i].end(); it++) {
						if (originalToNew[it->u] != -1 && originalToNew[it->v] != -1)
							ofs << originalToNew[it->u] << " " << originalToNew[it->v] << " " << it->weight << std::endl;
					}
				}
			}
			ofs.close();
			std::cout<<"write "<<output_file<<" successfully!"<<std::endl;

			//************output new graph file using original ids************
			string output_file_original(output_dir_name);
			output_file_original.append("_overlay_graph.original");
			ofstream ofs1(output_file_original);
			if(!ofs1.is_open()) std::cerr<<output_file_original<<" cannot be opened!"<<std::endl;
			ofs1<<numOfVertices<<" "<<numOfOverlayVertices<<std::endl;
			for (int i = 0; i < numOfVertices; i++) {
				if (es[i].size() != 0) {
					for (list<edge>::iterator it = es[i].begin(); it != es[i].end(); it++) {
						ofs1 << it->u << " " << it->v << " " << it->weight << std::endl;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
					}
				}
			}
			ofs1.close();
			std::cout<<"write "<<output_file_original<<" successfully!"<<std::endl;

			//ouput deleted vertices and its neighbors(left)
			//format: firstline-numOfVertices _numOfVertices
			//format:deleted_vertice degree
			//format:neighbor_0 neighbor_1 ...
			//format:weight_0 weight_1 ...
			// string output_region_file(output_dir_name);
			// output_region_file.append("_region_graph.txt");
			// ofstream ofsRegion(output_region_file);
			// bool isFirst;
			// if(!ofsRegion.is_open()) std::cerr<<output_region_file<<" cannot be opened!"<<std::endl;
			// ofsRegion<<numOfVertices<<" "<<_numOfVertices<<std::endl;
			// for(int i=0;i<numOfVertices;++i){
			// 	if(!isDeleted[i]) continue;
			// 	ofsRegion<<i<<" "<<es_original[i].size()<<std::endl;
			// 	if(es_original[i].size()!=0)
			// 	{
			// 		isFirst=true;
			// 		for (list<edge>::iterator it = es_original[i].begin(); it !=es_original[i].end(); it++){
			// 			if(isFirst) isFirst=false;
			// 			else ofsRegion<<" ";
			// 			ofsRegion<<it->v;
			// 		}
			// 		ofsRegion<<std::endl;
			// 		isFirst=true;
			// 		for (list<edge>::iterator it = es_original[i].begin(); it !=es_original[i].end(); it++){
			// 			if(isFirst) isFirst=false;
			// 			else ofsRegion<<" ";
			// 			ofsRegion<<it->weight;
			// 		}
			// 		ofsRegion<<std::endl;				
			// 	}
			// }
			// ofsRegion.close();
			// std::cout<<"write "<<output_region_file<<" successfully!"<<std::endl;		

			//************output deleted vertices with resume_graph************
			//format::firstline-numOfVertices _numOfVertices
			//format: new_id original_id deleted_degree deleted_vertice_0 deleted_vertice_1 ...
			string output_resume_file(output_dir_name);
			output_resume_file.append("_resume_graph.txt");
			ofstream ofsResume(output_resume_file);
			bool is_first;
			if(!ofsResume.is_open()) std::cerr<<output_resume_file<<" cannot be opened!"<<std::endl;
			ofsResume<<numOfVertices<<" "<<numOfOverlayVertices<<std::endl;
			vector<pair<int,int> > deletedEdges;
			for (int v = 0; v < numOfOverlayVertices; v++) {
				int i=newToOriginal[v];
				for (list<edge>::iterator it = es_original[i].begin(); it !=es_original[i].end(); it++){
					list<edge>::iterator iter;
					for (iter= es[i].begin(); iter !=es[i].end()&&iter->v!=it->v; iter++);
					if(iter==es[i].end()) deletedEdges.push_back(make_pair(it->v,it->weight));
				}
				int d_siz=deletedEdges.size();
				ofsResume<<v<<" "<<i<<" "<<d_siz<<std::endl;
				if(d_siz!=0)
				{
					is_first=true;
					for(int i=0;i<d_siz;++i){
						if(is_first) is_first=false;
						else ofsResume<<" ";
						ofsResume<<deletedEdges[i].first;					
					}
					ofsResume<<std::endl;
					is_first=true;
					for(int i=0;i<d_siz;++i){
						if(is_first) is_first=false;
						else ofsResume<<" ";
						ofsResume<<deletedEdges[i].second;					
					}
					ofsResume<<std::endl;
					vector<pair<int,int> >().swap(deletedEdges);
					deletedEdges.clear();
				}
			}
			ofsResume.close();
			std::cout<<"write "<<output_resume_file<<" successfully!"<<std::endl;

			//************output delete order************
			string output_delete_order_file(output_dir_name);
			output_delete_order_file.append("_delete_order.txt");
			ofstream ofsDorder(output_delete_order_file);				
			if(!ofsDorder.is_open()) std::cerr<<output_delete_order_file<<" cannot be opened!"<<std::endl;
			for(int i=0;i<_delete_inv.size();++i){
				ofsDorder<<_delete_inv[i]<<std::endl;
			}
			ofsDorder.close();
			std::cout<<"write "<<output_delete_order_file<<" successfully!"<<std::endl;

			std::cout<<"******************generate and output overlay graph finished!************"<<std::endl;
		}

	protected:
	/**
		void count_add_edges(int curr_id,int& add_link_cnt,list<edge>& add_edges_list,list<edge>& update_edges_list){
			add_link_cnt=0;
			edge tmp;
			int u,v,w;
			bool isLink;
			int neighbor_size=es[curr_id].size();
			//clear
			add_edges_list.clear();
			update_edges_list.clear();
			//not process degree=1
			if(neighbor_size>1){
				vector<pair<int,int> > neighbors;//store the neighbor and its distance
				add_link_cnt=(neighbor_size*(neighbor_size-1))/2;
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
				for (int i  = 0; i < neighbor_size-1; i++){
					for (int j = i+1; j < neighbor_size; j++){
						isLink=false;
						u=neighbors[i].first;
						v=neighbors[j].first;
						w=neighbors[i].second+neighbors[j].second;
						tmp.u=u;
						tmp.v=v;
						tmp.weight=w;
						for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
							if(it->v==v){
								--add_link_cnt;
								if(w<it->weight){
									update_edges_list.push_back(tmp);
								}
								isLink=true;
								break;
							}
						}
						if(!isLink) add_edges_list.push_back(tmp);
					}
				}
			}
		}
		int deleteNodeByFreqMinAddEdges(){
			vector<list<edge> > add_edges_vec(numOfVertices);
			vector<list<edge> > update_edges_vec(numOfVertices);
			benchmark::heap<2, int, int> lqueue(numOfVertices);
			int u,v,w;//temp variables u-start v-end
			edge temp;//edge temp variable
			int d_cnt=0;//current num of deleted nodes
			bool isOne=false;
			int add_link_cnt;
			//initialize the add edges
			for(int id=0;id<numOfVertices;++id){
				if(HFPointFlag[id]) continue;
				count_add_edges(id,add_link_cnt,add_edges_vec[id],update_edges_vec[id]);
				lqueue.update(id,add_link_cnt);
			}
			//start delete
			while(!lqueue.empty()){
				if(d_cnt>=_numOfDeletedVertices) break;
				//extract min and store info
				int curr_id,curr_add;
				lqueue.extract_min(curr_id,curr_add);
				isDeleted[curr_id]=true;
				d_cnt++;
				_delete_inv.push_back(curr_id);
				//delete curr_id to neighbors edges
				u=es[curr_id].front().v;
				es[curr_id].clear();
				//delete neighbors to curr_id edges
				for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
					if(it->v==curr_id){
						es[u].erase(it);
						break;
					}
				}
				//add needed edges
				//update lqueue
				if(es[curr_id].size()==1){//degree=1
					u=es[curr_id].front().v;
					//delete curr_id and clear all edges
					es[curr_id].clear();
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==curr_id){
							es[u].erase(it);
							if(!HFPointFlag[u]) lqueue.update(u,es[u].size());
							//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
							break;
						}
					}
				}else{//degree>1

				}
				vector<pair<int,int> > neighbors;//store the neighbor and its distance
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));


				//process neighbor pair()
				for(int i=0;i<neighbors.size()-1; i++){
					u=neighbors[i].first;
					for(int j=i+1;j<neighbors.size();j++){
						v=neighbors[j].first;
					}
				}
			}
		}
		**/

		/**
		 * @Author: wanjingyi
		 * @description: delete node by degree and query frequency(i.e., HFpoint is not deleted)
		 * @param {*}
		 * @return {real numOfDeleteVertices}
		*/  
  		int deleteNodeByFreqDegree(int bfsHop_num,int searchThreshold){
			std::cout<<"****************************deleteNode by freq start!***************************"<<std::endl;
			//variables
			vector<vector<int> > dist_table(thresholdDegree,vector<int>(numOfVertices,INF_WEIGHT));//store the distances
			vector<bool> vis(numOfVertices,false);
			vector<int> hops(numOfVertices,0);
			benchmark::heap<2, int, int> pqueue(numOfVertices);
			int numOfLfpoint=numOfVertices-numOfHfpoint;
			//cout<<"numOfLfpoint = "<<numOfLfpoint<<std::endl;
			benchmark::heap<2, int, int> lqueue(numOfVertices);
			int u,v,w;//temp variables u-start v-end
			edge temp;//edge temp variable
			int d_cnt=0;//current num of deleted nodes
			bool isOne=false;
			//********************first stage delete LFPoint by degree********************
			//get the initial degree
			for(int i=0;i<numOfVertices;++i){
				//judge isolated nodes
				int deg=es[i].size();
				if(deg==0){
					cout<<"erroe: isolated vertex "<<i<<std::endl;
					continue;
				}
				if(HFPointFlag[i]) continue;
				lqueue.update(i,deg);
			}
			//cout<<"numOfLfpoint="<<numOfLfpoint<<" lqueue.size()="<<lqueue.size()<<std::endl;
			//start delete vertices
			while(!lqueue.empty()){
				if(d_cnt>=_numOfDeletedVertices) break;
				int curr_id,curr_deg;
				lqueue.extract_min(curr_id,curr_deg);
				//cout<<"*****************"<<curr_id<<"-"<<curr_deg<<"*****************"<<std::endl;
				if(curr_deg>thresholdDegree){
					cout<<"break  "<<d_cnt<<":"<<curr_id<<"-"<<curr_deg<<" "<<std::endl;
					break;
				}
				if(curr_deg==0){//means left only one node
					isOne=true;
					break;
				}
				if(curr_deg==1){//degree 1
					u=es[curr_id].front().v;
					es[curr_id].clear();
					isDeleted[curr_id]=true;
					d_cnt++;
					_delete_inv.push_back(curr_id);
					v=curr_id;
					if(es[u].size()==1) cout<<"graph error:line segment "<<u<<"-"<<v<<std::endl;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==v){
							es[u].erase(it);
							if(!HFPointFlag[u]) lqueue.update(u,es[u].size());
							//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
							break;
						}
					}
				}else{//degree >=2
					//deal with each neighbor pair
					vector<pair<int,int> > neighbors;//store the neighbor and its distance
					for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++){
						neighbors.push_back(make_pair(it->v,it->weight));
					}
					int neighbor_cnt=neighbors.size();
					vector<vector<bool> > isLink(neighbor_cnt,vector<bool>(neighbor_cnt,false));//prune
					int add_edges_cnt=((neighbor_cnt-1)*neighbor_cnt)/2;//prune
					//delete curr_id and clear all edges
					es[curr_id].clear();
					isDeleted[curr_id]=true;
					d_cnt++;
					_delete_inv.push_back(curr_id);
					//delete neighbors to curr_id edges
					for(int i=0;i<neighbor_cnt;++i){
						u=neighbors[i].first;
						for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
							if(it->v== curr_id){
								es[u].erase(it);
								break;
							}
						}
					}
					//update isLink neighbor pair
					for(int i=0;i<neighbor_cnt-1; i++){
						u=neighbors[i].first;
						for(int j=i+1;j<neighbor_cnt;j++){
							v=neighbors[j].first;
							w=neighbors[i].second+neighbors[j].second;
							for (list<edge>::iterator iter_u = es[u].begin(); iter_u != es[u].end(); iter_u++){
								//if exist edge(u,v),update the weight(u,v) to minimun
								if(iter_u->v==v){
									--add_edges_cnt;//prune
									isLink[i][j]=true;//prune
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
									break;
								}
							}
						}
					}
					//if(add_edges_cnt<=searchThreshold) cout<<"("<<curr_id<<"-"<<curr_deg<<")" <<" pair_total_cnt="<<((neighbor_cnt-1)*neighbor_cnt)/2<<" add_cnt="<<add_edges_cnt<<std::endl;
					if(add_edges_cnt<=searchThreshold){//not need to bfs_search
						//process each pair: u-v menas neighbor pair
						for(int i=0;i<neighbor_cnt-1; i++){
							u=neighbors[i].first;
							for(int j=i+1;j<neighbor_cnt;j++){
								v=neighbors[j].first;
								w=neighbors[i].second+neighbors[j].second;
								//cout<<"u="<<u<<" v="<<v<<" w="<<w<<" dist_table="<<dist_table_u[v]<<" hop="<<hops[v]<<std::endl;
								if(!isLink[i][j]){
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
						//update queue
						for(int i=0;i<neighbor_cnt;++i){
							u=neighbors[i].first;
							if(!HFPointFlag[u]) lqueue.update(u,es[u].size());
						}
					}else{//need bfs search
						//process each pair: u-v menas neighbor pair
						for(int i=0;i<neighbor_cnt-1; i++){
							u=neighbors[i].first;
							vector<int>& dist_table_u=dist_table[i];
							bfs_sp_search(u,dist_table_u,bfsHop_num,pqueue,hops,vis);
							for(int j=i+1;j<neighbor_cnt;j++){
								v=neighbors[j].first;
								w=neighbors[i].second+neighbors[j].second;
								//cout<<"u="<<u<<" v="<<v<<" w="<<w<<" dist_table="<<dist_table_u[v]<<" hop="<<hops[v]<<std::endl;
								if(!isLink[i][j]&&w<dist_table_u[v]){
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
						//update queue
						for(int i=0;i<neighbor_cnt;++i){
							u=neighbors[i].first;
							if(!HFPointFlag[u]) lqueue.update(u,es[u].size());
							//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
							for(int j=0;j<numOfVertices;++j) dist_table[i][j]=INF_WEIGHT;
						}
					}
				}
			}//while ends
			cout<<"d_cnt="<<d_cnt<<std::endl;
			//********************second stage********************
			if(isOne){//only left one node(degree=0)
				//check
				if(numOfVertices-d_cnt!=1) cout<<"error:numOfVertices-d_cnt!=1!"<<std::endl;
				cout<<"****************************deleteNode by freq end!***************************"<<std::endl;
				return 1;
			}
			//continue to process and delete HFPoint
			cout<<"****************************deleteNode by freq end!***************************"<<std::endl;
			return d_cnt;
		}

		/**
		 * @Author: wanjingyi
		 * @description: delete node by degree and query frequency(i.e., HFpoint is not deleted)
		 * @param {*}
		 * @return {real numOfDeleteVertices}
		*/  
  		int deleteNodeByFreqDegree(int bfsHop_num){
			std::cout<<"****************************deleteNode by freq start!***************************"<<std::endl;
			//variables
			vector<vector<int> > dist_table(thresholdDegree,vector<int>(numOfVertices,INF_WEIGHT));//store the distances
			vector<bool> vis(numOfVertices,false);
			vector<int> hops(numOfVertices,0);
			benchmark::heap<2, int, int> pqueue(numOfVertices);
			int numOfLfpoint=numOfVertices-numOfHfpoint;
			//cout<<"numOfLfpoint = "<<numOfLfpoint<<std::endl;
			benchmark::heap<2, int, int> lqueue(numOfVertices);
			int u,v,w;//temp variables u-start v-end
			edge temp;//edge temp variable
			int d_cnt=0;//current num of deleted nodes
			bool isOne=false;
			//********************first stage delete LFPoint by degree********************
			//get the initial degree
			for(int i=0;i<numOfVertices;++i){
				//judge isolated nodes
				int deg=es[i].size();
				if(deg==0){
					cout<<"erroe: isolated vertex "<<i<<std::endl;
					continue;
				}
				if(HFPointFlag[i]) continue;
				lqueue.update(i,deg);
			}
			//cout<<"numOfLfpoint="<<numOfLfpoint<<" lqueue.size()="<<lqueue.size()<<std::endl;
			//start delete vertices
			while(!lqueue.empty()){
				if(d_cnt>=_numOfDeletedVertices) break;
				int curr_id,curr_deg;
				lqueue.extract_min(curr_id,curr_deg);
				//cout<<"*****************"<<curr_id<<"-"<<curr_deg<<"*****************"<<std::endl;
				if(curr_deg>thresholdDegree){
					cout<<"break  "<<d_cnt<<":"<<curr_id<<"-"<<curr_deg<<" "<<std::endl;
					break;
				}
				if(curr_deg==0){//means left only one node
					isOne=true;
					break;
				}
				if(curr_deg==1){//degree 1
					u=es[curr_id].front().v;
					es[curr_id].clear();
					isDeleted[curr_id]=true;
					d_cnt++;
					_delete_inv.push_back(curr_id);
					v=curr_id;
					if(es[u].size()==1) cout<<"graph error:line segment "<<u<<"-"<<v<<std::endl;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==v){
							es[u].erase(it);
							if(!HFPointFlag[u]) lqueue.update(u,es[u].size());
							//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
							break;
						}
					}
				}else{//degree >=2
					//deal with each neighbor pair
					vector<pair<int,int> > neighbors;//store the neighbor and its distance
					for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++){
						neighbors.push_back(make_pair(it->v,it->weight));
					}
					int neighbor_cnt=neighbors.size();
					//int add_edges=((neighbor_cnt-1)*neighbor_cnt)/2;//to be deleted
					//delete curr_id and clear all edges
					es[curr_id].clear();
					isDeleted[curr_id]=true;
					d_cnt++;
					_delete_inv.push_back(curr_id);
					//delete neighbors to curr_id edges
					for(int i=0;i<neighbor_cnt;++i){
						u=neighbors[i].first;
						for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++){
							if(it->v== curr_id){
								es[u].erase(it);
								break;
							}
						}
					}
					//update isLink neighbor pair
					for(int i=0;i<neighbor_cnt-1; i++){
						u=neighbors[i].first;
						for(int j=i+1;j<neighbor_cnt;j++){
							v=neighbors[j].first;
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
									break;
								}
							}
						}
					}
					//process each pair: u-v menas neighbor pair
					for(int i=0;i<neighbor_cnt-1; i++){
						u=neighbors[i].first;
						vector<int>& dist_table_u=dist_table[i];
						bfs_sp_search(u,dist_table_u,bfsHop_num,pqueue,hops,vis);
						for(int j=i+1;j<neighbor_cnt;j++){
							v=neighbors[j].first;
							w=neighbors[i].second+neighbors[j].second;
							//cout<<"u="<<u<<" v="<<v<<" w="<<w<<" dist_table="<<dist_table_u[v]<<" hop="<<hops[v]<<std::endl;
							if(w<dist_table_u[v]){
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
					//update queue
					for(int i=0;i<neighbor_cnt;++i){
						u=neighbors[i].first;
						if(!HFPointFlag[u]) lqueue.update(u,es[u].size());
						//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
						for(int j=0;j<numOfVertices;++j) dist_table[i][j]=INF_WEIGHT;
					}
				}
			}//while ends
			cout<<"d_cnt="<<d_cnt<<std::endl;
			//********************second stage********************
			if(isOne){//only left one node(degree=0)
				//check
				if(numOfVertices-d_cnt!=1) cout<<"error:numOfVertices-d_cnt!=1!"<<std::endl;
				cout<<"****************************deleteNode by freq end!***************************"<<std::endl;
				return 1;
			}
			//continue to process and delete HFPoint
			cout<<"****************************deleteNode by freq end!***************************"<<std::endl;
			return d_cnt;
		}

		/**
		 * @Author: wanjingyi
		 * @description: delete node by degree and query frequency(i.e., HFpoint is not deleted)
		 * @param {*}
		 * @return {real numOfDeleteVertices}
		*/  
  		int deleteNodeByFreqDegree(){
			//std::cout<<"****************************deleteNode by freq start!***************************"<<std::endl;
			int numOfLfpoint=numOfVertices-numOfHfpoint;
			//cout<<"numOfLfpoint = "<<numOfLfpoint<<std::endl;
			benchmark::heap<2, int, int> lqueue(numOfVertices);
			int u,v,w;//temp variables u-start v-end
			edge temp;//edge temp variable
			int d_cnt=0;//current num of deleted nodes
			bool isOne=false;
			//********************first stage delete LFPoint by degree********************
			//get the initial degree
			for(int i=0;i<numOfVertices;++i){
				//judge isolated nodes
				int deg=es[i].size();
				if(deg==0){
					cout<<"erroe: isolated vertex "<<i<<std::endl;
					continue;
				}
				if(HFPointFlag[i]) continue;
				lqueue.update(i,deg);
			}
			//cout<<"numOfLfpoint="<<numOfLfpoint<<" lqueue.size()="<<lqueue.size()<<std::endl;
			//start delete vertices
			while(!lqueue.empty()){
				//cout<<"lqueue size="<<lqueue.size()<<std::endl;//to be deleted
				if(d_cnt>=_numOfDeletedVertices) break;
				int curr_id,curr_deg;
				lqueue.extract_min(curr_id,curr_deg);
				if(curr_deg>thresholdDegree){
					cout<<"break  "<<d_cnt<<":"<<curr_id<<"-"<<curr_deg<<" "<<std::endl;
					break;
				}
				if(curr_deg==0){//means left only one node
					isOne=true;
					break;
				}
				if(curr_deg==1){//degree 1
					u=es[curr_id].front().v;
					es[curr_id].clear();
					isDeleted[curr_id]=true;
					d_cnt++;
					_delete_inv.push_back(curr_id);
					v=curr_id;
					if(es[u].size()==1) cout<<"graph error:line segment "<<u<<"-"<<v<<std::endl;
					for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
						if(it->v==v){
							es[u].erase(it);
							if(!HFPointFlag[u]) lqueue.update(u,es[u].size());
							//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
							break;
						}
					}
				}else{//degree >=2
					//deal with each neighbor pair
					vector<pair<int,int> > neighbors;//store the neighbor and its distance
					bool isLink=false;//indicate whether exist edge(u,v)
					for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++) neighbors.push_back(make_pair(it->v,it->weight));
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
								if(!HFPointFlag[u]) lqueue.update(u,es[u].size());
								//if(lqueue.contains(u)) lqueue.update(u,es[u].size());
								break;
							}
						}
					}
					//delete curr_id and clear all edges
					es[curr_id].clear();
					isDeleted[curr_id]=true;
					d_cnt++;
					_delete_inv.push_back(curr_id);
				}
			}//while ends
			cout<<"d_cnt="<<d_cnt<<std::endl;
			//********************second stage********************
			if(isOne){//only left one node(degree=0)
				//check
				if(numOfVertices-d_cnt!=1) cout<<"error:numOfVertices-d_cnt!=1!"<<std::endl;
				cout<<"****************************deleteNode by freq end!***************************"<<std::endl;
				return 1;
			}
			//continue to process and delete HFPoint
			cout<<"****************************deleteNode by freq end!***************************"<<std::endl;
			return d_cnt;
		}

		/*
		 *@description: initalize all variables
		 *@author: wanjingyi
		 *@date: 2021-01-21
		*/
		void initIterms(){
			std::cout<<"******************init varibales start!**************"<<std::endl;
			//initialize all variables
			isDeleted.resize(numOfVertices,0);
			isDeletedLFNode.resize(numOfVertices,0);
			level.resize(numOfVertices);
			HFPoint.reserve(numOfVertices);
			removedEdge.resize(numOfVertices);
			HFPointFlag.resize(numOfVertices,0);
			vector<bool> HFPointFlag(numOfVertices, 0);//indicate whether is HFpoint
			std::cout<<"******************init varibales finished!**************"<<std::endl;
		}
		/*
		*@description: _freq_rank and _freq_inv
		*@author: wanjingyi
		*@date: 2021-01-21
		*/
		void rankVerticesByFreq(){
			_freq_rank.resize(numOfVertices);
			_freq_inv.resize(numOfVertices);
			int r;
			for(r=0;r<HFPoint.size();++r){
				_freq_inv[r]=HFPoint[r];
				_freq_rank[HFPoint[r]]=r;
			}
			for(int v=0;v<numOfVertices;++v){
				if(HFPointFlag[v]) continue;
				_freq_inv[r]=v;
				_freq_rank[v]=r;
				++r;
			}
			if(r!=numOfVertices) std::cerr<<"order freq rank error!"<<std::endl;
		}

		/*
		 *@description: delete nodes with degree 1
		 *@author: gaoyongyong
		 *@return; current deleted node' neighbour
		 *@date: 2021-01-21
		*/
		int degreeOneDelete(vector<list<edge>>& es, vector<list<edge>>& removedEdge, vector<bool>& isDeleted, const int number) {
			int curNum = number;
			edge temp;
			temp.u = es[curNum].front().u;
			temp.v = es[curNum].front().v;
			temp.weight = es[curNum].front().weight;
			removedEdge[curNum].push_back(temp);//store the removed edge
			isDeleted[curNum] = 1;
			es[curNum].clear();
			curNum = temp.v;//deal the neighbor node's removedEdge
			for (list<edge>::iterator it = es[curNum].begin(); it != es[curNum].end(); it++) {
				if (it->v == temp.u) {
					es[curNum].erase(it);
					if (es[curNum].size() == 0)
						isDeleted[curNum] = 1;
					break;
				}
			}
			if (es[curNum].size() == 1)//if (u,v) has been removed, v's degree is 1 and continue to delete v
			degreeOneDelete(es, removedEdge, isDeleted, curNum);
			return curNum;
		}

	/*
	 *@description: delete vertices with degree 1 just in LFPoint
	 *@author: wanjingyi
	 *@return; current deleted node' neighbour
	 *@date: 2021-01-21
	*/
	int degreeOneDelete(vector<list<edge>>& es, vector<list<edge>>& removedEdge, vector<bool>& isDeleted, const int number, const vector<bool>& isDeleteLFPoint) {
		int curNum = number;
		edge temp;
		temp.u = es[curNum].front().u;
		temp.v = es[curNum].front().v;
		temp.weight = es[curNum].front().weight;
		removedEdge[curNum].push_back(temp);//store the removed edge
		isDeleted[curNum] = 1;
		es[curNum].clear();
		curNum = temp.v;//deal the neighbor node's removedEdge
		for (list<edge>::iterator it = es[curNum].begin(); it != es[curNum].end(); it++) {
			if (it->v == temp.u) {
				es[curNum].erase(it);
				if (es[curNum].size() == 0)
					isDeleted[curNum] = 1;
				break;
			}
		}
		if (es[curNum].size() == 1 && isDeleteLFPoint[curNum]) {//if (u,v) has been removed, v's degree is 1 and continue to delete v
			degreeOneDelete(es, removedEdge, isDeleted, curNum, isDeleteLFPoint);
		}
		return curNum;
	}

	/*
	 *@description: delete node with degree>=2
	 *@author: gaoyongyong
	 *@date: 2021-01-21
	*/
	vector<int> degreeDelete(vector<list<edge>>& es, vector<list<edge>>& removedEdge, vector<bool>& isDeleted, const int number, const int degree) {
		int edgeNum = degree * (degree - 1) / 2;//删除度为degree的点需要添加degree*(degree-1)/2条边，保证删除后图的正确性
		vector<edge> edgeProcess;//用于保存删除结点后，增加的边信息
		vector<edge> originalEdge(degree);//用于保存删除结点到其邻接点的所有边信息
		vector<int> neighbour(degree);//store the neighbors of curr deleted node
		edge temp;
		int u, v, w;
		int curNum = number;
		int flag = 0;
		int j = 0;
		//store the added edges after deletion
		for (list<edge>::iterator it = es[curNum].begin(); it != es[curNum].end(); it++, j++) {
			temp.u = it->u;
			temp.v = it->v;
			temp.weight = it->weight;
			removedEdge[curNum].push_back(temp);
			neighbour[j] = it->v;
			originalEdge[j].u = curNum;
			originalEdge[j].v = it->v;
			originalEdge[j].weight = it->weight;
		}
		//for every pair of neighbors(u,v) ,add a new edge
		for (int i = 0; i < degree; i++) {
			temp.u = originalEdge[i].v;
			temp.v = originalEdge[(i + 1) % degree].v;
			temp.weight = originalEdge[i].weight + originalEdge[(i + 1) % degree].weight;
			edgeProcess.push_back(temp);
		}
		//for every pair of nonneighbors(u,v),add a new edge
		for (int i = 0; i < degree - 2; i++) {
			int cease{ degree };//initial cases:the end position
			if (i == 0)//for the first point, only traverse end at the last second one
				cease = degree - 1;
			for (int j = i + 2; j < cease; j++) {//from the start point +2
				temp.u = originalEdge[i].v;
				temp.v = originalEdge[j].v;
				temp.weight = originalEdge[i].weight + originalEdge[j].weight;
				edgeProcess.push_back(temp);
			}
		}
		isDeleted[curNum] = 1;
		es[curNum].clear();
		// if(edgeNum!=edgeProcess.size()) //modified by wanjingyi
		// 	cout<<"erroer:edgeNum!=edgeProcess.size()!"<<std::endl;
		for (int i = 0; i < edgeNum; i++) {
			edgeDelete(es, removedEdge, isDeleted, edgeProcess[i].u, curNum);//删除curNum邻接点的边（u,curNum）
		}
		//process the degree with 2 and delete currNum's neighbors' edges（u,curNum）
		if (degree == 2) {
			edgeDelete(es, removedEdge, isDeleted, edgeProcess[1].u, curNum);
		}

		for (int i = 0; i < edgeNum; i++) {
			u = edgeProcess[i].u;
			v = edgeProcess[i].v;
			w = edgeProcess[i].weight;
			flag = 0;
			for (list<edge>::iterator it = removedEdge[u].begin(); it != removedEdge[u].end(); it++) {
				if (it->v == v) {
					if (it->weight < w) { 
						w = it->weight;
					}
				}
			}
			for (list<edge>::iterator it = removedEdge[v].begin(); it != removedEdge[v].end(); it++) {
				if (it->v == u) {
					if (it->weight < w) { 
						w = it->weight;
					}
				}
			}
			for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {
				if (it->v == v) {
					if (it->weight < w) {
						flag = 1;
						//diagonal.push_back(*it);
					}
					else
						es[u].erase(it);
					break;
				}
			}
			for (list<edge>::iterator it = es[v].begin(); it != es[v].end(); it++) {
				if (it->v == u) {
					if (it->weight < w) { 
						flag = 1;
						//diagonal.push_back(*it);
					}
					else
						es[v].erase(it);
					break;
				}
			}
			if (0 == flag) {
				temp.u = u;
				temp.v = v;
				temp.weight = w;
				es[u].push_back(temp);
				temp.u = v;
				temp.v = u;
				es[v].push_back(temp);
			}
		}

		//if the deleted neighbor's degree is 1 and process it
		for (int i = 0; i < degree; i++) {
			u = neighbour[i];
			if (es[u].size() == 1)
				degreeOneDelete(es, removedEdge, isDeleted, u);
		}
		
		//continue to process the deleted neighbor by ascending degree order
		for (int i = 2; i <= degree; i++) {
			for (int j = 0; j < degree; j++) {
				u = neighbour[j];
				if (es[u].size() == i) {
					vector<int> num = degreeDelete(es, removedEdge, isDeleted, u, i);
					for (int k = 0; k < num.size(); k++) {
						if (i < degree && es[num[k]].size() == i + 1)
							degreeDelete(es, removedEdge, isDeleted, num[k], i + 1);
					}
				}
			}
		}
		return neighbour;

	}

	/*
	 *@description: delete edges
	 *@author: gaoyongyong
	 *@date: 2021-01-21
	*/
	void edgeDelete(vector<list<edge>>& es, vector<list<edge>>& removedEdge, vector<bool>& isDeleted, const int u, const int curNum) {
		if (curNum == 4359||u==4359) {
			int num=curNum;
		}
		for (list<edge>::iterator it = es[u].begin(); it != es[u].end(); it++) {//删除结点u中边(u,curNum)
			if (it->v == curNum) {
				edge temp;
				temp.u = it->u;
				temp.v = it->v;
				temp.weight = it->weight;
				removedEdge[u].push_back(temp);
				es[u].erase(it);
				if (es[u].size() == 1)
					isDeleted[u] = 1;
				break;
			}
		}
	}

	/**
	 * @description: delete node by degree for directed graphs
	 * @param {*}
	 * @return {*}
	 * @author: Wan Jingyi
	 */ 
 	void deleteNodeByDegree_directed(){
		std::cout<<"******************directed delete node by degree start!************"<<std::endl;
		dtime=GetCurrentTimeSec();
		//initialize variables
		isDeleted.resize(numOfVertices,false);
		_delete_rank.resize(numOfVertices);
		_delete_inv.reserve(numOfVertices);
		numOfOriginalVertices=numOfVertices;
		numOfOverlayVertices=0;
		numOfDeletedVertices=0;
		_numOfEdges=0;
		_numOfEdges_r=0;
		_numOfEdges_f=0;
		benchmark::heap<2, int, int> dqueue(numOfVertices);
		int u,v,u_w,v_w,w;//temp variables u-start v-end
		edge temp;//edge temp variable
		int d_cnt=0;//current num of deleted nodes
		int curr_id,curr_deg,deg;
		bool isOne=false;

		//get the initial degree
		for(int i=0;i<numOfVertices;++i){
			if(es[i].size()==0&&r_es[i].size()==0){
				cout<<"erroe: isolated vertex "<<i<<std::endl;
				continue;
			}
			deg=es[i].size()*r_es[i].size();
			if(deg==1&&es[i].front().v==r_es[i].front().u) deg=0;
			//to be deleted
			std::cout<<i<<":es_size="<<es[i].size()<<" r_es_size="<<r_es[i].size()<<" deg="<<deg<<std::endl;
			dqueue.update(i,deg);
		}

		//start delete vertices
		while(!dqueue.empty()){
			dqueue.extract_min(curr_id,curr_deg);
			//to be deleted
			//std::cout<<d_cnt<<":"<<curr_id<<"-"<<curr_deg<<std::endl;
			if(curr_deg>thresholdDegree){
				cout<<"break  "<<d_cnt<<":"<<curr_id<<"-"<<curr_deg<<" "<<std::endl;
				break;
			}
			isDeleted[curr_id]=true;
			_delete_inv.push_back(curr_id);
			_delete_rank[curr_id]=d_cnt++;
			if(curr_deg==0){//line
				if(es[curr_id].size()==0&&r_es[curr_id].size()==0){//only left one vertice
					isOne=true;
				}
				else if(es[curr_id].size()==0&&r_es[curr_id].size()!=0){//all edges are in backward
					for(list<edge>::iterator it=r_es[curr_id].begin();it!=r_es[curr_id].end();++it){
						u=r_es[curr_id].front().u;
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
					u=r_es[curr_id].front().u;
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
			}else{//need to deal with each neighbor pair
				vector<pair<int,int> > neighbors;//store the forward neighbor and its distance
				vector<pair<int,int> > r_neighbors;//store the backward neighbor and its distance
				bool isLink=false;//indicate whether exist edge(u,v)
				for (list<edge>::iterator it = es[curr_id].begin(); it != es[curr_id].end(); it++){
					neighbors.push_back(make_pair(it->v,it->weight));
				}
				for (list<edge>::iterator it = r_es[curr_id].begin(); it != r_es[curr_id].end(); it++){
					r_neighbors.push_back(make_pair(it->u,it->weight));
				}
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
		}//while ends

		dtime=GetCurrentTimeSec()-dtime;
		numOfDeletedVertices=d_cnt;
		numOfOverlayVertices=numOfOriginalVertices-numOfDeletedVertices;
		std::cout<<"numOfDeletedVertices = "<<numOfDeletedVertices<<std::endl;
		std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<std::endl;
		std::cout<<"******************directed delete node by degree finished!************"<<std::endl;
		return;

	}

};

#endif