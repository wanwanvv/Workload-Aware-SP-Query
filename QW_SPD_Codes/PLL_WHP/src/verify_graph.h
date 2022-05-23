/*
 * @Descripttion: 
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-09-21 09:46:24
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-02-26 08:28:03
 */
#pragma once
#ifndef VERIFY_GRAPH_H
#define VERIFY_GRAPH_H

#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include "./graph.h"
#include "./paras.h"
#include "./heap.h"
#include "./utils.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define INF_WEIGHT SP_Constants::INF_WEIGHT
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
using namespace std;

/**
 * @Author: wanjingyi
 * @description: some functions for search
 */
class Search_graph{
    public:
        //variables
        typedef vector<NodeEdgeWeightPair> path;

        //constructions and deconstructions

        //******functions******

        /*
            *@description: compute all SP distances with omp
            *@author: wanjingyi
            *@date: 2020-12-24
        */
        void compute_all_SP_distances(WGraph &wgraph, vector<vector<EdgeWeight> >& SP_distances,NodeID n,int _num_threads=10){
            cout<<"********************compute all SP distances start!*******************"<<endl;
            int num_threads = _num_threads;
            cout<<"num_threads="<<num_threads<<endl;
            SP_distances.resize(n,vector<EdgeWeight>(n,INF_WEIGHT));
            //variables
            vector<benchmark::heap<2, EdgeWeight, NodeID> > pqueue(num_threads, benchmark::heap<2, EdgeWeight, NodeID>(n) );
            vector<vector<bool> > vis(num_threads, vector<bool>(n,false)); //vis for omp
            vector<vector<EdgeWeight> > distances(num_threads, vector<EdgeWeight>(n, INF_WEIGHT));
            cout<<"init iterms ......."<<endl;
            omp_set_num_threads(num_threads);//设置线程数
            #pragma omp parallel for schedule(dynamic)
                for (NodeID v = n - 1; v >=0; --v) {
                    //cout<<v<<" thread="<<omp_get_thread_num()<<endl;//to be deleted
                    source_Dijkstra(v,wgraph,vis[omp_get_thread_num()],pqueue[omp_get_thread_num()],distances[omp_get_thread_num()],SP_distances,n);
                }
            cout<<"********************compute all SP distances finished!*******************"<<endl;
        }

        /*
            *@description: compute all SP paths with omp
            *@author: wanjingyi
            *@date: 2020-12-24
        */
        void compute_all_SP_paths(WGraph &wgraph,vector<vector<EdgeWeight> >& SP_distances,vector<path>& SP_paths,vector<vector<int> >& SP_index,int& _numOfSPs,NodeID n,int _num_threads=10){
            cout<<"********************compute all SP paths start!*******************"<<endl;
            int num_threads = _num_threads;
            cout<<"num_threads="<<num_threads<<endl;
             _numOfSPs=0;
            SP_distances.resize(n,vector<EdgeWeight>(n,INF_WEIGHT));
            SP_paths.reserve(n*n);//reserve n*n paths
            SP_index.resize(n,vector<int>(n,-1));//-1 means no path
            //variables
            vector<benchmark::heap<2, EdgeWeight, NodeID> > pqueue(num_threads, benchmark::heap<2, EdgeWeight, NodeID>(n) );
            vector<vector<bool> > vis(num_threads, vector<bool>(n)); //vis for omp
            vector<vector<EdgeWeight> > distances(num_threads, vector<EdgeWeight>(n, INF_WEIGHT));
            vector<vector<path> >  paths(num_threads, vector<path>(n));
            omp_set_num_threads(num_threads);//设置线程数
            #pragma omp parallel for schedule(dynamic)
            for (NodeID v = n - 1; v >=0; --v) {
                source_Dijkstra(v,wgraph,vis[omp_get_thread_num()],pqueue[omp_get_thread_num()],distances[omp_get_thread_num()],paths[omp_get_thread_num()],SP_distances,SP_paths,SP_index);
            }
            _numOfSPs=SP_paths.size();
            std::cout<<"numofSPs = "<<_numOfSPs<<endl;
            cout<<"********************compute all SP paths finished!*******************"<<endl;
        }


        /**
         * @Author: wanjingyi
         * @description: dijkstra algorithm for SP distance
         * @param {graph class,start node index...}
         * @return {void}
         */
        void source_Dijkstra(NodeID source,WGraph &wgraph, vector<bool>& vis, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, vector<EdgeWeight>& distances,vector<vector<EdgeWeight> >& SP_distances,int n){
            pqueue.update(source,0);
            distances[source]=0;
            //paths[source].push_back(make_pair(source,0));
            NodeID u,v; EdgeWeight u_w,v_w,v_d;
            while(!pqueue.empty())
            {
                pqueue.extract_min(u,u_w);
                vis[u]=true;
                for (EdgeID eid = wgraph.vertices[u]; eid < wgraph.vertices[u + 1]; ++eid) {
                    v = wgraph.edges[eid].first;
                    v_w=wgraph.edges[eid].second;
                    v_d=distances[u]+v_w;
                    if(!vis[v]){
                        if(v_d<distances[v]){
                            distances[v]=v_d;
                            pqueue.update(v,v_d);
                        }
                    }
                }
            }
            vector<EdgeWeight> & SP_distances_source=SP_distances[source];
            //clear tmp list
            for(NodeID i=0;i<n;++i){
                SP_distances_source[i]=distances[i];
                vis[i]=false;
                distances[i]=INF_WEIGHT;
            }
            pqueue.clear();
        }

        /*
        *@description: Dijkstra algorithm for shortest paths
        *@author: wanjingyi
        *@date: 2020-12-24
        */
        void source_Dijkstra(NodeID source,WGraph &wgraph, vector<bool>& vis, benchmark::heap<2, EdgeWeight, NodeID>& pqueue, vector<EdgeWeight>& distances, vector<path>& paths,vector<vector<EdgeWeight> >& SP_distances,vector<path>& SP_paths,vector<vector<int> >& SP_index){
            pqueue.update(source,0);
            distances[source]=0;
            //paths[source].push_back(make_pair(source,0));
            NodeID u,v; EdgeWeight u_w,v_w,v_d;
            while(!pqueue.empty())
            {
                pqueue.extract_min(u,u_w);
                vis[u]=true;
                for (EdgeID eid = wgraph.vertices[u]; eid < wgraph.vertices[u + 1]; ++eid) {
                    v = wgraph.edges[eid].first;
                    v_w=wgraph.edges[eid].second;
                    v_d=distances[u]+v_w;
                    if(!vis[v]){
                        if(v_d<distances[v]){
                            distances[v]=v_d;
                            pqueue.update(v,v_d);
                            //store the path
                            if(u==source){
                                paths[v].push_back(make_pair(v,v_w));
                            }else{
                                path prevPath(paths[u]);
                                prevPath.push_back(make_pair(v,v_w));
                                paths[v]=prevPath;
                            }
                        }
                    }
                }
            }
            vector<EdgeWeight> & SP_distances_source=SP_distances[source];
            vector<int> &SP_index_source=SP_index[source];
            //clear tmp list
            for(NodeID i=0;i<numOfVertices;++i){
                if(!paths[i].empty()){
                    SP_paths.push_back(paths[i]);
                    SP_index_source[i]=SP_paths.size()-1;
                }
                SP_distances_source[i]=distances[i];
                vis[i]=false;
                distances[i]=INF_WEIGHT;
                paths[i].clear();
            }
            pqueue.clear();
        }
    
    protected:
        
};

/**
 * @Author: wanjingyi
 * @description: check shortest paths rightness
 */
class Compute_graph:public Search_graph{
public:
    //****************constructions************
    Compute_graph(WGraph& wgraph, char* outputDirName,int generatete_model,int num_of_sets,int num_of_data,int num_of_threads){
        if(generatete_model==0){

        }
    }

protected:
    void generate_different_distances_data_all(WGraph& wgraph,char* outputDirName,int k,int n,int num_of_threads){
        //comute all pair sp distances
        vector<vector<EdgeWeight> > SP_distances(numOfVertices,vector<EdgeWeight>(numOfVertices,INF_WEIGHT));
        compute_all_SP_distances(wgraph, SP_distances,numOfVertices,num_of_threads);
        //compute min an max distance
        EdgeWeight min_distance=INF_WEIGHT;
        EdgeWeight max_distance=0,curr_distance;
        NodeID s,t;
        for(s=0;s<numOfVertices;++s){
            if(DIRECTED_FLAG==true) t=s;
            else t=0;
            curr_distance=SP_distances[s][t];
            if(curr_distance!=0&&curr_distance!=INF_WEIGHT){
                min_distance=min(min_distance,curr_distance);
                max_distance=max(max_distance,curr_distance);
            }
        }
        double x=pow(10.0,log10(double(max_distance))/(double)k);
        cout.setf(ios::fixed);
        cout.precision(4);
        cout<<" min_distance="<<min_distance<<" max_distance="<<max_distance<<" x="<<x<<endl;
        //to be deleted
        double x1=pow(10.0,log10(double(1024))/(double)k);
        cout<<"x1="<<x1<<endl;
        //compute k range
        double cur=1.0;
        vector<double> range{cur}; 
        for(int i=0;i<k;++i){
            cur*=x;
            range.push_back(cur);
            cout<<cur<<" ";
        }
        cout<<endl;//to be deleted
        //classify the queries
        vector<vector<pair<NodeID,NodeID>>> query_data(k);
        for(s=0;s<numOfVertices;++s){
            if(DIRECTED_FLAG=false) t=s;
            else t=0;
            for(;t<numOfVertices;++t){
                curr_distance=SP_distances[s][t];
                if(curr_distance!=0&&curr_distance!=INF_WEIGHT){
                    int i=static_cast<int>(log(double(curr_distance))/(double)log(x));//floor
                    query_data[i].push_back(make_pair(s,t));
                }
            }
        }
        string output_prefix(outputDirName);
        for(int i=0;i<k;++i){
            string output_file=output_prefix+to_string(i)+".txt";
            ofstream ofs(output_file);
            if(!ofs.is_open()) cout<<"Cannot open "<<output_file<<endl;
            for(size_t j=0;j<query_data[i].size();++j){
                ofs<<query_data[i][j].first<<" "<<query_data[i][j].second<<endl;
            }
        }
    }

    void generate_different_distances_data(WGraph& wgraph,char* outputDirName,int k,int n,int num_of_threads){

    }

    pair<EdgeWeight,EdgeWeight> compute_max_and_min_distances(WGraph& wgraph){

    }

    void classify_st_pairs_by_distances(WGraph& wgraph,double x,vector<vector<pair<NodeID,NodeID>>>& query_data){

    }
};

#endif // !VERIFY_GRAPH_H

