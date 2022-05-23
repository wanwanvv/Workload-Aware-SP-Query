/*
 * @Description: 
 * @Author: wanjingyi
 * @Date: 2021-01-26 15:24:50
 * @LastEditTime: 2021-10-22 09:37:10
 */
#pragma once
#ifndef FREQUENCY_HIERARCHY_ORDERING_H
#define FREQUENCY_HIERARCHY_ORDERING_H

#include <algorithm>
#include <unordered_set>
#include <unordered_map>  
#include <cmath>
#include <time.h>
#include "./graph.h"
#include "./labels.h"
#include "./time_util.h" 
#include "./heap.h"
#include "./paras.h"
#include "./ordering.h"
#include "./construction.h"
#include "./utils.h"

#ifdef _WIN32
	#include<google/sparsehash/sparseconfig.h>
#endif
    #include<google/dense_hash_map>
	#include<google/dense_hash_set>

#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define MAX_VALUE 0x7FFFFFFF

bool cmp1(const pair<NodeID,int> a, const pair<NodeID,int> b) {
    return a.second>b.second;//自定义的比较函数
}

struct cmp_queue
    {template <typename T, typename U>
        bool operator()(T const &left, U const &right)
        {
        // 以 second 比较。输出结果为 Second 较大的在前 Second 相同时，先进入队列的元素在前。
            if (left.second < right.second)
                return true;
            return false;
        }
    };

/*
 *@description: the basis class of other synthesis ordering ways
 *@author: wanjingyi
 *@date: 2020-12-07
*/
class FOrdering:public Ordering{
    public:
        //*******************variables*******************
        NodeID numOfHFpoint;//High frequency point num
        vector<bool> HFPointIndex;//whether the point is high freq point 
        NodeID numOfOrderVertices;//num of vertices need to order
        //*******************constructions and deconstructions*******************
        FOrdering(){}
        ~FOrdering(){}

        //*******************functions*******************
        void save_rank(const char* order_dir_name) {
            string order_file(order_dir_name);
            order_file.append(".order");
            ofstream ofs(order_file);
            if(!ofs.is_open()) cout<<order_file<<" cannot be opend!"<<endl;
            for (int i = 0; i < numOfVertices; ++i) {
                ofs << inv[i] << endl;
            }
            ofs.close();
	    }

};

/*
 *@description: class used to order nodes with linear itemrs*coefficeient = weight sum
 *@author: wanjingyi
 *@date: 2020-12-13
*/
template<typename weightType>
class Linear_Ordering :public FOrdering{
    public:
        Processing::calcCoefficient<double> _calcCoef;
        HFLabel labels;
        vector<unsigned int> _freq;//get the query time of each node by index
        vector<unsigned int> _betwenness;//get the betwenness by node index
        vector<unsigned int> _coverage;//get the coverage by node index
        vector<unsigned int> _depth;//get the depth by node index
        vector<unsigned int> _degree;//get the degree by node index
        vector<weightType> _degree_weight;//degree tie breaking
        vector<weightType> _betwenness_weight;//betwenness tie breaking
        vector<weightType> _freq_weight;//freq tie breaking
        vector<NodeID> _freq_inv;//fetch the index by rank
        vector<NodeID> _freq_rank;//fetch the rank  by index
        vector<NodeID> _betwenness_inv;//fetch the index by rank
        vector<NodeID> _betwenness_rank;//fetch the rank  by index
        vector<NodeID> _coverage_inv;//fetch the index by rank
        vector<NodeID> _coverage_rank;//fetch the rank  by index
        vector<NodeID> _depth_inv;//fetch the index by rank
        vector<NodeID> _depth_rank;//fetch the rank  by index
        vector<NodeID> _degree_inv;//fetch the index by rank
        vector<NodeID> _degree_rank;//fetch the rank  by index
        long long total_sum_size=0,hf_sum_size=0;//total size variables
	    double total_ave_size=0,hf_ave_size=0;//average size variables
        ofstream debug_out;//output debug information
        unsigned int  max_freq=0;//max query time used to normalization
        unsigned int max_degree=0;//max degree of all nodes
        unsigned int max_coverage=0;//max coverage of all nodes
        unsigned int max_depth=0;//max depth of all nodes
        unsigned int max_betweennss=0;//max betweenness of all nodes
        unsigned int  min_freq=MAX_VALUE;//min query time used to normalization
        unsigned int min_degree=MAX_VALUE;//min degree of all nodes
        unsigned int min_coverage=MAX_VALUE;//min coverage of all nodes
        unsigned int min_depth=MAX_VALUE;//min depth of all nodes
        unsigned int min_betweenness=MAX_VALUE;//min betweenness of all nodes
        unsigned int  interval_freq=0;//max query time -min query time
        unsigned int interval_degree=0;//max_degree -min_degree
        unsigned int interval_coverage=0;//max_coverage-min_coverage
        unsigned int interval_depth=0;//max_depth-min_depth
        unsigned int interval_betweenness=0;//max_depth-min_depth
        vector<bool> isConstruct;//indicates whether the node has been constructed
        //***************************construction functions********************

        Linear_Ordering(){
        }
        ~Linear_Ordering(){

        }


        /*
        *@description: construction overlay graph used for undirected and weighted graph
        *@author: wanjingyi
        *@date: 2020-12-13
        */

        Linear_Ordering(WGraph& wgraph,Processing::calcCoefficient<double> calcCoef,const vector<unsigned int>& freq,const vector<unsigned int>& bet,const vector<unsigned int>& cov,int hfRate){
            //std::cout<<"Linear_Ordering undirected weighted graph...."<<endl;
            if(DIRECTED_FLAG==false) undirected_weighted_orderByValues_overlay(wgraph,calcCoef,freq,bet,cov,hfRate);
            else directed_weighted_orderByValues_overlay(wgraph,calcCoef,freq,bet,cov,hfRate);
            Relabel(wgraph);
        }

    protected: 
        /*
         *@description: order the nodes by detailed_values*Coefficient sum for overlay graph
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        // void undirected_weighted_orderByValues_overlay(WGraph& wgraph,Processing::calcCoefficient<double> calcCoef,const vector<unsigned int>& freq,const vector<unsigned int>& bet,const vector<unsigned int>& cov,int hfRate){
        //    initIterms_overlay();
        //     _calcCoef=calcCoef;//get the coefficient from command line
        //     if(_calcCoef.is_deg_mult) load_degree_overlay(wgraph);//degree
        //     if(_calcCoef.is_freq_mult) load_HFpoint_overlay(freq);
        //     if(_calcCoef.is_bet_mult) load_betwenness_overlay(bet);
        //     if(_calcCoef.is_cov_mult) load_coverage_overlay(cov);
        //     //to be deleted 
        //     //cout<<"_calcCoef.bet_mult="<<_calcCoef.bet_mult<<" _calcCoef.freq_mult="<<_calcCoef.freq_mult<<" _calcCoef.deg_mult="<<_calcCoef.deg_mult<<endl;
        //     calcWeightByValue();//cumpute the order
        // }
        void undirected_weighted_orderByValues_overlay(WGraph& wgraph,Processing::calcCoefficient<double> calcCoef,const vector<unsigned int>& freq,const vector<unsigned int>& bet,const vector<unsigned int>& cov,int hfRate){
           initIterms_overlay();
            _calcCoef=calcCoef;//get the coefficient from command line
            if(_calcCoef.is_freq_mult) load_HFpoint_overlay(freq);
            if(_calcCoef.is_bet_mult) load_betwenness_overlay(bet);
            if(_calcCoef.is_cov_mult) load_coverage_overlay(cov);
            //to be deleted 
            //cout<<"_calcCoef.bet_mult="<<_calcCoef.bet_mult<<" _calcCoef.freq_mult="<<_calcCoef.freq_mult<<" _calcCoef.deg_mult="<<_calcCoef.deg_mult<<endl;
            if(_calcCoef.is_deg_mult) calcWeightByValue_degree(wgraph);
            else calcWeightByValue();//cumpute the order
        }

        //Add for directed
        /*
         *@description: order the nodes by detailed_values*Coefficient sum for overlay graph
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void directed_weighted_orderByValues_overlay(WGraph& wgraph,Processing::calcCoefficient<double> calcCoef,const vector<unsigned int>& freq,const vector<unsigned int>& bet,const vector<unsigned int>& cov,int hfRate){
           initIterms_overlay();
            _calcCoef=calcCoef;//get the coefficient from command line
            if(_calcCoef.is_freq_mult) load_HFpoint_overlay(freq);
            if(_calcCoef.is_bet_mult) load_betwenness_overlay(bet);
            if(_calcCoef.is_cov_mult) load_coverage_overlay(cov);
            //to be deleted 
            //cout<<"_calcCoef.bet_mult="<<_calcCoef.bet_mult<<" _calcCoef.freq_mult="<<_calcCoef.freq_mult<<" _calcCoef.deg_mult="<<_calcCoef.deg_mult<<endl;
            if(_calcCoef.is_deg_mult) calcWeightByValue_degree_directed(wgraph);
            else calcWeightByValue();//cumpute the order
        }

        /*
         *@description: calculate the node weight by value
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void calcWeightByValue(){ 
            std::cout<<"*****************************calcWeightByValue begins!**************************"<<endl;
            vector<pair<weightType,NodeID> > orderWeight;
            srand(0x728f9d0);//srand seed
            orderWeight.reserve(numOfVertices);
            //std::cout<<" freq_mult="<<_calcCoef.freq_mult<<" deg_mult="<<_calcCoef.deg_mult<<" bet_mult="<<_calcCoef.bet_mult<<" cov_mult="<<_calcCoef.cov_mult<<" dep_mult="<<_calcCoef.dep_mult;
            for(NodeID v=0;v<numOfVertices;++v){
                weightType result=0;
                // result=(weightType)_freq[v]*_calcCoef.freq_mult+(weightType)_degree[v]*_calcCoef.deg_mult+
                // (weightType)_betwenness[v]*_calcCoef.bet_mult+(weightType)_coverage[v]*_calcCoef.cov_mult+
                // (weightType)_depth[v]*_calcCoef.dep_mult;//tie breaking +weightType(rand()) / RAND_MAX
                if(_calcCoef.is_deg_mult) result+=((weightType)_degree_weight[v]/(weightType)interval_degree)*_calcCoef.deg_mult;
                if(_calcCoef.is_freq_mult) result+=((weightType)_freq_weight[v]/(weightType)interval_freq)*_calcCoef.freq_mult;
                if(_calcCoef.is_bet_mult) result+=((weightType)_betwenness_weight[v]/(weightType)interval_betweenness)*_calcCoef.bet_mult;
                orderWeight.push_back(make_pair(result,v) );
            }
            sort(orderWeight.rbegin(),orderWeight.rend());
            //output orderWeight
            // string write_filename_prefix(orderWeightFileName);
			// string orderWeight_filename=write_filename_prefix.append(".orderWeight");
            // ofstream out(orderWeight_filename.c_str());//to be deleted
            // std::cout<<"orderWeight_filename="<<orderWeight_filename<<endl;
            for(NodeID i=0;i<numOfVertices;++i){
                NodeID v=orderWeight[i].second;
                inv[i]=v;
                rank[v]=i;
                //out<<i<<"-("<<v<<","<<orderWeight[i].first<<") "<<endl;//to be deleted
            }
           //out.close();//to be deleted
            //std::cout<<"*****************************calcWeightByValue finished!**************************"<<endl;
        }

        /*
         *@description: calculate the node weight by value
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void calcWeightByValue_degree(WGraph& wgraph){
            std::cout<<"*****************************calcWeightByValue_degree begins!**************************"<<endl;
            vector<pair<weightType,NodeID> > orderWeight;
            vector<NodeID> empyty_value_vertices;
            srand(0x728f9d0);//srand seed
            orderWeight.reserve(numOfVertices);
            for(NodeID v=0;v<numOfVertices;++v){
                weightType result=0;
                if(_calcCoef.is_freq_mult) result+=((weightType)_freq_weight[v]/(weightType)interval_freq)*_calcCoef.freq_mult;
                if(_calcCoef.is_bet_mult) result+=((weightType)_betwenness_weight[v]/(weightType)interval_betweenness)*_calcCoef.bet_mult;
                //if((weightType)_betwenness_weight[v]<1e-6) std::cout<<v<<":"<<(weightType)_betwenness_weight[v]<<endl;
                if(result<1e-6) empyty_value_vertices.push_back(v);
                else orderWeight.push_back(make_pair(result,v) );
            }
            sort(orderWeight.rbegin(),orderWeight.rend());
            NodeID r;
            for(r=0;r<orderWeight.size();++r){
                NodeID v=orderWeight[r].second;
                inv[r]=v;
                rank[v]=r;
            }
            //deal with the left vertices ranked by degree
            if(empyty_value_vertices.size()>0){
                cout<<"Rank left "<<numOfVertices-r<<" vertices by degree..."<<endl;
                vector<pair<weightType,NodeID> >().swap(orderWeight);
                orderWeight.clear();
                for(size_t i=0;i<empyty_value_vertices.size();++i) {
                    NodeID v=empyty_value_vertices[i];
                    orderWeight.push_back(make_pair((wgraph.vertices[v+1]-wgraph.vertices[v])+float(rand()) / RAND_MAX,v));
                }
                sort(orderWeight.rbegin(),orderWeight.rend());
                for(size_t i=0;i<orderWeight.size()&&r<numOfVertices;++r,++i){
                    NodeID v=orderWeight[i].second;
                    inv[r]=v;
                    rank[v]=r;
                }
            }
            //std::cout<<"*****************************calcWeightByValue finished!**************************"<<endl;
        }

        //Add for directed
        /*
         *@description: calculate the node weight by value
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void calcWeightByValue_degree_directed(WGraph& wgraph){
            std::cout<<"*****************************calcWeightByValue_degree directed begins!**************************"<<endl;
            vector<pair<weightType,NodeID> > orderWeight;
            vector<NodeID> empyty_value_vertices;
            srand(0x728f9d0);//srand seed
            orderWeight.reserve(numOfVertices);
            for(NodeID v=0;v<numOfVertices;++v){
                weightType result=0;
                if(_calcCoef.is_freq_mult) result+=((weightType)_freq_weight[v]/(weightType)interval_freq)*_calcCoef.freq_mult;
                if(_calcCoef.is_bet_mult) result+=((weightType)_betwenness_weight[v]/(weightType)interval_betweenness)*_calcCoef.bet_mult;
                if(result<1e-6) empyty_value_vertices.push_back(v);
                else orderWeight.push_back(make_pair(result,v) );
            }
            sort(orderWeight.rbegin(),orderWeight.rend());
            NodeID r;
            for(r=0;r<orderWeight.size();++r){
                NodeID v=orderWeight[r].second;
                inv[r]=v;
                rank[v]=r;
            }
            //deal with the left vertices ranked by degree
            if(empyty_value_vertices.size()>0){
                cout<<"Rank left "<<numOfVertices-r<<" vertices by degree..."<<endl;
                vector<pair<weightType,NodeID> >().swap(orderWeight);
                orderWeight.clear();
                for(size_t i=0;i<empyty_value_vertices.size();++i) {
                    NodeID v=empyty_value_vertices[i];
                    orderWeight.push_back(make_pair((wgraph.vertices[v + 1] - wgraph.vertices[v]) * (wgraph.r_vertices[v + 1] - wgraph.r_vertices[v]) +float(rand()) / RAND_MAX,v));
                }
                sort(orderWeight.rbegin(),orderWeight.rend());
                for(size_t i=0;i<orderWeight.size()&&r<numOfVertices;++r,++i){
                    NodeID v=orderWeight[i].second;
                    inv[r]=v;
                    rank[v]=r;
                }
            }
            //std::cout<<"*****************************calcWeightByValue finished!**************************"<<endl;
        }

        /*
        *@description: used to init all iterms' values
        *@author: wanjingyi
        *@date: 2020-12-13
        */
        bool initIterms(){
            //inv rank
            inv.resize(numOfVertices);
            rank.resize(numOfVertices);
            //isCnstruct
            isConstruct.resize(numOfVertices,false);
            //query frequency
            _freq.resize(numOfVertices,0);
            _freq_rank.resize(numOfVertices,0);
            _freq_inv.resize(numOfVertices);
            numOfHFpoint=0;
            HFPointIndex.resize(numOfVertices,0);
            //degree
            _degree.resize(numOfVertices,0);
            _degree_rank.resize(numOfVertices,0);
            _degree_inv.resize(numOfVertices);
            //betwenness
            _betwenness.resize(numOfVertices,0);
            _betwenness_inv.resize(numOfVertices);
            _betwenness_rank.resize(numOfVertices,0);
            //coverage
            _coverage.resize(numOfVertices,0);
            _coverage_inv.resize(numOfVertices);
            _coverage_rank.resize(numOfVertices,0);
            //depth
            _depth_rank.resize(numOfVertices,0);
            _depth.resize(numOfVertices,0);
            _depth_inv.resize(numOfVertices);
            //for tie breaking
            _degree_weight.resize(numOfVertices,0);
            _freq_weight.resize(numOfVertices,0);
            _betwenness_weight.resize(numOfVertices,0);
        }

        /*
        *@description: used to init all iterms' values
        *@author: wanjingyi
        *@date: 2020-12-13
        */
        bool initIterms_overlay(){
            //inv rank
            inv.resize(numOfVertices);
            rank.resize(numOfVertices);
            //for tie breaking
            _degree_weight.resize(numOfVertices,0);
            _freq_weight.resize(numOfVertices,0);
            _betwenness_weight.resize(numOfVertices,0);
        }

        /*
         *@description: tie-break betweenness
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void load_betweenness_overlay(const vector<unsigned int>& betwenness){
            NodeID id;NodeID bet;            
            srand(100);
            for(NodeID id=0;id<numOfVertices;++id){
                _betwenness_weight[id]=betwenness[id]+ float(rand()) / RAND_MAX;
                if(bet<min_betweenness) min_betweenness=bet;
                if(bet>max_betweennss) max_betweennss=bet;
            }
            interval_betweenness=max_betweennss-min_betweenness;
            std::cout<<"***********load betweenness overlay finished************"<<endl;
        }

        /*
         *@description: used to tie break the frequency values
         *@author: wanjingyi
         *@date: 2020-12-13
        */
        void load_HFpoint_overlay(const vector<unsigned int>& freq){//hfRate/1000
            srand(100);
            for(NodeID id=0;id<numOfVertices;++id){
                //weightType freq_w=(weightType)(freq[id]+ float(rand()) / RAND_MAX);
                weightType freq_w=static_cast<weightType>(freq[id]);
                _freq_weight[id]=freq_w;
                if(freq_w>max_freq) max_freq=freq_w;
                if(freq_w<min_freq) min_freq=freq_w;
			}
            interval_freq=max_freq-min_freq;
        }

        void load_betwenness_overlay(const vector<unsigned int>& bet){
            srand(100);
            for(NodeID id=0;id<numOfVertices;++id){
                //_betwenness_weight[id]=(weightType)(bet[id]+ float(rand()) / RAND_MAX);
                _betwenness_weight[id]=static_cast<weightType>(bet[id]);
                if(bet[id]<min_betweenness) min_betweenness=bet[id];
                if(bet[id]>max_betweennss) max_betweennss=bet[id];
			}
            interval_betweenness=max_betweennss-min_betweenness;
        }

        void load_coverage_overlay(const vector<unsigned int>& cov){
            return;
        }

        /*
         *@description: used to load frequency and hfpoint
         *@author: wanjingyi
         *@date: 2020-12-13
        */
        void load_HFpoint(char* load_filename,int hfRate){ //hfRate/1000
            //std::cout<<"***********load hfpoint begins************"<<endl;
            //std::cout<<" hfRate ="<<hfRate<<endl;
            if(hfRate==0) numOfHFpoint=numOfVertices;
            else numOfHFpoint = static_cast<int> ( (double)numOfVertices*hfRate/(double)HF_DIVIDION);
            if(numOfHFpoint<=0) cout<<"error:numOfHFpoint<=0"<<endl;
            //cout<<"numOfHFpoint  = "<<numOfHFpoint <<endl;
            ifstream in(load_filename);//input HFPoint file to ifstream
            if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
            vector<pair<float,NodeID> > queryFreq;//(freq,id)
            NodeID id,freq;
            srand(100);
            char line[24];
			long long  query_cnt=0;
            int hf_cnt=0;
            while (in.getline(line,sizeof(line)))
            {
                if(hf_cnt>=numOfHFpoint) break;
				stringstream ls(line);
                ls>>id>>freq;
                query_cnt+=freq;
                _freq[id]=freq;
               HFPointIndex[id]=true;
                weightType weight_tmp=freq+ float(rand()) / RAND_MAX;
                _freq_weight[id]=weight_tmp;
                queryFreq.push_back(make_pair(weight_tmp,id));
                ++hf_cnt;
            }
            if(hf_cnt<numOfHFpoint) numOfHFpoint=hf_cnt;
            for(NodeID v=0;v<numOfVertices;++v){
                if(_freq[v]==0){
                    queryFreq.push_back(make_pair(float(rand()) / RAND_MAX,v));
                }
            }
            //sort the node by query times
            //sort(queryFreq.rbegin(),queryFreq.rend());//descending order
            //initialize query freq information
            for(NodeID i=0;i<numOfVertices;++i){
                NodeID v=queryFreq[i].second;
                _freq_inv[i]=v;
                _freq_rank[v]=i;
            }
            max_freq=queryFreq[0].first;
            min_freq=0;
            interval_freq=max_freq-min_freq;
            /**
            std::cout<<"max_query_time = "<<max_freq<<endl;
            std::cout<<"min_query_time = "<<min_freq<<endl;
            std::cout<<"sum_query_time = "<<sum_freq<<endl;
            std::cout<<"interval_freq = "<<interval_freq<<endl;
            std::cout<<"***********load hfpoint finished************"<<endl;
            **/
        }

};

#endif
