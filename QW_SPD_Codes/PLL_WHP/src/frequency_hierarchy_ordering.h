/*
 * @Descripttion:  query workload aware ordering
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-10 21:28:01
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-11-03 16:32:00
 */
#pragma once
#ifndef FREQUENCY_HIERARCHY_ORDERING_H
#define FREQUENCY_HIERARCHY_ORDERING_H

#include<algorithm>
#include<unordered_set>
#include<time.h>
#include "graph.h"
#include "graph_search.h"
#include "labels.h"
#include "time_util.h" 
#include "heap.h"
#include "paras.h"
#include "ordering.h"
#include "construction.h"
#include<unordered_map>  
#include<cmath>
#ifdef _WIN32
	#include<google/sparsehash/sparseconfig.h>
#endif
    #include<google/dense_hash_map>
	#include<google/dense_hash_set>



#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define cover_value_type long long
#define countType unsigned int
#define MAX_VALUE 0x7FFFFFFF
//DEBUG
#define DEBUG_FLAG 1
char debugFileName[255] = "../dataset/manhatan/SHP/SOrder"; 

//pair compare 
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

class FOrdering:public Ordering{
    public:
        //*******************variables*******************
        NodeID _numOfHfpoint=0;//num of hfpoint
        vector<bool> HFPointFlag;//indicates whether node is hfpoint
        
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

template<typename weightType>
class Greedy_fordering : public FOrdering{
    public:
        //****************variables****************
        vector<int> _freq;//get the query time of each node by index
        vector<int> _betwenness;//get the betwenness by node index
        int  _max_freq;//max query time used to normalization
        int _max_betweenness;//max betweenness of all nodes
        int _min_betweenness;//min betweenness of all nodes
        long long _sum_freq;//total num of query time
        weightType  _total_freq_weight;
        vector<weightType> _freq_weight;//get the query time weight value of each node by index
        vector<weightType> _betwenness_weight;//get the betwenness weight value by node index
        //*****************constructions and deconstructions****************
        Greedy_fordering(WGraph& wgraph,char* bet_filename,char* point_freq_filename){ //undirected weighted graph
            undirected_weighted_order(wgraph, bet_filename,point_freq_filename);
        }
        Greedy_fordering(){}
        ~Greedy_fordering(){}

    protected:
        //************construction functions***************
        void undirected_weighted_order(WGraph& wgraph, char* bet_filename,char* point_freq_filename){
            load_and_preprocess_queryFreq(point_freq_filename);
            load_and_preprocess_betwenness(bet_filename);
            greedy_selection();
            Relabel(wgraph);
        }

        //************used functions**************
        void greedy_selection(){
            //initialize variables
            inv.resize(numOfVertices);
            rank.resize(numOfVertices);
            vector<bool> isSelected(numOfVertices,false);//inidicating whether has selected
            NodeID cnt=0;
            weightType delta_query_cost,min_delta_query_cost;
            NodeID min_id;
            bool isFirst;
            while (cnt<numOfVertices)
            {
                isFirst=true;
                for(NodeID v=0;v<numOfVertices;++v){
                    if(isSelected[v]) continue;
                    delta_query_cost=(_total_freq_weight-_freq_weight[v])*(AMPLIFY_FACTOR-_betwenness_weight[v]);
                    if(isFirst){
                        isFirst=false;
                        min_delta_query_cost=delta_query_cost;
                        min_id=v;
                    }else{
                        if(delta_query_cost<min_delta_query_cost){
                            min_delta_query_cost=delta_query_cost;
                            min_id=v;
                        }
                    }
                }
                //if(cnt<100) cout<<cnt<<" "<<min_id<<" "<<min_delta_query_cost<<endl;//to be deleted
                _total_freq_weight-=_freq_weight[min_id];
                rank[min_id]=cnt;
                inv[cnt]=min_id;
                isSelected[min_id]=true;
                ++cnt;
            }
            
        }

        void load_and_preprocess_betwenness(char* load_filename){
            //initialize
            _betwenness.resize(numOfVertices,0);
            _min_betweenness=INF_WEIGHT;
            _max_betweenness=0;
            _betwenness_weight.resize(numOfVertices,0);
            ifstream in(load_filename);//input coverage file to ifstream
            if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
            NodeID id;
            int bet;
            for(NodeID i=0;i<numOfVertices;++i){
                in>>id>>bet;
                _betwenness[id]=bet;
                
                if(bet<_min_betweenness) _min_betweenness=bet;
                if(bet>_max_betweenness) _max_betweenness=bet;
            }
            int _interval_bet=_max_betweenness-_min_betweenness;
            for(NodeID i=0;i<numOfVertices;++i){
                _betwenness_weight[i]=(weightType)_betwenness[i]*AMPLIFY_FACTOR/(weightType)_interval_bet;
            }
            in.close();
        }

        void load_and_preprocess_queryFreq(char* load_filename){
            //initialize
            _total_freq_weight=0;
            _sum_freq =0;
            _max_freq = 0;
            _freq.resize(numOfVertices,0);
            _freq_weight.resize(numOfVertices,0);
            //input HFPoint file to ifstream
            ifstream in(load_filename);
            if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
            char line[24];
            NodeID id;
            int freq;
            while (in.getline(line,sizeof(line)))
            {
                stringstream ls(line);
                ls>>id>>freq;
                _freq[id]=freq;
                _sum_freq+=freq;
                if(freq>_max_freq) _max_freq=freq;
            }
            in.close();
            for(NodeID v=0;v<numOfVertices;++v){
                weightType freq_tmp=(weightType)_freq[v]*AMPLIFY_FACTOR/(weightType)_max_freq;
                _freq_weight[v]=freq_tmp;
                _total_freq_weight+=freq_tmp;
            }
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
        //HFLabel labels;
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
        vector<NodeID> _degree_rank;//fetch the rank  by indexs
        ofstream debug_out;//output debug information
        unsigned int  total_freq=0;//total point freq time
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
        *@description: construction used for undirected and weighted graph
        *@author: wanjingyi
        *@date: 2020-12-13
        */
        Linear_Ordering(WGraph& wgraph,Processing::calcCoefficient<double> calcCoef,char* orderWeightFileName,char* queryFreqFileName,char* betwennessFileName,char* coverageFileName,int hfRate,int t_ordering_flag=0){
            //std::cout<<"Linear_Ordering undirected weighted graph...."<<endl;
            if(DIRECTED_FLAG==false) undirected_weighted_orderByValues(wgraph,calcCoef,orderWeightFileName,queryFreqFileName,betwennessFileName,coverageFileName,hfRate);
            else{
                directed_weighted_orderByValues(wgraph,calcCoef,orderWeightFileName,queryFreqFileName,betwennessFileName,coverageFileName,hfRate);
            }
            Relabel(wgraph);
        }

    protected:
        /*
         *@description: order the nodes by detailed_values*Coefficient sum
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void undirected_weighted_orderByValues(WGraph& wgraph,Processing::calcCoefficient<double> calcCoef,char* orderWeightFileName,char* queryFreqFileName,char* betwennessFileName,char* coverageFileName,int hfRate)
        {
            initIterms();//init all values
            _calcCoef=calcCoef;//get the coefficient from command line
            load_HFpoint(queryFreqFileName,hfRate);//frequency
            if(_calcCoef.is_bet_mult) getBetwennessFromFile(betwennessFileName);//betwenness
            if(_calcCoef.is_deg_mult) calcWeightByValue_degree(wgraph,orderWeightFileName);//cumpute the order
            else calcWeightByValue(orderWeightFileName);//cumpute the order
        }

        /*
         *@description: order the nodes by detailed_values*Coefficient sum
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void directed_weighted_orderByValues(WGraph& wgraph,Processing::calcCoefficient<double> calcCoef,char* orderWeightFileName,char* queryFreqFileName,char* betwennessFileName,char* coverageFileName,int hfRate)
        {
            initIterms();//init all values
            _calcCoef=calcCoef;//get the coefficient from command line
            load_HFpoint(queryFreqFileName,hfRate);//frequency
            if(_calcCoef.is_bet_mult) getBetwennessFromFile(betwennessFileName);//betwenness
            if(_calcCoef.is_deg_mult) calcWeightByValue_degree_directed(wgraph,orderWeightFileName);//cumpute the order
            else calcWeightByValue(orderWeightFileName);//cumpute the order
        }

        /*
         *@description: calculate the node weight by value
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void calcWeightByValue(char* orderWeightFileName){
            //std::cout<<"*****************************calcWeightByValue begins!**************************"<<endl;
            vector<pair<weightType,NodeID> > orderWeight;
            srand(0x728f9d0);//srand seed
            orderWeight.reserve(numOfVertices);
            ofstream weight_out;
            bool isOutputWeight=false;
            //std::cout<<" freq_mult="<<_calcCoef.freq_mult<<" deg_mult="<<_calcCoef.deg_mult<<" bet_mult="<<_calcCoef.bet_mult<<" cov_mult="<<_calcCoef.cov_mult<<" dep_mult="<<_calcCoef.dep_mult;
            for(NodeID v=0;v<numOfVertices;++v){
                weightType result=0;
                // result=(weightType)_freq[v]*_calcCoef.freq_mult+(weightType)_degree[v]*_calcCoef.deg_mult+
                // (weightType)_betwenness[v]*_calcCoef.bet_mult+(weightType)_coverage[v]*_calcCoef.cov_mult+
                // (weightType)_depth[v]*_calcCoef.dep_mult;//tie breaking +weightType(rand()) / RAND_MAX
                //if(_calcCoef.is_deg_mult) result+=((weightType)(_degree_weight[v]-min_degree)/(weightType)interval_degree)*_calcCoef.deg_mult;
                if(_calcCoef.is_freq_mult) result+=((weightType)(_freq_weight[v]-min_freq)/(weightType)interval_freq)*_calcCoef.freq_mult;
                if(_calcCoef.is_bet_mult) result+=((weightType)(_betwenness_weight[v]-min_betweenness)/(weightType)interval_betweenness)*_calcCoef.bet_mult;
                orderWeight.push_back(make_pair(result,v) );
            }
            sort(orderWeight.rbegin(),orderWeight.rend());
            //output orderWeight
            if(*orderWeightFileName!='\0'){
                string write_filename_prefix(orderWeightFileName);
                string orderWeight_filename=write_filename_prefix.append(".orderWeight");
                weight_out.open(orderWeight_filename.c_str());//to be deleted
                if(!weight_out.is_open()) std::cout<<"Cannot open orderWeight_filename="<<orderWeight_filename<<endl;
                else isOutputWeight=true;
            }

            for(NodeID i=0;i<numOfVertices;++i){
                NodeID v=orderWeight[i].second;
                inv[i]=v;
                rank[v]=i;
                if(isOutputWeight) weight_out<<i<<"-("<<v<<","<<orderWeight[i].first<<") "<<endl;//output order_weight
            }
            if(isOutputWeight) weight_out.close();//to be deleted
            //std::cout<<"*****************************calcWeightByValue finished!**************************"<<endl;
        }

        /*
         *@description: calculate the node weight by value
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void calcWeightByValue_degree(WGraph& wgraph,char* orderWeightFileName){
            //std::cout<<"*****************************calcWeightByValue begins!**************************"<<endl;
            vector<pair<weightType,NodeID> > orderWeight;
            vector<NodeID> empyty_value_vertices;
            srand(0x728f9d0);//srand seed
            orderWeight.reserve(numOfVertices);
            //std::cout<<" freq_mult="<<_calcCoef.freq_mult<<" deg_mult="<<_calcCoef.deg_mult<<" bet_mult="<<_calcCoef.bet_mult<<" cov_mult="<<_calcCoef.cov_mult<<" dep_mult="<<_calcCoef.dep_mult;
            ofstream weight_out;
            bool isOutputWeight=false;
            
            for(NodeID v=0;v<numOfVertices;++v){
                weightType result=0;
                // result=(weightType)_freq[v]*_calcCoef.freq_mult+(weightType)_degree[v]*_calcCoef.deg_mult+
                // (weightType)_betwenness[v]*_calcCoef.bet_mult+(weightType)_coverage[v]*_calcCoef.cov_mult+
                // (weightType)_depth[v]*_calcCoef.dep_mult;//tie breaking +weightType(rand()) / RAND_MAX
                //if(_calcCoef.is_deg_mult) result+=((weightType)(_degree_weight[v]-min_degree)/(weightType)interval_degree)*_calcCoef.deg_mult;
                if(_calcCoef.is_freq_mult) result+=((weightType)(_freq_weight[v]-min_freq)/(weightType)interval_freq)*_calcCoef.freq_mult;
                if(_calcCoef.is_bet_mult) result+=((weightType)(_betwenness_weight[v]-min_betweenness)/(weightType)interval_betweenness)*_calcCoef.bet_mult;
                if(result<1e-6) empyty_value_vertices.push_back(v);
                else orderWeight.push_back(make_pair(result,v) );
            }
            sort(orderWeight.rbegin(),orderWeight.rend());

            //output orderWeight
            if(*orderWeightFileName!='\0'){
                string write_filename_prefix(orderWeightFileName);
                string orderWeight_filename=write_filename_prefix.append(".orderWeight");
                weight_out.open(orderWeight_filename.c_str());//to be deleted
                if(!weight_out.is_open()) std::cout<<"Cannot open orderWeight_filename="<<orderWeight_filename<<endl;
                else isOutputWeight=true;
            }

            NodeID r;
            for(r=0;r<orderWeight.size();++r){
                NodeID v=orderWeight[r].second;
                inv[r]=v;
                rank[v]=r;
                if(isOutputWeight) weight_out<<r<<"-("<<v<<","<<orderWeight[r].first<<") "<<endl;//output order_weight
            }
           //out.close();//to be deleted
           //deal with the left vertices ranked by degree
            if(empyty_value_vertices.size()>0){
                cout<<"Rank left "<<numOfVertices-r<<" vertices by degree..."<<endl;
                if(isOutputWeight) weight_out<<"**************************************************"<<endl;
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
                    if(isOutputWeight) weight_out<<r<<"-("<<v<<","<<orderWeight[i].first<<") "<<endl;//output order_weight
                }
            }
            if(isOutputWeight) weight_out.close();
            //std::cout<<"*****************************calcWeightByValue finished!**************************"<<endl;
        }

        /*
         *@description: calculate the node weight by value
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void calcWeightByValue_degree_directed(WGraph& wgraph,char* orderWeightFileName){
            //std::cout<<"*****************************calcWeightByValue begins!**************************"<<endl;
            vector<pair<weightType,NodeID> > orderWeight;
            vector<NodeID> empyty_value_vertices;
            srand(0x728f9d0);//srand seed
            orderWeight.reserve(numOfVertices);
            //std::cout<<" freq_mult="<<_calcCoef.freq_mult<<" deg_mult="<<_calcCoef.deg_mult<<" bet_mult="<<_calcCoef.bet_mult<<" cov_mult="<<_calcCoef.cov_mult<<" dep_mult="<<_calcCoef.dep_mult;
            ofstream weight_out;
            bool isOutputWeight=false;
            
            for(NodeID v=0;v<numOfVertices;++v){
                weightType result=0;
                // result=(weightType)_freq[v]*_calcCoef.freq_mult+(weightType)_degree[v]*_calcCoef.deg_mult+
                // (weightType)_betwenness[v]*_calcCoef.bet_mult+(weightType)_coverage[v]*_calcCoef.cov_mult+
                // (weightType)_depth[v]*_calcCoef.dep_mult;//tie breaking +weightType(rand()) / RAND_MAX
                //if(_calcCoef.is_deg_mult) result+=((weightType)(_degree_weight[v]-min_degree)/(weightType)interval_degree)*_calcCoef.deg_mult;
                if(_calcCoef.is_freq_mult) result+=((weightType)(_freq_weight[v]-min_freq)/(weightType)interval_freq)*_calcCoef.freq_mult;
                if(_calcCoef.is_bet_mult) result+=((weightType)(_betwenness_weight[v]-min_betweenness)/(weightType)interval_betweenness)*_calcCoef.bet_mult;
                if(result<1e-6) empyty_value_vertices.push_back(v);
                else orderWeight.push_back(make_pair(result,v) );
            }
            sort(orderWeight.rbegin(),orderWeight.rend());

            //output orderWeight
            if(*orderWeightFileName!='\0'){
                string write_filename_prefix(orderWeightFileName);
                string orderWeight_filename=write_filename_prefix.append(".orderWeight");
                weight_out.open(orderWeight_filename.c_str());//to be deleted
                if(!weight_out.is_open()) std::cout<<"Cannot open orderWeight_filename="<<orderWeight_filename<<endl;
                else isOutputWeight=true;
            }

            NodeID r;
            for(r=0;r<orderWeight.size();++r){
                NodeID v=orderWeight[r].second;
                inv[r]=v;
                rank[v]=r;
                if(isOutputWeight) weight_out<<r<<"-("<<v<<","<<orderWeight[r].first<<") "<<endl;//output order_weight
            }
           //out.close();//to be deleted
           //deal with the left vertices ranked by degree
            if(empyty_value_vertices.size()>0){
                cout<<"Rank left "<<numOfVertices-r<<" vertices by degree..."<<endl;
                if(isOutputWeight) weight_out<<"**************************************************"<<endl;
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
                    if(isOutputWeight) weight_out<<r<<"-("<<v<<","<<orderWeight[i].first<<") "<<endl;//output order_weight
                }
            }
            if(isOutputWeight) weight_out.close();
            //std::cout<<"*****************************calcWeightByValue finished!**************************"<<endl;
        }

    /*
     *@description: calculate the node weight by query cost
     *@author: wanjingyi
     *@date: 2020-12-21
    */
    void calcWeightByCost(){

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
            _numOfHfpoint=0;
            HFPointFlag.resize(numOfVertices,0);
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
         *@description: used to read betwenness order from file
         *@author: wanjingyi
         *@date: 2020-12-13
        */
        void getBetwennessOrderFromFile(char* load_filename){
            std::cout<<"**********************getBetwennessOrderFromFile begins!**********************"<<endl;
            ifstream in(load_filename);//input betwenness file to ifstream
            if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
            NodeID id;
            for(NodeID i=0;i<numOfVertices;++i){
                in>>id;
                _betwenness_rank[id]=i;
                _betwenness_inv[i]=id;
            }
            in.close();
            std::cout<<"**********************getBetwennessOrderFromFile finished!**********************"<<endl;
        }

        /*
         *@description: read the node's betweenness from file
         *@author: wanjingyi
         *@date: 2020-12-16
        */
        void getBetwennessFromFile(char* load_filename){
            //std::cout<<"***********************getBetwennessFromFile begins!*************************"<<endl;
            ifstream in(load_filename);//input coverage file to ifstream
            if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
            NodeID id;NodeID bet;
            srand(100);
            for(NodeID i=0;i<numOfVertices;++i){
                in>>id>>bet;
                _betwenness[id]=bet;
                _betwenness_weight[id]=static_cast<weightType>(bet);
                if(bet<min_betweenness) min_betweenness=bet;
                if(bet>max_betweennss) max_betweennss=bet;
            }
            in.close();
            interval_betweenness=max_betweennss-min_betweenness;
            /**
            std::cout<<"max_dbetweennss = "<<max_betweennss<<endl;
            std::cout<<"min_betweenness = "<<min_betweenness<<endl;
            std::cout<<"interval_betweenness = "<<interval_betweenness<<endl;
            std::cout<<"***********************getBetwennessFromFile finished!!*************************"<<endl;
            **/
        }

        /*
         *@description: used to read coverage order from file
         *@author: wanjingyi
         *@date: 2020-12-13
        */
        void getCoverageOrderFromFile(char* load_filename){
            std::cout<<"***********************getCoverageOrderFromFile begins!*************************"<<endl;
            ifstream in(load_filename);//input coverage file to ifstream
            if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
            NodeID id;
            for(NodeID i=0;i<numOfVertices;++i){
                in>>id;
                _coverage_rank[id]=i;
                _coverage_inv[i]=id;
            }
            in.close();
            std::cout<<"***********************getCoverageOrderFromFile finished!!*************************"<<endl;
        }

        void getCoverageFromFile(char* load_filename){
            return;
        }

        void load_original_HFpoint(char* load_filename,int hfRate=50){//hfRate/1000

        }

        // void load_HFpoint(char* load_filename,int hfRate){ //hfRate/1000
        //     //std::cout<<"***********load hfpoint begins************"<<endl;
        //     //std::cout<<" hfRate ="<<hfRate<<endl;
        //     if(hfRate==0) _numOfHfpoint=numOfVertices;
        //     else _numOfHfpoint = static_cast<int> ( (double)numOfVertices*hfRate/(double)HF_DIVIDION);
        //     if(_numOfHfpoint<=0) cout<<"error:numOfHFpoint<=0"<<endl;
        //     //cout<<"numOfHFpoint  = "<<numOfHFpoint <<endl;
        //     ifstream in(load_filename);//input HFPoint file to ifstream
        //     if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
        //     vector<pair<float,NodeID> > queryFreq;//(freq,id)
        //     NodeID id,freq;
        //     srand(100);
        //     char line[24];
		// 	long long  query_cnt=0;
        //     int hf_cnt=0;
        //     while (in.getline(line,sizeof(line)))
        //     {
        //         if(hf_cnt>=_numOfHfpoint) break;
		// 		stringstream ls(line);
        //         ls>>id>>freq;
        //         query_cnt+=freq;
        //         _freq[id]=freq;
        //        HFPointFlag[id]=true;
        //         _freq_weight[id]=(weightType) freq;
        //         queryFreq.push_back(make_pair((weightType) freq,id));
        //         ++hf_cnt;
        //     }
        //     if(hf_cnt<_numOfHfpoint) _numOfHfpoint=hf_cnt;
        //     for(NodeID v=0;v<numOfVertices;++v){
        //         if(_freq[v]==0){
        //             queryFreq.push_back(make_pair(float(rand()) / RAND_MAX,v));
        //         }
        //     }
        //     //sort the node by query times
        //     //sort(queryFreq.rbegin(),queryFreq.rend());//descending order
        //     //initialize query freq information
        //     for(NodeID i=0;i<numOfVertices;++i){
        //         NodeID v=queryFreq[i].second;
        //         _freq_inv[i]=v;
        //         _freq_rank[v]=i;
        //     }
        //     max_freq=queryFreq[0].first;
        //     min_freq=0;
        //     interval_freq=max_freq-min_freq;
        //     /**
        //     std::cout<<"max_query_time = "<<max_freq<<endl;
        //     std::cout<<"min_query_time = "<<min_freq<<endl;
        //     std::cout<<"sum_query_time = "<<sum_freq<<endl;
        //     std::cout<<"interval_freq = "<<interval_freq<<endl;
        //     std::cout<<"***********load hfpoint finished************"<<endl;
        //     **/
        // }

        /*
         *@description: used to load frequency and hfpoint
         *@author: wanjingyi
         *@date: 2020-12-13
        */
        void load_HFpoint(char* HFPoint_file,int hfRate,bool isSort=false){ 
            _numOfHfpoint=0;
            if(is_directory(HFPoint_file)){
                string pointFreqPath(HFPoint_file);
                vector<string> freqid_files_tmp;
                get_filelist_from_dir(pointFreqPath,freqid_files_tmp,true);
                append_to_full_path(pointFreqPath,freqid_files_tmp);
                std::cout<<"freqid_files_tmp.size() = "<<freqid_files_tmp.size()<<std::endl;
                for(string file: freqid_files_tmp){
                    load_hfpoint_slice((char*)file.c_str(),hfRate);
                }
            }else{
                load_hfpoint_slice(HFPoint_file,hfRate);
            }
            for(NodeID id=0;id<numOfVertices;++id){
                _freq_weight[id]=static_cast<weightType> (_freq[id]);
                if(_freq[id]!=0){
                    HFPointFlag[id]=true;
                    ++_numOfHfpoint;                      
                    total_freq+=_freq[id];                                                                                                                                     
                }
                if(_freq[id]>max_freq) max_freq=_freq[id];
                if(_freq[id]<min_freq) min_freq=_freq[id];
            }
            interval_freq=max_freq-min_freq;
            std::cout<<"max_freq = "<<max_freq<<endl;
            std::cout<<"min_freq = "<<min_freq<<endl;
            std::cout<<"interval_freq = "<<interval_freq<<endl;
            if(isSort){
                vector<pair<float,NodeID> > queryFreq;//(freq,id)
                for(NodeID v=0;v<numOfVertices;++v){
                    queryFreq.push_back(make_pair(static_cast<weightType> (_freq[v]),v));
                }
                //sort the node by query times
                sort(queryFreq.rbegin(),queryFreq.rend());//descending order
                //initialize query freq information
                for(NodeID i=0;i<numOfVertices;++i){
                    NodeID v=queryFreq[i].second;
                    _freq_inv[i]=v;
                    _freq_rank[v]=i;
                }
            }
            std::cout<<"_numOfHfpoint="<<_numOfHfpoint<<std::endl;
            std::cout<<"*******************Load hfpoint successfully!*******************"<<std::endl;
        }

        /*
         *@description: used to load frequency and hfpoint
         *@author: wanjingyi
         *@date: 2020-12-13
        */
        void load_hfpoint_slice(char* load_filename,int hfRate){ //hfRate/1000
            int numOfHfpoint_slice;
            if(hfRate==0) numOfHfpoint_slice=numOfVertices;
            else numOfHfpoint_slice = static_cast<int> ( (double)numOfVertices*hfRate/(double)HF_DIVIDION);
            if(numOfHfpoint_slice<=0) cout<<"error:_numOfHfpoint<=0"<<endl;
            ifstream in(load_filename);//input HFPoint file to ifstream
            if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
            char line[24];
            int hf_cnt=0;
            NodeID id,freq;
            while (in.getline(line,sizeof(line)))
            {
                if(hf_cnt>=numOfHfpoint_slice) break;
				stringstream ls(line);
                ls>>id>>freq;
                _freq[id]+=freq;
                ++hf_cnt;
            }
            if(hf_cnt<_numOfHfpoint) numOfHfpoint_slice=hf_cnt;
            //std::cout<<"numOfHfpoint_slice = "<<numOfHfpoint_slice<<std::endl;
            return;
        }

};


#endif //FREQUENCY_HIERARCHY_ORDERING