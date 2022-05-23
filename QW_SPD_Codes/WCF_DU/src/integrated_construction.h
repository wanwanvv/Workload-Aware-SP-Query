/*
 * @Descripttion: used to labelings construction for integrated soluutions
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-04-13 10:28:56
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-11-07 09:42:37
 */
#ifndef  INTEGRATED_CONSTURCTION_H
#define INTEGRATED_CONSTURCTION_H

#include <vector>
#include <fstream>
#include "./paras.h"
#include "./utils.h"
#include "./time_util.h"
#include "./h2h_construction.h"
#include "./integrated_index.h"
#include "./ordering.h"
#include "./frequency_hierarchy_ordering.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define INF_WEIGHT SP_Constants::INF_WEIGHT

using namespace std;
using namespace time_util;

class Integrated_construction{
public:
    //**************variables**************
    WGraph overlay_wgraph;

    double beta=1;
    double time_build_index=0;
    double time_build_rmq=0;
    vector<bool> _isDeleted;
    vector<int> _newToOriginal;
    vector<int> _originalToNew;
    int numOfOverlayVertices=0;
    int numOfOverlayEdges=0;
    int numOfOriginalVertices=0;
    int numOfOriginalEdges=0;
    int numOfDeletedVertices=0;
    int numOfHfpoint=0;

    vector<unsigned int> betweenessValues;//get bet values for overlay left vertices
    vector<unsigned int> coverageValues;//get cov values for overlay left vertices
    vector<unsigned int> depthValues;//get dep values for overlay left vertices

    vector<bool> HFPointFlag;//hfpoint flag by node index
    vector<NodeID> HFPoint;//store the HFpoint
    vector<unsigned int> queryTime;//get query time for all vertices
    vector<unsigned int> queryTime_overlay;//get query time for overlay left vertices
    long long total_query_time;//total query time for all vertices
    int max_query_time;//max querytime for all vertices
    int min_query_time;//min querytime for all vertices

    //Add for directed
    vector<unsigned int> queryTime_s;//get query time for all start vertices
    vector<unsigned int> queryTime_t;//get query time for all termination vertices
    vector<unsigned int> queryTime_overlay_s;//get query time for overlay left start vertices
    vector<unsigned int> queryTime_overlay_t;//get query time for overlay left termination vertices

    //Performance_size_factors
    long long total_sum_size=0;
    double total_ave_size=0;
    long long hf_sum_size=0;
    double hf_ave_size=0;
    long long overlay_sum_size=0;
    double overlay_ave_size=0;
    long long tree_sum_size=0;
    double tree_ave_size=0;
    double normalization_query_cost=0;
    double standard_query_cost=0;

    //Performance_time_factors
    double _labeling_time=0;
    double _avg_labeling_time=0;
    double _ordering_time=0;
    double time_tree_index=0;
    double time_overlay_order=0;
    double time_overlay_index=0;
    double time_compute_rmq;
    
    //index
    Integrated_Index _integrated_index;

    //****************constructions and deconstructions*****************
    /**
     * @description: 
     * @param {int} hfRate
     * @return {*}
     * @author: Wan Jingyi
     */    
    Integrated_construction(char* graphFileName,char* labelFileName,char* labelSizeFileName,
                                                        char* outputDirName,char* analysisDirName,char* HFPoint_file,
                                                        int hfRate,int thresholdDegree,int bestDegree,
                                                        int index_model,int ordering_model,const Processing::calcCoefficient<double>& calcCoef,
                                                        int isExperiment,char* givenOrderFileName="",bool isSorted=false)
    {
        if(index_model==0){//the first integrated solution
            //load original graph
            WGraph wgraph;
            wgraph.load_graph(graphFileName);
            numOfOriginalVertices=numOfVertices;
            numOfOriginalEdges=numOfEdges;
            std::cout<<"numOfOriginalVertices = "<<numOfOriginalVertices<<" numOfOriginalEdges = "<<numOfOriginalEdges<<endl;
            //initialize variables
            _integrated_index.index_.resize(numOfOriginalVertices);
            _isDeleted.resize(numOfOriginalVertices,0);
            load_hfpoint(HFPoint_file,hfRate);
            //lower level h2h-index construction
            build_h2h_index(graphFileName,outputDirName,analysisDirName,thresholdDegree,isExperiment);
            _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
            std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" numOfOverlayEdges = "<<numOfOverlayEdges<<endl;
            load_hfpoint_overlay(numOfOverlayVertices);
            numOfVertices=numOfOverlayVertices;
            numOfEdges=numOfOverlayEdges;
            //build overlay index
            build_overlay_index(calcCoef,hfRate,ordering_model);
            numOfVertices=numOfOriginalVertices;
            numOfEdges=numOfOriginalEdges;
            _labeling_time=time_tree_index+time_overlay_index;
            _ordering_time=time_overlay_order;
            _avg_labeling_time=_labeling_time/(double)numOfVertices;
            std::cout<<"time_tree_index = "<<time_tree_index*1e6<<std::endl;
            std::cout<<"time_overlay_index = "<<time_overlay_index*1e6<<" time_overlay_order = "<<time_overlay_order*1e6<<std::endl;
            std::cout<<"_labeling_time = "<<_labeling_time*1e6<<" _avg_labeling_time = "<<_avg_labeling_time*1e6<<std::endl;
            //ouput isDeletedFile
            save_isDeleted(outputDirName);
            //save indexFile to binary file
            _integrated_index.save_labels(labelFileName);
            //save label size file
            _integrated_index.save_label_size(labelSizeFileName);
            if(isExperiment==0){
                //output analyze
                save_analysis_result(analysisDirName);
                //save order
                _integrated_index.save_rank(analysisDirName);
                //save label size by order
                _integrated_index.save_label_size_byOrder(analysisDirName);
                //output index list txt file for debug
                _integrated_index.write_labels(analysisDirName);
            }else if(isExperiment==1){
                //appendt analysis
                append_analysis_result(analysisDirName);
            }else if(isExperiment==2){//isExperiment==2 update model
                //append to file
                ofstream ofs(analysisDirName,ios::app|ios::out);//append way
                if(!ofs.is_open()) cout<<"Cannot open "<<analysisDirName<<endl;
                ofs.setf(ios::fixed);
                ofs.precision(4);
                ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
                ofs.close();
            }else{
                
            }
            
        }else if(index_model==1){//the second integrated solution
            //load original graph
            WGraph wgraph;
            wgraph.load_graph(graphFileName);
            numOfOriginalVertices=numOfVertices;
            numOfOriginalEdges=numOfEdges;
            std::cout<<"numOfOriginalVertices = "<<numOfOriginalVertices<<" numOfOriginalEdges = "<<numOfOriginalEdges<<endl;
            //initialize variables
            _integrated_index.index_.resize(numOfOriginalVertices);
            _isDeleted.resize(numOfOriginalVertices,0);
            load_hfpoint(HFPoint_file,hfRate);
            //deleted, overlay label and h2h-index construction
            delete_and_build_forest(graphFileName,outputDirName,analysisDirName,calcCoef,hfRate,ordering_model,thresholdDegree,bestDegree,isExperiment,givenOrderFileName,isSorted);
            _labeling_time=time_tree_index+time_overlay_index;
            _ordering_time=time_overlay_order;
            _avg_labeling_time=_labeling_time/(double)numOfVertices;
            std::cout<<"time_tree_index = "<<time_tree_index*1e6<<std::endl;
            std::cout<<"time_overlay_index = "<<time_overlay_index*1e6<<" time_overlay_order = "<<time_overlay_order*1e6<<std::endl;
            std::cout<<"_labeling_time = "<<_labeling_time*1e6<<" _avg_labeling_time = "<<_avg_labeling_time*1e6<<std::endl;
            //ouput isDeletedFile
            save_isDeleted(outputDirName);
            //save indexFile to binary file
            _integrated_index.save_labels_1(labelFileName);
            //save label size file
            if(*labelSizeFileName!='\0') _integrated_index.save_label_size(labelSizeFileName);
            if(isExperiment==0){
                //output analysis
                save_analysis_result(analysisDirName);
                //save different labels for debug
                _integrated_index.save_labels_differ_analysis_1(analysisDirName);
                //save order
                _integrated_index.save_rank(analysisDirName);
                //save label size by order
                _integrated_index.save_label_size_byOrder(analysisDirName);
                //output index list txt file for debug
                _integrated_index.write_labels_1(analysisDirName);
                //output overllay_graph
                overlay_wgraph.save_overlay_graph(outputDirName,_newToOriginal);
            }else if(isExperiment==1){
                //appendt analysis
                append_analysis_result(analysisDirName);
            }else if(isExperiment==2){//update model
                //append to file
                ofstream ofs(analysisDirName,ios::app|ios::out);//append way
                if(!ofs.is_open()) cout<<"Cannot open "<<analysisDirName<<endl;
                ofs.setf(ios::fixed);
                ofs.precision(4);
                ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
                ofs.close();
            }else{
                
            }
        }else if(index_model==2){//the third integrated solution
            //load original graph
            WGraph wgraph;
            wgraph.load_graph(graphFileName);
            numOfOriginalVertices=numOfVertices;
            numOfOriginalEdges=numOfEdges;
            std::cout<<"numOfOriginalVertices="<<numOfOriginalVertices<<" numOfOriginalEdges="<<numOfOriginalEdges<<endl;
            //initialize variables
            _integrated_index.index_.resize(numOfOriginalVertices);
            _isDeleted.resize(numOfOriginalVertices,0);
            load_hfpoint(HFPoint_file,hfRate);
            //deleted, overlay label and h2h-index construction
            delete_and_pushPLL_forest(graphFileName,outputDirName,analysisDirName,calcCoef,hfRate,ordering_model,thresholdDegree,bestDegree,isExperiment,givenOrderFileName);
            _labeling_time=time_tree_index+time_overlay_index;
            _ordering_time=time_overlay_order;
            _avg_labeling_time=_labeling_time/(double)numOfVertices;
            std::cout<<"time_tree_index = "<<time_tree_index*1e6<<std::endl;
            std::cout<<"time_overlay_index = "<<time_overlay_index*1e6<<" time_overlay_order = "<<time_overlay_order*1e6<<std::endl;
            std::cout<<"_labeling_time = "<<_labeling_time*1e6<<" _avg_labeling_time = "<<_avg_labeling_time*1e6<<std::endl;
            //ouput isDeletedFile
            save_isDeleted(outputDirName);
            //save indexFile to binary file
            _integrated_index.save_labels_2(labelFileName);
            //save label size file
            if(*labelSizeFileName!='\0') _integrated_index.save_label_size(labelSizeFileName);
            if(isExperiment==0){
                 _integrated_index.save_labels_differ_analysis_2(analysisDirName);
                //output analysis
                save_analysis_result(analysisDirName);
                //save order
                _integrated_index.save_rank(analysisDirName);
                //save label size by order
                _integrated_index.save_label_size_byOrder(analysisDirName);
                //output index list txt file for debug
                _integrated_index.write_labels_1(analysisDirName);
                //output overllay_graph
                overlay_wgraph.save_overlay_graph(outputDirName,_newToOriginal);
            }else if(isExperiment==1){
                //appendt analysis
                append_analysis_result(analysisDirName);
            }else if(isExperiment==2){
                ofstream ofs(analysisDirName,ios::app|ios::out);//append way
                if(!ofs.is_open()) cout<<"Cannot open "<<analysisDirName<<endl;
                ofs.setf(ios::fixed);
                ofs.precision(4);
                ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
                ofs.close();
            }else{//update model

            }
        }

    }

    /**
     * @description: integrated-index construction for directed graph
     * @param {*}
     * @return {*}
     * @author: Wan Jingyi
     */    
    Integrated_construction(char* graphFileName,char* labelFileName,char* labelSizeFileName,
                                                        char* outputDirName,char* analysisDirName,char* HFPoint_file,
                                                        char* queryPair_file,int hfRate,int thresholdDegree,
                                                        int bestDegree,int index_model,int ordering_model,
                                                        const Processing::calcCoefficient<double>& calcCoef,int isExperiment,int graph_type,
                                                        char* givenOrderFileName="",bool isSorted=false)
    {
        if(index_model==0){
            std::cout<<"The first integrated-index solution..."<<std::endl;
        }else if(index_model==1){
            WGraph wgraph;
            wgraph.load_graph(graphFileName,graph_type);
            numOfOriginalVertices=numOfVertices;
            numOfOriginalEdges=numOfEdges;
            //initialize variables
            _integrated_index.bindex_.resize(numOfOriginalVertices);
            _isDeleted.resize(numOfOriginalVertices,0);
            load_hfpoint(HFPoint_file,hfRate);
            //deleted, overlay label and h2h-index construction
            delete_and_build_forest_directed(graphFileName,outputDirName,analysisDirName,calcCoef,hfRate,ordering_model,thresholdDegree,bestDegree,isExperiment,graph_type,givenOrderFileName,isSorted);
            //test rightness
            // Check_overlay_graph check_overlay_graph;
            // check_overlay_graph.check_left_vertices_SP(wgraph,overlay_wgraph,_newToOriginal,outputDirName,numOfOriginalVertices, numOfOverlayVertices,5);
            _labeling_time=time_tree_index+time_overlay_index;
            _ordering_time=time_overlay_order;
            _avg_labeling_time=_labeling_time/(double)numOfVertices;
            std::cout<<"time_tree_index = "<<time_tree_index*1e6<<std::endl;
            std::cout<<"time_overlay_index = "<<time_overlay_index*1e6<<" time_overlay_order = "<<time_overlay_order*1e6<<std::endl;
            std::cout<<"_labeling_time = "<<_labeling_time*1e6<<" _avg_labeling_time = "<<_avg_labeling_time*1e6<<std::endl;
            //ouput isDeletedFile
            save_isDeleted(outputDirName);
            //save indexFile to binary file
            _integrated_index.save_labels_1_directed(labelFileName);
            //save label size file
            if(*labelSizeFileName!='\0') _integrated_index.save_label_size_directed(labelSizeFileName);
            if(isExperiment==0){
                //output analysis
                save_analysis_result_directed(analysisDirName,queryPair_file);
                // //save different labels for debug
                // _integrated_index.save_labels_differ_analysis_1(analysisDirName);
                //save order
                _integrated_index.save_rank_directed(analysisDirName);
                // //save label size by order
                // _integrated_index.save_label_size_byOrder(analysisDirName);
                //output index list txt file for debug
                _integrated_index.write_labels_1_directed(analysisDirName);
                //output overllay_graph
                overlay_wgraph.save_overlay_graph(outputDirName,_newToOriginal,graph_type);
            }else if(isExperiment==1){
                //appendt analysis
                append_analysis_result_directed(analysisDirName,queryPair_file);
            }else if(isExperiment==2){
                ofstream ofs(analysisDirName,ios::app|ios::out);//append way
                if(!ofs.is_open()) cout<<"Cannot open "<<analysisDirName<<endl;
                ofs.setf(ios::fixed);
                ofs.precision(4);
                ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
                ofs.close();
            }else{//update model

            }
        }else if(index_model==2){
            //load original graph
            WGraph wgraph;
            wgraph.load_graph(graphFileName,graph_type);
            numOfOriginalVertices=numOfVertices;
            numOfOriginalEdges=numOfEdges;
            //initialize variables
            _integrated_index.bindex_.resize(numOfOriginalVertices);
            _isDeleted.resize(numOfOriginalVertices,0);
            load_hfpoint(HFPoint_file,hfRate);
            //deleted, overlay label and h2h-index construction
            delete_and_pushPLL_forest_directed(graphFileName,outputDirName,analysisDirName,calcCoef,hfRate,ordering_model,thresholdDegree,bestDegree,isExperiment,graph_type,givenOrderFileName);
            _labeling_time=time_tree_index+time_overlay_index;
            _ordering_time=time_overlay_order;
            _avg_labeling_time=_labeling_time/(double)numOfVertices;
            std::cout<<"time_tree_index = "<<time_tree_index*1e6<<std::endl;
            std::cout<<"time_overlay_index = "<<time_overlay_index*1e6<<" time_overlay_order = "<<time_overlay_order*1e6<<std::endl;
            std::cout<<"_labeling_time = "<<_labeling_time*1e6<<" _avg_labeling_time = "<<_avg_labeling_time*1e6<<std::endl;
            //ouput isDeletedFile
            save_isDeleted(outputDirName);
            //save indexFile to binary file
            _integrated_index.save_labels_2_directed(labelFileName);
            //save label size file
            if(*labelSizeFileName!='\0') _integrated_index.save_label_size_directed(labelSizeFileName);
            if(isExperiment==0){
                //output analysis
                save_analysis_result_directed(analysisDirName,queryPair_file);
                // //save different labels for debug
                // _integrated_index.save_labels_differ_analysis_1(analysisDirName);
                //save order
                _integrated_index.save_rank_directed(analysisDirName);
                // //save label size by order
                // _integrated_index.save_label_size_byOrder(analysisDirName);
                //output index list txt file for debug
                _integrated_index.write_labels_1_directed(analysisDirName);
                //output overllay_graph
                overlay_wgraph.save_overlay_graph(outputDirName,_newToOriginal,graph_type);
            }else if(isExperiment==1){
                //appendt analysis
                append_analysis_result_directed(analysisDirName,queryPair_file);
            }else if(isExperiment==2){
                ofstream ofs(analysisDirName,ios::app|ios::out);//append way
                if(!ofs.is_open()) cout<<"Cannot open "<<analysisDirName<<endl;
                ofs.setf(ios::fixed);
                ofs.precision(4);
                ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
                ofs.close();
            }else{//update model

            }
        }
    }

    /**
     * @description: Multiple threads version
     * @param {int} hfRate
     * @return {*}
     * @author: Wan Jingyi
     */    
    Integrated_construction(char* graphFileName,char* labelFileName,char* labelSizeFileName,
                                                        char* outputDirName,char* analysisDirName,char* HFPoint_file,
                                                        int hfRate,int thresholdDegree,int bestDegree,
                                                        int index_model,int ordering_model,const Processing::calcCoefficient<double>& calcCoef,
                                                        int isExperiment,int num_threads,char* givenOrderFileName="",bool isSorted=false)
    {
        if(index_model==0){//the first integrated solution
            std::cout<<"the first solution doesn't implement parallel computation..."<<std::endl;
        }else if(index_model==1){//the second integrated solution
            //load original graph
            WGraph wgraph;
            wgraph.load_graph(graphFileName);
            numOfOriginalVertices=numOfVertices;
            numOfOriginalEdges=numOfEdges;
            std::cout<<"numOfOriginalVertices = "<<numOfOriginalVertices<<" numOfOriginalEdges = "<<numOfOriginalEdges<<endl;
            //initialize variables
            _integrated_index.index_.resize(numOfOriginalVertices);
            _isDeleted.resize(numOfOriginalVertices,0);
            load_hfpoint(HFPoint_file,hfRate);
            //deleted, overlay label and h2h-index construction
            delete_and_build_forest(graphFileName,outputDirName,analysisDirName,givenOrderFileName,calcCoef,hfRate,ordering_model,thresholdDegree,bestDegree,isExperiment,num_threads,isSorted);
            _labeling_time=time_tree_index+time_overlay_index;
            _ordering_time=time_overlay_order;
            _avg_labeling_time=_labeling_time/(double)numOfVertices;
            std::cout<<"time_tree_index = "<<time_tree_index*1e6<<std::endl;
            std::cout<<"time_overlay_index = "<<time_overlay_index*1e6<<" time_overlay_order = "<<time_overlay_order*1e6<<std::endl;
            std::cout<<"_labeling_time = "<<_labeling_time*1e6<<" _avg_labeling_time = "<<_avg_labeling_time*1e6<<std::endl;
            //ouput isDeletedFile
            save_isDeleted(outputDirName);
            //save indexFile to binary file
            _integrated_index.save_labels_1(labelFileName);
            //save label size file
            if(*labelSizeFileName!='\0') _integrated_index.save_label_size(labelSizeFileName);
            if(isExperiment==0){
                //output analysis
                save_analysis_result(analysisDirName);
                //save different labels for debug
                _integrated_index.save_labels_differ_analysis_1(analysisDirName);
                //save order
                _integrated_index.save_rank(analysisDirName);
                //save label size by order
                 _integrated_index.save_label_size_byOrder(analysisDirName);
                //output index list txt file for debug
                _integrated_index.write_labels_1(analysisDirName);
            }else if(isExperiment==1){
                //appendt analysis
                append_analysis_result(analysisDirName);
            }else if(isExperiment==2){
                ofstream ofs(analysisDirName,ios::app|ios::out);//append way
                if(!ofs.is_open()) cout<<"Cannot open "<<analysisDirName<<endl;
                ofs.setf(ios::fixed);
                ofs.precision(4);
                ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
                ofs.close();
            }else{//update model
                
            }

        }else if(index_model==2){//the third integrated solution
            //load original graph
            WGraph wgraph;
            wgraph.load_graph(graphFileName);
            numOfOriginalVertices=numOfVertices;
            numOfOriginalEdges=numOfEdges;
            std::cout<<"numOfOriginalVertices="<<numOfOriginalVertices<<" numOfOriginalEdges="<<numOfOriginalEdges<<endl;
            //initialize variables
            _integrated_index.index_.resize(numOfOriginalVertices);
            _isDeleted.resize(numOfOriginalVertices,0);
            load_hfpoint(HFPoint_file,hfRate);
            //deleted, overlay label and h2h-index construction
            delete_and_pushPLL_forest(graphFileName,outputDirName,analysisDirName,calcCoef,hfRate,ordering_model,thresholdDegree,bestDegree,isExperiment,num_threads,givenOrderFileName);
            _labeling_time=time_tree_index+time_overlay_index;
            _ordering_time=time_overlay_order;
            _avg_labeling_time=_labeling_time/(double)numOfVertices;
            std::cout<<"time_tree_index = "<<time_tree_index*1e6<<std::endl;
            std::cout<<"time_overlay_index = "<<time_overlay_index*1e6<<" time_overlay_order = "<<time_overlay_order*1e6<<std::endl;
            std::cout<<"_labeling_time = "<<_labeling_time*1e6<<" _avg_labeling_time = "<<_avg_labeling_time*1e6<<std::endl;
            //ouput isDeletedFile
            save_isDeleted(outputDirName);
            //save indexFile to binary file
            _integrated_index.save_labels_2(labelFileName);
            //save label size file
            if(*labelSizeFileName!='\0') _integrated_index.save_label_size(labelSizeFileName);
            if(isExperiment==0){
                 _integrated_index.save_labels_differ_analysis_2(analysisDirName);
                //output analysis
                save_analysis_result(analysisDirName);
                //save order
                _integrated_index.save_rank(analysisDirName);
                //save label size by order
                _integrated_index.save_label_size_byOrder(analysisDirName);
                //output index list txt file for debug
                _integrated_index.write_labels_1(analysisDirName);
            }else if(isExperiment==1){
                //appendt analysis
                append_analysis_result(analysisDirName);
            }else if(isExperiment==2){
                ofstream ofs(analysisDirName,ios::app|ios::out);//append way
                if(!ofs.is_open()) cout<<"Cannot open "<<analysisDirName<<endl;
                ofs.setf(ios::fixed);
                ofs.precision(4);
                ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
                ofs.close();
            }else{//update model
                
            }
        }

    }

    /**
     * @description: Multiple threads version,add for directed
     * @param {int} hfRate
     * @return {*}
     * @author: Wan Jingyi
     */    
    Integrated_construction(char* graphFileName,char* labelFileName,char* labelSizeFileName,
                                                        char* outputDirName,char* analysisDirName,char* HFPoint_file,
                                                        char* queryPair_file,int hfRate,int thresholdDegree,
                                                        int bestDegree,int index_model,int ordering_model,
                                                        const Processing::calcCoefficient<double>& calcCoef,int isExperiment,int num_threads,
                                                        int graph_type,char* givenOrderFileName="",bool isSorted=false)
    {
        if(index_model==0){
            std::cout<<"Multiplt threads: the first integrated-index solution..."<<std::endl;
        }else if(index_model==1){
            std::cout<<"Multiplt threads: the second integrated-index solution..."<<std::endl;
            WGraph wgraph;
            wgraph.load_graph(graphFileName,graph_type);
            numOfOriginalVertices=numOfVertices;
            numOfOriginalEdges=numOfEdges;
            //initialize variables
            _integrated_index.bindex_.resize(numOfOriginalVertices);
            _isDeleted.resize(numOfOriginalVertices,0);
            load_hfpoint(HFPoint_file,hfRate);
            //deleted, overlay label and h2h-index construction
            delete_and_build_forest_directed(graphFileName,outputDirName,analysisDirName,givenOrderFileName,calcCoef,hfRate,ordering_model,thresholdDegree,bestDegree,isExperiment,graph_type,num_threads,isSorted);
            //test rightness
            // Check_overlay_graph check_overlay_graph;
            // check_overlay_graph.check_left_vertices_SP(wgraph,overlay_wgraph,_newToOriginal,outputDirName,numOfOriginalVertices, numOfOverlayVertices,5);
            _labeling_time=time_tree_index+time_overlay_index;
            _ordering_time=time_overlay_order;
            _avg_labeling_time=_labeling_time/(double)numOfVertices;
            std::cout<<"time_tree_index = "<<time_tree_index*1e6<<std::endl;
            std::cout<<"time_overlay_index = "<<time_overlay_index*1e6<<" time_overlay_order = "<<time_overlay_order*1e6<<std::endl;
            std::cout<<"_labeling_time = "<<_labeling_time*1e6<<" _avg_labeling_time = "<<_avg_labeling_time*1e6<<std::endl;
            //ouput isDeletedFile
            save_isDeleted(outputDirName);
            //save indexFile to binary file
            _integrated_index.save_labels_1_directed(labelFileName);
            //save label size file
            if(*labelSizeFileName!='\0') _integrated_index.save_label_size_directed(labelSizeFileName);
            if(isExperiment==0){
                //output analysis
                save_analysis_result_directed(analysisDirName,queryPair_file);
                // //save different labels for debug
                // _integrated_index.save_labels_differ_analysis_1(analysisDirName);
                //save order
                _integrated_index.save_rank_directed(analysisDirName);
                // //save label size by order
                // _integrated_index.save_label_size_byOrder(analysisDirName);
                //output index list txt file for debug
                _integrated_index.write_labels_1_directed(analysisDirName);
            }else if(isExperiment==1){
                //appendt analysis
                append_analysis_result_directed(analysisDirName,queryPair_file);
            }else if(isExperiment==2){
                ofstream ofs(analysisDirName,ios::app|ios::out);//append way
                if(!ofs.is_open()) cout<<"Cannot open "<<analysisDirName<<endl;
                ofs.setf(ios::fixed);
                ofs.precision(4);
                ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
                ofs.close();
            }else{//update model
                
            }
        }else if(index_model==2){
            std::cout<<"Multiplt threads: the third integrated-index solution..."<<std::endl;
            //load original graph
            WGraph wgraph;
            wgraph.load_graph(graphFileName,graph_type);
            numOfOriginalVertices=numOfVertices;
            numOfOriginalEdges=numOfEdges;
            //initialize variables
            _integrated_index.bindex_.resize(numOfOriginalVertices);
            _isDeleted.resize(numOfOriginalVertices,0);
            load_hfpoint(HFPoint_file,hfRate);
            //deleted, overlay label and h2h-index construction
            delete_and_pushPLL_forest_directed(graphFileName,outputDirName,analysisDirName,calcCoef,hfRate,ordering_model,thresholdDegree,bestDegree,isExperiment,num_threads,graph_type,givenOrderFileName);
            _labeling_time=time_tree_index+time_overlay_index;
            _ordering_time=time_overlay_order;
            _avg_labeling_time=_labeling_time/(double)numOfVertices;
            std::cout<<"time_tree_index = "<<time_tree_index*1e6<<std::endl;
            std::cout<<"time_overlay_index = "<<time_overlay_index*1e6<<" time_overlay_order = "<<time_overlay_order*1e6<<std::endl;
            std::cout<<"_labeling_time = "<<_labeling_time*1e6<<" _avg_labeling_time = "<<_avg_labeling_time*1e6<<std::endl;
            //ouput isDeletedFile
            save_isDeleted(outputDirName);
            //save indexFile to binary file
            _integrated_index.save_labels_2_directed(labelFileName);
            //save label size file
            if(*labelSizeFileName!='\0') _integrated_index.save_label_size_directed(labelSizeFileName);
            if(isExperiment==0){
                //output analysis
                save_analysis_result_directed(analysisDirName,queryPair_file);
                // //save different labels for debug
                // _integrated_index.save_labels_differ_analysis_1(analysisDirName);
                //save order
                _integrated_index.save_rank_directed(analysisDirName);
                // //save label size by order
                // _integrated_index.save_label_size_byOrder(analysisDirName);
                //output index list txt file for debug
                _integrated_index.write_labels_1_directed(analysisDirName);
            }else if(isExperiment==1){
                //appendt analysis
                append_analysis_result_directed(analysisDirName,queryPair_file);
            }else if(isExperiment==2){
                ofstream ofs(analysisDirName,ios::app|ios::out);//append way
                if(!ofs.is_open()) cout<<"Cannot open "<<analysisDirName<<endl;
                ofs.setf(ios::fixed);
                ofs.precision(4);
                ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
                ofs.close();
            }else{//update model
                
            }
        }
    }

    ~Integrated_construction(){}

    void get_time(double& labeling_time,double& ordering_time,double& rmq_time){
        labeling_time=_labeling_time;
        ordering_time=_ordering_time;
        rmq_time=time_compute_rmq;
    }

    // void get_label_size(vector<double>& node_size){
    //     for (NodeID v = 0; v < numOfVertices; ++v) {
    //         node_size[v]=_integrated_index.index_[v].siz;
    //     }
    //     return;
    // }

    void get_label_size(vector<double>& node_size){
        for (NodeID v = 0; v < numOfVertices; ++v) {
            if(DIRECTED_FLAG==true){
                node_size[v]=_integrated_index.bindex_[v].siz+_integrated_index.bindex_[v].r_siz;
            }else{
                node_size[v]=_integrated_index.index_[v].siz;
            }
        }
        return;
    }

protected:
    //*****************build functions**************
    /**
     * @description: single thread used to build forest for the third solution
     * @param {const} Processing
     * @param {int} hfRate
     * @param {int} ordering_model
     * @param {int} thresholdDegree
     * @param {int} bestDegree
     * @param {int} isExperiment
     * @return {*}
     * @author: Wan Jingyi
     */    
    void delete_and_pushPLL_forest(char* graphFileName,char* outputDirName,char* analysisDirName ,const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,int thresholdDegree,int bestDegree,int isExperiment,char* givenOrderFileName=""){
        H2H_construction h2h_construction(overlay_wgraph,graphFileName,HFPointFlag,
                                                                                    _isDeleted,_newToOriginal,_originalToNew,
                                                                                    numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,
                                                                                    thresholdDegree,bestDegree);
        _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
        std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" numOfOverlayEdges = "<<numOfOverlayEdges<<endl;
        //construction overlay labels
        load_hfpoint_overlay(numOfVertices);
        //build overlay index
        build_overlay_index(calcCoef,hfRate,ordering_model,givenOrderFileName);
        numOfVertices=numOfOriginalVertices;
        numOfEdges=numOfOriginalEdges;
        //push forest(including merge to integrated index)
        h2h_construction.push_forest(_integrated_index,time_tree_index,time_compute_rmq);
        //save lca variables
        h2h_construction.saveLcaVariables_forest(outputDirName);
        if(isExperiment==0){
            //output forest order
            h2h_construction.outputOrder(analysisDirName);
            //output structure information for debug
            h2h_construction.outputForestStructure(analysisDirName);
            //ouput analysis file
            h2h_construction.outputAnalysis_forest(analysisDirName);
        }else if(isExperiment==1){
            //ouput analysis file
            h2h_construction.appendAnalysis_forest(analysisDirName);
        }else{

        }

    }

    /**
     * @description: single thread used to build forest for the third solution,add for directed
     * @param {const} Processing
     * @param {int} hfRate
     * @param {int} ordering_model
     * @param {int} thresholdDegree
     * @param {int} bestDegree
     * @param {int} isExperiment
     * @return {*}
     * @author: Wan Jingyi
     */    
    void delete_and_pushPLL_forest_directed(char* graphFileName,char* outputDirName,char* analysisDirName ,const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,int thresholdDegree,int bestDegree,int isExperiment,int graphType,char* givenOrderFileName=""){
        H2H_construction h2h_construction(overlay_wgraph,graphFileName,HFPointFlag,
                                                                                    _isDeleted,_newToOriginal,_originalToNew,
                                                                                    numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,
                                                                                    graphType,thresholdDegree,bestDegree);
        _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
        std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" numOfOverlayEdges = "<<numOfOverlayEdges<<endl;
        //construction overlay labels
        load_hfpoint_overlay(numOfVertices);
        //build overlay index
        build_overlay_index_directed(calcCoef,hfRate,ordering_model,givenOrderFileName,outputDirName);
        numOfVertices=numOfOriginalVertices;
        numOfEdges=numOfOriginalEdges;
        //push forest(including merge to integrated index)
        h2h_construction.push_forest_directed(_integrated_index,time_tree_index,time_compute_rmq);
        //save lca variables
        h2h_construction.saveLcaVariables_forest(outputDirName);
        if(isExperiment==0){
            //output forest order
            h2h_construction.outputOrder(analysisDirName);
            //output structure information for debug
            h2h_construction.outputForestStructure(analysisDirName);
            //ouput analysis file
            h2h_construction.outputAnalysis_forest(analysisDirName);
        }else if(isExperiment==1){
            //ouput analysis file
            h2h_construction.appendAnalysis_forest(analysisDirName);
        }else{

        }

    }

    //multiThreads version
    void delete_and_pushPLL_forest(char* graphFileName,char* outputDirName,char* analysisDirName ,const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,int thresholdDegree,int bestDegree,int isExperiment,int num_threads,char* givenOrderFileName=""){
        H2H_construction h2h_construction(overlay_wgraph,graphFileName,HFPointFlag,
                                                                                    _isDeleted,_newToOriginal,_originalToNew,
                                                                                    numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,
                                                                                    thresholdDegree,bestDegree);
        //h2h_construction.outputForestStructure(analysisDirName);//to be deleted
        _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
        std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" numOfOverlayEdges = "<<numOfOverlayEdges<<endl;
        //construction overlay labels
        load_hfpoint_overlay(numOfVertices);
        //build overlay index
        build_overlay_index(calcCoef,hfRate,ordering_model,givenOrderFileName);
        numOfVertices=numOfOriginalVertices;
        numOfEdges=numOfOriginalEdges;
        //push forest(including merge to integrated index)
        h2h_construction.push_forest(_integrated_index,time_tree_index,time_compute_rmq,num_threads);
        //save lca variables
        h2h_construction.saveLcaVariables_forest(outputDirName);
        if(isExperiment==0){
            //output forest order
            h2h_construction.outputOrder(analysisDirName);
            //output structure information for debug
            h2h_construction.outputForestStructure(analysisDirName);
            //ouput analysis file
            h2h_construction.outputAnalysis_forest(analysisDirName);
        }else if(isExperiment==1){
            //ouput analysis file
            h2h_construction.appendAnalysis_forest(analysisDirName);
        }else{
        }
    }

    //multiThreads version for the third solution, add for directed
    void delete_and_pushPLL_forest_directed(char* graphFileName,char* outputDirName,char* analysisDirName ,const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,int thresholdDegree,int bestDegree,int isExperiment,int num_threads,int graphType,char* givenOrderFileName=""){
        H2H_construction h2h_construction(overlay_wgraph,graphFileName,HFPointFlag,
                                                                                    _isDeleted,_newToOriginal,_originalToNew,
                                                                                    numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,
                                                                                    graphType,thresholdDegree,bestDegree);
        //h2h_construction.outputForestStructure(analysisDirName);//to be deleted
        _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
        std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" numOfOverlayEdges = "<<numOfOverlayEdges<<endl;
        //construction overlay labels
        load_hfpoint_overlay(numOfVertices);
        //build overlay index
        build_overlay_index_directed(calcCoef,hfRate,ordering_model,givenOrderFileName,outputDirName);
        numOfVertices=numOfOriginalVertices;
        numOfEdges=numOfOriginalEdges;
        //push forest(including merge to integrated index)
        h2h_construction.push_forest_directed(_integrated_index,time_tree_index,time_compute_rmq,num_threads);
        //save lca variables
        h2h_construction.saveLcaVariables_forest(outputDirName);
        if(isExperiment==0){
            //output forest order
            h2h_construction.outputOrder(analysisDirName);
            //output structure information for debug
            h2h_construction.outputForestStructure(analysisDirName);
            //ouput analysis file
            h2h_construction.outputAnalysis_forest(analysisDirName);
        }else if(isExperiment==1){
            //ouput analysis file
            h2h_construction.appendAnalysis_forest(analysisDirName);
        }else{

        }

    }

    //Single thread for the second solution (undirected)
    void delete_and_build_forest(char* graphFileName,char* outputDirName,char* analysisDirName,const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,int thresholdDegree,int bestDegree,int isExperiment,char* givenOrderFileName="",bool isSorted=false){
        H2H_construction h2h_construction(graphFileName,overlay_wgraph,HFPointFlag,
                                                                                    _isDeleted,_newToOriginal,_originalToNew,
                                                                                    numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,
                                                                                    thresholdDegree,queryTime,bestDegree,
                                                                                    isSorted);
        _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
        std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" numOfOverlayEdges = "<<numOfOverlayEdges<<endl;
        //construction overlay labels
        load_hfpoint_overlay(numOfVertices);
        //build overlay index
        build_overlay_index(calcCoef,hfRate,ordering_model,givenOrderFileName);
        numOfVertices=numOfOriginalVertices;
        numOfEdges=numOfOriginalEdges;
        //compute forest(including merge to integrated index)
        h2h_construction.compute_forest(_integrated_index,time_tree_index,time_compute_rmq);
        //save lca variables
        h2h_construction.saveLcaVariables_forest(outputDirName);
        if(isExperiment==0){
            //output forest order
            h2h_construction.outputOrder(analysisDirName);
            //output structure information for debug
            h2h_construction.outputForestStructure(analysisDirName);
            //ouput analysis file
            h2h_construction.outputAnalysis_forest(analysisDirName);
        }else if(isExperiment==1){
            //ouput analysis file
            h2h_construction.appendAnalysis_forest(analysisDirName);
        }else{

        }
    }

    //Add for directed single thread
    void delete_and_build_forest_directed(char* graphFileName,char* outputDirName,char* analysisDirName,const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,int thresholdDegree,int bestDegree,int isExperiment,int graphType,char* givenOrderFileName="",bool isSorted=false)
    {
      
        H2H_construction h2h_construction(graphFileName,overlay_wgraph,HFPointFlag,
                                                                                    _isDeleted,_newToOriginal,_originalToNew,
                                                                                    numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,
                                                                                    thresholdDegree,queryTime,bestDegree,
                                                                                    graphType,isSorted);
        _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
        //build overlay index
        load_hfpoint_overlay(numOfOverlayVertices);
        build_overlay_index_directed(calcCoef,hfRate,ordering_model,givenOrderFileName,outputDirName);
        numOfVertices=numOfOriginalVertices;
        numOfEdges=numOfOriginalEdges;
        //compute forest(including merge to integrated index)
        h2h_construction.compute_forest_directed(_integrated_index,time_tree_index,time_compute_rmq);
        //save lca variables
        h2h_construction.saveLcaVariables_forest(outputDirName);
        if(isExperiment==0){
            //output forest order
            h2h_construction.outputOrder(analysisDirName);
            //output structure information for debug
            h2h_construction.outputForestStructure(analysisDirName);
            //ouput analysis file
            h2h_construction.outputAnalysis_forest(analysisDirName);
        }else if(isExperiment==1){
            //ouput analysis file
            h2h_construction.appendAnalysis_forest(analysisDirName);
        }
    }

    //multiple threads version
    void delete_and_build_forest(char* graphFileName,char* outputDirName,char* analysisDirName,char* givenOrderFileName,const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,int thresholdDegree,int bestDegree,int isExperiment,int num_threads,bool isSorted=false){
        H2H_construction h2h_construction(graphFileName,overlay_wgraph,HFPointFlag,
                                                                                    _isDeleted,_newToOriginal,_originalToNew,
                                                                                    numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,
                                                                                    thresholdDegree,queryTime,bestDegree,
                                                                                    isSorted);
        //h2h_construction.outputForestStructure(analysisDirName);//to be deleted
        _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
        std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" numOfOverlayEdges = "<<numOfOverlayEdges<<endl;
        //construction overlay labels
        load_hfpoint_overlay(numOfVertices);
        //build overlay index
        build_overlay_index(calcCoef,hfRate,ordering_model,givenOrderFileName);
        numOfVertices=numOfOriginalVertices;
        numOfEdges=numOfOriginalEdges;
        //compute forest(including merge to integrated index)
        h2h_construction.compute_forest(_integrated_index,time_tree_index,time_compute_rmq,num_threads);
        //save lca variables
        h2h_construction.saveLcaVariables_forest(outputDirName);
        if(isExperiment==0){
            //output forest order
            h2h_construction.outputOrder(analysisDirName);
            //output structure information for debug
            h2h_construction.outputForestStructure(analysisDirName);
            //ouput analysis file
            h2h_construction.outputAnalysis_forest(analysisDirName);
        }else if(isExperiment==1){
            //ouput analysis file
            h2h_construction.appendAnalysis_forest(analysisDirName);
        }else{

        }
    }

    /**
     * @description: Multiple threads for the second solution,add for directed
     * @param {const} Processing
     * @param {int} hfRate
     * @param {int} ordering_model
     * @param {int} thresholdDegree
     * @param {int} bestDegree
     * @param {int} isExperiment
     * @param {int} num_threads
     * @return {*}
     * @author: Wan Jingyi
     */    
    void delete_and_build_forest_directed(char* graphFileName,char* outputDirName,char* analysisDirName,char* givenOrderFileName,const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,int thresholdDegree,int bestDegree,int isExperiment,int graphType,int num_threads,bool isSorted=false){
        H2H_construction h2h_construction(graphFileName,overlay_wgraph,HFPointFlag,
                                                                                    _isDeleted,_newToOriginal,_originalToNew,
                                                                                    numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,
                                                                                    thresholdDegree,queryTime,bestDegree,
                                                                                    graphType,isSorted);
        _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
        //construction overlay labels
        load_hfpoint_overlay(numOfVertices);
        //build overlay index
        build_overlay_index_directed(calcCoef,hfRate,ordering_model,givenOrderFileName,outputDirName);
        numOfVertices=numOfOriginalVertices;
        numOfEdges=numOfOriginalEdges;
        //compute forest(including merge to integrated index)
        h2h_construction.compute_forest_directed(_integrated_index,time_tree_index,time_compute_rmq,num_threads);
        //save lca variables
        h2h_construction.saveLcaVariables_forest(outputDirName);
        if(isExperiment==0){
            //output forest order
            h2h_construction.outputOrder(analysisDirName);
            //output structure information for debug
            h2h_construction.outputForestStructure(analysisDirName);
            //ouput analysis file
            h2h_construction.outputAnalysis_forest(analysisDirName);
        }else if(isExperiment==1){
            //ouput analysis file
            h2h_construction.appendAnalysis_forest(analysisDirName);
        }
    }

    /**
     * @description: build h2h-index for the  first solution
     * @param {int} thresholdDegree
     * @param {bool} isExperiment
     * @return {*}
     * @author: Wan Jingyi
     */    
    void build_h2h_index(char* graphFileName,char* outputDirName,char* analysisDirName,int thresholdDegree,bool isExperiment){
        H2H_construction h2h_construction(_integrated_index.index_,overlay_wgraph,graphFileName,
                                                                                    HFPointFlag,time_tree_index,time_compute_rmq,
                                                                                    _isDeleted,_newToOriginal,_originalToNew,
                                                                                    numOfOverlayVertices,numOfDeletedVertices,numOfOverlayEdges,
                                                                                    thresholdDegree);
        _integrated_index._numOfOverlayVertices=numOfOverlayVertices;
        std::cout<<"numOfOverlayVertices = "<<numOfOverlayVertices<<" numOfOverlayEdges = "<<numOfOverlayEdges<<endl;
        //save  lca-variables(".lca") to binary file
        h2h_construction.saveLcaVariables(outputDirName);
        if(isExperiment){

        }else{
            //output tree order
            //h2h_construction.outputOrder(analysisDirName);
            //output information for debug
            //h2h_construction.outputTreeStructure(analysisDirName);
            //analyze tree
            //h2h_construction._dp_tree.analysis_and_save(analysisDirName,queryTime,HFPointFlag,numOfHfpoint,max_query_time,total_query_time);
            h2h_construction._dp_tree.compute_tree_performance();
            //ouput analysis file
            h2h_construction.outputAnalysis(analysisDirName);
        }
        std::cout<<"********************Build h2h index successfully!**************"<<std::endl;
        return;
    }

    void build_overlay_index(const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,char* givenOrderFileName=""){
        if(ordering_model==0){//lhp
            double  _order_time,_hub_time;
            Betweenness_Ordering betweenness_ordering(16, beta, overlay_wgraph, numOfVertices,_order_time,_hub_time,betweenessValues);
            time_overlay_order=GetCurrentTimeSec();
            Linear_Ordering<double> linear_ordering(overlay_wgraph,calcCoef,queryTime_overlay,betweenessValues,coverageValues,hfRate);
            time_overlay_order=GetCurrentTimeSec()-time_overlay_order;
            time_overlay_order+=_order_time;
            time_overlay_index=GetCurrentTimeSec();
            PL_W pl_w(overlay_wgraph, linear_ordering); 
            time_overlay_index=GetCurrentTimeSec()-time_overlay_index;
            //merge to integrated_index
            merge_overlay_labels(pl_w.labels.index_,linear_ordering.rank);
        }else if(ordering_model==1){//dhp
            time_overlay_order=GetCurrentTimeSec();
            Degree_Ordering degree_order(overlay_wgraph);
            time_overlay_order=GetCurrentTimeSec()-time_overlay_order;
            time_overlay_index=GetCurrentTimeSec();
            PL_W pl_w(overlay_wgraph, degree_order);
            time_overlay_index=GetCurrentTimeSec()-time_overlay_index;
            //merge to integrated_index
            merge_overlay_labels(pl_w.labels.index_,degree_order.rank);
        }else if(ordering_model==2){//given order
            time_overlay_order=GetCurrentTimeSec();
            Given_Ordering given_order(givenOrderFileName, overlay_wgraph,_originalToNew);
            time_overlay_order=GetCurrentTimeSec()-time_overlay_order;
            time_overlay_index=GetCurrentTimeSec();
            PL_W pl_w(overlay_wgraph, given_order);
            time_overlay_index=GetCurrentTimeSec()-time_overlay_index;
            //merge to integrated_index
            merge_overlay_labels(pl_w.labels.index_,given_order.rank);
        }
        //cout<<"numOfVertices="<<numOfVertices<<endl;//for debug
        std::cout<<"********************Build overlay index successfully!**************"<<std::endl;
        return;
    }


    void build_overlay_index_directed(const Processing::calcCoefficient<double>& calcCoef,int hfRate,int ordering_model,char* givenOrderFileName="",char* outputDirName=""){
        if(ordering_model==0){//lhp
            double  _order_time,_hub_time;
            Betweenness_Ordering betweenness_ordering(16, beta, overlay_wgraph, numOfVertices,_order_time,_hub_time,betweenessValues);
            time_overlay_order=GetCurrentTimeSec();
            Linear_Ordering<double> linear_ordering(overlay_wgraph,calcCoef,queryTime_overlay,betweenessValues,coverageValues,hfRate);
            time_overlay_order=GetCurrentTimeSec()-time_overlay_order;
            time_overlay_order+=_order_time;
            time_overlay_index=GetCurrentTimeSec();
            PL_W pl_w(overlay_wgraph, linear_ordering, DIRECTED_FLAG);
            time_overlay_index=GetCurrentTimeSec()-time_overlay_index;
            //test rightness (to be deleted)
            // string label_file(outputDirName);
            // label_file.append("_overlay.label");
            // pl_w.dlabels.save_labels(label_file.c_str());
            //merge to integrated_index
            merge_overlay_labels_directed(pl_w.dlabels.index_,pl_w.dlabels.bindex_,linear_ordering.rank);
        }else if(ordering_model==1){//dhp
            time_overlay_order=GetCurrentTimeSec();
            Degree_Ordering degree_order(overlay_wgraph);
            time_overlay_order=GetCurrentTimeSec()-time_overlay_order;
            time_overlay_index=GetCurrentTimeSec();
            PL_W pl_w(overlay_wgraph, degree_order, DIRECTED_FLAG);
            time_overlay_index=GetCurrentTimeSec()-time_overlay_index;
            //merge to integrated_index
            merge_overlay_labels_directed(pl_w.dlabels.index_,pl_w.dlabels.bindex_,degree_order.rank);
        }else if(ordering_model==2){//given order
            time_overlay_order=GetCurrentTimeSec();
            Given_Ordering given_order(givenOrderFileName, overlay_wgraph,_originalToNew);
            time_overlay_order=GetCurrentTimeSec()-time_overlay_order;
            time_overlay_index=GetCurrentTimeSec();
            PL_W pl_w(overlay_wgraph, given_order, DIRECTED_FLAG);
            time_overlay_index=GetCurrentTimeSec()-time_overlay_index;
            //merge to integrated_index
            merge_overlay_labels_directed(pl_w.dlabels.index_,pl_w.dlabels.bindex_,given_order.rank);
        }
        std::cout<<"********************Build overlay index successfully!**************"<<std::endl;
        return;
    }

    //*****************useful functions*******************
    void load_hfpoint_slice(char* HFPoint_file,int hfRate){
        int numOfHfpoint_slice;
        if(hfRate==0){
            numOfHfpoint_slice=numOfOriginalVertices;
        }else{
            numOfHfpoint_slice= static_cast<int> ( (double)numOfOriginalVertices*hfRate/(double)HF_DIVIDION);
        }
        if(numOfHfpoint_slice<=0) cout<<"error:numOfHfpoint_slice<=0"<<std::endl;
        std::ifstream in(HFPoint_file);//input HFPoint file to ifstream
        if(!in.is_open()) {std::cerr<<"Cannot open "<<HFPoint_file<<std::endl;}
        int id;int t;
        int cnt_hf=0;
        char line[24];
        while (in.getline(line,sizeof(line))){
            if(cnt_hf>=numOfHfpoint_slice) break;
            std::stringstream hp(line);
            hp>>id>>t;		
            queryTime[id]+=t;
            ++cnt_hf;
        }
        if(cnt_hf<numOfHfpoint_slice){
            numOfHfpoint_slice=cnt_hf;
        }    
        //std::cout<<"numOfHfpoint_slice = "<<numOfHfpoint_slice<<std::endl;
        return;
    }

    void load_hfpoint(char* HFPoint_file,int hfRate){
        numOfHfpoint= 0;//first line is the number of HFpoints
        HFPointFlag.resize(numOfOriginalVertices,0);
        queryTime.resize(numOfOriginalVertices,0);
        min_query_time=INF_WEIGHT;
        max_query_time=0;
        total_query_time=0;
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
        for(NodeID id=0;id<numOfOriginalVertices;++id){
            if(queryTime[id]!=0){
                HFPointFlag[id]=true;
                ++numOfHfpoint;
                total_query_time+=queryTime[id];                                                                                                                                                                      
                if(queryTime[id]>max_query_time) max_query_time=queryTime[id];
                if(queryTime[id]<min_query_time) min_query_time=queryTime[id];
                HFPoint.push_back(id);
            }
        }
        std::cout<<"numOfHfpoint = "<<numOfHfpoint<<std::endl;
        std::cout<<"max_query_time = "<<max_query_time<<std::endl;
        std::cout<<"min_query_time = "<<min_query_time<<std::endl;
        std::cout<<"toal_query_time = "<<total_query_time<<std::endl;
        std::cout<<"*******************Load hfpoint successfully!*******************"<<std::endl;
    }

    void load_hfpoint_overlay(int num){
        queryTime_overlay.resize(num,0);
        for(NodeID i=0;i<numOfOriginalVertices;++i){
            int v=_originalToNew[i];
            if(v!=-1) queryTime_overlay[v]=queryTime[i];
        }
        std::cout<<"**************Load hfpoint overlay successfully!*******************"<<std::endl;
    }

    void merge_overlay_labels(const vector<index_t1>& overlay_index,const vector<NodeID>& overlay_rank){
        vector<integrated_index_t>& index_=_integrated_index.index_;
        for(NodeID id=0;id<numOfVertices;++id){
            int k=overlay_index[id].spt_v.size();
            NodeID v=_newToOriginal[id];
            if(_isDeleted[v]) continue;
            index_[v].spt_v.resize(k);
            index_[v].spt_d.resize(k);
            for(NodeID i=0;i<k-1;++i){
                NodeID spt_v_t=overlay_index[id].spt_v[i];
                index_[v].spt_v[i]=spt_v_t;
                index_[v].spt_d[i]=overlay_index[id].spt_d[i];
                if(spt_v_t>_integrated_index._max_overlay_hub) _integrated_index._max_overlay_hub=spt_v_t;
            }
            index_[v].spt_v[k-1]=numOfOriginalVertices;
            index_[v].spt_d[k-1]=INF_WEIGHT;
            //size and level
            index_[v].siz=k-1;
            index_[v].level=overlay_rank[id];
        }
        cout<<"_max_overlay_hub="<<_integrated_index._max_overlay_hub<<endl;
        std::cout<<"*******************Merge overlay labels successfully!*******************"<<std::endl;
    }

    void merge_overlay_labels_directed(const vector<index_t1>& overlay_index,const vector<index_t1>& overlay_bindex,const vector<NodeID>& overlay_rank){
        vector<integrated_bindex_t>& bindex_=_integrated_index.bindex_;
        for(NodeID id=0;id<numOfVertices;++id){
            NodeID v=_newToOriginal[id];
            if(_isDeleted[v]) continue;
            //*********************forward*********************
            int k=overlay_index[id].spt_v.size();
            bindex_[v].spt_v.resize(k);
            bindex_[v].spt_d.resize(k);
            for(NodeID i=0;i<k-1;++i){
                NodeID spt_v_t=overlay_index[id].spt_v[i];
                bindex_[v].spt_v[i]=spt_v_t;
                bindex_[v].spt_d[i]=overlay_index[id].spt_d[i];
                if(spt_v_t>_integrated_index._max_overlay_hub) _integrated_index._max_overlay_hub=spt_v_t;
            }
            bindex_[v].spt_v[k-1]=numOfOriginalVertices;
            bindex_[v].spt_d[k-1]=INF_WEIGHT;
            //size
            bindex_[v].siz=k-1;
            //*********************backward*********************
            k=overlay_bindex[id].spt_v.size();
            bindex_[v].r_spt_v.resize(k);
            bindex_[v].r_spt_d.resize(k);
            for(NodeID i=0;i<k-1;++i){
                NodeID spt_v_t=overlay_bindex[id].spt_v[i];
                bindex_[v].r_spt_v[i]=spt_v_t;
                bindex_[v].r_spt_d[i]=overlay_bindex[id].spt_d[i];
                if(spt_v_t>_integrated_index._max_overlay_hub_r) _integrated_index._max_overlay_hub_r=spt_v_t;
            }
            bindex_[v].r_spt_v[k-1]=numOfOriginalVertices;
            bindex_[v].r_spt_d[k-1]=INF_WEIGHT;
            //r_size 
            bindex_[v].r_siz=k-1;
            //level
            bindex_[v].level=overlay_rank[id];
        }
        cout<<"_max_overlay_hub="<<_integrated_index._max_overlay_hub<<endl;
        cout<<"_max_overlay_hub_r="<<_integrated_index._max_overlay_hub_r<<endl;
        std::cout<<"*******************Merge overlay labels successfully!*******************"<<std::endl;
    }

    void save_isDeleted(char* outputDirName){
        string outputIsDeleted_file(outputDirName);
        outputIsDeleted_file.append(".isDeleted");
        ofstream ofsIsDeleted(outputIsDeleted_file);
        if(!ofsIsDeleted.is_open()) std::cerr<<outputIsDeleted_file<<" cannot be opened!"<<std::endl;
        ofsIsDeleted<<numOfOriginalVertices<<" "<<numOfOverlayVertices<<std::endl;
        for (int i = 0; i < numOfOriginalVertices; i++) {
            ofsIsDeleted <<i<<" "<< _isDeleted[i] << std::endl;
        }
        ofsIsDeleted.close();
        std::cout<<"*******************write isDeleted successfully!*******************"<<std::endl;
    }

    //all analysis
    void save_analysis_result(char* analysisDirName){
        //variables 
		double total_sum_size=0,total_ave_size=0;
		double hf_sum_size=0.0,hf_ave_size=0.0;
        double overlay_sum_size=0.0,overlay_ave_size=0.0,tree_sum_size=0.0,tree_ave_size=0.0;
		double total_performance=0,total_performance_s=0,ratio;
		for(int v=0;v<numOfOriginalVertices;++v){
			double label_size=_integrated_index.index_[v].siz;
			total_sum_size+=label_size;
			if(HFPointFlag[v]) hf_sum_size+=label_size;
			ratio=(double)queryTime[v]/((double)(max_query_time*DIVISION_FACTOR));
			total_performance+=ratio*label_size;
            if(_isDeleted[v]) tree_sum_size+=label_size;
            else overlay_sum_size+=label_size;
		}
		total_ave_size= (double) total_sum_size/(double) numOfOriginalVertices;
		hf_ave_size= hf_sum_size/(double) numOfHfpoint;
        overlay_ave_size=overlay_sum_size/(double)numOfOverlayVertices;
        tree_ave_size=tree_sum_size/(double)numOfDeletedVertices;
		total_performance_s=total_performance*((double)max_query_time);
		//output result
		cout<<"numOfOriginalVertices="<<numOfOriginalVertices<<" total_ave_size="<<total_ave_size<<endl;
		cout<<"numOfHFpoint="<<numOfHfpoint<<" hf_sum_size="<<hf_sum_size<<" hf_ave_size="<<hf_ave_size<<endl;
        cout<<"numOfOverlayVertices="<<numOfOverlayVertices<<" overlay_sum_size="<<overlay_sum_size<<" overlay_ave_size="<<overlay_ave_size<<endl;
        cout<<"numOfDeletedVertices="<<numOfDeletedVertices<<" tree_sum_size="<<tree_sum_size<<" tree_ave_size="<<tree_ave_size<<endl;
		cout<<"nomalization performance_result="<<total_performance<<" standard performance_result="<<total_performance_s<<endl;
        //write to file
        string analysisFileName(analysisDirName);
        analysisFileName.append(".analysis");
        ofstream ofs(analysisFileName);
        if(!ofs.is_open()) cout<<"Cannot open "<<analysisFileName<<endl;
        ofs.precision(4);
        ofs<<"_labeling_time="<<_labeling_time*1e6<<" _ordering_time="<<_ordering_time*1e6<<" time_compute_rmq="<<time_compute_rmq*1e6<<endl;
		ofs<<"numOfHFpoint="<<numOfHfpoint<<" hf_sum_size="<<hf_sum_size<<" hf_ave_size="<<hf_ave_size<<endl;
        ofs<<"numOfOverlayVertices="<<numOfOverlayVertices<<" overlay_sum_size="<<overlay_sum_size<<" overlay_ave_size="<<overlay_ave_size<<endl;
        ofs<<"numOfDeletedVertices="<<numOfDeletedVertices<<" tree_sum_size="<<tree_sum_size<<" tree_ave_size="<<tree_ave_size<<endl;
		ofs<<"nomalization performance_result="<<total_performance<<" standard performance_result="<<total_performance_s<<endl;
        ofs.close();
        std::cout<<"Save performance all successfully!"<<std::endl;//for debug
    }

    //all analysis for experiment model
    void append_analysis_result(char* analysisFileName){
        //variables 
		double total_sum_size=0,total_ave_size=0;
		double hf_sum_size=0.0,hf_ave_size=0.0;
        double overlay_sum_size=0.0,overlay_ave_size=0.0,tree_sum_size=0.0,tree_ave_size=0.0;
		double total_performance=0,total_performance_s=0,ratio;
		for(int v=0;v<numOfOriginalVertices;++v){
			double label_size=_integrated_index.index_[v].siz;
			total_sum_size+=label_size;
			if(HFPointFlag[v]) hf_sum_size+=label_size;
			ratio=(double)queryTime[v]/((double)(max_query_time*DIVISION_FACTOR));
			total_performance+=ratio*label_size;
            if(_isDeleted[v]) tree_sum_size+=label_size;
            else overlay_sum_size+=label_size;
		}
		total_ave_size= (double) total_sum_size/(double) numOfOriginalVertices;
		hf_ave_size= hf_sum_size/(double) numOfHfpoint;
        overlay_ave_size=overlay_sum_size/(double)numOfOverlayVertices;
        tree_ave_size=tree_sum_size/(double)numOfDeletedVertices;
		total_performance_s=total_performance*((double)max_query_time);
		//append to file
		ofstream ofs(analysisFileName,ios::app|ios::out);//append way
		if(!ofs.is_open()) cout<<"Cannot open "<<analysisFileName<<endl;
		ofs.setf(ios::fixed);
		ofs.precision(4);
        ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
		ofs<<total_sum_size<<" "<<total_ave_size<<" "<<hf_sum_size<<" "<<hf_ave_size<<" ";
        ofs<<overlay_sum_size<<" "<<overlay_ave_size<<" "<<tree_sum_size<<" "<<tree_ave_size<<" ";
        ofs<<total_performance_s<<" ";
		ofs.close();
        std::cout<<"Append performance oll successfully!"<<std::endl;//for debug
    }

    //all analysis
    void save_analysis_result_directed(char* analysisDirName,char* queryPairFileName){
        //load query time s and t
        load_workload_directed(queryPairFileName,queryTime_s,queryTime_t,numOfVertices);
        //variables 
		double total_sum_size=0,total_ave_size=0;
		double hf_sum_size=0.0,hf_ave_size=0.0;
        double overlay_sum_size=0.0,overlay_ave_size=0.0,tree_sum_size=0.0,tree_ave_size=0.0;
		double total_performance=0,ratio_s,ratio_t;
        double isize=0,r_isize=0;
		for(int v=0;v<numOfOriginalVertices;++v){
			isize=_integrated_index.bindex_[v].siz;
			r_isize=_integrated_index.bindex_[v].r_siz;
            total_sum_size+=isize+r_isize;
			if(HFPointFlag[v]) hf_sum_size+=isize+r_isize;
			ratio_s=(double)queryTime_s[v]/(double)DIVISION_FACTOR;
            ratio_t=(double)queryTime_t[v]/(double)DIVISION_FACTOR;
			total_performance+=ratio_s*isize+ratio_t*r_isize;
            if(_isDeleted[v]) tree_sum_size+=isize+r_isize;
            else overlay_sum_size+=isize+r_isize;
		}
		total_ave_size= (double) total_sum_size/(double) numOfOriginalVertices;
		hf_ave_size= hf_sum_size/(double) numOfHfpoint;
        overlay_ave_size=overlay_sum_size/(double)numOfOverlayVertices;
        tree_ave_size=tree_sum_size/(double)numOfDeletedVertices;
		//output result
		cout<<"numOfOriginalVertices="<<numOfOriginalVertices<<" total_ave_size="<<total_ave_size<<endl;
		cout<<"numOfHFpoint="<<numOfHfpoint<<" hf_sum_size="<<hf_sum_size<<" hf_ave_size="<<hf_ave_size<<endl;
        cout<<"numOfOverlayVertices="<<numOfOverlayVertices<<" overlay_sum_size="<<overlay_sum_size<<" overlay_ave_size="<<overlay_ave_size<<endl;
        cout<<"numOfDeletedVertices="<<numOfDeletedVertices<<" tree_sum_size="<<tree_sum_size<<" tree_ave_size="<<tree_ave_size<<endl;
		cout<<"standard performance_result="<<total_performance<< endl;
        //write to file
        string analysisFileName(analysisDirName);
        analysisFileName.append(".analysis");
        ofstream ofs(analysisFileName);
        if(!ofs.is_open()) cout<<"Cannot open "<<analysisFileName<<endl;
        ofs.precision(4);
        ofs<<"_labeling_time="<<_labeling_time*1e6<<" _ordering_time="<<_ordering_time*1e6<<" time_compute_rmq="<<time_compute_rmq*1e6<<endl;
		ofs<<"numOfHFpoint="<<numOfHfpoint<<" hf_sum_size="<<hf_sum_size<<" hf_ave_size="<<hf_ave_size<<endl;
        ofs<<"numOfOverlayVertices="<<numOfOverlayVertices<<" overlay_sum_size="<<overlay_sum_size<<" overlay_ave_size="<<overlay_ave_size<<endl;
        ofs<<"numOfDeletedVertices="<<numOfDeletedVertices<<" tree_sum_size="<<tree_sum_size<<" tree_ave_size="<<tree_ave_size<<endl;
		ofs<<"standard performance_result="<<total_performance<<endl;
        ofs.close();
        std::cout<<"Save performance all successfully!"<<std::endl;//for debug
    }

    //all analysis for experiment model
    void append_analysis_result_directed(char* analysisFileName,char* queryPairFileName ){
        //load query time s and t
        load_workload_directed(queryPairFileName,queryTime_s,queryTime_t,numOfVertices);
        //variables 
		double total_sum_size=0,total_ave_size=0;
		double hf_sum_size=0.0,hf_ave_size=0.0;
        double overlay_sum_size=0.0,overlay_ave_size=0.0,tree_sum_size=0.0,tree_ave_size=0.0;
		double total_performance=0,ratio_s,ratio_t;
        double isize=0,r_isize=0;
		for(int v=0;v<numOfOriginalVertices;++v){
			isize=_integrated_index.bindex_[v].siz;
			r_isize=_integrated_index.bindex_[v].r_siz;
            total_sum_size+=isize+r_isize;
			if(HFPointFlag[v]) hf_sum_size+=isize+r_isize;
			ratio_s=(double)queryTime_s[v]/(double)DIVISION_FACTOR;
            ratio_t=(double)queryTime_t[v]/(double)DIVISION_FACTOR;
			total_performance+=ratio_s*isize+ratio_t*r_isize;
            if(_isDeleted[v]) tree_sum_size+=isize+r_isize;
            else overlay_sum_size+=isize+r_isize;
		}
		total_ave_size= (double) total_sum_size/(double) numOfOriginalVertices;
		hf_ave_size= hf_sum_size/(double) numOfHfpoint;
        overlay_ave_size=overlay_sum_size/(double)numOfOverlayVertices;
        tree_ave_size=tree_sum_size/(double)numOfDeletedVertices;
		//append to file
		ofstream ofs(analysisFileName,ios::app|ios::out);//append way
		if(!ofs.is_open()) cout<<"Cannot open "<<analysisFileName<<endl;
		ofs.setf(ios::fixed);
		ofs.precision(4);
        ofs<<_labeling_time*1e6<<" "<<_ordering_time*1e6<<" "<<time_compute_rmq*1e6<<" ";
		ofs<<total_sum_size<<" "<<total_ave_size<<" "<<hf_sum_size<<" "<<hf_ave_size<<" ";
        ofs<<overlay_sum_size<<" "<<overlay_ave_size<<" "<<tree_sum_size<<" "<<tree_ave_size<<" ";
        ofs<<total_performance<<" ";
		ofs.close();
        std::cout<<"Append performance oll successfully!"<<std::endl;//for debug
    }

};

#endif