/*
 * @Descripttion: Construct pll in parallel paradigm
 * @version: 
 * @Author: wanjingyi
 * @Date: 2021-02-22 16:37:30
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-09-12 17:50:24
 */
#ifndef _PARALLEL_INDEXING_LABELING_H
#define _PARALLEL_INDEXING_LABELING_H

#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "../command.h"
#include "../src/paras.h"
#include "../src/graph.h"
#include "../src/time_util.h"
#include "../src/ordering.h"
#include "../src/parallel_construction.h"

using namespace time_util; 
using namespace PADO;

namespace command{
    class ParallelIndexingLabeling:public Command{
        public:
             void exit_with_help(){
                printf("Usage:\n");
                printf("-------------------------------------------------------------------\n");
                printf("bin/wcf_du_run -d [directedGraphFlag] -w [weightedGraphFlag] -g [graphFileName] [-o [orderingModel]] [-b [batchSize]] -a [outputDirName] [-n [numOfThreads]] \n");
                printf("Parameter explanation:\n");
                printf("-------------------------------------------------------------------\n");
                printf("[orderingModel] = 0: given order, need to input order file;1: degree order (default)\n");
                printf("[batchSize]: batch size for parallel, default is 512 \n");
                printf("[numOfThreads]: num of threads, default is 5\n");
                exit(1);
            }

            int main(int argc, char *argv[])
            {
                char graphFileName[255] = "";
                char inputOrderFileName[255] = "";
                char outputDirName[255] = "";
               // char queryFreqFileName[256]=""; //point-pair query time filename(node s,node t,queryTime)
                int t_directed_flag = 0;
                int t_weighted_flag = 0;
                int t_ordering_flag = 1;
                int t_vectorized_flag = 0;
                int t_multithread_flag = 0;
                int num_of_threads = 5;
                int batchSize=512;
        
                for(int i = 2; i < argc; i++){
                    if(argv[i][0] != '-') break;
                    if(++i >= argc) exit_with_help();
                    switch (argv[i-1][1]){
                        case 'g':
                            strcpy(graphFileName,argv[i]);
                            cout<<"graphFileName = "<<graphFileName<<endl;
                            break;
                        case 'd':
                            t_directed_flag = atoi(argv[i]);
                            break;
                        case 'w':
                            t_weighted_flag = atoi(argv[i]);
                            break;
                        case 'o':
                            t_ordering_flag = atoi(argv[i]);
                            break;
                        case 'b':
                            batchSize = atoi(argv[i]);
                            cout<<"batchSize="<<batchSize<<endl;
                            break;
                        case 'a':
                            strcpy(outputDirName,argv[i]);
                            cout<<"outputDirName = "<<outputDirName<<endl;
                            break;
                        case 'v':
                            t_vectorized_flag = atoi(argv[i]);
                            break;
                        case 'f':
                            strcpy(inputOrderFileName,argv[i]);
                            std::cout<<"inputOrderFileName="<<inputOrderFileName<<std::endl;
                            break;
                        case 'm':
                            t_multithread_flag= atoi(argv[i]);
                            break;
                        case 'n':
                            num_of_threads = atoi(argv[i]);
                            break;
                        default:
                            exit_with_help();           
                            break;
                    }
                }      

                if(*outputDirName=='\0'){
                    cout<<"outputDirName cannot be null!!!"<<std::endl;
                    exit_with_help();
                }

                if (t_directed_flag == 1) DIRECTED_FLAG = true;
                if (t_weighted_flag == 1) WEIGHTED_FLAG = true;    

                Graph graph;
                WGraph wgraph;
                if(WEIGHTED_FLAG == true){
                    wgraph.load_graph(graphFileName);
                }else{
                    graph.load_graph(graphFileName);
                }
                cout << numOfVertices << " nodes and " << numOfEdges << " arcs " << endl;             

                if (numOfVertices == 0 || numOfEdges == 0){
                    cout << "Corruptted graph file" << endl;
                    return 0;
                }

                if (DIRECTED_FLAG == true){
                    cout<<"DIRECTED_FLAG=true ";
                    if( WEIGHTED_FLAG == true){
                        cout<<"WEIGHTED_FLAG=true ";
                    }else{
                        cout<<"WEIGHTED_FLAG=false ";
                    }
                }else{//Undirected
                    cout<<"DIRECTED_FLAG=false ";
                    if( WEIGHTED_FLAG == true){//weighted
                        cout<<"WEIGHTED_FLAG=true ";
                        if(t_vectorized_flag==0){// No vectorization
                            cout<<"vectorization=false ";
                            if(t_multithread_flag==0){// Single Thread
                                cout<<"multithread=false "<<endl;
                                if(t_ordering_flag==0){//given order
                                    cout<<"Generate given order start!"<<endl;
                                    Given_Ordering given_ordering(inputOrderFileName,wgraph,true);
                                    cout<<"Generate given order successfully!"<<endl;
                                    //double _labeling_time=GetCurrentTimeSec();
                                    WeightedVertexCentricPLL *VCPLL = new WeightedVertexCentricPLL(wgraph,batchSize); // 512 is the batch size
                                    //_labeling_time=GetCurrentTimeSec()-_labeling_time;
                                    VCPLL->reorder_labels(given_ordering.rank);
                                    //cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                                    //double ave_labeling_time=_labeling_time/(double) numOfVertices;
                                    //cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                                    //save labels
                                    given_ordering.save_rank(outputDirName);//no backprefix for ordering class
                                    cout<<"save rank succesfully!"<<endl;
                                    VCPLL->save_labels(outputDirName);
                                    cout<<"save labels succesfully!"<<endl;
                                    VCPLL->save_label_size(outputDirName,given_ordering.inv);
                                    cout<<"save_label_size succesfully!"<<endl;
                                    VCPLL->write_labels(outputDirName,given_ordering.inv);
                                    cout<<"write_labels succesfully!"<<endl;
                                    delete VCPLL;
                                    return 0;
                                }else if (t_ordering_flag==1){//degree order
                                    cout<<"Generate degree order start!"<<endl;
                                    Degree_Ordering degree_order(wgraph,true);
                                    //to be deleted
                                    // cout<<"wgraph:"<<endl;
                                    // for(NodeID r=0;r<numOfVertices;++r){
                                    //     cout<<r<<":";
                                    //     for (EdgeID eid = wgraph.vertices[r]; eid < wgraph.vertices[r + 1]; ++eid) {
                                    //         cout<<wgraph.edges[eid].first<<"-"<<wgraph.edges[eid].second<<" ";
                                    //     }
                                    //     cout<<endl;
                                    // }
                                    cout<<"Generate degree order successfully!"<<endl;
                                    WeightedVertexCentricPLL *VCPLL = new WeightedVertexCentricPLL(wgraph,batchSize); // 512 is the batch size
                                    VCPLL->reorder_labels(degree_order.rank);
                                    //save labels
                                    degree_order.save_rank(outputDirName);//no backprefix for ordering class
                                    cout<<"save rank succesfully!"<<endl;
                                    VCPLL->save_labels(outputDirName);
                                    cout<<"save labels succesfully!"<<endl;
                                    VCPLL->save_label_size(outputDirName,degree_order.inv);
                                    cout<<"save_label_size succesfully!"<<endl;
                                    VCPLL->write_labels(outputDirName,degree_order.inv);
                                    cout<<"write_labels succesfully!"<<endl;
                                    delete VCPLL;
                                    return 0;
                                }else if(t_ordering_flag==2){//betweenness
                                    
                                }
                            }else{// Multithread
                                cout<<"multithread=true "<<endl;
                                if(t_ordering_flag==0){//given order
                                    cout<<"Generate given order start!"<<endl;
                                    Given_Ordering given_ordering(inputOrderFileName,wgraph,true);
                                    cout<<"Generate given order successfully!"<<endl;
                                    //omp_set_num_threads(num_of_threads);
                                    WeightedParaVertexCentricPLL *VCPLL =new WeightedParaVertexCentricPLL(wgraph,batchSize,num_of_threads);
                                    VCPLL->reorder_labels(given_ordering.rank);
                                    //save labels
                                    given_ordering.save_rank(outputDirName);//no backprefix for ordering class
                                    cout<<"save rank succesfully!"<<endl;
                                    VCPLL->save_labels(outputDirName);
                                    cout<<"save labels succesfully!"<<endl;
                                    VCPLL->save_label_size(outputDirName,given_ordering.inv);
                                    cout<<"save_label_size succesfully!"<<endl;
                                    VCPLL->write_labels(outputDirName,given_ordering.inv);
                                    cout<<"write_labels succesfully!"<<endl;
                                    delete VCPLL;
                                    return 0;
                                }else if (t_ordering_flag==1){//degree order
                                    cout<<"Generate degree order start!"<<endl;
                                    Degree_Ordering degree_order(wgraph,true);
                                    cout<<"Generate degree order successfully!"<<endl;
                                    //omp_set_num_threads(num_of_threads);
                                    WeightedParaVertexCentricPLL *VCPLL =new WeightedParaVertexCentricPLL(wgraph,batchSize,num_of_threads);
                                    VCPLL->reorder_labels(degree_order.rank);
                                    //save labels
                                    degree_order.save_rank(outputDirName);//no backprefix for ordering class
                                    cout<<"save rank succesfully!"<<endl;
                                    VCPLL->save_labels(outputDirName);
                                    cout<<"save labels succesfully!"<<endl;
                                    VCPLL->save_label_size(outputDirName,degree_order.inv);
                                    cout<<"save_label_size succesfully!"<<endl;
                                    VCPLL->write_labels(outputDirName,degree_order.inv);
                                    cout<<"write_labels succesfully!"<<endl;
                                    delete VCPLL;
                                    return 0;
                                }
                            }
                        }else{//vectorization
                            cout<<"vectorization=true ";
                            if(t_multithread_flag==0){// Single Thread
                                cout<<"multithread=false "<<endl;

                            }else{// Multithread
                                cout<<"multithread=true "<<endl;
                                //     if(t_ordering_flag==0){//given order
                                //     cout<<"Generate given order start!"<<endl;
                                //     Given_Ordering given_ordering(inputOrderFileName,wgraph,true);
                                //     cout<<"Generate given order successfully!"<<endl;
                                //     double _labeling_time=GetCurrentTimeSec();
                                //     //omp_set_num_threads(num_of_threads);
                                //     WeightedParaVertexCentricPLLVec *VCPLL =new WeightedParaVertexCentricPLLVec(wgraph,batchSize,num_of_threads);
                                //     _labeling_time=GetCurrentTimeSec()-_labeling_time;
                                //     VCPLL->reorder_labels(given_ordering.rank);
                                //     cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                                //     double ave_labeling_time=_labeling_time/(double) numOfVertices;
                                //     cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                                //     //save labels
                                //     given_ordering.save_rank(outputDirName);//no backprefix for ordering class
                                //     cout<<"save rank succesfully!"<<endl;
                                //     VCPLL->save_labels(outputDirName);
                                //     cout<<"save labels succesfully!"<<endl;
                                //     VCPLL->save_label_size(outputDirName,given_ordering.inv);
                                //     cout<<"save_label_size succesfully!"<<endl;
                                //     VCPLL->write_labels(outputDirName,given_ordering.inv);
                                //     cout<<"write_labels succesfully!"<<endl;
                                //     delete VCPLL;
                                //     return 0;
                                // }else if (t_ordering_flag==1){//degree order
                                //     cout<<"Generate degree order start!"<<endl;
                                //     Degree_Ordering degree_order(wgraph,true);
                                //     cout<<"Generate degree order successfully!"<<endl;
                                //     double _labeling_time=GetCurrentTimeSec();
                                //     //omp_set_num_threads(num_of_threads);
                                //     WeightedParaVertexCentricPLL *VCPLL =new WeightedParaVertexCentricPLL(wgraph,batchSize,num_of_threads);
                                //     _labeling_time=GetCurrentTimeSec()-_labeling_time;
                                //     VCPLL->reorder_labels(degree_order.rank);
                                //     cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                                //     double ave_labeling_time=_labeling_time/(double) numOfVertices;
                                //     cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                                //     //save labels
                                //     degree_order.save_rank(outputDirName);//no backprefix for ordering class
                                //     cout<<"save rank succesfully!"<<endl;
                                //     VCPLL->save_labels(outputDirName);
                                //     cout<<"save labels succesfully!"<<endl;
                                //     VCPLL->save_label_size(outputDirName,degree_order.inv);
                                //     cout<<"save_label_size succesfully!"<<endl;
                                //     VCPLL->write_labels(outputDirName,degree_order.inv);
                                //     cout<<"write_labels succesfully!"<<endl;
                                //     delete VCPLL;
                                //     return 0;
                                // }
                            }
                        }
                    }else{
                        cout<<"WEIGHTED_FLAG=false ";
                    }
                }

                return 0;
            }//main ends
    };
}

#endif // !_PARALLEL_INDEXING_LABELING_H