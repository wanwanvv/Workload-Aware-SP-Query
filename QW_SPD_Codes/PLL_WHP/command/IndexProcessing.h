/*
 * @Descripttion: index construction
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-10 11:27:10
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-02-13 09:54:31
 */

#ifndef _COMMAND_INDEXING_PROCESSING
#define _COMMAND_INDEXING_PROCESSING


#include "../src/graph.h"
#include "../src/ordering.h"
#include "../src/construction.h"
#include "../src/labels.h"
#include "../src/time_util.h"
#include "../command.h"
#include "../src/coverage_ordering.h"
#include <cstring>

using namespace time_util; 

namespace command{
    class IndexProcessing: public Command{
    public:
        void exit_with_help(){
            printf("Usage:\n");
            printf("\tbin/pll_whp_run -x -d [directedGraphFlag] -w [weightedGraphFlag] -s [experimentFlag] -m [orderingSchemes] [-a [s_beta]] \n \\ 
            -g [graphFileName] -e [outputLabelFileName] -i [outputLabelSizeFileName] [-f [outputResultDirName]] [-q [queryFreqFileName]] [-r [queryPairFileName]]\n \\ 
            [-b [outputBetweennessFile]] [-c [graphType]]\n");
            printf("-------------------------------------------------------------------\n");
            printf("Parameter explanation:\n");
            printf("\t-d: 0 or 1, for undirected and directed graphs, default is 0\n");
            printf("\t-w: 0 or 1, for unweighted and weighted graphs, default is 0\n");
            printf("\t-s: default label\n \t\t\t1: path label\n \t\t\t2: bp label\n \t\t\t3: HLC label\n \t\t\t4: HLCM label\n");
            printf("\t-m: 0-DHP \n \t\t\t1-BHP\n \t\t\t2-SHP\n");
            printf("\t-a: [Optional]trade-off parameter for BHP, default is 1\n");
            printf("\t-b: optional, output the betweeness values\n");
            printf("-------------------------------------------------------------------\n");
            exit(1);
        }
        int main(int argc, char *argv[])
        {
            char graphFileName[255] = "";
            char labelFileName[255] = "";
            char labelSizeFileName[255] = "";
            char outputBetFileName[255] = "";
            int t_directed_flag = 0;
            int t_weighted_flag = 0;
            int t_ordering_flag = 0;
            int t_experiment_flag=0;//1-experiment 0-debug model
            double beta = 1;
            char outputDirName[255] = ""; 
            char givenOrderFileName[255] = ""; //given order filename for pll
            char queryFreqFileName[256]=""; //point query time filename(node t,freq)
            char queryPairFileName[256]=""; //pair query time filename(s,t,freq)
			bool isOutputAnalysis=false;//whether write relevant files
            bool isOutputLabelSize=false;//whether write labelsize to file
            bool isOutputHF=false;//whether write analysis size to file
            bool isOutputHF_directed=false;//whether write analysis size to file
            int graphType=0; //0-full connected 1-u is start and t is destination
            int hfRate=0; //high frequency point rate 1/1000
            int k_deg=0,k_freq=0,k_cov=0,k_dep=0,k_bet=0;//coefficients(int)
            bool is_deg=false, is_freq=false, is_cov=false, is_dep=false,is_bet=false;
            
            // if(argc < 14) // modified by wanjingyi
            //     exit_with_help();
            
            for(int i = 2; i < argc; i++){
                if(argv[i][0] != '-') break;
                if(++i >= argc)
                    exit_with_help();
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
                    case 'c':
                        graphType = atoi(argv[i]);
                        break;
                    case 's':
                        t_experiment_flag = atoi(argv[i]);
                        std::cout<<"t_experiment_flag = "<<t_experiment_flag<<endl;
                        break;
                    case 'm':
                        t_ordering_flag = atoi(argv[i]);
                        break;
                    case 'e':
                        strcpy(labelFileName, argv[i]);
                        cout<<"labelFileName = "<<labelFileName<<endl;
                        break;
                    case 'b':
                        strcpy(outputBetFileName, argv[i]);
                        cout<<"outputBetFileName = "<<outputBetFileName<<endl;
                        break;
                    case 'a':
                        beta = atof(argv[i]);
                        break;
                    case 'f':
                        strcpy(outputDirName,argv[i]);
                        if(*outputDirName!='\0') isOutputAnalysis=true;
                        cout<<"outputDirName = "<<outputDirName<<endl;
                        break;
                    case 'h':
                        strcpy(givenOrderFileName,argv[i]);
                        cout<<"givenOrderFileName = "<<givenOrderFileName<<endl;
                        break;
                    case 'q':
                        strcpy(queryFreqFileName,argv[i]);
                        if(*queryFreqFileName!='\0'){
                            isOutputHF=true;
                            std::cout<<"queryFreqFileName = "<<queryFreqFileName<<endl;
                        }else cerr<<"queryFreqFileName cann't be null!"<<endl;     
                        break;
                    case 'r':
                        strcpy(queryPairFileName, argv[i]);
                        if(*queryPairFileName!='\0'){
                            isOutputHF_directed=true;
                            std::cout<<"queryPairFileName = "<<queryPairFileName<<endl;
                        }else cerr<<"queryPairFileName cann't be null!"<<endl;          
                        break;
                    case 'i':
                        strcpy(labelSizeFileName,argv[i]);
                        cout<<"labelSizeFileName = "<<labelSizeFileName<<endl;
                        if(*outputDirName!='\0') isOutputLabelSize=true;
                        break;
                    default:
                        cout<<"default!"<<std::endl;
                        exit_with_help();
                        break;
                }
            }
            
            if (t_directed_flag == 1)
                DIRECTED_FLAG = true;
            if (t_weighted_flag == 1)
                WEIGHTED_FLAG = true;
            
            if(t_directed_flag != 1 && t_directed_flag != 0)
                exit_with_help();
            if(t_weighted_flag != 1 && t_weighted_flag != 0)
                exit_with_help();
            if(beta<=0)
                exit_with_help();
            if(*labelFileName=='\0'){
                cout<<"labelFileName cannot be null!"<<std::endl;
                exit_with_help();
            }
            if(t_ordering_flag==3&&*givenOrderFileName=='\0') exit_with_help();
            
            Graph graph;
            WGraph wgraph;
            if(WEIGHTED_FLAG == true){
                wgraph.load_graph(graphFileName,graphType);
            }else{
                graph.load_graph(graphFileName);
            }
            cout << numOfVertices << " nodes and " << numOfEdges << " arcs " << endl;
            
            if (numOfVertices == 0 || numOfEdges == 0){
                cout << "Corruptted graph file" << endl;
                return 0;
            }
            
            // Indexing
            if (DIRECTED_FLAG == true){
                if( WEIGHTED_FLAG == true){
                    if(t_ordering_flag == 0){
                        cout<<"using Degree_Ordering......"<<endl;
                        double _labeling_time,_hub_time,_ordering_time;
                        _ordering_time = GetCurrentTimeSec();
                        Degree_Ordering degree_order(wgraph);
                        _ordering_time=GetCurrentTimeSec()-_ordering_time;
                        _hub_time=GetCurrentTimeSec();
                        PL_W pl_w(wgraph, degree_order, DIRECTED_FLAG);
                        _hub_time = GetCurrentTimeSec() - _hub_time;
                        _labeling_time=_ordering_time+_hub_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                        cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            pl_w.dlabels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) pl_w.dlabels.save_label_size(labelSizeFileName,degree_order.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                degree_order.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                pl_w.dlabels.write_labels(outputDirName,degree_order.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF_directed){
                                    pl_w.dlabels.load_hfpoint_and_qt(queryPairFileName,hfRate);
                                    pl_w.dlabels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            pl_w.dlabels.save_labels(labelFileName);
                            if(isOutputLabelSize) pl_w.dlabels.save_label_size(labelSizeFileName);
                            //pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            //pl_w.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                            if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }
                            cout<<"save experiment result successfully!!"<<endl;
                        }
                        return 0;
                    }else if (t_ordering_flag == 1){
                        cout<<"using Betweenness_Ordering......"<<endl;
                        double _labeling_time,_hub_time,_ordering_time;
                        Betweenness_Ordering betweenness_ordering(16, beta, wgraph, numOfVertices,_ordering_time,_hub_time,outputBetFileName);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        _labeling_time=_ordering_time+_hub_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                        cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            betweenness_ordering.dlabels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) betweenness_ordering.dlabels.save_label_size(labelSizeFileName,betweenness_ordering.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                betweenness_ordering.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                betweenness_ordering.dlabels.write_labels(outputDirName,betweenness_ordering.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF_directed){
                                    betweenness_ordering.dlabels.load_hfpoint_and_qt(queryPairFileName,hfRate);
                                    betweenness_ordering.dlabels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            betweenness_ordering.dlabels.save_labels(labelFileName);
                            if(isOutputLabelSize) betweenness_ordering.dlabels.save_label_size(labelSizeFileName);
                            //betweenness_ordering.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            //betweenness_ordering.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                                if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }                               
                            cout<<"save experiment result successfully!!"<<endl;
                        }
                        return 0;
                    }else if (t_ordering_flag == 2){ // t_ordering_flag == 2
                        double _labeling_time = GetCurrentTimeSec();
                        Coverage_Ordering coverage_ordering(wgraph, DIRECTED_FLAG, true);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout<<" average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            coverage_ordering.dlabels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) coverage_ordering.dlabels.save_label_size(labelSizeFileName, coverage_ordering.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                coverage_ordering.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                coverage_ordering.dlabels.write_labels(outputDirName , coverage_ordering.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF_directed){
                                    coverage_ordering.dlabels.load_hfpoint_and_qt(queryPairFileName,hfRate);
                                    coverage_ordering.dlabels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            coverage_ordering.labels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) coverage_ordering.labels.save_label_size(labelSizeFileName, coverage_ordering.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << "#" <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }
                            cout<<"save experiment result successfully!!"<<endl;
                        }
                        return 0;
                    }else if(t_ordering_flag==3){
                        cout<<"using Given_Ordering......"<<endl;
                        double _labeling_time,_hub_time,_ordering_time;
                        _ordering_time = GetCurrentTimeSec();
                        Given_Ordering given_order(givenOrderFileName,wgraph);
                        _ordering_time=GetCurrentTimeSec()-_ordering_time;
                        _hub_time=GetCurrentTimeSec();
                        PL_W pl_w(wgraph, given_order, DIRECTED_FLAG);
                        _hub_time = GetCurrentTimeSec() - _hub_time;
                        _labeling_time=_ordering_time+_hub_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                        cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            pl_w.dlabels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) pl_w.dlabels.save_label_size(labelSizeFileName,given_order.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                given_order.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                pl_w.dlabels.write_labels(outputDirName,given_order.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF_directed){
                                    pl_w.dlabels.load_hfpoint_and_qt(queryPairFileName,hfRate);
                                    pl_w.dlabels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            pl_w.dlabels.save_labels(labelFileName);
                            if(isOutputLabelSize) pl_w.dlabels.save_label_size(labelSizeFileName);
                            //pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            //pl_w.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                            if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }
                            cout<<"save experiment result successfully!!"<<endl;
                        }
                        return 0;
                    }else if(t_ordering_flag==4){//random order
                        cout<<"using Random_Ordering......"<<endl;
                        double _labeling_time,_hub_time,_ordering_time;
                        _ordering_time = GetCurrentTimeSec();
                        Random_Ordering random_order(wgraph);
                        _ordering_time=GetCurrentTimeSec()-_ordering_time;
                        _hub_time=GetCurrentTimeSec();
                         PL_W pl_w(wgraph, random_order, DIRECTED_FLAG);
                        _hub_time = GetCurrentTimeSec() - _hub_time;
                        _labeling_time=_ordering_time+_hub_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                        cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            pl_w.dlabels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) pl_w.dlabels.save_label_size(labelSizeFileName,random_order.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                random_order.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                pl_w.dlabels.write_labels(outputDirName,random_order.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF_directed){
                                    pl_w.dlabels.load_hfpoint_and_qt(queryPairFileName,hfRate);
                                    pl_w.dlabels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            pl_w.dlabels.save_labels(labelFileName);
                            if(isOutputLabelSize) pl_w.dlabels.save_label_size(labelSizeFileName);
                            //pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            //pl_w.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                            if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }
                            cout<<"save experiment result successfully!!"<<endl;
                        }
                        return 0;
                    }else if(t_ordering_flag==5){//frequency order
                        cout<<"using Frequency_Ordering......"<<endl;
                        double _labeling_time,_hub_time,_ordering_time;
                        _ordering_time = GetCurrentTimeSec();
                        Frequency_Ordering frequency_order(queryFreqFileName,wgraph);
                        _ordering_time=GetCurrentTimeSec()-_ordering_time;
                        _hub_time=GetCurrentTimeSec();
                         PL_W pl_w(wgraph, frequency_order, DIRECTED_FLAG);
                        _hub_time = GetCurrentTimeSec() - _hub_time;
                        _labeling_time=_ordering_time+_hub_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                        cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            pl_w.dlabels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) pl_w.dlabels.save_label_size(labelSizeFileName,frequency_order.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                frequency_order.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                pl_w.dlabels.write_labels(outputDirName,frequency_order.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF_directed){
                                    pl_w.dlabels.load_hfpoint_and_qt(queryPairFileName,hfRate);
                                    pl_w.dlabels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            pl_w.dlabels.save_labels(labelFileName);
                            if(isOutputLabelSize) pl_w.dlabels.save_label_size(labelSizeFileName);
                            //pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            //pl_w.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                            if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }
                            cout<<"save experiment result successfully!!"<<endl;
                        }
                        return 0;
                    }
                }else{
                    //***********************directed unweighted*********************//
                    if(t_ordering_flag == 0){
                        double _labeling_time = GetCurrentTimeSec();
                        Degree_Ordering degree_order(graph);
                        PL pl(graph, degree_order, DIRECTED_FLAG);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout<<" average indexing time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        pl.dlabels.save_labels(labelFileName);
                        return 0;
                    }else if (t_ordering_flag == 1){
                        double _labeling_time = GetCurrentTimeSec();
                        Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout<<" average indexing time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        string orderFileName(labelFileName);
                        betweenness_ordering.save_rank(orderFileName.c_str());
                        betweenness_ordering.dlabels.save_labels(labelFileName);
                        return 0;
                    }else if(t_ordering_flag == 2){ 
                        double _labeling_time = GetCurrentTimeSec();
                        Coverage_Ordering coverage_ordering(graph, DIRECTED_FLAG);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout<<" average indexing time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        coverage_ordering.dlabels.save_labels(labelFileName);
                        return 0;
                    }else{//t_ordering_flag == 3
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(givenOrderFileName,graph);
                        PL pl(graph, given_order, DIRECTED_FLAG);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout<<" average indexing time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        pl.dlabels.save_labels(labelFileName);
                        return 0;
                    }
                    //***********************directed unweighted*********************//
                }
            }else{ //undirected
                if( WEIGHTED_FLAG == true){
                    if(t_ordering_flag == 0){
                        cout<<"using Degree_Ordering......"<<endl;
                        double _labeling_time,_hub_time,_ordering_time;
                        _ordering_time = GetCurrentTimeSec();
                        Degree_Ordering degree_order(wgraph);
                        _ordering_time=GetCurrentTimeSec()-_ordering_time;
                        _hub_time=GetCurrentTimeSec();
                        PL_W pl_w(wgraph, degree_order);
                        _hub_time = GetCurrentTimeSec() - _hub_time;
                        _labeling_time=_ordering_time+_hub_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                        cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            pl_w.labels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) pl_w.labels.save_label_size(labelSizeFileName,degree_order.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                degree_order.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                pl_w.labels.write_labels(outputDirName,degree_order.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF){
                                    pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                                    pl_w.labels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            pl_w.labels.save_labels(labelFileName);
                            if(isOutputLabelSize) pl_w.labels.save_label_size(labelSizeFileName);
                            //pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            //pl_w.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                            if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }
                            cout<<"save experiment result successfully!!"<<endl;
                        }

                        return 0;
                    }else if (t_ordering_flag == 1){
                        cout<<"using Betweenness_Ordering......"<<endl;
                        double _labeling_time,_hub_time,_ordering_time;
                        Betweenness_Ordering betweenness_ordering(16, beta, wgraph, numOfVertices,_ordering_time,_hub_time,outputBetFileName);
                        _labeling_time=_ordering_time+_hub_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                        cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            betweenness_ordering.labels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) betweenness_ordering.labels.save_label_size(labelSizeFileName,betweenness_ordering.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                betweenness_ordering.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                betweenness_ordering.labels.write_labels(outputDirName,betweenness_ordering.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF){
                                    betweenness_ordering.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                                    betweenness_ordering.labels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            betweenness_ordering.labels.save_labels(labelFileName);
                            if(isOutputLabelSize) betweenness_ordering.labels.save_label_size(labelSizeFileName);
                            //betweenness_ordering.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            //betweenness_ordering.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                                if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }                               
                            cout<<"save experiment result successfully!!"<<endl;
                        }
                        return 0;
                    }else if(t_ordering_flag == 2){ // t_ordering_flag == 2
                        cout<<"using Coverage_Ordering......"<<endl;
                        double _labeling_time = GetCurrentTimeSec();
                        Coverage_Ordering coverage_ordering(wgraph);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            coverage_ordering.labels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) coverage_ordering.labels.save_label_size(labelSizeFileName, coverage_ordering.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                coverage_ordering.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                coverage_ordering.labels.write_labels(outputDirName , coverage_ordering.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF){
                                    coverage_ordering.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                                    coverage_ordering.labels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            coverage_ordering.labels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) coverage_ordering.labels.save_label_size(labelSizeFileName, coverage_ordering.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << "#" <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }
                            cout<<"save experiment result successfully!!"<<endl;
                        }
                        return 0;
                    }else if(t_ordering_flag == 3){
                            cout<<"using Given_Ordering......"<<endl;
                            double _labeling_time,_hub_time,_ordering_time;
                            _ordering_time = GetCurrentTimeSec();
                            Given_Ordering given_order(givenOrderFileName, wgraph);
                            _ordering_time=GetCurrentTimeSec()-_ordering_time;
                            _hub_time=GetCurrentTimeSec();
                    	    PL_W pl_w(wgraph, given_order);
                            _hub_time = GetCurrentTimeSec() - _hub_time;
                            _labeling_time = _ordering_time+_hub_time;
                            double ave_labeling_time=_labeling_time/(double) numOfVertices;
                            cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                            cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                            cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                            cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                            if(t_experiment_flag==0){
                                pl_w.labels.save_labels(labelFileName);
                                cout<<"**********************save_labels successfully!**********************"<<endl;
                                pl_w.labels.save_label_size(labelSizeFileName,given_order.inv);
                                cout<<"**********************save_label_size succesfully!**********************"<<endl;
                                if(isOutputAnalysis){
                                    pl_w.labels.write_labels(outputDirName,given_order.inv);
                                    cout<<"**********************write_labels succesfully!**********************"<<endl;
                                    pl_w.labels.save_anaylysis_size(queryFreqFileName,outputDirName,hfRate);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }else if(t_experiment_flag==1){
                                pl_w.labels.save_labels(labelFileName);
                                pl_w.labels.save_label_size(labelSizeFileName);
                                //pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                                //pl_w.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                                 if(isOutputAnalysis){
                                    ofstream ofs(outputDirName,ios::app|ios::out);
                                    if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                    ofs.setf(ios::fixed);
                                    ofs.precision(4);
                                    ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                    ofs.close();
                                }                               
                            }
                        return 0;
                    }else if(t_ordering_flag==4){//random order
                        cout<<"using Random_Ordering......"<<endl;
                        double _labeling_time,_hub_time,_ordering_time;
                        _ordering_time = GetCurrentTimeSec();
                        Random_Ordering random_order(wgraph);
                        _ordering_time=GetCurrentTimeSec()-_ordering_time;
                        _hub_time=GetCurrentTimeSec();
                        PL_W pl_w(wgraph, random_order);
                        _hub_time = GetCurrentTimeSec() - _hub_time;
                        _labeling_time=_ordering_time+_hub_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                        cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            pl_w.labels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) pl_w.labels.save_label_size(labelSizeFileName,random_order.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                random_order.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                pl_w.labels.write_labels(outputDirName,random_order.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF){
                                    pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                                    pl_w.labels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            pl_w.labels.save_labels(labelFileName);
                            if(isOutputLabelSize) pl_w.labels.save_label_size(labelSizeFileName);
                            //pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            //pl_w.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                            if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }
                            cout<<"save experiment result successfully!!"<<endl;
                        }   
                        return 0;
                    }else if(t_ordering_flag==5){//frequency order
                        cout<<"using Frequency_Ordering......"<<endl;
                        double _labeling_time,_hub_time,_ordering_time;
                        _ordering_time = GetCurrentTimeSec();
                        Frequency_Ordering frequency_order(queryFreqFileName,wgraph);
                        _ordering_time=GetCurrentTimeSec()-_ordering_time;
                        _hub_time=GetCurrentTimeSec();
                        PL_W pl_w(wgraph, frequency_order);
                        _hub_time = GetCurrentTimeSec() - _hub_time;
                        _labeling_time=_ordering_time+_hub_time;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                        cout << "pushing time:" << _hub_time *1e6 <<  " microseconds" << endl;
                        cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        if(t_experiment_flag==0){
                            pl_w.labels.save_labels(labelFileName);
                            cout<<"**********************save_labels successfully!**********************"<<endl;
                            if(isOutputLabelSize) pl_w.labels.save_label_size(labelSizeFileName,frequency_order.inv);
                            cout<<"**********************save_label_size succesfully!**********************"<<endl;
                            if(isOutputAnalysis){
                                frequency_order.save_rank(outputDirName);
                                cout<<"**********************save_rank succesfully!**********************"<<endl;
                                pl_w.labels.write_labels(outputDirName,frequency_order.inv);
                                cout<<"**********************write_labels succesfully!**********************"<<endl;
                                if(isOutputHF){
                                    pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                                    pl_w.labels.save_anaylysis_size(outputDirName);
                                    cout<<"**********************save_analysis_size succesfully!**********************"<<endl;
                                }
                            }
                        }else if(t_experiment_flag==1){
                            pl_w.labels.save_labels(labelFileName);
                            if(isOutputLabelSize) pl_w.labels.save_label_size(labelSizeFileName);
                            //pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            //pl_w.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                            if(isOutputAnalysis){
                                ofstream ofs(outputDirName,ios::app|ios::out);
                                if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                                ofs.setf(ios::fixed);
                                ofs.precision(4);
                                ofs  << _ordering_time *1e6 <<  " " <<_labeling_time *1e6 <<" ";
                                ofs.close();
                            }
                            cout<<"save experiment result successfully!!"<<endl;
                        }
                        return 0;
                    }
                }else{
                     //unweighted
                    if(t_ordering_flag == 0){
                        double _labeling_time = GetCurrentTimeSec();
                        Degree_Ordering degree_order(graph);
                        PL pl(graph, degree_order);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout<<" average indexing time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        
                        pl.labels.save_labels(labelFileName);
                        /*********save label size written by wanjingyi****/
                        if(isOutputLabelSize) pl.labels.save_label_size(labelSizeFileName,degree_order.rank);
                        cout<<"save_label_size succesfully!"<<endl;
                        /*********save label size written by wanjingyi****/

                        /************write_labels to file modified by wanjingyi***********/
                        if(isOutputAnalysis) pl.labels.write_labels(outputDirName,degree_order.inv);
                        cout<<"write_labels succesfully!"<<endl;
                        /************write_labels to file modified by wanjingyi***********/
                        return 0;
                    }else if (t_ordering_flag == 1){
                        double _labeling_time = GetCurrentTimeSec();
                        Betweenness_Ordering betweenness_ordering(16, beta, graph, numOfVertices);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout<<" average indexing time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;

                        betweenness_ordering.labels.save_labels(labelFileName);
                        return 0;
                    }else if( t_ordering_flag == 2){
                        double _labeling_time = GetCurrentTimeSec();
                        vector<NodeID> border(numOfVertices);
                        Coverage_Ordering coverage_ordering(graph, border);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout<<" average indexing time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        
                        coverage_ordering.labels.save_labels(labelFileName);
                        return 0;
                    }else if( t_ordering_flag == 3){
                        double _labeling_time = GetCurrentTimeSec();
                        Given_Ordering given_order(givenOrderFileName,graph);
                        PL pl(graph, given_order);
                        _labeling_time = GetCurrentTimeSec() - _labeling_time;
                        cout << "Indexing time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                        double ave_labeling_time=_labeling_time/(double) numOfVertices;
                        cout<<" average indexing time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                        pl.dlabels.save_labels(labelFileName);
                    }
                }
            }
        }
    };
}
#endif 
