/*
 * @Descripttion:  whp-index construction
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-10 21:34:21
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-02-27 12:24:01
 */

#ifndef _COMMAND_CONSTRUCTION_LABELING
#define _COMMAND_CONSTRUCTION_LABELING

#include <cstring>
#include "../command.h"
#include "../src/time_util.h"
#include "../src/paras.h"
#include "../src/frequency_hierarchy_ordering.h"

namespace command{
    class ConstructLabel:public Command{
        public:
            void exit_with_help(){
                printf("Usage:\n");
                printf("\tbin/pll_whp_run -p -d [directedGraphFlag] -w [weightedGraphFlag] -s [experimentFlag] -g [graphFileName] -e [exportLabelFileName]\n \\ 
                 [-q [queryFreqFileName]] [-i [labelSizeFileName]] [-f [outputResultDirName]] [-b [betwennessFileName]] [-k [freqParameter]] [-l [betParameter]] \n \\ 
                 [-c [graphType]] [-o [orderWeightFileName]] [-r [queryPairFileName]] \n");
                printf("-------------------------------------------------------------------\n");
                printf("Parameter explanation:\n");
                printf("\t-d: 0 or 1, for undirected and directed graphs, default is 0\n");
                printf("\t-w: 0 or 1, for unweighted and weighted graphs, default is 0\n");
                printf("\t-s: default label\n \t\t\t1: path label\n \t\t\t2: bp label\n \t\t\t3: HLC label\n \t\t\t4: HLCM label\n");
                printf("\t-k: [0,100], value of parameter beta\n");
                printf("\t-l: 100-freqParameter, weight value of topological importance\n");
                printf("\t-b: betweenness value of each vertex computed once before\n");
            }

            int main(int argc, char* argv[]){
                char graphFileName[255] = "";
                char labelFileName[255] = "";
                char labelSizeFileName[255]="";
                int t_directed_flag = 0;
                int t_weighted_flag = 0;
                int t_ordering_flag = 0; 
                int t_indexing_flag=0; //indexing model opetions:
                int t_experiment_flag=0;//0-debug model 1-only output 
                char queryFreqFileName[255] = ""; //query frequency file name
                char queryPairFileName[255] = ""; //query frequency file name
                char outputDirName[255] = ""; 
                char betwennessFileName[255]=""; //abetwenness file computed before
                char coverageFileName[255]=""; //coverageFileName computed before
                char orderWeightFileName[255]=""; //output prder weight filename
                bool isOutputAnalysis=false;//whether write relevant files to file
                bool isOutputLabelSize=false;//whether write label size to file
                bool isOutputHF=false;//whether write analysis size to file
                bool isOutputHF_directed=false;//whether write analysis size to file (for directed)
                int hfRate=0;//default is 50/1000
                int graphType=0; //0-full connected 1-u is start and t is destination
                int k_deg=0,k_freq=0,k_cov=0,k_dep=0,k_bet=0;//coefficients(int)
                bool is_deg=false, is_freq=false, is_cov=false, is_dep=false,is_bet=false;

                if(argc<16) exit_with_help();

                for(int i = 2; i < argc; i++){
                    if(argv[i][0] != '-') break;
                    if(++i >= argc)
                        exit_with_help();
                    switch (argv[i-1][1]){
                        case 'd':
                            t_directed_flag = atoi(argv[i]);
                            break;
                        case 'w':
                            t_weighted_flag = atoi(argv[i]);
                            break;
                        case 's':
                            t_experiment_flag = atoi(argv[i]);
                            std::cout<<"t_experiment_flag = "<<t_experiment_flag<<endl;
                            break;
                        case 'c':
                            graphType = atoi(argv[i]);
                            break;
                        case 'g':
                            strcpy(graphFileName,argv[i]);
                            std::cout<<"graphFileName = "<<graphFileName<<endl;
                            break;
                        case 'e':
                            strcpy(labelFileName, argv[i]);
                            std::cout<<"labelFileName = "<<labelFileName<<endl;
                            break;
                        case 'o':
                            strcpy(orderWeightFileName, argv[i]);
                            std::cout<<"orderWeightFileName = "<<orderWeightFileName<<endl;
                            break;
                        case 'q':
                            strcpy(queryFreqFileName, argv[i]);
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
                            if(*labelSizeFileName!='\0'){
                                isOutputLabelSize=true;
                                cout<<"labelSizeFileName = "<<labelSizeFileName<<endl;
                            }
                            break;
                        case 'f':
                            strcpy(outputDirName, argv[i]);
                            if(*outputDirName!='\0'){
                                isOutputAnalysis=true;
                                std::cout<<"outputDirName = "<<outputDirName<<endl;
                            }else cerr<<"outputDirName cannot be null!"<<endl;          
                            break;
                        case 'b':
                            strcpy(betwennessFileName, argv[i]);
                            std::cout<<"betwennessFileName="<<betwennessFileName<<endl;
                            break;
                        case 'j':
                            k_deg = atoi(argv[i]);
                            if(k_deg!=0) is_deg=true;
                            break;
                        case 'k':
                            k_freq = atoi(argv[i]);
                            if(k_freq!=0) is_freq=true;
                            break;
                        case 'l':
                            k_bet = atoi(argv[i]);
                            if(k_bet!=0) is_bet=true;
                            break;
                        case 'u':
                            k_cov = atoi(argv[i]);
                            if(k_cov!=0) is_cov=true;
                            break;
                        case 'v':
                            k_dep = atoi(argv[i]);
                            if(k_dep!=0) is_dep=true;
                            break;
                        default:
                            exit_with_help();
                            break;
                    }
                }

            if(!is_deg&&!is_dep&&!is_freq&&!is_cov){
                cout<<"All weight parameters are 0!"<<std::endl;
                exit_with_help();
            }
            if(is_bet&&*betwennessFileName=='\0'){
                cout<<"betwennessFileName cannot be null!"<<std::endl;
                exit_with_help();  
            }
            if(*labelFileName=='\0'){
                cout<<"labelFileName cannot be null!"<<std::endl;
                exit_with_help();
            }

            if (t_directed_flag == 1)
                DIRECTED_FLAG = true;
            if (t_weighted_flag == 1)
                WEIGHTED_FLAG = true;

            WGraph wgraph;
            Graph graph;

            if(WEIGHTED_FLAG){
                wgraph.load_graph(graphFileName,graphType);
            }
            else graph.load_graph(graphFileName);
            std::cout << numOfVertices << " nodes and " << numOfEdges << " arcs " << endl;
            
            if (numOfVertices == 0 || numOfEdges == 0){
                std::cout << "Corruptted graph file" << endl;
                return 0;
            }

            //indexing

            std::cout<<"t_special_flag==0 default label"<<endl;
            if(DIRECTED_FLAG==true){
                std::cout<<"DIRECTED_FLAG==true"<<endl;
                if(WEIGHTED_FLAG==true){
                    std::cout<<"WEIGHTED_FLAG==true"<<endl;
                    std::cout<<"************2-precompute items weight*coefficient selection push algorithm*******"<<endl;
                    double _labeling_time, _hub_time,_ordering_time;//output time
                    _ordering_time=GetCurrentTimeSec();
                    Processing::calcCoefficient<double> calcCoef((double)k_deg/100,(double)k_freq/100,(double)k_cov/100,(double)k_dep/100,(double)k_bet/100,is_deg, is_freq, is_cov, is_dep,is_bet);
                    Linear_Ordering<double> linear_ordering(wgraph,calcCoef,orderWeightFileName,queryFreqFileName,betwennessFileName,coverageFileName,hfRate,t_ordering_flag);
                    _ordering_time=GetCurrentTimeSec()-_ordering_time;
                    _hub_time = GetCurrentTimeSec();
                    PL_W pl_w(wgraph, linear_ordering, DIRECTED_FLAG);
                    _hub_time = GetCurrentTimeSec()-_hub_time;
                    _labeling_time=_ordering_time+_hub_time;
                    double ave_labeling_time=_labeling_time/(double) numOfVertices;
                    cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                    cout << "pushing time:" <<_hub_time *1e6 <<  " microseconds" << endl;
                    cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                    cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                    //save labels to binary file
                    pl_w.dlabels.save_labels(labelFileName);
                    cout<<"***************save_labels successffully!(binary)*******************"<<endl;
                    //save label size to txt file
                    if(isOutputLabelSize) pl_w.dlabels.save_label_size(labelSizeFileName);
                    cout<<"**********************save_label_size succesfully!**********************"<<endl;
                    if(t_experiment_flag==0){
                        //output analysis file
                        if(isOutputAnalysis){
                            //save rank to txt file
                            linear_ordering.save_rank(outputDirName);
                            cout<<"******************save_rank successffully!!****************"<<endl;
                            pl_w.dlabels.write_labels(outputDirName,linear_ordering.inv);
                            cout<<"**********************write_labels succesfully!**********************"<<endl;
                            if(isOutputHF_directed){
                                pl_w.dlabels.load_hfpoint_and_qt(queryPairFileName,hfRate);
                                pl_w.dlabels.save_anaylysis_size(outputDirName);                
                                cout<<"******************save analysis successfully!!****************"<<endl;           
                            }
                        }
                    }else if(t_experiment_flag==1){
                        if(isOutputAnalysis&&isOutputHF_directed){
                            pl_w.dlabels.load_hfpoint_and_qt(queryPairFileName,hfRate);
                            pl_w.dlabels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                            cout<<"save experiment result successfully!!"<<endl;       
                        }

                    }else {//update model
                        if(isOutputAnalysis){
                            ofstream ofs(outputDirName,ios::app|ios::out);//append way
                            if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                            ofs.setf(ios::fixed);
                            ofs.precision(4);
                            ofs<<_ordering_time* 1e6 <<" "<<_labeling_time* 1e6 <<" ";
                            ofs.close();             
                        }
                    }
                    return 0;

                }else{
                    std::cout<<"WEIGHTED_FLAG==false"<<endl;
                }
                
            }else{
                if(WEIGHTED_FLAG==true){
                    std::cout<<"WEIGHTED_FLAG==true"<<endl;
  
                    std::cout<<"************2-precompute items weight*coefficient selection push algorithm*******"<<endl;
                    //compute the coefficients  of order
                    double _labeling_time, _hub_time,_ordering_time;//output time
                    _ordering_time=GetCurrentTimeSec();
                    Processing::calcCoefficient<double> calcCoef((double)k_deg/100,(double)k_freq/100,(double)k_cov/100,(double)k_dep/100,(double)k_bet/100,is_deg, is_freq, is_cov, is_dep,is_bet);
                    Linear_Ordering<double> linear_ordering(wgraph,calcCoef,orderWeightFileName,queryFreqFileName,betwennessFileName,coverageFileName,hfRate,t_ordering_flag);
                    _ordering_time=GetCurrentTimeSec()-_ordering_time;
                    _hub_time = GetCurrentTimeSec();
                    PL_W pl_w(true,wgraph, linear_ordering);
                    _hub_time = GetCurrentTimeSec()-_hub_time;
                    _labeling_time=_ordering_time+_hub_time;
                    double ave_labeling_time=_labeling_time/(double) numOfVertices;
                    cout << "ordering time:" << _ordering_time *1e6 <<  " microseconds" << endl;
                    cout << "pushing time:" <<_hub_time *1e6 <<  " microseconds" << endl;
                    cout << "labeling time:" << _labeling_time *1e6 <<  " microseconds" << endl;
                    cout<<"average labeling time:"<<ave_labeling_time*1e6 <<  " microseconds" << endl;
                    //save labels to binary file
                    pl_w.labels.save_labels(labelFileName);
                    cout<<"***************save_labels successffully!(binary)*******************"<<endl;
                    //save labels to txt file
                    if(isOutputLabelSize) pl_w.labels.save_label_size(labelSizeFileName);
                    cout<<"**********************save_label_size succesfully!**********************"<<endl;
                    if(t_experiment_flag==0){
                        //output analysis file
                        if(isOutputAnalysis){
                            //save rank to txt file
                            linear_ordering.save_rank(outputDirName);
                            cout<<"******************save_rank successffully!!****************"<<endl;
                            pl_w.labels.write_labels(outputDirName,linear_ordering.inv);
                            cout<<"**********************write_labels succesfully!**********************"<<endl;
                            if(isOutputHF){
                                pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                                pl_w.labels.save_anaylysis_size(outputDirName);                
                                cout<<"******************save analysis successfully!!****************"<<endl;           
                            }
                        }
                    }else if(t_experiment_flag==1){
                        if(isOutputAnalysis&&isOutputHF){
                            pl_w.labels.load_hfpoint_and_qt(queryFreqFileName,hfRate);
                            pl_w.labels.append_experiment_result(outputDirName,_labeling_time,_ordering_time);
                            cout<<"save experiment result successfully!!"<<endl;       
                        }
                    }else{//update model
                        if(isOutputAnalysis){
                            ofstream ofs(outputDirName,ios::app|ios::out);//append way
                            if(!ofs.is_open()) cout<<"Cannot open "<<outputDirName<<endl;
                            ofs.setf(ios::fixed);
                            ofs.precision(4);
                            ofs<<_ordering_time* 1e6 <<" "<<_labeling_time* 1e6 <<" ";
                            ofs.close();
                        }
                    }
                    return 0;

                }else{
                    std::cout<<"WEIGHTED_FLAG==false"<<endl;
                }
            }
        }
            
    };
}

#endif