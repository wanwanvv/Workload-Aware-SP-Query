/*
 * @Descripttion: command header file used to labels construction for integrated solutions
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-04-13 09:21:59
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-11-04 15:12:09
 */
#pragma once
#ifndef _COMMAND_INTEGRATED_LABELING
#define _COMMAND_INTEGRATED_LABELING

#include <iostream>
#include <fstream>
#include <cstring>
#include "../src/time_util.h"
#include "../command.h"
#include "../src/integrated_construction.h"

using namespace time_util; 
namespace command{
    class IntegrationIndexing:public Command{
        public:
            void exit_with_help(){
                printf("Usage:\n");
                printf("wcf_du_run -z -d [directedGraphFlag] -w [weightedGraphFlag] -o [coreOrderingModel] -s [experimentFlag] \n -i [indexingModel] -m [removeMaxDegree] -k [freqParameter] -l [betParameter] [-n [bestRemoveDegree]] \n-c [numOfThreads] -g [graphFileName] -q [queryFreqFileName] -t [interFileOutputPrefix] -a [outputResultPrefix] -e [outputLabelFileName] -b [outputLabelSizeFileName] [-f [inputOrderFileName]]\n");
                printf("-------------------------------------------------------------------\n");
                printf("Parameter explanation:\n");
                printf("[coreOrderingModel] = 0: WHP; 1: DHP; 2: Given order\n");
                printf("[indexingModel] = 0: naive WCF-Index (outdated); 1: WCF-Index; 2: WCF_variant-Index\n");
                printf("[bestRemoveDegree] : <=removeMaxDegree, inter overlay graph for minimum indexing\n");
                printf("[freqParameter] = [0,100], value of parameter beta\n");
                printf("[betParameter] = 100-freqParameter, weight value of topological importance\n");
                printf("[numOfThreads]: num of threads, default is 1\n");
                printf("[removeMaxDegree]: the max degree of nodes to be removed (>=1)\n");
                printf("-------------------------------------------------------------------\n");
                exit(1);
            }

            int main(int argc, char* argv[]){
                char graphFileName[255] = "";
                char labelFileName[255] = "";
                char labelSizeFileName[255] = "";
                char outputDirName[255] = "";
                char analysisDirName[255] = "";
                char queryFreqFileName[255] = "";
                char queryPairFileName[255] = "";
                char givenOrderFileName[255] = "";
                int t_directed_flag = 0;
                int t_weighted_flag = 0;
                int t_ordering_flag=0;//0-lhp,1-dhp for overlay 2-hop label
                int t_indexing_flag=0;//0-first solution,1-second solution
                int t_experiment_flag=0;//0-debug model,1-expriment model
                bool isOutputAnalysis=false;
                int hfRate=0;
                int maxDegree=0;
                int bestDegree=-1;
                int graphType=0;
                int k_deg=0,k_freq=0,k_cov=0,k_dep=0,k_bet=0;//coefficients(int)
                bool is_deg=false, is_freq=false, is_cov=false, is_dep=false,is_bet=false;
                bool is_multiThreads=false;
                int num_threads=5;
                int is_sorted=0; //whether sort the border vertices by frequency
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
                        case 'i':
                            t_indexing_flag = atoi(argv[i]);
                            break;
                        case 'm':
                            maxDegree = atoi(argv[i]);
                            break;
                        case 'n':
                            bestDegree = atoi(argv[i]);
                            break;
                        case 'c':
                           num_threads = atoi(argv[i]);
                           if(num_threads>1) is_multiThreads=true;
                           std::cout<<"is_multiThreads="<<is_multiThreads<<std::endl;
                            break;
                        case 'o':
                            t_ordering_flag = atoi(argv[i]);
                            break;
                        case 's':
                            t_experiment_flag = atoi(argv[i]);
                            break;
                        case 'g':
                            strcpy(graphFileName,argv[i]);
                            std::cout<<"graphFileName="<<graphFileName<<std::endl;
                            break;
                        case 't':
                            strcpy(outputDirName,argv[i]);
                            std::cout<<"outputDirName="<<outputDirName<<std::endl;
                            break;
                        case 'a':
                            strcpy(analysisDirName, argv[i]);
                            isOutputAnalysis=true;
                            std::cout<<"analysisDirName="<<analysisDirName<<std::endl;
                            break;
                        case 'q':
                            strcpy(queryFreqFileName,argv[i]);
                            std::cout<<"queryFreqFileName="<<queryFreqFileName<<std::endl;
                            break;
                        case 'f':
                            strcpy(givenOrderFileName,argv[i]);
                            std::cout<<"givenOrderFileName="<<givenOrderFileName<<std::endl;
                            break;
                        case 'e':
                            strcpy(labelFileName,argv[i]);
                            std::cout<<"labelFileName="<<labelFileName<<std::endl;
                            break;
                        case 'b':
                            strcpy(labelSizeFileName,argv[i]);
                            std::cout<<"labelSizeFileName="<<labelSizeFileName<<std::endl;
                            break;
                        case 'r':
                            hfRate = atoi(argv[i]);
                            cout<<"hfRate = "<<hfRate<<endl;
                            break;
                        case 'u':
                            graphType = atoi(argv[i]);
                            break;
                        case 'p':
                            is_sorted = atoi(argv[i]);
                            break;
                        case 'j':
                            k_deg = atoi(argv[i]);
                            if(k_deg!=0){
                                is_deg=true;
                                cout<<"k_deg = "<<k_deg<<endl;
                            }
                            break;
                        case 'k':
                            k_freq = atoi(argv[i]);
                            if(k_freq!=0){
                                is_freq=true;
                                cout<<"k_freq = "<<k_freq<<endl;
                            }
                            break;
                        case 'l':
                            k_bet = atoi(argv[i]);
                            if(k_bet!=0){
                                is_bet=true;
                                cout<<"k_bet = "<<k_bet<<endl;
                            }
                            break;
                        case 'v':
                            strcpy(queryPairFileName,argv[i]);
                            std::cout<<"queryPairFileName="<<queryPairFileName<<std::endl;
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

                if(t_directed_flag){
                    std::cout<<"directed ";
                    if(t_weighted_flag){
                        std::cout<<"weighted graph!";
                        //coeffient struct
                        Processing::calcCoefficient<double> calcCoef((double)k_deg/100,(double)k_freq/100,(double)k_cov/100,(double)k_dep/100,(double)k_bet/100,is_deg, is_freq, is_cov, is_dep,is_bet);
                        if(bestDegree==-1) bestDegree=maxDegree;
                        if(is_multiThreads){
                            std::cout<<"Multiple threads..."<<std::endl;
                            Integrated_construction integrated_construction(graphFileName,labelFileName,labelSizeFileName,outputDirName,analysisDirName,queryFreqFileName,queryPairFileName,hfRate,maxDegree,bestDegree,t_indexing_flag,t_ordering_flag,calcCoef,t_experiment_flag,num_threads,graphType,givenOrderFileName,is_sorted);
                        }else{
                            std::cout<<"Single threads..."<<std::endl;
                            Integrated_construction integrated_construction(graphFileName,labelFileName,labelSizeFileName,outputDirName,analysisDirName,queryFreqFileName,queryPairFileName,hfRate,maxDegree,bestDegree,t_indexing_flag,t_ordering_flag,calcCoef,t_experiment_flag,graphType,givenOrderFileName,is_sorted);
                        }
                    }else{
                        std::cout<<"unweighted  graph!"<<std::endl;
                    }
                }else{
                    std::cout<<"undirected ";
                    if(t_weighted_flag){
                        std::cout<<"weighted graph!";
                        //coeffient struct
                        Processing::calcCoefficient<double> calcCoef((double)k_deg/100,(double)k_freq/100,(double)k_cov/100,(double)k_dep/100,(double)k_bet/100,is_deg, is_freq, is_cov, is_dep,is_bet);
                        //bestDegree to save overlay graph
                        if(bestDegree==-1) bestDegree=maxDegree;
                        if(is_multiThreads) Integrated_construction integrated_construction(graphFileName,labelFileName,labelSizeFileName,outputDirName,analysisDirName,queryFreqFileName,hfRate,maxDegree,bestDegree,t_indexing_flag,t_ordering_flag,calcCoef,t_experiment_flag,num_threads,givenOrderFileName,is_sorted);
                        else Integrated_construction integrated_construction(graphFileName,labelFileName,labelSizeFileName,outputDirName,analysisDirName,queryFreqFileName,hfRate,maxDegree,bestDegree,t_indexing_flag,t_ordering_flag,calcCoef,t_experiment_flag,givenOrderFileName,is_sorted);
                    }else{std::cout<<"unweighted  graph!"<<std::endl;}
                }

                return 0;
            }
    };
}

#endif