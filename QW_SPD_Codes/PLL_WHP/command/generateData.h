/*
 * @Descripttion: 
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-09-21 09:39:01
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2022-02-17 19:48:40
 */
#pragma once
#ifndef _COMMAND_VERIFY_H
#define _COMMAND_VERIFY_H

#include <iostream>
#include <fstream>
#include <cstring>
#include "../src/graph_search.h"
#include "../src/verify_graph.h"
#include "../command.h"

namespace command{
    class GenerateData: public Command{
        public:
            void exit_with_help(){
                printf("Usage:\n");
                printf("\t... \n");//./cache_run -c -d 0 -w 1 -g  ../dataset/test.graph -t ../debug/test -l ../dataset/query.log -b 10
                printf("Parameter explanation:\n");
                printf("Examples:\n");
            }

            int main(int argc,char* argv[]){
                char graphFileName[255] = "";
                char outputDirName[255] = "";
                int t_directed_flag = 0;
                int t_weighted_flag = 0;
                int num_of_threads= 5;
                int num_of_sets=10;
                int num_of_data=100000;
                int t_model_flag=0;
                int graphType=0; //0-full connected 1-u is start and t is destination

                for(int i = 2; i < argc; i++){
                    if(argv[i][0] != '-') break;
                    if(++i >= argc) exit_with_help();
                    switch (argv[i-1][1]){
                        case 'd':
                            t_directed_flag = atoi(argv[i]);
                            break;
                        case 'w':
                            t_weighted_flag = atoi(argv[i]);
                            break;
                        case 'm':
                            t_model_flag = atoi(argv[i]);
                            break;
                        case 'c':
                            graphType = atoi(argv[i]);
                            break;
                        case 'g':
                            strcpy(graphFileName,argv[i]);
                            cout<<"graphFileName = "<<graphFileName<<endl;
                            break;               
                        case 'a':
                            strcpy(outputDirName,argv[i]);
                            cout<<"outputDirName = "<<outputDirName<<endl;
                            break;
                        case 'k':
                            num_of_sets= atoi(argv[i]);
                            cout<<"num_of_sets = "<<num_of_sets<<endl;
                            break;   
                        case 't':
                            num_of_threads = atoi(argv[i]);
                            cout<<"num_of_threads = "<<num_of_threads<<endl;
                            break;   
                        case 'n':
                            num_of_data = atoi(argv[i]);
                            cout<<"num_of_data = "<<num_of_data<<endl;
                            break;   
                        default:
                            exit_with_help();
                            break;
                    }
                }

                if (t_directed_flag == 1) DIRECTED_FLAG = true;
                if (t_weighted_flag == 1) WEIGHTED_FLAG = true;


                if(DIRECTED_FLAG==true){
                    if(WEIGHTED_FLAG == true){
                        cout<<"directed weighted graph "<<endl;
                        WGraph wgraph;
                        wgraph.load_graph(graphFileName,graphType);
                        Compute_graph compute_graph(wgraph,outputDirName,t_model_flag,num_of_sets,num_of_data,num_of_threads);
                        return 0;
                    }else{}
                }else{
                    if(WEIGHTED_FLAG == true){
                        cout<<"undirected weighted graph "<<endl;
                        WGraph wgraph;
                        wgraph.load_graph(graphFileName,graphType);
                        Compute_graph compute_graph(wgraph,outputDirName,t_model_flag,num_of_sets,num_of_data,num_of_threads);
                        return 0;
                    }else{}
                }

                return 0;
            }
    };
}

#endif