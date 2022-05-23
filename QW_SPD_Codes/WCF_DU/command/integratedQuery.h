/*
 * @Descripttion: command header file used to query distances for integrated solutions
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-04-13 09:22:15
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-11-03 14:49:28
 */
#pragma once
#ifndef _COMMAND_INTEGRATED_QUERY
#define _COMMAND_INTEGRATED_QUERY

#include <iostream>
#include <fstream>
#include <cstring>
#include "../src/time_util.h"
#include "../command.h"
#include "../src/integrated_index.h"

using namespace time_util; 
namespace command{
    class IntegrationQuery:public Command{
        public:
            void exit_with_help(){
                printf("Usage:\n");
                printf("\tmy_sspexp_run -l -d [directedGraphFlag] -w [weightedGraphFlag] -s [specialFlag] -m [indexingSchemes] -o [OrderingSchemes] -g [graphFileName] \n -e [exportLabelFileName] [-q [query_freq_file]] [-h [HFPoint_file]] [-p label_list_file] [-i [label_size_file]] [-f size_analysis_file\n [-b betweenness_filename] [-c coverage_filename] [-r high_frequency rate default-5%%] [-j -k -l -u -v coeffient of params(degree,frequency,betwenness,coverage,depth;0~10)]\n");
                printf("-------------------------------------------------------------------\n");
            }

            int main(int argc, char* argv[]){
                int numQuery = 0;
                int warmup = numQuery / 2;
                char graphFileName[255] = "";
                char indexFileName[255] = "";
                char  variableFileName[255] = "";
                char outputAnaDirName[255]=""; 
                char qdisFileName[255]=""; 
                char queryDataFileName[255]="";
                char isDeletedFileName[255]="";
                char labelSizeFileName[255]="";
                int t_directed_flag = 0;
                int t_weighted_flag = 0;
                int t_special_flag = 0;
                int t_batch_flag = 0;
                int dataType=0;//0-random 1-realdata 2-simulation
                int queryModel=0;//0-first solution,1-second solution
                bool isOutputAnaly=false;
                bool isOutputDis=false;
                bool isOutputQc=false;//indicate whether output query cost(true if input label size filename)
                int multiQuery = 1;
                int isShuffle=0;
                int hfRate=0;
                double time_rmq_1=0;

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
                            t_special_flag = atoi(argv[i]);
                            break;
                        case 'm':
                            queryModel = atoi(argv[i]);
                            break;
                        case 'p':
                            dataType = atoi(argv[i]);
                            break;
                        case 'b':
                            t_batch_flag = atoi(argv[i]);
                            break;
                        case 'x':
                           isShuffle = atoi(argv[i]);
                           cout<<"isShuffle="<<isShuffle<<std::endl;
                            break;
                        case 'n':
                            numQuery = atoi(argv[i]);
                            warmup = numQuery / 2;
                            break;
                        case 'u':
                            multiQuery = atoi(argv[i]);
                            std::cout<<"multiQuery = "<<multiQuery<<std::endl;
                            break;
                        case 'g':
                            strcpy(graphFileName,argv[i]);
                            std::cout<<"graphFileName="<<graphFileName<<std::endl;
                            break;
                        case 'y':
                            strcpy(isDeletedFileName, argv[i]);
                            cout<<"isDeletedFileName = "<<isDeletedFileName<<std::endl;
                            break;
                        case 'q':
                            strcpy(queryDataFileName,argv[i]); ///modified by wanjingyi
                            cout<<"queryDataFileName = "<<queryDataFileName<<std::endl;
                            break;
                        case 'a':
                            strcpy(outputAnaDirName,argv[i]); 
                            isOutputAnaly=true;
                            cout<<"outputAnaDirName = "<<outputAnaDirName<<std::endl;
                            break;
                        case 'f':
                            strcpy(qdisFileName,argv[i]); 
                            isOutputDis=true;
                            cout<<"qdisFileName = "<<qdisFileName<<std::endl;
                            break;
                        case 'i':
                            strcpy(indexFileName, argv[i]);
                            cout<<"indexFileName = "<<indexFileName<<std::endl;
                            break;
                        case 'v':
                            strcpy(variableFileName, argv[i]);
                            cout<<"variableFileName = "<<variableFileName<<endl;
                            break;
                        case 'e':
                            strcpy(labelSizeFileName, argv[i]);
                            if(*labelSizeFileName!='\0'){
                                std::cout<<"labelSizeFileName = "<<labelSizeFileName<<endl;
                                isOutputQc=true;
                            }
                            break;
                        default:
                            exit_with_help();
                            break;
                    }
                }

                if (t_directed_flag == 1)
                    DIRECTED_FLAG = true;
                if (t_weighted_flag == 1)
                    WEIGHTED_FLAG = true;

                if(dataType==1&&(*queryDataFileName=='\0')){
                    std::cout<<"queryDataFileName cannot be null!"<<std::endl;
                    exit_with_help();
                }
                if(t_batch_flag==1&&!is_directory(queryDataFileName)){
                    std::cout<<"Please input a query directory!"<<std::endl;
                    exit_with_help();
                }
                
                srand(31101982);//random seed

                if(t_special_flag == 0){ // default labels
                    cout<<"t_special_flag=0-default labels ";
                    if (DIRECTED_FLAG == true){//directed 
                        cout<<"DIRECTED_FLAG=true ";
                        if(WEIGHTED_FLAG){//weighted
                            cout<<"WEIGHTED_FLAG=true ";
                            if(queryModel==0){//the first integrated solution
                                cout<<"queryModel=0-the first integrated solution";
                            }else if(queryModel==1){
                                cout<<"queryModel=1-the second integrated solution";
                                //*****************Load index*****************
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels_1_directed(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                time_rmq_1=GetCurrentTimeSec();
                                integrated_index.load_lca_variables(variableFileName);
                                time_rmq_1=GetCurrentTimeSec()-time_rmq_1;
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                if(isOutputQc) integrated_index.load_label_size_directed(labelSizeFileName);
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag == 0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        //load query files and get its size
                                        string queryDataPath(queryDataFileName);
                                        vector<string> query_data_files;
                                        get_filelist_from_dir(queryDataPath,query_data_files,true);
                                        append_to_full_path(queryDataPath,query_data_files);
                                        int numOfQueryDataFile=query_data_files.size();
                                        //initialize variables
                                        double qtime,avg_qtime;
                                        vector<double> qtime_slices;
                                        vector<double> avg_qtime_slices;
                                        vector<double> query_cost_slices;
                                        vector<pair<int, int> > queries_data;
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            char* query_data_file=(char*)query_data_files[i].c_str();
                                            //load real query data
                                            if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                            else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                            warmup=numQuery/2;
                                            vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                            if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                            vector<pair<int, int> > queries(queries_data);
                                            queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                            vector<pair<int, int> > ().swap(queries_data);
                                            queries_data.clear();
                                            //warm up
                                            for(int i = 0; i < warmup; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	
                                                integrated_index.query_integrated_1_directed(s,t);
                                            }
                                            qtime=GetCurrentTimeSec();
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	 
                                                distances[i-warmup]=integrated_index.query_integrated_1_directed(s,t);
                                            }
                                            qtime = GetCurrentTimeSec() - qtime;
                                            avg_qtime = qtime/ (double)numQuery;
                                            qtime_slices.push_back(qtime);
                                            avg_qtime_slices.push_back(avg_qtime);
                                            if(isOutputQc){
                                                double qc=integrated_index.get_query_cost_directed(query_data_file);
                                                query_cost_slices.push_back(qc);
                                            }
                                        }
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            double query_cost_all_slices=0;
                                            double qtime_all_slices=0;
                                            double avg_qtime_all_slices=0;
                                            //output query time
                                            for(int i=0;i<numOfQueryDataFile;++i){
                                                if(isOutputQc){
                                                    ofs_time<<query_cost_slices[i]<<" ";
                                                    query_cost_all_slices+=query_cost_slices[i];
                                                }
                                                qtime_all_slices+=qtime_slices[i];
                                                avg_qtime_all_slices+=avg_qtime_slices[i];
                                                ofs_time<<qtime_slices[i] * 1e6<<" "<<avg_qtime_slices[i] * 1e6 <<endl;
                                            }
                                            double query_cost_avg_slice=query_cost_all_slices/(double)numOfQueryDataFile;
                                            double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfQueryDataFile;
                                            ofs_time<<"numOfQueryDataFile="<<numOfQueryDataFile<<endl;
                                            if(isOutputQc) ofs_time<<"query_cost_all_slices="<<query_cost_all_slices<<" query_cost_avg_slice="<<query_cost_avg_slice<<endl;
                                            ofs_time  <<"qtime_all_slices="<<qtime_all_slices * 1e6 << " " <<"avg_qtime_avg_slice="<<avg_qtime_avg_slice * 1e6 <<std::endl;
                                            ofs_time.close();
                                        }
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }else if(queryModel==2){
                                cout<<"queryModel=2-the third integrated solution";
                                //*****************Load index*****************
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels_2_directed(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                time_rmq_1=GetCurrentTimeSec();
                                integrated_index.load_lca_variables(variableFileName);
                                time_rmq_1=GetCurrentTimeSec()-time_rmq_1;
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                if(isOutputQc) integrated_index.load_label_size_directed(labelSizeFileName);
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        //load query files and get its size
                                        string queryDataPath(queryDataFileName);
                                        vector<string> query_data_files;
                                        get_filelist_from_dir(queryDataPath,query_data_files,true);
                                        append_to_full_path(queryDataPath,query_data_files);
                                        int numOfQueryDataFile=query_data_files.size();
                                        //cout<<"numOfQueryDataFile = "<<numOfQueryDataFile<<endl;
                                        //initialize variables
                                        double qtime,avg_qtime;
                                        vector<double> qtime_slices;
                                        vector<double> avg_qtime_slices;
                                        vector<double> query_cost_slices;
                                        vector<pair<int, int> > queries_data;
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            char* query_data_file=(char*)query_data_files[i].c_str();
                                            //load real query data
                                            if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                            else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                            warmup=numQuery/2;
                                            vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                            if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                            vector<pair<int, int> > queries(queries_data);
                                            queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                            vector<pair<int, int> > ().swap(queries_data);
                                            queries_data.clear();
                                            //warm up
                                            for(int i = 0; i < warmup; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	
                                                integrated_index.query_integrated_2_directed(s,t);
                                            }
                                            qtime=GetCurrentTimeSec();
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	 
                                                distances[i-warmup]=integrated_index.query_integrated_2_directed(s,t);
                                            }
                                            qtime = GetCurrentTimeSec() - qtime;
                                            avg_qtime = qtime/ (double)numQuery;
                                            qtime_slices.push_back(qtime);
                                            avg_qtime_slices.push_back(avg_qtime);
                                            if(isOutputQc){
                                                double qc=integrated_index.get_query_cost_directed(query_data_file);
                                                query_cost_slices.push_back(qc);
                                            }
                                        }
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            double query_cost_all_slices=0;
                                            double qtime_all_slices=0;
                                            double avg_qtime_all_slices=0;
                                            //output query time
                                            for(int i=0;i<numOfQueryDataFile;++i){
                                                if(isOutputQc){
                                                    ofs_time<<query_cost_slices[i]<<" ";
                                                    query_cost_all_slices+=query_cost_slices[i];
                                                }
                                                qtime_all_slices+=qtime_slices[i];
                                                avg_qtime_all_slices+=avg_qtime_slices[i];
                                                ofs_time<<qtime_slices[i] * 1e6<<" "<<avg_qtime_slices[i] * 1e6 <<endl;
                                            }
                                            double query_cost_avg_slice=query_cost_all_slices/(double)numOfQueryDataFile;
                                            double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfQueryDataFile;
                                            ofs_time<<"numOfQueryDataFile="<<numOfQueryDataFile<<endl;
                                            if(isOutputQc) ofs_time<<"query_cost_all_slices="<<query_cost_all_slices<<" query_cost_avg_slice="<<query_cost_avg_slice<<endl;
                                            ofs_time  <<"qtime_all_slices="<<qtime_all_slices * 1e6 << " " <<"avg_qtime_avg_slice="<<avg_qtime_avg_slice * 1e6 <<std::endl;
                                            ofs_time.close();
                                        }
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }
                        }else{
                            cout<<"WEIGHTED_FLAG=false "<<endl;
                        }
                    }else{//undirected
                        cout<<"DIRECTED_FLAG=false ";
                        if(WEIGHTED_FLAG){//weighted
                            cout<<"WEIGHTED_FLAG=true ";
                            if(queryModel==0){//the first integrated solution
                                cout<<"queryModel=0-defualt ";
                                //*****************Load index*****************
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                time_rmq_1=GetCurrentTimeSec();
                                integrated_index.load_lca_variables(variableFileName);
                                time_rmq_1=GetCurrentTimeSec()-time_rmq_1;
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }else if(queryModel==1){//the second integrated solution
                                cout<<"queryModel=1-the second integrated solution";
                                //*****************Load index*****************
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels_1(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                time_rmq_1=GetCurrentTimeSec();
                                integrated_index.load_lca_variables(variableFileName);
                                time_rmq_1=GetCurrentTimeSec()-time_rmq_1;
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                if(isOutputAnaly) integrated_index.load_label_size(labelSizeFileName);
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated_1(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated_1(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        //load query files and get its size
                                        string queryDataPath(queryDataFileName);
                                        vector<string> query_data_files;
                                        get_filelist_from_dir(queryDataPath,query_data_files,true);
                                        append_to_full_path(queryDataPath,query_data_files);
                                        int numOfQueryDataFile=query_data_files.size();
                                        //initialize variables
                                        double qtime,avg_qtime;
                                        vector<double> qtime_slices;
                                        vector<double> avg_qtime_slices;
                                        vector<double> query_cost_slices;
                                        vector<pair<int, int> > queries_data;
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            char* query_data_file=(char*)query_data_files[i].c_str();
                                            //load real query data
                                            if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                            else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                            warmup=numQuery/2;
                                            vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                            if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                            vector<pair<int, int> > queries(queries_data);
                                            queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                            vector<pair<int, int> > ().swap(queries_data);
                                            queries_data.clear();
                                            //warm up
                                            for(int i = 0; i < warmup; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	
                                                integrated_index.query_integrated_1(s,t);
                                            }
                                            qtime=GetCurrentTimeSec();
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	 
                                                distances[i-warmup]=integrated_index.query_integrated_1(s,t);
                                            }
                                            qtime = GetCurrentTimeSec() - qtime;
                                            avg_qtime = qtime/ (double)numQuery;
                                            qtime_slices.push_back(qtime);
                                            avg_qtime_slices.push_back(avg_qtime);
                                            if(isOutputQc){
                                                double qc=integrated_index.get_query_cost(query_data_file);
                                                query_cost_slices.push_back(qc);
                                            }
                                        }
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            double query_cost_all_slices=0;
                                            double qtime_all_slices=0;
                                            double avg_qtime_all_slices=0;
                                            //output query time
                                            for(int i=0;i<numOfQueryDataFile;++i){
                                                if(isOutputQc){
                                                    ofs_time<<query_cost_slices[i]<<" ";
                                                    query_cost_all_slices+=query_cost_slices[i];
                                                }
                                                qtime_all_slices+=qtime_slices[i];
                                                avg_qtime_all_slices+=avg_qtime_slices[i];
                                                ofs_time<<qtime_slices[i] * 1e6<<" "<<avg_qtime_slices[i] * 1e6 <<endl;
                                            }
                                            double query_cost_avg_slice=query_cost_all_slices/(double)numOfQueryDataFile;
                                            double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfQueryDataFile;
                                            ofs_time<<"numOfQueryDataFile="<<numOfQueryDataFile<<endl;
                                            if(isOutputQc) ofs_time<<"query_cost_all_slices="<<query_cost_all_slices<<" query_cost_avg_slice="<<query_cost_avg_slice<<endl;
                                            ofs_time  <<"qtime_all_slices="<<qtime_all_slices * 1e6 << " " <<"avg_qtime_avg_slice="<<avg_qtime_avg_slice * 1e6 <<std::endl;
                                            ofs_time.close();
                                        }
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated_1(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }else if(queryModel==2){
                                cout<<"queryModel 2-the third integrated solution";
                                //*****************Load index*****************
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels_2(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                time_rmq_1=GetCurrentTimeSec();
                                integrated_index.load_lca_variables(variableFileName);
                                time_rmq_1=GetCurrentTimeSec()-time_rmq_1;
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                if(isOutputQc) integrated_index.load_label_size(labelSizeFileName);
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated_2(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated_2(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        //load query files and get its size
                                        string queryDataPath(queryDataFileName);
                                        vector<string> query_data_files;
                                        get_filelist_from_dir(queryDataPath,query_data_files,true);
                                        append_to_full_path(queryDataPath,query_data_files);
                                        int numOfQueryDataFile=query_data_files.size();
                                        //initialize variables
                                        double qtime,avg_qtime;
                                        vector<double> qtime_slices;
                                        vector<double> avg_qtime_slices;
                                        vector<double> query_cost_slices;
                                        vector<pair<int, int> > queries_data;
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            char* query_data_file=(char*)query_data_files[i].c_str();
                                            //load real query data
                                            if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                            else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                            warmup=numQuery/2;
                                            vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                            if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                            vector<pair<int, int> > queries(queries_data);
                                            queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                            vector<pair<int, int> > ().swap(queries_data);
                                            queries_data.clear();
                                            //warm up
                                            for(int i = 0; i < warmup; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	
                                                integrated_index.query_integrated_2(s,t);
                                            }
                                            qtime=GetCurrentTimeSec();
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	 
                                                distances[i-warmup]=integrated_index.query_integrated_2(s,t);
                                            }
                                            qtime = GetCurrentTimeSec() - qtime;
                                            avg_qtime = qtime/ (double)numQuery;
                                            qtime_slices.push_back(qtime);
                                            avg_qtime_slices.push_back(avg_qtime);
                                            if(isOutputQc){
                                                double qc=integrated_index.get_query_cost(query_data_file);
                                                query_cost_slices.push_back(qc);
                                            }
                                        }
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            double query_cost_all_slices=0;
                                            double qtime_all_slices=0;
                                            double avg_qtime_all_slices=0;
                                            //output query time
                                            for(int i=0;i<numOfQueryDataFile;++i){
                                                if(isOutputQc){
                                                    ofs_time<<query_cost_slices[i]<<" ";
                                                    query_cost_all_slices+=query_cost_slices[i];
                                                }
                                                qtime_all_slices+=qtime_slices[i];
                                                avg_qtime_all_slices+=avg_qtime_slices[i];
                                                ofs_time<<qtime_slices[i] * 1e6<<" "<<avg_qtime_slices[i] * 1e6 <<endl;
                                            }
                                            double query_cost_avg_slice=query_cost_all_slices/(double)numOfQueryDataFile;
                                            double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfQueryDataFile;
                                            ofs_time<<"numOfQueryDataFile="<<numOfQueryDataFile<<endl;
                                            if(isOutputQc) ofs_time<<"query_cost_all_slices="<<query_cost_all_slices<<" query_cost_avg_slice="<<query_cost_avg_slice<<endl;
                                            ofs_time  <<"qtime_all_slices="<<qtime_all_slices * 1e6 << " " <<"avg_qtime_avg_slice="<<avg_qtime_avg_slice * 1e6 <<std::endl;
                                            ofs_time.close();
                                        }
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        cout<<"generate query data successfully!"<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2(s,t);
                                        }
                                        cout<<"warm up successfully!"<<endl;
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated_2(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Build rmq time:" << time_rmq_1 * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            string analy_filename_prefix(outputAnaDirName);
                                            string analy_filename=analy_filename_prefix.append(".qtime");
                                            ofstream ofs_time(analy_filename);
                                            if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                            ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<" time_rmq_1="<<time_rmq_1 * 1e6<<std::endl;
                                            ofs_time.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }
                        }else{//unweighted
                        }
                    }
                }else if(t_special_flag==1){
                    cout<<"t_special_flag=1-experiment ";
                    if (DIRECTED_FLAG == true){//directed 
                        cout<<"DIRECTED_FLAG=true ";
                        if(WEIGHTED_FLAG){//weighted
                            cout<<"WEIGHTED_FLAG=true ";
                            if(queryModel==0){//the first integrated solution
                                cout<<"queryModel=0-defualt ";
                            }else if(queryModel==1){
                                cout<<"queryModel=1-the second integrated solution";
                                //*****************Load index*****************
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels_1_directed(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                time_rmq_1=GetCurrentTimeSec();
                                integrated_index.load_lca_variables(variableFileName);
                                time_rmq_1=GetCurrentTimeSec()-time_rmq_1;
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                if(isOutputQc) integrated_index.load_label_size_directed(labelSizeFileName);
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        //*****************Load index*****************
                                        //load query files and get its size
                                        string queryDataPath(queryDataFileName);
                                        vector<string> query_data_files;
                                        get_filelist_from_dir(queryDataPath,query_data_files,true);
                                        append_to_full_path(queryDataPath,query_data_files);
                                        int numOfQueryDataFile=query_data_files.size();
                                        //initialize variables
                                        double qtime,avg_qtime;
                                        vector<double> qtime_slices;
                                        vector<double> avg_qtime_slices;
                                        vector<double> query_cost_slices;
                                        vector<pair<int, int> > queries_data;
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            char* query_data_file=(char*)query_data_files[i].c_str();
                                            //load real query data
                                            if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                            else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                            warmup=numQuery/2;
                                            vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                            if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                            vector<pair<int, int> > queries(queries_data);
                                            queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                            vector<pair<int, int> > ().swap(queries_data);
                                            queries_data.clear();
                                            //warm up
                                            for(int i = 0; i < warmup; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	
                                                integrated_index.query_integrated_1_directed(s,t);
                                            }
                                            qtime=GetCurrentTimeSec();
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	 
                                                distances[i-warmup]=integrated_index.query_integrated_1_directed(s,t);
                                            }
                                            qtime = GetCurrentTimeSec() - qtime;
                                            avg_qtime = qtime/ (double)numQuery;
                                            qtime_slices.push_back(qtime);
                                            avg_qtime_slices.push_back(avg_qtime);
                                            if(isOutputQc){
                                                double qc=integrated_index.get_query_cost_directed(query_data_file);
                                                query_cost_slices.push_back(qc);
                                            }
                                        }    
                                        //output query time
                                        if(isOutputAnaly){
                                            double query_cost_all_slices=0;
                                            double qtime_all_slices=0;
                                            double avg_qtime_all_slices=0;
                                            for(int i=0;i<numOfQueryDataFile;++i){
                                                if(isOutputQc){
                                                    query_cost_all_slices+=query_cost_slices[i];
                                                }
                                                qtime_all_slices+=qtime_slices[i];
                                                avg_qtime_all_slices+=avg_qtime_slices[i];
                                            }
                                            double query_cost_avg_slice=query_cost_all_slices/(double)numOfQueryDataFile;
                                            double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfQueryDataFile;
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs.setf(ios::fixed);
                                            ofs.precision(4);
                                            ofs<<numOfQueryDataFile<<" "<<endl;
                                            if(isOutputQc) ofs<<query_cost_all_slices<<" "<<query_cost_avg_slice<<" ";
                                            ofs  << qtime_all_slices * 1e6 <<  " " <<avg_qtime_avg_slice * 1e6 <<" " ;
                                            ofs.close();
                                        }
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated_1_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }else if(queryModel==2){
                                cout<<"queryModel=2-the third integrated solution ";
                                //*****************Load index*****************
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels_2_directed(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                integrated_index.load_lca_variables(variableFileName);
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                if(isOutputQc) integrated_index.load_label_size_directed(labelSizeFileName);
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        //load query files and get its size
                                        string queryDataPath(queryDataFileName);
                                        vector<string> query_data_files;
                                        get_filelist_from_dir(queryDataPath,query_data_files,true);
                                        append_to_full_path(queryDataPath,query_data_files);
                                        int numOfQueryDataFile=query_data_files.size();
                                        //initialize variables
                                        double qtime,avg_qtime;
                                        vector<double> qtime_slices;
                                        vector<double> avg_qtime_slices;
                                        vector<double> query_cost_slices;
                                        vector<pair<int, int> > queries_data;
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            char* query_data_file=(char*)query_data_files[i].c_str();
                                            //load real query data
                                            if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                            else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                            warmup=numQuery/2;
                                            vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                            if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                            vector<pair<int, int> > queries(queries_data);
                                            queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                            vector<pair<int, int> > ().swap(queries_data);
                                            queries_data.clear();
                                            //warm up
                                            for(int i = 0; i < warmup; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	
                                                integrated_index.query_integrated_2_directed(s,t);
                                            }
                                            qtime=GetCurrentTimeSec();
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	 
                                                distances[i-warmup]=integrated_index.query_integrated_2_directed(s,t);
                                            }
                                            qtime = GetCurrentTimeSec() - qtime;
                                            avg_qtime = qtime/ (double)numQuery;
                                            qtime_slices.push_back(qtime);
                                            avg_qtime_slices.push_back(avg_qtime);
                                            if(isOutputQc){
                                                double qc=integrated_index.get_query_cost_directed(query_data_file);
                                                query_cost_slices.push_back(qc);
                                            }
                                        }    
                                        //output query time
                                        if(isOutputAnaly){
                                            double query_cost_all_slices=0;
                                            double qtime_all_slices=0;
                                            double avg_qtime_all_slices=0;
                                            for(int i=0;i<numOfQueryDataFile;++i){
                                                if(isOutputQc){
                                                    query_cost_all_slices+=query_cost_slices[i];
                                                }
                                                qtime_all_slices+=qtime_slices[i];
                                                avg_qtime_all_slices+=avg_qtime_slices[i];
                                            }
                                            double query_cost_avg_slice=query_cost_all_slices/(double)numOfQueryDataFile;
                                            double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfQueryDataFile;
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs.setf(ios::fixed);
                                            ofs.precision(4);
                                            ofs<<numOfQueryDataFile<<" "<<endl;
                                            if(isOutputQc) ofs<<query_cost_all_slices<<" "<<query_cost_avg_slice<<" ";
                                            ofs  << qtime_all_slices * 1e6 <<  " " <<avg_qtime_avg_slice * 1e6 <<" " ;
                                            ofs.close();
                                        }
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated_2_directed(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }
                        }else{
                            cout<<"WEIGHTED_FLAG=false ";
                        }  
                    }else{//undirected
                        cout<<"DIRECTED_FLAG=false ";
                        if(WEIGHTED_FLAG){//weighted
                            cout<<"WEIGHTED_FLAG=true ";
                            if(queryModel==0){//the first integrated solution
                                cout<<"queryModel=0-defualt ";
                                //*****************Load index*****************
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                time_rmq_1=GetCurrentTimeSec();
                                integrated_index.load_lca_variables(variableFileName);
                                time_rmq_1=GetCurrentTimeSec()-time_rmq_1;
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                         //output query time
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                         //output query time
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                         //output query time
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }else if(queryModel==1){//the second integrated solution
                                cout<<"queryModel=1-the second integrated solution ";
                                //*****************Load index*****************
                                cout<<"t_batch_flag==0-oneQuery "<<endl;
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels_1(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                integrated_index.load_lca_variables(variableFileName);
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                if(isOutputQc) integrated_index.load_label_size(labelSizeFileName);
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag == 0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated_1(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        //*****************Load index*****************
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated_1(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        string queryDataPath(queryDataFileName);
                                        vector<string> query_data_files;
                                        get_filelist_from_dir(queryDataPath,query_data_files,true);
                                        append_to_full_path(queryDataPath,query_data_files);
                                        int numOfQueryDataFile=query_data_files.size();
                                        //initialize variables
                                        double qtime,avg_qtime;
                                        vector<double> qtime_slices;
                                        vector<double> avg_qtime_slices;
                                        vector<double> query_cost_slices;
                                        vector<pair<int, int> > queries_data;
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            char* query_data_file=(char*)query_data_files[i].c_str();
                                            //load real query data
                                            if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                            else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                            warmup=numQuery/2;
                                            vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                            if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                            vector<pair<int, int> > queries(queries_data);
                                            queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                            vector<pair<int, int> > ().swap(queries_data);
                                            queries_data.clear();
                                            //warm up
                                            for(int i = 0; i < warmup; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	
                                                integrated_index.query_integrated_1(s,t);
                                            }
                                            qtime=GetCurrentTimeSec();
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	 
                                                distances[i-warmup]=integrated_index.query_integrated_1(s,t);
                                            }
                                            qtime = GetCurrentTimeSec() - qtime;
                                            avg_qtime = qtime/ (double)numQuery;
                                            qtime_slices.push_back(qtime);
                                            avg_qtime_slices.push_back(avg_qtime);
                                            if(isOutputQc){
                                                double qc=integrated_index.get_query_cost(query_data_file);
                                                query_cost_slices.push_back(qc);
                                            }
                                        }    
                                        //output query time
                                        if(isOutputAnaly){
                                            double query_cost_all_slices=0;
                                            double qtime_all_slices=0;
                                            double avg_qtime_all_slices=0;
                                            for(int i=0;i<numOfQueryDataFile;++i){
                                                if(isOutputQc){
                                                    query_cost_all_slices+=query_cost_slices[i];
                                                }
                                                qtime_all_slices+=qtime_slices[i];
                                                avg_qtime_all_slices+=avg_qtime_slices[i];
                                            }
                                            double query_cost_avg_slice=query_cost_all_slices/(double)numOfQueryDataFile;
                                            double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfQueryDataFile;
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs.setf(ios::fixed);
                                            ofs.precision(4);
                                            ofs<<numOfQueryDataFile<<" "<<endl;
                                            if(isOutputQc) ofs<<query_cost_all_slices<<" "<<query_cost_avg_slice<<" ";
                                            ofs  << qtime_all_slices * 1e6 <<  " " <<avg_qtime_avg_slice * 1e6 <<" " ;
                                            ofs.close();
                                        }
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_1(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated_1(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }else if(queryModel==2){//the third integrated solution
                                cout<<"queryModel=2-the third integrated solution ";
                                //*****************Load index*****************
                                Integrated_Index integrated_index;
                                integrated_index.load_is_deleted(isDeletedFileName);
                                std::cout<<"Load is_deleted successfully!"<<std::endl;
                                integrated_index.load_labels_2(indexFileName);
                                std::cout<<"Load labels successfully!"<<std::endl;
                                integrated_index.load_lca_variables(variableFileName);
                                std::cout<<"Load lca_variables successfully!"<<std::endl;
                                if(isOutputQc) integrated_index.load_label_size(labelSizeFileName);
                                //*****************Load index*****************
                                if(dataType==0){//random
                                    cout<<"dataType=0-random ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries(numQuery + warmup);
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        //generate random query data
                                        for(int i = 0; i < numQuery + warmup; ++i){
                                            int s = rand()%numOfVertices;
                                            int t = rand()%numOfVertices;
                                            if(s == t){
                                                i--;
                                                continue;
                                            }		
                                            queries[i]= make_pair(s, t);
                                        }
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup + numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                           distances[i-warmup]=integrated_index.query_integrated_2(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(qdisFileName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup + numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }else if(dataType==1){//realdata
                                    cout<<"dataType=1-realdata ";
                                    if(t_batch_flag==0){
                                        cout<<"t_batch_flag==0-oneQuery "<<endl;
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries_data;
                                        //load real query data
                                        if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=integrated_index.query_integrated_2(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        string queryDataPath(queryDataFileName);
                                        vector<string> query_data_files;
                                        get_filelist_from_dir(queryDataPath,query_data_files,true);
                                        append_to_full_path(queryDataPath,query_data_files);
                                        int numOfQueryDataFile=query_data_files.size();
                                        //initialize variables
                                        double qtime,avg_qtime;
                                        vector<double> qtime_slices;
                                        vector<double> avg_qtime_slices;
                                        vector<double> query_cost_slices;
                                        vector<pair<int, int> > queries_data;
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            char* query_data_file=(char*)query_data_files[i].c_str();
                                            //load real query data
                                            if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                            else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                            warmup=numQuery/2;
                                            vector<int> distances(numQuery,INF_WEIGHT);//output distances
                                            if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());
                                            vector<pair<int, int> > queries(queries_data);
                                            queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                            vector<pair<int, int> > ().swap(queries_data);
                                            queries_data.clear();
                                            //warm up
                                            for(int i = 0; i < warmup; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	
                                                integrated_index.query_integrated_2(s,t);
                                            }
                                            qtime=GetCurrentTimeSec();
                                            for(int i = warmup; i < warmup+numQuery; ++i){
                                                int s = queries[i].first;
                                                int t = queries[i].second;	 
                                                distances[i-warmup]=integrated_index.query_integrated_2(s,t);
                                            }
                                            qtime = GetCurrentTimeSec() - qtime;
                                            avg_qtime = qtime/ (double)numQuery;
                                            qtime_slices.push_back(qtime);
                                            avg_qtime_slices.push_back(avg_qtime);
                                            if(isOutputQc){
                                                double qc=integrated_index.get_query_cost(query_data_file);
                                                query_cost_slices.push_back(qc);
                                            }
                                        }    
                                        //output query time
                                        if(isOutputAnaly){
                                            double query_cost_all_slices=0;
                                            double qtime_all_slices=0;
                                            double avg_qtime_all_slices=0;
                                            for(int i=0;i<numOfQueryDataFile;++i){
                                                if(isOutputQc){
                                                    query_cost_all_slices+=query_cost_slices[i];
                                                }
                                                qtime_all_slices+=qtime_slices[i];
                                                avg_qtime_all_slices+=avg_qtime_slices[i];
                                            }
                                            double query_cost_avg_slice=query_cost_all_slices/(double)numOfQueryDataFile;
                                            double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfQueryDataFile;
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs.setf(ios::fixed);
                                            ofs.precision(4);
                                            ofs<<numOfQueryDataFile<<" "<<endl;
                                            if(isOutputQc) ofs<<query_cost_all_slices<<" "<<query_cost_avg_slice<<" ";
                                            ofs  << qtime_all_slices * 1e6 <<  " " <<avg_qtime_avg_slice * 1e6 <<" " ;
                                            ofs.close();
                                        }
                                        return 0;
                                    }
                                }else if(dataType==2){//simulation
                                    cout<<"dataType=2-simulation ";
                                    if(t_batch_flag==0){
                                        //variables
                                        double qtime,avg_qtime;
                                        vector<pair<int, int> > queries;
                                        //generate simulation query data
                                        for(int s = 0; s < numOfVertices; ++s){
                                            for(int t=s;t<numOfVertices;++t)
                                            queries.push_back(make_pair(s,t));
                                        }
                                        numQuery=queries.size();
                                        warmup=numQuery/2;
                                        vector<int> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries.begin(), queries.end());
                                        cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            integrated_index.query_integrated_2(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = 0; i < numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            distances[i]=integrated_index.query_integrated_2(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << endl;
                                        cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
                                        if(isOutputAnaly){
                                            ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                            if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                            ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" ";
                                            ofs.close();
                                        }
                                        if(isOutputDis){
                                            string dis_output_filename(outputAnaDirName);
                                            dis_output_filename.append(".distances");
                                            ofstream ofs(dis_output_filename);
                                            if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<endl;
                                            for(int i = 0; i < numQuery; ++i){
                                                ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<endl;
                                            }
                                            ofs.close();
                                        }
                                        return 0;
                                    }else if(t_batch_flag==1){
                                        cout<<"t_batch_flag==1-batchQuery "<<endl;
                                        return 0;
                                    }
                                }
                            }
                        }
                    }
                }

                return 0;
            }
    };
}

#endif