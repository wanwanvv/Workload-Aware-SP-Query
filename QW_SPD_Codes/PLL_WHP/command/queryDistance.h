/*
 * @Descripttion: query distances
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-10 12:22:22
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-11-07 11:38:53
 */
#ifndef _COMMAND_QUERY_DISTANCE
#define _COMMAND_QUERY_DISTANCE

#include <cstring>
#include <string>
#include <iostream>
#include <unistd.h>
#include <vector>
#include "../src/graph.h"
#include "../src/labels.h"
#include "../src/time_util.h"
#include "../command.h"
#include "../src/utils.h"

#define QUERY_NUM 10000000
using namespace time_util; 
unsigned int sp[QUERY_NUM];
unsigned int tp[QUERY_NUM];

namespace command{
    class QueryProcessing: public Command{
        public:
            void exit_with_help(){
                printf("Usage:\n");
                printf("bin/pll_whp_run -t -d [directedGraphFlag] -w [weightedGraphFlag] -s [experimentFlag] -p [dataType] \n \
                -l [inputLabelFileName] [-t [queryPairFileName]] [-i [queryFreqFileName]] [-a [outputResultDirName]] \n \
                [-f [outputDistanceFileName]] [-u [multiQuery]] [-n [numOfQuery]] [-x [shuffleFlag]] \n");
                printf("-------------------------------------------------------------------\n");
                printf("Parameter explanation:\n");
                printf("-d: 0 or 1, for undirected and directed graphs, default is 0 \n");
                printf("-w: 0 or 1, for unweighted and weighted graphs, default is 0 \n");
                printf("-s: 0-default model, 1-experiment model \n");
                printf("-u: mutiple value of input query pairs, default is 1 \n");
                printf("-p: 0-generate query data random,1-user-defined input query data,2-gen3rate all s-t pair data for verification\n");
                printf("-------------------------------------------------------------------\n");
                
                exit(1);
            }
            
            int main(int argc, char* argv[]){
                char labelFileName[255] = "";
                int t_directed_flag = 0;
                int t_weighted_flag = 0;
                int t_special_flag = 0;
                int numQuery = 0;
                int warmup = 0;
                int t_batch_flag = 0; //0-one time slice quey 1-batch time slice query, default is 0
                char outputAnaDirName[255]=""; 
                char qdisFileName[255]=""; 
                char queryDataFileName[255]="";
                char labelSizeFileName[255] = "";
                int dataType=0;//0-random 1-realdata 2-simulation
                bool isOutputDis=false;//indicating whether to output the result 
                bool isOutputAnaly=false;
                bool isOutputQc=false;//indicate whether output query cost(true if input label size filename)
                int multiQuery = 1;//multi times
                int isShuffle=0;
                
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
                        case 's':
                            t_special_flag = atoi(argv[i]);
                            break;
                        case 'b':
                            t_batch_flag = atoi(argv[i]);
                            break;
                        case 'l':
                            strcpy(labelFileName, argv[i]);
                            cout<<"labelFileName = "<<labelFileName<<std::endl;
                            break;
                        case 'a':
                            strcpy(outputAnaDirName,argv[i]); 
                            if(*outputAnaDirName!='\0') isOutputAnaly=true;
                            cout<<"outputAnaDirName = "<<outputAnaDirName<<std::endl;
                            break;
                        case 'f':
                            strcpy(qdisFileName,argv[i]); 
                            isOutputDis=true;
                            cout<<"qdisFileName = "<<qdisFileName<<std::endl;
                            break;
                        case 't':
                            strcpy(queryDataFileName,argv[i]); ///modified by wanjingyi
                            if(*queryDataFileName=='\0') cerr<<"queryDataFileName cannot be null!"<<std::endl;
                            cout<<"queryDataFileName = "<<queryDataFileName<<std::endl;
                            break;
                        case 'z':
                            strcpy(labelSizeFileName, argv[i]);
                            if(*labelSizeFileName!='\0'){
                                std::cout<<"labelSizeFileName = "<<labelSizeFileName<<endl;
                                isOutputQc=true;
                            }
                            break;
                        case 'p':
                            dataType = atoi(argv[i]);
                            break;
                        case 'u':
                            multiQuery = atoi(argv[i]);
                            std::cout<<"multiQuery = "<<multiQuery<<std::endl;
                            break;
                        case 'x':
                           isShuffle = atoi(argv[i]);
                           cout<<"isShuffle="<<isShuffle<<std::endl;
                            break;
                        case 'n':
                            numQuery = atoi(argv[i]);
                            warmup = numQuery / 2;
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

                if(t_special_flag==0){ // default debug model
                    cout<<"t_special_flag=0-debug model";
                    if (DIRECTED_FLAG == true){//directed 
                        cout<<"DIRECTED_FLAG=true ";
                        if(WEIGHTED_FLAG){//weighted
                            cout<<"WEIGHTED_FLAG=true ";
                            if(dataType==0){//random
                                cout<<"dataType=0-random"<<std::endl;
                                //variables
                                double qtime,avg_qtime;
                                vector<pair<int, int> > queries(numQuery + warmup);
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                //load labels
                                HFDLabel lab;
                                lab.index_p=NULL;
                                lab.bindex_p=NULL;
                                lab.load_labels(labelFileName);
                                cout<<"load labels successfully!"<<std::endl;
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
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                cout<<"generate query data successfully!"<<std::endl;
                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    lab.query_p(s,t);
                                }
                                cout<<"warm up successfully!"<<std::endl;
                                qtime=GetCurrentTimeSec();
                                for(int i = warmup; i < warmup + numQuery; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    distances[i-warmup]=lab.query_p(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                //output query time
                                if(isOutputAnaly){
                                    string analy_filename_prefix(outputAnaDirName);
                                    string analy_filename=analy_filename_prefix.append(".qtime");
                                    ofstream ofs_time(analy_filename);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                    ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<std::endl;
                                    ofs_time.close();
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = warmup; i < warmup + numQuery; ++i){
                                        ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<std::endl;
                                    }
                                    ofs.close();
                                }
                                return 0;
                            }else if (dataType==1)//realData
                            {
                                cout<<"dataType=1-realData ";
                                //load labels
                                HFDLabel lab;
                                lab.index_p=NULL;
                                lab.bindex_p=NULL;
                                lab.load_labels(labelFileName);
                                if(isOutputQc) lab.load_label_size(labelSizeFileName);
                                if(t_batch_flag==1){
                                    cout<<"t_batch_flag=1-batchQuery "<<endl;
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
                                        if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<i<<" numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            lab.query_p(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i=warmup;i<warmup+numQuery;++i)
                                        {
                                            // int s = *sp_tmp++;
                                            // int t = *tp_tmp++;
                                            int s = queries[i].first;
                                            int t = queries[i].second;
                                            distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                            //lab.query(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        qtime_slices.push_back(qtime);
                                        avg_qtime_slices.push_back(avg_qtime);
                                        if(isOutputQc){
                                            double qc=lab.get_query_cost(query_data_file);
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
                                        ofs_time  <<"qtime_all_slices="<<qtime_all_slices * 1e6 << " avg_qtime_avg_slice="<<avg_qtime_avg_slice * 1e6 <<std::endl;
                                        ofs_time.close();
                                        
                                    }
                                }else{
                                    cout<<"t_batch_flag=0-oneQuery "<<endl;
                                    //variables
                                    double qtime,avg_qtime;
                                    vector<pair<int, int> > queries_data;
                                    //read real data from file
                                    //numQuery=load_query_pair_time(pointFreqFileName,queryDataFileName,queries,hfRate,numOfHFpoint,multiQuery);
                                    if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                    else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                    warmup=numQuery/2;
                                    vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                    if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                    vector<pair<int, int> > queries(queries_data);
                                    queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                    vector<pair<int, int> > ().swap(queries_data);
                                    queries_data.clear();
                                    cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                    cout<<"generate query data successfully!"<<std::endl;
                                    //warm up
                                    for(int i = 0; i < warmup; ++i){
                                        int s = queries[i].first;
                                        int t = queries[i].second;	
                                        lab.query_p(s,t);
                                    }
                                    cout<<"warm up successfully!"<<std::endl;
                                    qtime=GetCurrentTimeSec();
                                    for(int i=warmup;i<warmup+numQuery;++i)
                                    {
                                        // int s = *sp_tmp++;
                                        // int t = *tp_tmp++;
                                        int s = queries[i].first;
                                        int t = queries[i].second;
                                        distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                        //lab.query(s,t);
                                    }
                                    qtime = GetCurrentTimeSec() - qtime;
                                    avg_qtime = qtime/ (double)numQuery;
                                    cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                    cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                    //output query time
                                    if(isOutputAnaly){
                                        if((*queryDataFileName!='\0')&&isOutputQc){
                                            lab.load_hfpoint_and_qt(queryDataFileName,0);
                                            lab.append_experiment_result_pll(outputAnaDirName);
                                        }
                                        string analy_filename_prefix(outputAnaDirName);
                                        string analy_filename=analy_filename_prefix.append(".qtime");
                                        ofstream ofs_time(analy_filename);
                                        if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                        ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<std::endl;
                                        ofs_time.close();
                                        
                                    }
                                    if(isOutputDis){
                                        string dis_output_filename(qdisFileName);
                                        dis_output_filename.append(".distances");
                                        ofstream ofs(dis_output_filename);
                                        if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<std::endl;
                                        }
                                        ofs.close();
                                    }
                                }
                                return 0;

                            }else if (dataType==2)//simulation
                            {
                                cout<<"dataType=2-simulation"<<std::endl;
                                //load labels
                                HFDLabel lab;
                                lab.index_p=NULL;
                                lab.bindex_p=NULL;
                                lab.load_labels(labelFileName);
                                cout<<"load labels successfully!"<<std::endl;
                                //variables
                                double qtime,avg_qtime;
                                vector<pair<NodeID, NodeID> > queries;
                                //generate simulation query data
                                for(NodeID s = 0; s < numOfVertices; ++s){
                                    for(NodeID t=s;t<numOfVertices;++t)
                                    queries.push_back(make_pair(s,t));
                                }
                                numQuery=queries.size();
                                warmup=numQuery/2;
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                if(isShuffle) random_shuffle(queries.begin(), queries.end());
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                cout<<"generate query data successfully!"<<std::endl;
                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    lab.query_p(s,t);
                                }
                                cout<<"warm up successfully!"<<std::endl;
                                qtime=GetCurrentTimeSec();
                                for(int i = 0; i < numQuery; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    distances[i]=lab.query_p(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                //output query time
                                if(isOutputAnaly){
                                    string analy_filename_prefix(outputAnaDirName);
                                    string analy_filename=analy_filename_prefix.append(".qtime");
                                    ofstream ofs_time(analy_filename);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                    ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<std::endl;
                                    ofs_time.close();
                                    
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = 0; i < numQuery; ++i){
                                        ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<std::endl;
                                    }
                                    ofs.close();
                                }

                            }
                        }else{
                            cout<<"WEIGHTED_FLAG=false";
                        }
                    }else{//undirected
                        cout<<"DIRECTED_FLAG=false ";
                        if(WEIGHTED_FLAG){//weighted
                            cout<<"WEIGHTED_FLAG=true ";
                            if(dataType==0){//random
                                cout<<"dataType=0-random"<<std::endl;;
                                //variables
                                double qtime,avg_qtime;
                                vector<pair<int, int> > queries(numQuery + warmup);
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                //load labels
                                HFLabel lab;
                                lab.index_p=NULL;
                                lab.load_labels(labelFileName);
                                cout<<"load labels successfully!"<<std::endl;
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
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                cout<<"generate query data successfully!"<<std::endl;
                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    lab.query_p(s,t);
                                }
                                cout<<"warm up successfully!"<<std::endl;
                                qtime=GetCurrentTimeSec();
                                for(int i = warmup; i < warmup + numQuery; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    distances[i-warmup]=lab.query_p(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                //output query time
                                if(isOutputAnaly){
                                    string analy_filename_prefix(outputAnaDirName);
                                    string analy_filename=analy_filename_prefix.append(".qtime");
                                    ofstream ofs_time(analy_filename);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                    ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<std::endl;
                                    ofs_time.close();
                                    
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = warmup; i < warmup + numQuery; ++i){
                                        ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<std::endl;
                                    }
                                    ofs.close();
                                }
                                return 0;
                            }else if (dataType==1)//realData
                            {
                                cout<<"dataType=1-realData ";   
                                //load labels
                                HFLabel lab;
                                lab.index_p=NULL;
                                lab.load_labels(labelFileName);
                                if(isOutputQc) lab.load_label_size(labelSizeFileName);
                                if(t_batch_flag==1){
                                    cout<<"t_batch_flag=1-batchQuery "<<endl;
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
                                        if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<i<<" numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            lab.query_p(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i=warmup;i<warmup+numQuery;++i)
                                        {
                                            // int s = *sp_tmp++;
                                            // int t = *tp_tmp++;
                                            int s = queries[i].first;
                                            int t = queries[i].second;
                                            distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                            //lab.query(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        qtime_slices.push_back(qtime);
                                        avg_qtime_slices.push_back(avg_qtime);
                                        if(isOutputQc){
                                            double qc=lab.get_query_cost(query_data_file);
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
                                }else{
                                    cout<<"t_batch_flag=0-oneQuery"<<endl;
                                    //variables
                                    double qtime,avg_qtime;
                                    vector<pair<int, int> > queries_data;
                                    //read real data from file
                                    //numQuery=load_query_pair_time(pointFreqFileName,queryDataFileName,queries,hfRate,numOfHFpoint,multiQuery);
                                    if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                    else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                    warmup=numQuery/2;
                                    vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                    if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                    vector<pair<int, int> > queries(queries_data);
                                    queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                    vector<pair<int, int> > ().swap(queries_data);
                                    queries_data.clear();
                                    cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                    cout<<"generate query data successfully!"<<std::endl;
                                    //warm up
                                    for(int i = 0; i < warmup; ++i){
                                        int s = queries[i].first;
                                        int t = queries[i].second;	
                                        lab.query_p(s,t);
                                    }
                                    cout<<"warm up successfully!"<<std::endl;
                                    qtime=GetCurrentTimeSec();
                                    for(int i=warmup;i<warmup+numQuery;++i)
                                    {
                                        // int s = *sp_tmp++;
                                        // int t = *tp_tmp++;
                                        int s = queries[i].first;
                                        int t = queries[i].second;
                                        distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                        //lab.query(s,t);
                                    }
                                    qtime = GetCurrentTimeSec() - qtime;
                                    avg_qtime = qtime/ (double)numQuery;
                                    cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                    cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                    //output query time
                                    if(isOutputAnaly){
                                        if(isOutputQc){
                                            //lab.load_label_size(labelSizeFileName);
                                            lab.load_hfpoint_and_qt(queryDataFileName,0);
                                            lab.append_experiment_result_pll(outputAnaDirName);
                                        }
                                        string analy_filename_prefix(outputAnaDirName);
                                        string analy_filename=analy_filename_prefix.append(".qtime");
                                        ofstream ofs_time(analy_filename);
                                        if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                        ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<std::endl;
                                        ofs_time.close();
                                        
                                    }
                                    if(isOutputDis){
                                        string dis_output_filename(qdisFileName);
                                        dis_output_filename.append(".distances");
                                        ofstream ofs(dis_output_filename);
                                        if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<std::endl;
                                        }
                                        ofs.close();
                                    }
                                }
                                return 0;

                            }else if (dataType==2)//simulation
                            {
                                cout<<"dataType=2-simulation"<<std::endl;
                                //load labels
                                HFLabel lab;
                                lab.index_p=NULL;
                                lab.load_labels(labelFileName);
                                cout<<"load labels successfully!"<<std::endl;
                                //variables
                                double qtime,avg_qtime;
                                vector<pair<NodeID, NodeID> > queries;
                                //generate simulation query data
                                for(NodeID s = 0; s < numOfVertices; ++s){
                                    for(NodeID t=s;t<numOfVertices;++t)
                                    queries.push_back(make_pair(s,t));
                                }
                                numQuery=queries.size();
                                warmup=numQuery/2;
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                if(isShuffle) random_shuffle(queries.begin(), queries.end());
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                cout<<"generate query data successfully!"<<std::endl;
                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    lab.query_p(s,t);
                                }
                                cout<<"warm up successfully!"<<std::endl;
                                qtime=GetCurrentTimeSec();
                                for(int i = 0; i < numQuery; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    distances[i]=lab.query_p(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                //output query time
                                if(isOutputAnaly){
                                    string analy_filename_prefix(outputAnaDirName);
                                    string analy_filename=analy_filename_prefix.append(".qtime");
                                    ofstream ofs_time(analy_filename);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                    ofs_time  <<"qtime="<<qtime * 1e6 << " " <<"avg_qtime="<<avg_qtime * 1e6 <<std::endl;
                                    ofs_time.close();
                                    
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = 0; i < numQuery; ++i){
                                        ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<std::endl;
                                    }
                                    ofs.close();
                                }

                            }
                        }else{
                            cout<<"WEIGHTED_FLAG=false ";
                        }
                    }
                    return 0;
                }else if(t_special_flag==1){ //experiment model
                    cout<<"t_special_flag=1-experiment model ";
                    if (DIRECTED_FLAG == true){//directed 
                        cout<<"DIRECTED_FLAG=true ";
                        if(WEIGHTED_FLAG){//weighted
                            cout<<"WEIGHTED_FLAG=true ";
                            if(dataType==0){//random
                                cout<<"dataType=0-random"<<std::endl;;
                                //variables
                                double qtime,avg_qtime;
                                vector<pair<int, int> > queries(numQuery + warmup);
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                //load labels
                                HFDLabel lab;
                                lab.index_p=NULL;
                                lab.bindex_p=NULL;
                                lab.load_labels(labelFileName);
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
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    lab.query_p(s,t);
                                }
                                qtime=GetCurrentTimeSec();
                                for(int i = warmup; i < warmup + numQuery; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    distances[i-warmup]=lab.query_p(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                //output query time
                                if(isOutputAnaly){
                                    ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                    if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                    ofs.setf(ios::fixed);
                                    ofs.precision(4);
                                    ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" " ;
                                    ofs.close();
                                    
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = warmup; i < warmup + numQuery; ++i){
                                        ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<std::endl;
                                    }
                                    ofs.close();
                                }
                                return 0;
                            }else if (dataType==1)//realData
                            {
                                cout<<"dataType=1-realData ";
                                //load labels
                                HFDLabel lab;
                                lab.index_p=NULL;
                                lab.bindex_p=NULL;
                                lab.load_labels(labelFileName);
                                //load label size
                                if(isOutputQc) lab.load_label_size(labelSizeFileName);
                                if(t_batch_flag==1){
                                    cout<<"t_batch_flag=1-batchQuery "<<endl;
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
                                        if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<i<<" numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            lab.query_p(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i=warmup;i<warmup+numQuery;++i)
                                        {
                                            // int s = *sp_tmp++;
                                            // int t = *tp_tmp++;
                                            int s = queries[i].first;
                                            int t = queries[i].second;
                                            distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                            //lab.query(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        qtime_slices.push_back(qtime);
                                        avg_qtime_slices.push_back(avg_qtime);
                                        if(isOutputQc){
                                            double qc=lab.get_query_cost(query_data_file);
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
                                }else{
                                    cout<<"t_batch_flag=0-oneQuery "<<endl;
                                    //variables
                                    double qtime,avg_qtime;
                                    vector<pair<int, int> > queries_data;
                                    //read real data from file
                                    //numQuery=load_query_pair_time(pointFreqFileName,queryDataFileName,queries,hfRate,numOfHFpoint,multiQuery);
                                    if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                    else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                    warmup=numQuery/2;
                                    vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                    if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                    vector<pair<int, int> > queries(queries_data);
                                    queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                    vector<pair<int, int> > ().swap(queries_data);
                                    queries_data.clear();
                                    cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                    //warm up
                                    for(int i = 0; i < warmup; ++i){
                                        int s = queries[i].first;
                                        int t = queries[i].second;	
                                        lab.query_p(s,t);
                                    }
                                    qtime=GetCurrentTimeSec();
                                    for(int i=warmup;i<warmup+numQuery;++i)
                                    {
                                        // int s = *sp_tmp++;
                                        // int t = *tp_tmp++;
                                        int s = queries[i].first;
                                        int t = queries[i].second;
                                        distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                        //lab.query(s,t);
                                    }
                                    qtime = GetCurrentTimeSec() - qtime;
                                    avg_qtime = qtime/ (double)numQuery;
                                    cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                    cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                    //output query time
                                    if(isOutputAnaly){
                                        if((*queryDataFileName!='\0')&&isOutputAnaly){
                                            lab.load_hfpoint_and_qt(queryDataFileName,0);
                                            lab.append_experiment_result_pll(outputAnaDirName);
                                        }
                                        ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                        if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                        ofs.setf(ios::fixed);
                                        ofs.precision(4);
                                        ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" " ;
                                        ofs.close();
                                        
                                    }
                                    if(isOutputDis){
                                        string dis_output_filename(qdisFileName);
                                        dis_output_filename.append(".distances");
                                        ofstream ofs(dis_output_filename);
                                        if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<std::endl;
                                        }
                                        ofs.close();
                                    }
                                }
                                return 0;

                            }else if (dataType==2)//simulation
                            {
                                cout<<"dataType=2-simulation"<<std::endl;
                                //load labels
                                HFDLabel lab;
                                lab.index_p=NULL;
                                lab.bindex_p=NULL;
                                lab.load_labels(labelFileName);
                                //variables
                                double qtime,avg_qtime;
                                vector<pair<NodeID, NodeID> > queries;
                                //generate simulation query data
                                for(NodeID s = 0; s < numOfVertices; ++s){
                                    for(NodeID t=s;t<numOfVertices;++t)
                                    queries.push_back(make_pair(s,t));
                                }
                                numQuery=queries.size();
                                warmup=numQuery/2;
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                if(isShuffle) random_shuffle(queries.begin(), queries.end());
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    lab.query_p(s,t);
                                }
                                qtime=GetCurrentTimeSec();
                                for(int i = 0; i < numQuery; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    distances[i]=lab.query_p(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                //output query time
                                if(isOutputAnaly){
                                    ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                    if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                    ofs.setf(ios::fixed);
                                    ofs.precision(4);
                                    ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" " ;
                                    ofs.close();
                                    
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = 0; i < numQuery; ++i){
                                        ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<std::endl;
                                    }
                                    ofs.close();
                                }

                            }
                        }else{
                            cout<<"WEIGHTED_FLAG=false";
                        }
                    }else{//undirected
                        cout<<"DIRECTED_FLAG=false ";
                        if(WEIGHTED_FLAG){//weighted
                            cout<<"WEIGHTED_FLAG=true ";
                            if(dataType==0){//random
                                cout<<"dataType=0-random"<<std::endl;
                                //variables
                                double qtime,avg_qtime;
                                vector<pair<int, int> > queries(numQuery + warmup);
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                //load labels
                                HFLabel lab;
                                lab.index_p=NULL;
                                lab.load_labels(labelFileName);
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
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    lab.query_p(s,t);
                                }
                                qtime=GetCurrentTimeSec();
                                for(int i = warmup; i < warmup + numQuery; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    distances[i-warmup]=lab.query_p(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                //output query time
                                if(isOutputAnaly){
                                    ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                    if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                    ofs.setf(ios::fixed);
                                    ofs.precision(4);
                                    ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" " ;
                                    ofs.close();
                                    
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = warmup; i < warmup + numQuery; ++i){
                                        ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<std::endl;
                                    }
                                    ofs.close();
                                }
                                return 0;

                            }else if (dataType==1)//realData
                            {
                                cout<<"dataType=1-realData"; 
                                //load labels
                                HFLabel lab;
                                lab.index_p=NULL;
                                lab.load_labels(labelFileName);
                                if(isOutputQc) lab.load_label_size(labelSizeFileName);
                                if(t_batch_flag==1){
                                    cout<<"t_batch_flag=1-batchQuery "<<endl;
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
                                        if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<i<<" numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            lab.query_p(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i=warmup;i<warmup+numQuery;++i)
                                        {
                                            // int s = *sp_tmp++;
                                            // int t = *tp_tmp++;
                                            int s = queries[i].first;
                                            int t = queries[i].second;
                                            distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                            //lab.query(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        qtime_slices.push_back(qtime);
                                        avg_qtime_slices.push_back(avg_qtime);
                                        if(isOutputQc){
                                            double qc=lab.get_query_cost(query_data_file);
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
                                }else{                                
                                    cout<<"t_batch_flag=0-oneQuery"<<endl;
                                    //variables
                                    double qtime,avg_qtime;
                                    vector<pair<int, int> > queries_data;
                                    int numOfHFpoint;
                                    //read real data from file
                                    //numQuery=load_query_pair_time(pointFreqFileName,queryDataFileName,queries,hfRate,numOfHFpoint,multiQuery);
                                    if(numQuery==0) numQuery=load_query_pair_time(queryDataFileName,queries_data,multiQuery); 
                                    else numQuery=load_query_pair_time(queryDataFileName,queries_data,numQuery,true); 
                                    warmup=numQuery/2;
                                    vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                    if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                    vector<pair<int, int> > queries(queries_data);
                                    queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                    vector<pair<int, int> > ().swap(queries_data);
                                    queries_data.clear();
                                    cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                    //warm up
                                    for(int i = 0; i < warmup; ++i){
                                        int s = queries[i].first;
                                        int t = queries[i].second;	
                                        lab.query_p(s,t);
                                    }
                                    qtime=GetCurrentTimeSec();
                                    for(int i=warmup;i<warmup+numQuery;++i)
                                    {
                                        // int s = *sp_tmp++;
                                        // int t = *tp_tmp++;
                                        int s = queries[i].first;
                                        int t = queries[i].second;
                                        distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                        //lab.query(s,t);
                                    }
                                    qtime = GetCurrentTimeSec() - qtime;
                                    avg_qtime = qtime/ (double)numQuery;
                                    cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                    cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                    //output query time
                                    if(isOutputAnaly){
                                        if(isOutputQc){
                                            //lab.load_label_size(labelSizeFileName);
                                            lab.load_hfpoint_and_qt(queryDataFileName,0);
                                            lab.append_experiment_result_pll(outputAnaDirName);
                                        }
                                        ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                        if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                        ofs.setf(ios::fixed);
                                        ofs.precision(4);
                                        ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" " ;
                                        ofs.close();
                                        
                                    }
                                    if(isOutputDis){
                                        string dis_output_filename(qdisFileName);
                                        dis_output_filename.append(".distances");
                                        ofstream ofs(dis_output_filename);
                                        if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i-warmup]<<std::endl;
                                        }
                                        ofs.close();
                                    }
                                }
                                return 0;

                            }else if (dataType==2)//simulation
                            {
                                cout<<"dataType=2-simulation"<<std::endl;

                                //load labels
                                HFLabel lab;
                                lab.index_p=NULL;
                                lab.load_labels(labelFileName);
                                //variables
                                double qtime,avg_qtime;
                                vector<pair<NodeID, NodeID> > queries;
                                //generate simulation query data
                                for(NodeID s = 0; s < numOfVertices; ++s){
                                    for(NodeID t=s;t<numOfVertices;++t)
                                    queries.push_back(make_pair(s,t));
                                }
                                numQuery=queries.size();
                                warmup=numQuery/2;
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                if(isShuffle) random_shuffle(queries.begin(), queries.end());
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    lab.query_p(s,t);
                                }
                                qtime=GetCurrentTimeSec();
                                for(int i = 0; i < numQuery; ++i){
                                    int s = queries[i].first;
                                    int t = queries[i].second;	
                                    distances[i]=lab.query_p(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                if(isOutputAnaly){
                                    ofstream ofs(outputAnaDirName,ios::app|ios::out);
                                    if(!ofs.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                    ofs.setf(ios::fixed);
                                    ofs.precision(4);
                                    ofs  << qtime * 1e6 <<  " " <<avg_qtime * 1e6 <<" " ;
                                    ofs.close();
                                    
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = 0; i < numQuery; ++i){
                                        ofs<<queries[i].first<<" "<<queries[i].second<<" "<<distances[i]<<std::endl;
                                    }
                                    ofs.close();
                                }
                                return 0;

                        }else{
                            cout<<"WEIGHTED_FLAG=false ";
                        }
                    }
                    return 0;
                }
            }else if(t_special_flag==2){//for spatial experiment model
                cout<<"t_special_flag=2-spatial experiment model";
                if (DIRECTED_FLAG == true){//directed
                    cout<<"DIRECTED_FLAG=true ";
                    if(WEIGHTED_FLAG){//weighted
                        cout<<"WEIGHTED_FLAG=true ";
                        if(dataType==0){//random
                            cout<<"dataType=0-random"<<std::endl;
                        }else if(dataType==1){
                            cout<<"dataType=1-realData ";
                            //load labels
                            HFDLabel lab;
                            lab.index_p=NULL;
                            lab.bindex_p=NULL;
                            lab.load_labels(labelFileName);
                            if(isOutputQc) lab.load_label_size(labelSizeFileName);
                            if(t_batch_flag==1){
                                cout<<"t_batch_flag=1-batchQuery "<<endl;
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
                                    if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                    else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                    warmup=numQuery/2;
                                    vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                    if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                    vector<pair<int, int> > queries(queries_data);
                                    queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                    vector<pair<int, int> > ().swap(queries_data);
                                    queries_data.clear();
                                    cout<<i<<" numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                    //warm up
                                    for(int i = 0; i < warmup; ++i){
                                        int s = queries[i].first;
                                        int t = queries[i].second;	
                                        lab.query_p(s,t);
                                    }
                                    qtime=GetCurrentTimeSec();
                                    for(int i=warmup;i<warmup+numQuery;++i)
                                    {
                                        // int s = *sp_tmp++;
                                        // int t = *tp_tmp++;
                                        int s = queries[i].first;
                                        int t = queries[i].second;
                                        distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                        //lab.query(s,t);
                                    }
                                    qtime = GetCurrentTimeSec() - qtime;
                                    avg_qtime = qtime/ (double)numQuery;
                                    qtime_slices.push_back(qtime);
                                    avg_qtime_slices.push_back(avg_qtime);
                                    if(isOutputQc){
                                        double qc=lab.get_query_cost(query_data_file);
                                        query_cost_slices.push_back(qc);
                                    }
                                }
                                if(isOutputAnaly){
                                    string analy_filename_prefix(outputAnaDirName);
                                    string analy_filename=analy_filename_prefix.append(".qtime");
                                    ofstream ofs_time(analy_filename);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                    //output query time
                                    for(int i=0;i<numOfQueryDataFile;++i){
                                        ofs_time<<avg_qtime_slices[i] * 1e6 <<" ";
                                    }
                                    ofs_time<<endl;
                                    ofs_time.close();
                                }
                            }else{
                                cout<<"t_batch_flag=0-oneQuery "<<endl;
                            }
                            return 0;
                        }else{//simulation
                            cout<<"dataType=2-simulation"<<std::endl;
                        }
                    }else{
                        cout<<"WEIGHTED_FLAG=false ";
                    }
                }else{
                    cout<<"DIRECTED_FLAG=false ";
                    if(WEIGHTED_FLAG){//weighted
                        cout<<"WEIGHTED_FLAG=true ";
                        if(dataType==0){//random
                            cout<<"dataType=0-random"<<std::endl;
                        }else if(dataType==1){
                            cout<<"dataType=1-realData ";
                                //load labels
                                HFLabel lab;
                                lab.index_p=NULL;
                                lab.load_labels(labelFileName);
                                if(isOutputQc) lab.load_label_size(labelSizeFileName);
                                if(t_batch_flag==1){
                                    cout<<"t_batch_flag=1-batchQuery "<<endl;
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
                                        if(numQuery==0) numQuery=load_query_pair_time(query_data_file,queries_data,multiQuery); 
                                        else numQuery=load_query_pair_time(query_data_file,queries_data,numQuery,true); 
                                        warmup=numQuery/2;
                                        vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                        if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                        vector<pair<int, int> > queries(queries_data);
                                        queries.insert(queries.begin(),queries_data.begin(),queries_data.begin()+warmup);
                                        vector<pair<int, int> > ().swap(queries_data);
                                        queries_data.clear();
                                        cout<<i<<" numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                        //warm up
                                        for(int i = 0; i < warmup; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	
                                            lab.query_p(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i=warmup;i<warmup+numQuery;++i)
                                        {
                                            // int s = *sp_tmp++;
                                            // int t = *tp_tmp++;
                                            int s = queries[i].first;
                                            int t = queries[i].second;
                                            distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
                                            //lab.query(s,t);
                                        }
                                        qtime = GetCurrentTimeSec() - qtime;
                                        avg_qtime = qtime/ (double)numQuery;
                                        qtime_slices.push_back(qtime);
                                        avg_qtime_slices.push_back(avg_qtime);
                                        if(isOutputQc){
                                            double qc=lab.get_query_cost(query_data_file);
                                            query_cost_slices.push_back(qc);
                                        }
                                    }
                                    if(isOutputAnaly){
                                        string analy_filename_prefix(outputAnaDirName);
                                        string analy_filename=analy_filename_prefix.append(".qtime");
                                        ofstream ofs_time(analy_filename);
                                        if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                        //output query time
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            ofs_time<<avg_qtime_slices[i] * 1e6 <<" ";
                                        }
                                        ofs_time<<std::endl;
                                        ofs_time.close();
                                    }
                                }else{
                                    cout<<"t_batch_flag=0-oneQuery"<<endl;
                                }
                            return 0;
                        }else{//simulation
                            cout<<"dataType=2-simulation"<<std::endl;
                        }
                    }else{
                        cout<<"WEIGHTED_FLAG=false ";
                    }
                }
            }
            return 0;
        }
    };
}

#endif