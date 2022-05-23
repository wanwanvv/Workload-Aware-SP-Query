/*
 * @Descripttion: Query shortest path distance
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-10-14 08:39:46
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-02-27 16:28:32
 */
#ifndef _COMMAND_QUERY_DISTANCE
#define _COMMAND_QUERY_DISTANCE

#include <cstring>
#include <string>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <unordered_map>

#include "../src/labels.h"
#include "../command.h"
#include "../src/time_util.h"

#define QUERY_NUM 10000000
using namespace time_util; 

namespace command{
    class QueryProcessing: public Command{
        public:
            void exit_with_help(){
                printf("Usage:\n");
                printf("bin/pll_whp_run -q -d [directedGraphFlag] -w [weightedGraphFlag] -s [experimentFlag] -p [dataType] -l [inputLabelFileName] \n \\
                [-t [queryPairFileName]] [-i [queryFreqFileName]] [-a [outputResultDirName]] [-f [outputDistanceFileName]] [-z [inputLabelSizeFileName]] \n \\
                 [-u [multiQuery]] [-n [numOfQuery]] [-x [shuffleFlag]] \n");
                printf("-------------------------------------------------------------------\n");
                printf("Parameter explanation:\n");
                printf("\t-d: 0 or 1, for undirected and directed graphs, default is 0 \n");
                printf("\t-w: 0 or 1, for unweighted and weighted graphs, default is 0 \n");
                printf("\t-s: 0-default model, 1-experiment model \n");
                printf("\t-u: mutiple value of input query pairs, default is 1 \n");
                printf("\t-p: 0-generate query data random,1-user-defined input query data,2-gen3rate all s-t pair data for verification\n");
                printf("\t-n: num of s-t query pairs, default is 0\n");
                printf("\t-x:  0 or 1, shuffle the order of input query pairs, default is 0\n");
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
                 int outputModel=0;//0-ouput atotal result, 1-output each time slice result, default is 0, default is 0
                
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
                        case 'o':
                            outputModel = atoi(argv[i]);
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
                                }else if(t_batch_flag==0){
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
                                        ofs<<numOfQueryDataFile<<" ";
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
                                        if((*queryDataFileName!='\0')&&isOutputQc){
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
                                        ofs<<numOfQueryDataFile<<" ";
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
                                    ofstream ofs_time(outputAnaDirName,ios::app|ios::out);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                    //output query time
                                    for(int i=0;i<numOfQueryDataFile;++i){
                                        ofs_time<<avg_qtime_slices[i] * 1e6 <<" ";
                                    }
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
                                        ofstream ofs_time(outputAnaDirName,ios::app|ios::out);
                                        if(!ofs_time.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                        //output query time
                                        for(int i=0;i<numOfQueryDataFile;++i){
                                            ofs_time<<avg_qtime_slices[i] * 1e6 <<" ";
                                        }
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
            }else if(t_special_flag==3){//spatial temporal model(multiple labels and query dirs)
                cout<<"t_special_flag=3-spatial temporal experiment model ";
                if (DIRECTED_FLAG == true){//directed
                    cout<<"DIRECTED_FLAG=true ";
                    if(WEIGHTED_FLAG){//weighted
                        cout<<"WEIGHTED_FLAG=true ";
                        if(dataType==0){//random
                            cout<<"dataType=0-random"<<std::endl;
                        }else if(dataType==1){
                            cout<<"dataType=1-realData ";
                            if(t_batch_flag==1){
                                cout<<"t_batch_flag=1-batchQuery "<<endl;
                                //init prefix
                                string queryDataFileName_prefix=string(queryDataFileName);
                                if(queryDataFileName_prefix[queryDataFileName_prefix.size()-1]!='/') queryDataFileName_prefix.append("/");
                                string labelFileName_prefix=string(labelFileName);
                                if(labelFileName_prefix[labelFileName_prefix.size()-1]!='/') labelFileName_prefix.append("/");
                                string labelSizeFileName_prefix=string(labelSizeFileName);
                                if(labelSizeFileName_prefix[labelSizeFileName_prefix.size()-1]!='/') labelSizeFileName_prefix.append("/");
                                //get num of time intervals
                                vector<string> label_files;
                                get_filelist_from_dir(labelFileName_prefix,label_files);
                                int time_interval_cnt=label_files.size();
                                std::cout<<"time_interval_cnt="<<time_interval_cnt<<std::endl;
                                //initialize variables
                                double qtime,avg_qtime;
                                vector<double> qtime_slices;
                                vector<double> avg_qtime_slices;
                                vector<double> query_cost_slices;
                                vector<string> all_query_data_files;
                                for(int i=0;i<time_interval_cnt;++i){
                                    std::cout<<"***********************Time interval index="<<i<<" query start!************************"<<std::endl;
                                    //each time interval variable
                                    string queryDataFileName_t=queryDataFileName_prefix+"TIP-"+to_string(i);
                                    string labelFileName_t=labelFileName_prefix+to_string(i)+".label";
                                    string labelSizeFileName_t=labelSizeFileName_prefix+to_string(i)+".size";
                                    //load label
                                    HFDLabel lab;
                                    lab.index_p=NULL;
                                    lab.bindex_p=NULL;
                                    lab.load_labels((char*)labelFileName_t.c_str());
                                    if(isOutputQc) lab.load_label_size((char*)labelSizeFileName_t.c_str());
                                    //load query files and get its size
                                    vector<pair<int, int> > queries_data;
                                    string queryDataPath(queryDataFileName_t);
                                    vector<string> query_data_files;
                                    vector<string> full_query_data_files;
                                    get_filelist_from_dir(queryDataPath,query_data_files,true);
                                    append_to_full_path(queryDataPath,query_data_files,full_query_data_files);
                                    int numOfQueryDataFile=query_data_files.size();
                                    for(int i=0;i<numOfQueryDataFile;++i){
                                        char* query_data_file=(char*)full_query_data_files[i].c_str();
                                        all_query_data_files.push_back(query_data_files[i]);
                                        std::cout<<query_data_file<<std::endl;
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
                                            lab.query_p(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
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
                                    std::cout<<"***********************Time interval index="<<i<<" query successfully!************************"<<std::endl;
                                }
                                int numOfTotalQueryDataFile=qtime_slices.size();
                                std::cout<<"numOfTotalQueryDataFile="<<numOfTotalQueryDataFile<<std::endl;
                                if(isOutputAnaly){
                                    ofstream ofs_time(outputAnaDirName,ios::app|ios::out);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                    if(outputModel==0){
                                        double query_cost_all_slices=0;
                                        double qtime_all_slices=0;
                                        double avg_qtime_all_slices=0;
                                        for(int i=0;i<numOfTotalQueryDataFile;++i){
                                            if(isOutputQc){
                                                query_cost_all_slices+=query_cost_slices[i];
                                            }
                                            qtime_all_slices+=qtime_slices[i];
                                            avg_qtime_all_slices+=avg_qtime_slices[i];
                                        }
                                        double query_cost_avg_slice=query_cost_all_slices/(double)numOfTotalQueryDataFile;
                                        double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfTotalQueryDataFile;
                                        //output query time
                                        if(isOutputQc) ofs_time<<query_cost_all_slices<<" "<<query_cost_avg_slice<<" ";
                                        ofs_time  << qtime_all_slices * 1e6 <<  " " <<avg_qtime_avg_slice * 1e6 <<" " ;
                                        
                                    }else if(outputModel==1){
                                        for(int i=0;i<numOfTotalQueryDataFile;++i){
                                            int day_index=get_day_index(all_query_data_files[i]);
                                            int cut_index=get_cut_index(all_query_data_files[i]);
                                            ofs_time<<day_index<<" "<<cut_index<<" ";
                                            if(isOutputQc) ofs_time<<query_cost_slices[i]<<" ";
                                            ofs_time  << avg_qtime_slices[i] * 1e6 << std::endl;
                                        }
                                    }
                                    ofs_time.close();
                                }
                                return 0;
                            }else{
                                cout<<"t_batch_flag=0-oneQuery "<<endl;
                            }
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
                            if(t_batch_flag==1){
                                cout<<"t_batch_flag=1-batchQuery "<<endl;
                                //init prefix
                                string queryDataFileName_prefix=string(queryDataFileName);
                                if(queryDataFileName_prefix[queryDataFileName_prefix.size()-1]!='/') queryDataFileName_prefix.append("/");
                                string labelFileName_prefix=string(labelFileName);
                                if(labelFileName_prefix[labelFileName_prefix.size()-1]!='/') labelFileName_prefix.append("/");
                                string labelSizeFileName_prefix=string(labelSizeFileName);
                                if(labelSizeFileName_prefix[labelSizeFileName_prefix.size()-1]!='/') labelSizeFileName_prefix.append("/");
                                //get num of time intervals
                                vector<string> label_files;
                                get_filelist_from_dir(labelFileName_prefix,label_files);
                                int time_interval_cnt=label_files.size();
                                std::cout<<"time_interval_cnt="<<time_interval_cnt<<std::endl;
                                //initialize variables
                                double qtime,avg_qtime;
                                vector<double> qtime_slices;
                                vector<double> avg_qtime_slices;
                                vector<double> query_cost_slices;
                                vector<string> all_query_data_files;
                                for(int i=0;i<time_interval_cnt;++i){
                                    std::cout<<"***********************Time interval index="<<i<<" query start!************************"<<std::endl;
                                    //each time interval variable
                                    string queryDataFileName_t=queryDataFileName_prefix+"TIP-"+to_string(i);
                                    string labelFileName_t=labelFileName_prefix+to_string(i)+".label";
                                    string labelSizeFileName_t=labelSizeFileName_prefix+to_string(i)+".size";
                                    //load label
                                    HFLabel lab;
                                    lab.index_p=NULL;
                                    lab.load_labels((char*)labelFileName_t.c_str());
                                    if(isOutputQc) lab.load_label_size((char*)labelSizeFileName_t.c_str());
                                    //load query files and get its size
                                    vector<pair<int, int> > queries_data;
                                    string queryDataPath(queryDataFileName_t);
                                    vector<string> query_data_files;
                                    vector<string> full_query_data_files;
                                    get_filelist_from_dir(queryDataPath,query_data_files,true);
                                    append_to_full_path(queryDataPath,query_data_files,full_query_data_files);
                                    int numOfQueryDataFile=query_data_files.size();
                                    for(int i=0;i<numOfQueryDataFile;++i){
                                        char* query_data_file=(char*)full_query_data_files[i].c_str();
                                        all_query_data_files.push_back(query_data_files[i]);
                                        std::cout<<query_data_file<<std::endl;
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
                                            lab.query_p(s,t);
                                        }
                                        qtime=GetCurrentTimeSec();
                                        for(int i = warmup; i < warmup+numQuery; ++i){
                                            int s = queries[i].first;
                                            int t = queries[i].second;	 
                                            distances[i-warmup]=lab.query_p(s,t);//modified by wanjingyi
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
                                    std::cout<<"***********************Time interval index="<<i<<" query successfully!************************"<<std::endl;
                                }
                                int numOfTotalQueryDataFile=qtime_slices.size();
                                std::cout<<"numOfTotalQueryDataFile="<<numOfTotalQueryDataFile<<std::endl;
                                if(isOutputAnaly){
                                    ofstream ofs_time(outputAnaDirName,ios::app|ios::out);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<outputAnaDirName<<endl;
                                    if(outputModel==0){
                                        double query_cost_all_slices=0;
                                        double qtime_all_slices=0;
                                        double avg_qtime_all_slices=0;
                                        for(int i=0;i<numOfTotalQueryDataFile;++i){
                                            if(isOutputQc){
                                                query_cost_all_slices+=query_cost_slices[i];
                                            }
                                            qtime_all_slices+=qtime_slices[i];
                                            avg_qtime_all_slices+=avg_qtime_slices[i];
                                        }
                                        double query_cost_avg_slice=query_cost_all_slices/(double)numOfTotalQueryDataFile;
                                        double avg_qtime_avg_slice=avg_qtime_all_slices/(double)numOfTotalQueryDataFile;
                                        //output query time
                                        if(isOutputQc) ofs_time<<query_cost_all_slices<<" "<<query_cost_avg_slice<<" ";
                                        ofs_time  << qtime_all_slices * 1e6 <<  " " <<avg_qtime_avg_slice * 1e6 <<" " ;
                                        
                                    }else if(outputModel==1){
                                        for(int i=0;i<numOfTotalQueryDataFile;++i){
                                            int day_index=get_day_index(all_query_data_files[i]);
                                            int cut_index=get_cut_index(all_query_data_files[i]);
                                            ofs_time<<day_index<<" "<<cut_index<<" ";
                                            if(isOutputQc) ofs_time<<query_cost_slices[i]<<" ";
                                            ofs_time  << avg_qtime_slices[i] * 1e6 << std::endl;
                                        }
                                    }
                                    ofs_time.close();
                                }
                                return 0;
                            }else{
                                cout<<"t_batch_flag=0-oneQuery "<<endl;
                            }
                        }else{//simulation
                            cout<<"dataType=2-simulation"<<std::endl;
                        }
                    }else{
                        cout<<"WEIGHTED_FLAG=false ";
                    }
                }
            }else if(t_special_flag==4){//spatial temporal model for k label sets(multiple labels and querydata file)
                cout<<"t_special_flag=4-spatial temporal k label sets experiment model ";
                if (DIRECTED_FLAG == true){//directed
                    cout<<"DIRECTED_FLAG=true ";
                    if(WEIGHTED_FLAG){
                        cout<<"WEIGHTED_FLAG=true ";
                        if(dataType==0){//random
                            cout<<"dataType=0-random"<<std::endl;
                        }else if(dataType==1){
                            if(t_batch_flag==0){
                                cout<<"t_batch_flag=0-Query for k label sets"<<endl;
                                //load label size for k label sets
                                string labelSizeFilePath(labelSizeFileName);
                                vector<string> label_size_files;
                                get_filelist_from_dir_string(labelSizeFilePath,label_size_files,true);
                                append_to_full_path(labelSizeFilePath,label_size_files);
                                int numOfLabels=label_size_files.size();
                                vector<HFDLabel> labs(numOfLabels);//store k labels
                                cout<<"numOfLabels="<<numOfLabels<<endl;//debug for deleted
                                for(int i=0;i<numOfLabels;++i){
                                    cout<<i<<":"<<label_size_files[i]<<endl;//debug for deleted
                                    labs[i].load_label_size((char*)label_size_files[i].c_str());
                                }
                                //generate query data,compute  min(|L(s)|+|L(t)|) and compute query cost
                                double query_cost=0;
                                vector<vector<int>> queries_data; //{s,t,min_index}
                                vector<vector<bool>> flags(numOfLabels,vector<bool>(numOfVertices,0));//indicate i labels load which vertices
                                ifstream in(queryDataFileName);
                                NodeID s,t,freq;
                                while (in>>s>>t>>freq)
                                {
                                    vector<NodeID> k_size;
                                    for(int i=0;i<numOfLabels;++i) k_size.push_back(labs[i].spt_v_num_f[s]+labs[i].spt_v_num_r[t]);
                                    int min_index=min_element(k_size.begin(),k_size.end())-k_size.begin();
                                    flags[min_index][s]=true;
                                    flags[min_index][t]=true;
                                    query_cost+=(double)k_size[min_index]*(freq/(double)DIVISION_FACTOR);
                                    for(int i=0;i<freq*multiQuery;++i) queries_data.push_back({s,t,min_index});
                                }
                                int numQuery_real=queries_data.size();
                                if(numQuery!=0){
                                    while(numQuery_real<numQuery-(int)queries_data.size()){
                                        for(int i=0;i<numQuery_real; ++i)
                                        queries_data.push_back(queries_data[i]);
                                    }
                                }
                                numQuery=queries_data.size();

                                //load k labels
                                string labelFilePath(labelFileName);
                                vector<string> label_files;
                                get_filelist_from_dir_string(labelFilePath,label_files,true);
                                append_to_full_path(labelFilePath,label_files);
                                for(int i=0;i<numOfLabels;++i){
                                    cout<<i<<":"<<label_files[i]<<endl;//debug for deleted
                                    labs[i].load_labels((char*)label_files[i].c_str(),flags[i]);
                                }
                                //generate query data
                                warmup=numQuery/2;
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                cout<<"generate query data successfully!"<<std::endl;
                                double qtime,avg_qtime;

                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries_data[i][0];
                                    int t = queries_data[i][1];	
                                    labs[queries_data[i][2]].query_p(s,t);
                                }
                                cout<<"warm up successfully!"<<std::endl;
                                qtime=GetCurrentTimeSec();

                                //query stage
                                for(int i=0;i<numQuery;++i)
                                {
                                    // int s = *sp_tmp++;
                                    // int t = *tp_tmp++;
                                    int s = queries_data[i][0];
                                    int t = queries_data[i][1];
                                    distances[i]=labs[queries_data[i][2]].query_p(s,t);//modified by wanjingyi
                                    //lab.query(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                
                                //output query time
                                if(isOutputAnaly){
                                    string analy_filename(outputAnaDirName);
                                    ofstream ofs_time(analy_filename,ios::app|ios::out);
                                    ofs_time.setf(ios::fixed);
                                    ofs_time.precision(4);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                    ofs_time  <<qtime * 1e6 << " "<<avg_qtime * 1e6 <<" "<<query_cost<<" ";
                                    ofs_time.close();
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = 0; i < numQuery; ++i){
                                        ofs<<queries_data[i][0]<<" "<<queries_data[i][1]<<" "<<distances[i]<<std::endl;
                                    }
                                    ofs.close();
                                }
                            }else{
                                cout<<"t_batch_flag=1-Query for k label sets"<<endl;
                            }
                        }else{//simulation
                            cout<<"dataType=2-simulation"<<std::endl;
                        }
                    }else{
                        cout<<"WEIGHTED_FLAG=false ";
                    }
                }else{
                    cout<<"DIRECTED_FLAG=false ";
                    if(WEIGHTED_FLAG){
                        cout<<"WEIGHTED_FLAG=true ";
                        if(dataType==0){//random
                            cout<<"dataType=0-random"<<std::endl;
                        }else if(dataType==1){
                            if(t_batch_flag==0){
                                cout<<"t_batch_flag=0-Query for k label sets"<<endl;
                                //load label size for k label sets
                                string labelSizeFilePath(labelSizeFileName);
                                vector<string> label_size_files;
                                get_filelist_from_dir_string(labelSizeFilePath,label_size_files,true);
                                append_to_full_path(labelSizeFilePath,label_size_files);
                                int numOfLabels=label_size_files.size();
                                vector<HFLabel> labs(numOfLabels);//store k labels
                                cout<<"numOfLabels="<<numOfLabels<<endl;//debug for deleted
                                for(int i=0;i<numOfLabels;++i){
                                    cout<<i<<":"<<label_size_files[i]<<endl;//debug for deleted
                                    labs[i].load_label_size((char*)label_size_files[i].c_str());
                                }
                                //generate query data,compute  min(|L(s)|+|L(t)|) and compute query cost
                                double query_cost=0;
                                vector<vector<int>> queries_data; //{s,t,min_index}
                                vector<vector<bool>> flags(numOfLabels,vector<bool>(numOfVertices,0));//indicate i labels load which vertices
                                ifstream in(queryDataFileName);
                                NodeID s,t,freq;
                                while (in>>s>>t>>freq)
                                {
                                    vector<NodeID> k_size;
                                    for(int i=0;i<numOfLabels;++i) k_size.push_back(labs[i].spt_v_num[s]+labs[i].spt_v_num[t]);
                                    int min_index=min_element(k_size.begin(),k_size.end())-k_size.begin();
                                    flags[min_index][s]=true;
                                    flags[min_index][t]=true;
                                    query_cost+=(double)k_size[min_index]*(freq/(double)DIVISION_FACTOR);
                                    for(int i=0;i<freq*multiQuery;++i) queries_data.push_back({s,t,min_index});
                                }
                                int numQuery_real=queries_data.size();
                                cout<<"numQuery="<<numQuery<<" queries_data.size()="<<numQuery_real<<endl;
                                if(numQuery!=0){
                                    while(numQuery_real<numQuery-(int)queries_data.size()){
                                        for(int i=0;i<numQuery_real; ++i)
                                        queries_data.push_back(queries_data[i]);
                                    }
                                }
                                numQuery=queries_data.size();
                                
                                //load k labels
                                string labelFilePath(labelFileName);
                                vector<string> label_files;
                                get_filelist_from_dir_string(labelFilePath,label_files,true);
                                append_to_full_path(labelFilePath,label_files);
                                for(int i=0;i<numOfLabels;++i){
                                    cout<<i<<":"<<label_files[i]<<endl;//debug for deleted
                                    labs[i].load_labels((char*)label_files[i].c_str(),flags[i]);
                                }
                                //generate query data
                                warmup=numQuery/2;
                                vector<EdgeWeight> distances(numQuery,INF_WEIGHT);
                                if(isShuffle) std::random_shuffle(queries_data.begin(), queries_data.end());  
                                cout<<"numQuery="<<numQuery<<" warmup="<<warmup<<std::endl;
                                cout<<"generate query data successfully!"<<std::endl;
                                double qtime,avg_qtime;

                                //warm up
                                for(int i = 0; i < warmup; ++i){
                                    int s = queries_data[i][0];
                                    int t = queries_data[i][1];	
                                    labs[queries_data[i][2]].query_p(s,t);
                                }
                                cout<<"warm up successfully!"<<std::endl;
                                qtime=GetCurrentTimeSec();

                                //query stage
                                for(int i=0;i<numQuery;++i)
                                {
                                    // int s = *sp_tmp++;
                                    // int t = *tp_tmp++;
                                    int s = queries_data[i][0];
                                    int t = queries_data[i][1];
                                    distances[i]=labs[queries_data[i][2]].query_p(s,t);//modified by wanjingyi
                                    //lab.query(s,t);
                                }
                                qtime = GetCurrentTimeSec() - qtime;
                                avg_qtime = qtime/ (double)numQuery;
                                cout << "Total query time:" << qtime * 1e6 <<  " microseconds" << std::endl;
                                cout << "Average query time:" << avg_qtime * 1e6 <<  " microseconds" << std::endl;
                                
                                //output query time
                                if(isOutputAnaly){
                                    string analy_filename(outputAnaDirName);
                                    ofstream ofs_time(analy_filename,ios::app|ios::out);
                                    ofs_time.setf(ios::fixed);
                                    ofs_time.precision(4);
                                    if(!ofs_time.is_open()) cout<<"Cannot open "<<analy_filename<<endl;
                                    ofs_time  <<qtime * 1e6 << " " <<avg_qtime * 1e6 <<" "<<query_cost<<" ";
                                    ofs_time.close();
                                }
                                if(isOutputDis){
                                    string dis_output_filename(qdisFileName);
                                    dis_output_filename.append(".distances");
                                    ofstream ofs(dis_output_filename);
                                    if(!ofs.is_open()) cout<<"Cannot open"<<dis_output_filename<<std::endl;
                                    for(int i = 0; i < numQuery; ++i){
                                        ofs<<queries_data[i][0]<<" "<<queries_data[i][1]<<" "<<distances[i]<<std::endl;
                                    }
                                    ofs.close();
                                }
                            }else{
                                cout<<"t_batch_flag=1-batchQuery for k label sets"<<endl;
                            }
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