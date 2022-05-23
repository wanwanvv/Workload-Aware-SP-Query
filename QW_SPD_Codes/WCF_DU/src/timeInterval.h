/*
 * @Descripttion: 
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-10-29 14:35:51
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-11-11 11:53:08
 */
#pragma once
#ifndef _TIME_INTERVAL_H
#define _TIME_INTERVAL_H

#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include "./integrated_construction.h"
#include "./integrated_construction.h"
#include "./paras.h"
#include "./time_util.h"
#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG

using namespace std;

/**
 * @description: the class is used to construct index after time interval partitioning
 * @param {*}
 * @return {*}
 * @author: Wan Jingyi
 */
class TimeInterval{
    public:
        //******************class variables*******************
        int _time_slice_cnt=0; //num of total time slices
        int _time_interval_cnt=0; //num of total time intervals after partitioning
        int _day_slice_cnt;//num of time slices in a day
        vector<string> _file_names; //store the all query file names
        vector<int> _partition_indexes;//store the update time slices offset index
        string _slice_result_file;//output result of each slice
        string _total_result_file;//output reuslt for total
        string _point_freq_dir_prefix;//"/" at end
        string _new_point_freq_dir_prefix;//"/" at end
        string _query_data_dir_prefix;//"/" at end
        string _new_query_data_dir_prefix;//"/" at end
        vector<Processing::calcCoefficient<double> > _paras;//store the paras for wcf-indexing
        vector<int> _maxDegree_paras;
        double _up_threshold=0;

        vector<int> _idQueryTime;//store the point-freq
        vector<double> _node_size;//store the node size
        vector<double> _query_cost;//store all vertices' query cost by time slice index

        double _total_labeling_time=0;
        double _total_ordering_time=0;
        double _total_compute_rmq_time=0;
        //*****************construction functions**************
        //used for RL-TIP algorithm
        TimeInterval(char* graphFileName,char* pointFreqDirName,char* updateIndexFileName,
                                    char* parasFileName,char* analyDirName,char* labelsDirName,
                                    char* labelSizeDirName,char* outputDirName,char* givenOrderFileName,
                                    int indexing_model,int ordering_model,int threads_num,
                                    int numOfDaySlices,int graphType,char* queryDataDirName){
            load_tip_information(numOfDaySlices,pointFreqDirName,updateIndexFileName);
            make_and_copy_partition_tmp_dirs();
            construct_labels_rl_tip(graphFileName,parasFileName,analyDirName,labelsDirName,labelSizeDirName,outputDirName,givenOrderFileName,indexing_model,ordering_model,threads_num,graphType);
            clear_partition_tmp_dirs();
            if(*queryDataDirName!='\0') partition_query_wokloads(queryDataDirName);
        }

        //used for GTP Algorithm
        TimeInterval(char* graphFileName,char* pointFreqDirName,char* parasFileName,
                                    char* analyDirName,char* labelsDirName,char* labelSizeDirName,
                                    char* outputDirName,char* givenOrderFileName,int indexing_model,
                                    int ordering_model,int threads_num,int graphType,
                                    double up_threshold,char* queryDataDirName)
        {
            load_gtp_information(pointFreqDirName,up_threshold);
            construct_labels_gtp(graphFileName,parasFileName,analyDirName,labelsDirName,labelSizeDirName,outputDirName,givenOrderFileName,indexing_model,ordering_model,threads_num,graphType,up_threshold);
            if(*queryDataDirName!='\0') partition_query_wokloads(queryDataDirName);
        }

        TimeInterval(){}

    protected:

        //*************************util functions for GTP*******************
        void partition_query_wokloads(char* queryDataDirName){
            string freq_str(queryDataDirName);
            _query_data_dir_prefix=freq_str;
            if(_query_data_dir_prefix[_query_data_dir_prefix.size()-1]!='/') _query_data_dir_prefix.append("/");
            _new_query_data_dir_prefix=_query_data_dir_prefix+"TIP-";
            string command_mkdir,command_copy,command_chmod;//commands
            int i=0,j=0;
            for(i=0;i<_time_interval_cnt;++i){
                //make point freq dir
                string new_query_data_dir=_new_query_data_dir_prefix+to_string(i)+"/";
                command_mkdir="mkdir -p "+new_query_data_dir;
                system(command_mkdir.c_str());
                for(;j<_partition_indexes[i];++j){
                    //copy point freq file
                    string curr_query_data_file=_query_data_dir_prefix+_file_names[j];
                    string new_query_data_file=new_query_data_dir+_file_names[j];
                    command_copy="cp "+curr_query_data_file+" "+ new_query_data_file;
                    system(command_copy.c_str());
                }
                //change the privilege
                command_chmod="chmod -R 777 "+new_query_data_dir;
                system(command_chmod.c_str());
            }
            std::cout<<"Partition_query_wokloads successfully!"<<std::endl;
        }
        /*
         *@description: judege whether the timeslicing need to be updated
         *@author: wanjingyi
         *@date: 2021-01-05
        */
        bool isSameCluster(double cost_t1,double cost_t2){
            double dis=(cost_t2-cost_t1)/cost_t1;
            //std::cout<<"dis="<<dis<<endl;//to be deleted
            return dis<=_up_threshold;
        }

        double get_time_slice_query_cost(char* load_filename){
            //initialize iterms
            vector<int> ().swap(_idQueryTime);
            _idQueryTime.resize(numOfVertices,0);
            double performance_result=0;//total performance function 
            ifstream in(load_filename);//input HFPoint file to ifstream
            if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
            size_t i;
            unsigned int query_time=0;
            NodeID q_id;
            for(i=0;in>>q_id>>query_time;++i){
                _idQueryTime[q_id]=query_time;
            }
            in.close();
			for (NodeID v = 0; v < numOfVertices; ++v) 
			{
				double isize = _node_size[v] ;
				//compute the query performance function
				double ratio=(double)_idQueryTime[v]/(double)DIVISION_FACTOR;
				performance_result+=ratio*(double)isize;
			}
            return performance_result;
        }
    
        //*************************util functions for RL-TIP*******************
        void load_tip_information(int numOfDaySlices,char* pointFreqDirName,char* updateIndexFileName){
            //init input data dir variables
            string freq_str(pointFreqDirName);
            _point_freq_dir_prefix=freq_str;
            if(_point_freq_dir_prefix[_point_freq_dir_prefix.size()-1]!='/') _point_freq_dir_prefix.append("/");
            _new_point_freq_dir_prefix=_point_freq_dir_prefix+"TIP-";
            //initialize num of time slices in a day
            _day_slice_cnt=numOfDaySlices;
            //get total query data files
            string point_freq_path(pointFreqDirName);
            get_filelist_from_dir(point_freq_path,_file_names,true);
            _time_slice_cnt=_file_names.size();
            std::cout<<"_day_slice_cnt="<<_day_slice_cnt<<" _time_slice_cnt="<<_time_slice_cnt<<std::endl;
            //init day and cut range variables
            load_partition_indexes(updateIndexFileName);
            //_update_indexs.push_back(TIME_SLICE_CNT);
            cout<<"_time_interval_cnt="<<_time_interval_cnt<<endl;
            cout<<"Load_tip_information successfully!"<<endl;
            return;
        }

        void load_gtp_information(char* pointFreqDirName,double up_threshold){
            //init input data dir variables
            _up_threshold=up_threshold;
            string freq_str(pointFreqDirName);
            _point_freq_dir_prefix=freq_str;
            if(_point_freq_dir_prefix[_point_freq_dir_prefix.size()-1]!='/') _point_freq_dir_prefix.append("/");
            //get total query data files
            string point_freq_path(pointFreqDirName);
            get_filelist_from_dir(point_freq_path,_file_names,true);
            _time_slice_cnt=_file_names.size();
            std::cout<<" _time_slice_cnt="<<_time_slice_cnt<<std::endl;
            cout<<"Load_gtp_information successfully!"<<endl;
            return;
        }

        void load_partition_indexes(char* updateIndexFileName){
            //get the first file and last file
            int day_start_index=get_day_index(_file_names[0]);
            int cut_start_index=get_cut_index(_file_names[0]);
            int day_end_index=get_day_index(_file_names[_time_slice_cnt-1]);
            int cut_end_index=get_cut_index(_file_names[_time_slice_cnt-1]);
            //debug
            //std::cout<<"day_start_index="<<day_start_index<<" cut_start_index="<<cut_start_index<<std::endl;
            //std::cout<<"day_end_index="<<day_end_index<<" cut_end_index="<<cut_end_index<<std::endl;
            ifstream in(updateIndexFileName);
            if(!in.is_open()) std::cout<<"Cannot open "<<updateIndexFileName<<"!"<<std::endl;
            int day_index,cut_index;
            vector<pair<int,int> > partition_indexes;
            for (int i = 0; in >> day_index>>cut_index;i++) {
                partition_indexes.push_back(make_pair(day_index,cut_index));
            }
            in.close();
            partition_indexes.push_back(make_pair(day_end_index,cut_end_index+1));
            _time_interval_cnt=partition_indexes.size();
            //debug
            //cout<<"partition_indexes:";
            //for(size_t i=0;i<partition_indexes.size();++i) std::cout<<"("<<partition_indexes[i].first<<","<<partition_indexes[i].second<<") ";
            //cout<<endl;
            std::cout<<"_time_interval_cnt="<<_time_interval_cnt<<std::endl;
            //partition the time slices to time interval
            int last_day_index=day_start_index,curr_day_index;
            int last_cut_index=cut_start_index,curr_cut_index;
            int last_partition_index=0,partition_len;
            for(int i=0;i<_time_interval_cnt;++i){
                curr_day_index=partition_indexes[i].first;
                curr_cut_index=partition_indexes[i].second;
                partition_len=(curr_day_index-last_day_index)*_day_slice_cnt+curr_cut_index-last_cut_index;
                last_partition_index+=partition_len;
                _partition_indexes.push_back(last_partition_index);
                //debug
                //std::cout<<i<<":partition_len="<<partition_len<<" last_day_index="<<last_day_index<<" curr_day_index="<<curr_day_index<<" curr_cut_index="<<curr_cut_index<<" last_cut_index="<<last_cut_index<<endl;
                last_day_index=curr_day_index;
                last_cut_index=curr_cut_index;
            }
            //debug
            //std::cout<<"_partition_indexes:";
            //for(int i=0;i<_time_interval_cnt;++i){
            //    std::cout<<_partition_indexes[i]<<" ";
            //}
            //std::cout<<std::endl;
            //check
            if(_partition_indexes.size()!=_time_interval_cnt||_partition_indexes[_time_interval_cnt-1]!=_time_slice_cnt) std::cout<<"Error:  _partition_indexes compute failed!"<<std::endl;
            return;
        }

        void make_and_copy_partition_tmp_dirs(){
            string command_mkdir,command_copy,command_chmod;//commands
            int i=0,j=0;
            for(i=0;i<_time_interval_cnt;++i){
                //make point freq dir
                string new_point_freq_dir=_new_point_freq_dir_prefix+to_string(i)+"/";
                command_mkdir="mkdir -p "+new_point_freq_dir;
                system(command_mkdir.c_str());
                for(;j<_partition_indexes[i];++j){
                    //copy point freq file
                    string curr_point_freq_file=_point_freq_dir_prefix+_file_names[j];
                    string new_point_freq_file=new_point_freq_dir+_file_names[j];
                    command_copy="cp "+curr_point_freq_file+" "+ new_point_freq_file;
                    system(command_copy.c_str());
                }
                //change the privilege
                command_chmod="chmod -R 777 "+new_point_freq_dir;
                system(command_chmod.c_str());
            }
            std::cout<<"Make_and_copy_partition_tmp_dirs successfully!"<<std::endl;
            return;
        }

        void clear_partition_tmp_dirs(){
            string command_clear;
            for(int i=0;i<_time_interval_cnt;++i){
                string new_point_freq_dir=_new_point_freq_dir_prefix+to_string(i)+"/";
                command_clear="rm -rf "+new_point_freq_dir;
                system(command_clear.c_str());
            }
            std::cout<<"Clear_partition_tmp_dirs successfully!"<<std::endl;
            return;
        }

        /**
         * @description: load wcf-index parameters from file
         * @param {*}
         * @return {*}
         * @author: Wan Jingyi
         */        
        void load_parameters(char* parasFileName,int num){
            int maxDegree=0,k_freq=0,k_bet=0,k_deg=0,k_cov=0,k_dep=0;
            bool is_deg=false, is_freq=false, is_cov=false, is_dep=false,is_bet=false;
            if(*parasFileName!='\0'){
                ifstream in(parasFileName);//input HFPoint file to ifstream
                if(!in.is_open()) {cerr<<"Cannot open "<<parasFileName<<endl;}
                size_t i;
                for(i=0;in>>maxDegree>>k_freq>>k_deg;++i){
                    k_bet=100-k_freq;
                    is_freq=is_bet=true;
                    if(k_deg!=0) is_deg=true;
                    Processing::calcCoefficient<double> calcCoef((double)k_deg/100,(double)k_freq/100,(double)k_cov/100,(double)k_dep/100,(double)k_bet/100,is_deg, is_freq, is_cov, is_dep,is_bet);
                    _paras.push_back(calcCoef);
                    _maxDegree_paras.push_back(maxDegree);
                }
                in.close();
            }else{
                maxDegree=10;
                k_freq=10;
                k_bet=90;
                k_deg=100;
                is_deg=is_freq=is_bet=true;
                Processing::calcCoefficient<double> calcCoef((double)k_deg/100,(double)k_freq/100,(double)k_cov/100,(double)k_dep/100,(double)k_bet/100,is_deg, is_freq, is_cov, is_dep,is_bet);
                for(int i=0;i<num;++i){
                    _paras.push_back(calcCoef);
                    _maxDegree_paras.push_back(maxDegree);
                }
            }
            if(_paras.size()!=num) std::cout<<"Error:_paras.size()!=num ! _paras.size()="<<_paras.size()<<std::endl;
            std::cout<<"Load_parameters successfully!"<<std::endl;
        }
        
        void build_wcf_index(char* graphFileName,char* labelFileName,char* labelSizeFileName,char* outputDirName,char* analysisDirName,char* queryFreqDirName,char* givenOrderFileName,const Processing::calcCoefficient<double>& calcCoef,int maxDegree,int indexing_model,int ordering_model,int threads_num,int graphType,bool is_load_size=false,int is_sorted=0,int hfRate=0,int experiment_model=3){
            bool is_multiThreads=false;
            int bestDegree=maxDegree;
            if(threads_num!=0) is_multiThreads=true;
            char* queryPairDirName="";
            double _labeling_time,_ordering_time,_compute_rmq_time=0;
            if(DIRECTED_FLAG==true){
                if(is_multiThreads){
                    Integrated_construction integrated_construction(graphFileName,labelFileName,labelSizeFileName,outputDirName,analysisDirName,queryFreqDirName,queryPairDirName,hfRate,maxDegree,bestDegree,indexing_model,ordering_model,calcCoef,experiment_model,threads_num,graphType,givenOrderFileName,is_sorted);
                    integrated_construction.get_time(_labeling_time,_ordering_time,_compute_rmq_time);
                    if(is_load_size) {
                        vector<double> ().swap(_node_size);
                        _node_size.clear();
                        _node_size.resize(numOfVertices,0);
                        integrated_construction.get_label_size(_node_size);
                    }
                }else{
                    Integrated_construction integrated_construction(graphFileName,labelFileName,labelSizeFileName,outputDirName,analysisDirName,queryFreqDirName,queryPairDirName,hfRate,maxDegree,bestDegree,indexing_model,ordering_model,calcCoef,experiment_model,graphType,givenOrderFileName,is_sorted);
                    integrated_construction.get_time(_labeling_time,_ordering_time,_compute_rmq_time);
                    if(is_load_size) {
                        vector<double> ().swap(_node_size);
                        _node_size.clear();
                        _node_size.resize(numOfVertices,0);
                        integrated_construction.get_label_size(_node_size);
                    }
                }
            }else{
                if(is_multiThreads){
                    Integrated_construction integrated_construction(graphFileName,labelFileName,labelSizeFileName,outputDirName,analysisDirName,queryFreqDirName,hfRate,maxDegree,bestDegree,indexing_model,ordering_model,calcCoef,experiment_model,threads_num,givenOrderFileName,is_sorted);
                    integrated_construction.get_time(_labeling_time,_ordering_time,_compute_rmq_time);
                    if(is_load_size) {
                        vector<double> ().swap(_node_size);
                        _node_size.clear();
                        _node_size.resize(numOfVertices,0);
                        integrated_construction.get_label_size(_node_size);
                    }
                }else{
                    Integrated_construction integrated_construction(graphFileName,labelFileName,labelSizeFileName,outputDirName,analysisDirName,queryFreqDirName,hfRate,maxDegree,bestDegree,indexing_model,ordering_model,calcCoef,experiment_model,givenOrderFileName,is_sorted);
                    integrated_construction.get_time(_labeling_time,_ordering_time,_compute_rmq_time);
                    if(is_load_size) {
                        vector<double> ().swap(_node_size);
                        _node_size.clear();
                        _node_size.resize(numOfVertices,0);
                        integrated_construction.get_label_size(_node_size);
                    }
                }
            }
            _total_labeling_time+=_labeling_time;
            _total_ordering_time+=_ordering_time;
            _total_compute_rmq_time+=_compute_rmq_time;
            return;
        }

        void construct_labels_rl_tip(char* graphFileName,char* parasFileName,char* analyDirName,
                                                    char* labelsDirName,char* labelSizeDirName,char* outputDirName,
                                                    char* givenOrderFileName,int indexing_model,int ordering_model,
                                                    int threads_num,int graphType){
            //prepare parameters for wcf-index construction of each time interval
            load_parameters(parasFileName,_time_interval_cnt);
            string analy_dir_name_prefix(analyDirName);
            if(analy_dir_name_prefix[analy_dir_name_prefix.size()-1]!='/') analy_dir_name_prefix.append("/");
            _slice_result_file=analy_dir_name_prefix+"slice.performance";
            _total_result_file=analy_dir_name_prefix+"total.performance";
            if (access(_slice_result_file.c_str(), 0) == 0)//文件存在
            {
                if (remove(_slice_result_file.c_str()) == 0)
                {
                    //printf("Remove successfully");
                }
                else
                {
                    //printf("Remove failed!");
                }
            }
            if (access(_total_result_file.c_str(), 0) == 0)//文件存在
            {
                if (remove(_total_result_file.c_str()) == 0)
                {
                    //printf("Remove successfully");
                }
                else
                {
                    //printf("Remove failed!");
                }
            }
            string label_file_prefix(labelsDirName);
            if(label_file_prefix[label_file_prefix.size()-1]!='/') label_file_prefix.append("/");
            string label_size_prefix(labelSizeDirName);
            if(label_size_prefix[label_size_prefix.size()-1]!='/') label_size_prefix.append("/");
            string inter_file_prefix(outputDirName);
            if(inter_file_prefix[inter_file_prefix.size()-1]!='/') inter_file_prefix.append("/");
            for(int i=0;i<_time_interval_cnt;++i){
                string point_freq_dir=_new_point_freq_dir_prefix+to_string(i)+"/";
                string label_file_name=label_file_prefix+to_string(i)+".label";
                string label_size_name=label_size_prefix+to_string(i)+".size";
                string inter_file_name=inter_file_prefix+to_string(i);
                build_wcf_index(graphFileName,(char*)label_file_name.c_str(),(char*) label_size_name.c_str(),(char*) inter_file_name.c_str(),(char*) _slice_result_file.c_str(),(char*) point_freq_dir.c_str(),givenOrderFileName,_paras[i],_maxDegree_paras[i],indexing_model,ordering_model,threads_num,graphType);
            }
            //output result in append way
            ofstream ofs(analyDirName,ios::app|ios::out);
			if(!ofs.is_open()) cout<<"Cannot open "<<analyDirName<<endl;
			ofs.setf(ios::fixed);
			ofs.precision(4);
            ofs<<_time_interval_cnt<<" "<<_total_labeling_time*1e6<<" "<<_total_ordering_time*1e6<<" "<<_total_compute_rmq_time*1e6<<" ";
            ofs.close();
            std::cout<<"Construct_labels_rl_tip successfully!"<<std::endl;
            return;
        }

        void construct_labels_gtp(char* graphFileName,char* parasFileName,char* analyDirName,
                                                    char* labelsDirName,char* labelSizeDirName,char* outputDirName,
                                                    char* givenOrderFileName,int indexing_model,int ordering_model,
                                                    int threads_num,int graphType,double up_threshold)
        {
            //initialize variables
            _query_cost.resize(_time_slice_cnt,0);
            double total_qc_time=0;
            int cnt_qc_compute;
            load_parameters(parasFileName,_time_slice_cnt);
            string analy_dir_name_prefix(analyDirName);
            if(analy_dir_name_prefix[analy_dir_name_prefix.size()-1]!='/') analy_dir_name_prefix.append("/");
            _slice_result_file=analy_dir_name_prefix+"slice.performance";
            _total_result_file=analy_dir_name_prefix+"total.performance";
            if (access(_slice_result_file.c_str(), 0) == 0)//文件存在
            {
                if (remove(_slice_result_file.c_str()) == 0)
                {
                    //printf("Remove successfully");
                }
                else
                {
                    //printf("Remove failed!");
                }
            }
            if (access(_total_result_file.c_str(), 0) == 0)//文件存在
            {
                if (remove(_total_result_file.c_str()) == 0)
                {
                    //printf("Remove successfully");
                }
                else
                {
                    //printf("Remove failed!");
                }
            }
            string label_file_prefix(labelsDirName);
            if(label_file_prefix[label_file_prefix.size()-1]!='/') label_file_prefix.append("/");
            string label_size_prefix(labelSizeDirName);
            if(label_size_prefix[label_size_prefix.size()-1]!='/') label_size_prefix.append("/");
            string inter_file_prefix(outputDirName);
            if(inter_file_prefix[inter_file_prefix.size()-1]!='/') inter_file_prefix.append("/");
            string command_removedir,command_makedir,command_chmod;
            if (access(label_file_prefix.c_str(), 0) == 0){
                command_removedir="rm -rf "+ label_file_prefix;
                system(command_removedir.c_str());
                command_makedir="mkdir "+label_file_prefix;
                system(command_makedir.c_str());
                command_chmod="chmod -R 777 "+label_file_prefix;
                system(command_chmod.c_str());
            }
            if (access(label_size_prefix.c_str(), 0) == 0){
                command_removedir="rm -rf "+ label_size_prefix;
                system(command_removedir.c_str());
                command_makedir="mkdir "+label_size_prefix;
                system(command_makedir.c_str());
                command_chmod="chmod -R 777 "+label_size_prefix;
                system(command_chmod.c_str());
            }
            if (access(inter_file_prefix.c_str(), 0) == 0){
                command_removedir="rm -rf "+ inter_file_prefix;
                system(command_removedir.c_str());
                command_makedir="mkdir "+inter_file_prefix;
                system(command_makedir.c_str());
                command_chmod="chmod -R 777 "+inter_file_prefix;
                system(command_chmod.c_str());
            }
            int curr_index=0,k=0;
            bool isUpdate=false;
            _time_interval_cnt=0;
            //greedy update
            while(curr_index<_time_slice_cnt){
                //construct label
                string point_freq_file=_point_freq_dir_prefix+_file_names[curr_index];
                string label_file_name=label_file_prefix+to_string(_time_interval_cnt)+".label";
                string label_size_name=label_size_prefix+to_string(_time_interval_cnt)+".size";
                string inter_file_name=inter_file_prefix+to_string(_time_interval_cnt);
                build_wcf_index(graphFileName,(char*)label_file_name.c_str(),(char*) label_size_name.c_str(),(char*) inter_file_name.c_str(),(char*) _slice_result_file.c_str(),(char*) point_freq_file.c_str(),givenOrderFileName,_paras[curr_index],_maxDegree_paras[curr_index],indexing_model,ordering_model,threads_num,graphType,true);
                //compute query cost
                double qc_time=time_util::GetCurrentTimeSec();
                double c_query_cost=get_time_slice_query_cost((char*) point_freq_file.c_str());
                qc_time=time_util::GetCurrentTimeSec()-qc_time;
                total_qc_time+=qc_time;
                _query_cost[curr_index]=c_query_cost;
                isUpdate=false;
                for(k=curr_index+1;k<_time_slice_cnt;++k)
                {
                    //compute query cost
                    string k_point_freq_file=_point_freq_dir_prefix+_file_names[k];
                    double qc_time=time_util::GetCurrentTimeSec();
                    double k_query_cost=get_time_slice_query_cost((char*) k_point_freq_file.c_str());
                    qc_time=time_util::GetCurrentTimeSec()-qc_time;
                    ++cnt_qc_compute;
                    total_qc_time+=qc_time;
                    //judge whether update
                    if(isSameCluster(c_query_cost,k_query_cost)){
                        _query_cost[k]=k_query_cost;
                    }else{
                        isUpdate=true;
                    }
                    if(isUpdate){
                        ++_time_interval_cnt;
                        _partition_indexes.push_back(k);
                        break;
                    }
                }
                curr_index=k;
            }
            _partition_indexes.push_back(_time_slice_cnt);
            ++_time_interval_cnt;
            //debug
            cout<<"_partition_indexes: size="<<_partition_indexes.size()<<endl;
            for(size_t i=0;i<_partition_indexes.size();++i) cout<<_partition_indexes[i]<<" ";
            cout<<endl;
            //save result to file
            ofstream ofs(analyDirName,ios::app|ios::out);
			if(!ofs.is_open()) cout<<"Cannot open "<<analyDirName<<endl;
			ofs.setf(ios::fixed);
			ofs.precision(4);
            ofs<<_time_interval_cnt<<" "<<_total_labeling_time*1e6<<" "<<_total_ordering_time*1e6<<" "<<_total_compute_rmq_time*1e6<<" ";
            ofs.close();
            std::cout<<"_time_interval_cnt="<<_time_interval_cnt<<std::endl;
            std::cout<<"Construct_labels_gtp successfully!"<<std::endl;
        }

};

#endif