/*
 * @Descripttion: 
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-10-15 09:17:30
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-02-24 16:19:55
 */
#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <stdint.h>
#include <limits.h>
#include <sys/stat.h>
#include "./paras.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges

	/**
	 * @description: sort the files by day index and cut index
	 * @param {*}
	 * @return {*}
	 * @author: Wan Jingyi
	 */
	struct CmpByDayIndex
	{  
		bool operator()(const std::string& k1, const std::string& k2) 
		{ 
			int day_id_1=std::stoi(k1.substr(k1.find("y")+1,k1.find("-")-k1.find("y")-1));
			int cut_id_1=std::stoi(k1.substr(k1.find("t")+1,k1.find(".")-k1.find("t")-1));
			int day_id_2=std::stoi(k2.substr(k2.find("y")+1,k2.find("-")-k2.find("y")-1));
			int cut_id_2=std::stoi(k2.substr(k2.find("t")+1,k2.find(".")-k2.find("t")-1));
			if (day_id_1== day_id_2)
				return cut_id_1 < cut_id_2;
			else{
				return day_id_1 < day_id_2;
			}
		} 
	};

    /*
        *@description: get all files of directory
        *@author: wanjingyi
        *@date: 2020-12-29
    */
    void get_filelist_from_dir(std::string _path, std::vector<std::string>& _files,bool isSorted=false)
    {
        DIR* dir;	
        dir = opendir(_path.c_str());
        struct dirent* ptr;
        std::vector<std::string> file;
        while((ptr = readdir(dir)) != NULL)
        {
            if(ptr->d_name[0] == '.') {continue;}
            file.push_back(ptr->d_name);
        }
        closedir(dir);
        //if(isSorted) std::sort(file.begin(),file.end());
        if(isSorted) std::sort(file.begin(),file.end(),CmpByDayIndex());
        _files = file;
    }

    /*
        *@description: get all files of directory
        *@author: wanjingyi
        *@date: 2020-12-29
    */
    void get_filelist_from_dir_string(std::string _path, std::vector<std::string>& _files,bool isSorted=false)
    {
        DIR* dir;	
        dir = opendir(_path.c_str());
        struct dirent* ptr;
        std::vector<std::string> file;
        while((ptr = readdir(dir)) != NULL)
        {
            if(ptr->d_name[0] == '.') {continue;}
            file.push_back(ptr->d_name);
        }
        closedir(dir);
        if(isSorted) std::sort(file.begin(),file.end());
        _files = file;
    }

	int get_day_index(std::string file){
		int day_id=std::stoi(file.substr(file.find("y")+1,file.find("-")-file.find("y")-1));
		return day_id;
	}

	int get_cut_index(std::string file){
		int cut_id=std::stoi(file.substr(file.find("t")+1,file.find(".")-file.find("t")-1));
		return cut_id;
	}

    /*
     *@description: according the files get the full paths and sort
     *@author: wanjingyi
     *@date: 2021-01-06
    */
    void append_to_full_path(const std::string prefix,std::vector<std::string>& files)
    {
        for(int i=0;i<files.size();++i)
        {
            std::string freqidFilename=files[i];//pure filename
		    std::string freqidFilename_full_prefix(prefix);
			if(freqidFilename_full_prefix[freqidFilename_full_prefix.size()-1]!='/') freqidFilename_full_prefix.append("/");
		    std::string freqidFilename_full_path=freqidFilename_full_prefix.append(freqidFilename);//full path name
            files[i]=freqidFilename_full_path;
        }
    }
   
    void append_to_full_path(std::string prefix,std::vector<std::string> files,std::vector<std::string>& full_files)
    {
        for(int i=0;i<files.size();++i)
        {
            std::string freqidFilename=files[i];//pure filename
		    std::string freqidFilename_full_prefix(prefix);
			if(freqidFilename_full_prefix[freqidFilename_full_prefix.size()-1]!='/') freqidFilename_full_prefix.append("/");
		    std::string freqidFilename_full_path=freqidFilename_full_prefix.append(freqidFilename);//full path name
            full_files.push_back(freqidFilename_full_path);
        }
    }

	void load_workload_hf_freq(char* load_filename, std::vector<int>& freq){
		std::ifstream in(load_filename);//input query file to ifstream
		if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
		int s,t,f;
		while(in>>s>>t>>f){
			freq[s]+=f;
			freq[t]+=f;
		}
		in.close();
		return;
	}

        /*
         *@description: read s-t query freq data from file
         *@author: wanjingyi
         *@date: 2021-01-07
        */
		int load_query_pair_time(char* load_filename,std::vector<std::pair<int, int> >& queries,int multi=1){
			int cnt=0;//count the total times
			std::ifstream in(load_filename);//input query file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int s,t; 
			int query_time=0; 
			char line[24];
			int query_cnt=0;
            //int same_cnt=0;//to be deleted
			 while (in.getline(line,sizeof(line)))
			 {
				 std::stringstream ls(line);
				 ls>>s>>t>>query_time;
                //  std::cout<<s<<" "<<t<<" "<<query_time<<endl;
                //  if(s==t){
                //      same_cnt+=query_time;
                //      continue;
                //  }
                 for(int i=0;i<query_time*multi;++i) queries.push_back(std::make_pair(s,t));
				 cnt+=query_time*multi;
			 }
			 in.close();
			 //cout<<"s-t same num="<<same_cnt<<" ratio="<<(double)same_cnt/(double)cnt<<endl;//to be deleted
			 return cnt;
		}

        /*
         *@description: read s-t query freq data from file
         *@author: wanjingyi
         *@date: 2021-01-07
        */
		int load_query_pair_time(char* load_filename,std::vector<std::pair<int, int> >& queries,int numQuery,bool is_num_query){
			std::ifstream in(load_filename);//input query file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int s,t; 
			int query_time=0; 
			char line[24];
			int query_cnt=0;
			std::vector<std::pair<int, int> > queries_real;
            //int same_cnt=0;//to be deleted
			while (in.getline(line,sizeof(line)))
			{
				std::stringstream ls(line);
				ls>>s>>t>>query_time;
				for(int i=0;i<query_time;++i) queries_real.push_back(std::make_pair(s,t));
			}
			in.close();
			int cnt=numQuery;
			// while(cnt>0){
			// 	if(cnt>=queries_real.size()){
			// 		queries.insert(queries.end(),queries_real.begin(),queries_real.end());
			// 		cnt-=queries_real.size();
			// 	}else{
			// 		queries.insert(queries.end(),queries_real.begin(),queries_real.begin()+cnt);
			// 		cnt-=cnt;
			// 	}
			// }
			while(cnt>=queries_real.size()){
				queries.insert(queries.end(),queries_real.begin(),queries_real.end());
				cnt-=queries_real.size();
			}
			if(queries.size()!=numQuery) std::cout<<"Error: queries.size()!=numQuery!"<<" numQuey="<<numQuery<<" queries.size()="<<queries.size()<<std::endl;
			return queries.size();
		}

		int load_query_pair_time_filter(char* load_filename,std::vector<std::pair<int, int> >& queries,int numQuery,bool* isDeleted){
			std::ifstream in(load_filename);//input query file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int s,t; 
			int query_time=0; 
			std::vector<std::pair<int, int> > queries_real;
			while (in>>s>>t>>query_time)
			{
				if(isDeleted[s]||isDeleted[t]) continue;
				for(int i=0;i<query_time;++i) queries_real.push_back(std::make_pair(s,t));
			}
			in.close();
			int cnt=numQuery;
			while(cnt>=queries_real.size()){
				queries.insert(queries.end(),queries_real.begin(),queries_real.end());
				cnt-=queries_real.size();
			}
			if(cnt>0) queries.insert(queries.end(),queries_real.begin(),queries_real.begin()+cnt);
			if(queries.size()!=numQuery) std::cout<<"Error: queries.size()!=numQuery!"<<" numQuey="<<numQuery<<" queries.size()="<<queries.size()<<std::endl;
			return queries.size();
		}

	/**
	 * @description: judege whether is a directory
	 * @param {*}
	 * @return {*} true-directory, false-file
	 * @author: Wan Jingyi
	 */ 
 	bool is_directory(char* path){
		struct stat s;
		if (stat(path,&s)==0){
			if(s.st_mode & S_IFDIR){
				//std::cout<<"it's a directory"<<std::endl;
				return true;
			}else if (s.st_mode & S_IFREG){
				//std::cout<<"it's a file"<<std::endl;
				return false;
			}else{
				std::cout<<"not file not directory"<<std::endl;
				exit(-1);
			}
		}else{
			std::cout<<"error:"<<path<<" doesn't exist"<<std::endl;
			exit(-1);
		}
	}

#endif