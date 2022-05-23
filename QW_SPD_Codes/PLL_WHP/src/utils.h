/*
 * @Descripttion:  some util functions
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-10 21:24:40
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-02-13 20:24:07
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

//******************parallel functions****************

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

	int get_day_index(std::string file){
		int day_id=std::stoi(file.substr(file.find("y")+1,file.find("-")-file.find("y")-1));
		return day_id;
	}

	int get_cut_index(std::string file){
		int cut_id=std::stoi(file.substr(file.find("t")+1,file.find(".")-file.find("t")-1));
		return cut_id;
	}

// Compare and Swap
template <typename V_T>
bool CAS(V_T *ptr, V_T old_val, V_T new_val)
//inline bool CAS(void *ptr, V_T old_val, V_T new_val)
{
//	return __atomic_compare_exchange(ptr, &old_val, &new_val, false, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST);
	if (1 == sizeof(V_T)) {
		return __atomic_compare_exchange(reinterpret_cast<uint8_t *>(ptr), reinterpret_cast<uint8_t *>(&old_val),
										 reinterpret_cast<uint8_t *>(&new_val), false, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST);
	} else if (2 == sizeof(V_T)) {
		return __atomic_compare_exchange(reinterpret_cast<uint16_t *>(ptr), reinterpret_cast<uint16_t *>(&old_val),
										 reinterpret_cast<uint16_t *>(&new_val), false, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST);
	} else if (4 == sizeof(V_T)) {
		return __atomic_compare_exchange(reinterpret_cast<uint32_t *>(ptr), reinterpret_cast<uint32_t *>(&old_val),
										 reinterpret_cast<uint32_t *>(&new_val), false, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST);
	} else if (8 == sizeof(V_T)) {
		return __atomic_compare_exchange(reinterpret_cast<uint64_t *>(ptr), reinterpret_cast<uint64_t *>(&old_val),
										 reinterpret_cast<uint64_t *>(&new_val), false, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST);
	} else {
		printf("CAS cannot support the type.\n");
		exit(EXIT_FAILURE);
	}
}



        /*
         *@description: read s-t query freq data from file
         *@author: wanjingyi
         *@date: 2021-01-07
        */
		int load_query_pair_time(char* query_freq_filename,char* query_pair_filename,std::vector<std::pair<int, int> >& queries,int hf_rate,int &numOfHFpoint,int multi=1){
			//tmps
            char line[24];
            int id,i=0,s,t,query_time=0; 
            //load hf point
            numOfHFpoint=0;
            std::vector<bool> HFinGraphIndex(numOfVertices,false);
            if(hf_rate==0) numOfHFpoint=numOfVertices;
			else numOfHFpoint=static_cast<int>((double)(numOfVertices*hf_rate)/(double)HF_DIVIDION);
			std::cout<<"initial numOfHFpoint="<<numOfHFpoint<<" initial hfRate="<<hf_rate<<std::endl;
			std::ifstream in_freq(query_freq_filename);//input HFPoint file to ifstream
			if(!in_freq.is_open()) {std::cerr<<"Cannot open "<<query_freq_filename<<std::endl;}
			//read each line representing HFpoint to vector 
			while (in_freq.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>id>>query_time;		
				if(i>=numOfHFpoint) break;
                HFinGraphIndex[id]=true;
				i++;
			}
			if(i<numOfHFpoint){
				numOfHFpoint=i;
			}
			hf_rate=(double)numOfHFpoint/(double)numOfVertices;
            std::cout<<"real numOfHFpoint="<<numOfHFpoint<<" real hfRate="<<hf_rate<<std::endl;
            //load query pair
            int cnt=0;//count the total times
			std::ifstream in_pair(query_pair_filename);//input query file to ifstream
			if(!in_pair.is_open()) {std::cerr<<"Cannot open "<<query_pair_filename<<std::endl;}
            while (in_pair.getline(line,sizeof(line)))
            {
                std::stringstream ls(line);
                ls>>s>>t>>query_time;
                if(HFinGraphIndex[s]&&HFinGraphIndex[t]){
                    for(int i=0;i<query_time*multi;++i) queries.push_back(std::make_pair(s,t));
                    cnt+=query_time*multi;
                }
            }
            in_pair.close();
            return cnt;
		}

        /*
         *@description: load_query_num for unsigned int pointer
         *@author: wanjingyi
         *@date: 2021-01-14
        */
		int load_query_pair_time(char* query_freq_filename,char* query_pair_filename,unsigned int* sp,unsigned int* tp,int hf_rate,int& numOfHFpoint,int multi=1){
			//tmps
            char line[24];
            int id,i=0,s,t,query_time=0; 
            //load hf point
            numOfHFpoint=0;
            std::vector<bool> HFinGraphIndex(numOfVertices,false);
            if(hf_rate==0) numOfHFpoint=numOfVertices;
			else numOfHFpoint=static_cast<int>((double)(numOfVertices*hf_rate)/(double)HF_DIVIDION);
			std::cout<<"initial numOfHFpoint="<<numOfHFpoint<<" initial hfRate="<<hf_rate<<std::endl;
			std::ifstream in_freq(query_freq_filename);//input HFPoint file to ifstream
			if(!in_freq.is_open()) {std::cerr<<"Cannot open "<<query_freq_filename<<std::endl;}
			//read each line representing HFpoint to vector 
			while (in_freq.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>id>>query_time;		
				if(i>=numOfHFpoint) break;
                HFinGraphIndex[id]=true;
				i++;
			}
			if(i<numOfHFpoint){
				numOfHFpoint=i;
			}
			hf_rate=(double)numOfHFpoint/(double)numOfVertices;
            std::cout<<"real numOfHFpoint="<<numOfHFpoint<<" real hfRate="<<hf_rate<<std::endl;
            int cnt=0;//count the total times
			std::ifstream in_pair(query_pair_filename);//input query file to ifstream
			if(!in_pair.is_open()) {std::cerr<<"Cannot open "<<query_pair_filename<<std::endl;}
            //int same_cnt=0;//to be deleted
            while (in_pair.getline(line,sizeof(line)))
            {
                std::stringstream ls(line);
                ls>>s>>t>>query_time;
                if(!HFinGraphIndex[s]||!HFinGraphIndex[t]) continue;
                for(int i=0;i<query_time*multi;++i){
                    *sp++=s;
                    *tp++=t;
                }
				cnt+=query_time*multi;
            }
            in_pair.close();
            return cnt;
		}

//******************parallel functions****************

    //compare string by ascending number value order
    struct CmpByKeyNum 
    {  
        bool operator()(const std::string& k1, const std::string& k2) 
        {  
            if (k1.length() != k2.length())
                return k1.length() < k2.length();
            else
                return k1 < k2;
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

    void get_filelist_from_dir(std::string _path, std::vector<std::string>& _files,bool isSorted,int index)
    {
        DIR* dir;	
        dir = opendir(_path.c_str());
        struct dirent* ptr;
        std::vector<std::string> file;
        while((ptr = readdir(dir)) != NULL)
        {
            if(ptr->d_name[0] == '.') {continue;}
			std::string file_name=ptr->d_name;
			std::string id=file_name.substr(file_name.find("y")+1,file_name.find("-")-file_name.find("y")-1);
            if(id==std::to_string(index)) file.push_back(file_name);
        }
        closedir(dir);
        //if(isSorted) std::sort(file.begin(),file.end());
        if(isSorted) std::sort(file.begin(),file.end(),CmpByKeyNum());
        _files = file;
    }

    void get_filelist_from_dir(std::string _path, std::vector<std::string>& _files,std::string match_prefix,bool isSorted)
    {
        DIR* dir;	
        dir = opendir(_path.c_str());
        struct dirent* ptr;
        std::vector<std::string> file;
        while((ptr = readdir(dir)) != NULL)
        {
            if(ptr->d_name[0] == '.') {continue;}
			std::string file_name=ptr->d_name;
			std::string prefix=file_name.substr(0,file_name.find("-"));
            if(prefix.compare(match_prefix)==0) file.push_back(ptr->d_name);
        }
        closedir(dir);
        //if(isSorted) std::sort(file.begin(),file.end());
        if(isSorted) std::sort(file.begin(),file.end(),CmpByKeyNum());
        _files = file;
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

 /**
  * @Author: wanjingyi
  * @description: read vertices order from file for(getBetweennessOrderFromFIle getCoverageOrderFromFIle...)
  * @param {*}
  */
	void get_order_from_file(char* load_filename,std::vector<int>& rank,std::vector<int>& inv,int n){
            // rank.resize(numOfVertices);
			// inv.resize(numOfVertices);
			std::ifstream in(load_filename);//input betwenness file to ifstream
            if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
            int id;
            for(int i=0;i<n;++i){
                in>>id;
                rank[id]=i;
                inv[i]=id;
				//to be deleted
            }
            in.close();
            std::cout<<"get order from "<<load_filename<<" finished!"<<std::endl;
	}

        /*
         *@description: read s-t query freq data from file
         *@author: wanjingyi
         *@date: 2021-01-07
        */
		int load_query_pair_time_overlay(char* load_filename,char* newToOriginalFileName,std::vector<std::pair<int, int> >& queries,int multi=1){
			//load new to original
			int numOfOriginalVertices,numOfOverlayVertices;
			std::vector<int> newToOriginal;
			std::vector<int> originalToNew;
			std::ifstream in(newToOriginalFileName);
			if(!in.is_open()) std::cerr<<newToOriginalFileName<<" cannot be opened!"<<std::endl;
            in>>numOfOriginalVertices>>numOfOverlayVertices;
			std::cout<<"numOfOriginalVertices="<<numOfOriginalVertices<<" numOfOverlayVertices="<<numOfOverlayVertices<<std::endl;
			originalToNew.resize(numOfOriginalVertices, -1);
			newToOriginal.resize(numOfOverlayVertices);
			int u,v;
            int i;
			//build the map relationship 
			for (i = 0; in >> u >> v;i++) {
				newToOriginal[u]=v;
				originalToNew[v]=u;
			}
			if(i!=numOfOverlayVertices) std::cerr<<"readlines error:i!=numOfOverlayVertices!"<<std::endl;
			in.close();
			int cnt=0;//count the total times
			std::ifstream in_query(load_filename);//input query file to ifstream
			if(!in_query.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int s,t,s_new,t_new; 
			int query_time=0; 
			char line[24];
			int query_cnt=0;
			 while (in_query.getline(line,sizeof(line)))
			 {
				std::stringstream ls(line);
				ls>>s>>t>>query_time;
				s_new=originalToNew[s];
				t_new=originalToNew[t];
				if(s_new!=-1&&t_new!=-1){
					for(int i=0;i<query_time*multi;++i) queries.push_back(std::make_pair(s_new,t_new));
					cnt+=query_time*multi;
				}
			 }
			 in_query.close();
			 return cnt;
		}		

        /*
         *@description: read s-t query freq data from file
         *@author: wanjingyi
         *@date: 2021-01-07
        */
		int load_query_pair_time_overlay(char* load_filename,char* newToOriginalFileName,std::vector<std::pair<int, int> >& queries,int numQuery,bool isNumQuery){
			//load new to original
			int numOfOriginalVertices,numOfOverlayVertices;
			std::vector<int> newToOriginal;
			std::vector<int> originalToNew;
			std::ifstream in(newToOriginalFileName);
			if(!in.is_open()) std::cerr<<newToOriginalFileName<<" cannot be opened!"<<std::endl;
            in>>numOfOriginalVertices>>numOfOverlayVertices;
			std::cout<<"numOfOriginalVertices="<<numOfOriginalVertices<<" numOfOverlayVertices="<<numOfOverlayVertices<<std::endl;
			originalToNew.resize(numOfOriginalVertices, -1);
			newToOriginal.resize(numOfOverlayVertices);
			int u,v;
            int i;
			//build the map relationship 
			for (i = 0; in >> u >> v;i++) {
				newToOriginal[u]=v;
				originalToNew[v]=u;
			}
			if(i!=numOfOverlayVertices) std::cerr<<"readlines error:i!=numOfOverlayVertices!"<<std::endl;
			in.close();
			std::ifstream in_query(load_filename);//input query file to ifstream
			if(!in_query.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int s,t,s_new,t_new; 
			int query_time=0; 
			char line[24];
			int query_cnt=0;
			std::vector<std::pair<int, int> > queries_real;
			while (in_query.getline(line,sizeof(line)))
			{
				std::stringstream ls(line);
				ls>>s>>t>>query_time;
				s_new=originalToNew[s];
				t_new=originalToNew[t];
				if(s_new!=-1&&t_new!=-1){
					for(int i=0;i<query_time;++i) queries_real.push_back(std::make_pair(s_new,t_new));
				}
			}
			in_query.close();
			int cnt=numQuery;
			while(cnt>0){
				if(cnt>=queries_real.size()){
					queries.insert(queries.end(),queries_real.begin(),queries_real.end());
					cnt-=queries_real.size();
				}else{
					queries.insert(queries.end(),queries_real.begin(),queries_real.begin()+cnt);
					cnt-=cnt;
				}
			}
			if(queries.size()!=numQuery) std::cout<<"Error: queries.size()!=numQuery!"<<" numQuey="<<numQuery<<" queries.size()="<<queries.size()<<std::endl;
			return queries.size();
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
                //  std::cout<<s<<" "<<t<<" "<<query_time<<std::endl;
                //  if(s==t){
                //      same_cnt+=query_time;
                //      continue;
                //  }
                 for(int i=0;i<query_time*multi;++i) queries.push_back(std::make_pair(s,t));
				 cnt+=query_time*multi;
			 }
			 in.close();
			 //cout<<"s-t same num="<<same_cnt<<" ratio="<<(double)same_cnt/(double)cnt<<std::endl;//to be deleted
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
			while(cnt>0){
				if(cnt>=queries_real.size()){
					queries.insert(queries.end(),queries_real.begin(),queries_real.end());
					cnt-=queries_real.size();
				}else{
					queries.insert(queries.end(),queries_real.begin(),queries_real.begin()+cnt);
					cnt-=cnt;
				}
			}
			if(queries.size()!=numQuery) std::cout<<"Error: queries.size()!=numQuery!"<<" numQuey="<<numQuery<<" queries.size()="<<queries.size()<<std::endl;
			return queries.size();
		}

        /*
         *@description: load_query_num for unsigned int pointer
         *@author: wanjingyi
         *@date: 2021-01-14
        */
		int load_query_pair_time(char* load_filename,unsigned int* sp,unsigned int* tp,int multi=1){
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
                //  if(s==t){
                //      same_cnt+=query_time;
                //      continue;
                //  }
                 for(int i=0;i<query_time*multi;++i){
                     *sp++=s;
                     *tp++=t;
                 }
				 cnt+=query_time*multi;
			 }
			 in.close();
			 //cout<<"s-t same num="<<same_cnt<<" ratio="<<(double)same_cnt/(double)cnt<<std::endl;//to be deleted
			 return cnt;
		}

        /*
		 *@description: load hfpoint
		 *@author: wanjingyi
		 *@date: 2021-01-16
		*/
		void load_hfpoint(char* load_filename,int& numOfHfpoint){//5%%
			std::cout<<"initial numOfHfpoint = "<<numOfHfpoint<<std::endl;
			std::ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int id;int t,i=0;
			char line[24];
			//read each line representing HFpoint to vector 
			while (in.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>id>>t;		
				if(i>=numOfHfpoint)	break;
				i++;
			}
			if(i<numOfHfpoint){
				numOfHfpoint=i;
				std::cout<<"real numOfHfpoint = "<<numOfHfpoint<<std::endl;
			}
		}		

        /*
		 *@description: load hfpoint
		 *@author: wanjingyi
		 *@date: 2021-01-16
		*/
		void  load_hfpoint(char* load_filename,int& numOfHfpoint,std::vector<int>& HFPoint,std::vector<bool>& is_hfpoint){//5%%
			std::cout<<"initial numOfHfpoint = "<<numOfHfpoint<<std::endl;
			std::ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int id;int t,i=0;
			char line[24];
			//read each line representing HFpoint to vector 
			while (in.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>id>>t;		
				if(i>=numOfHfpoint)	break;
				HFPoint.push_back(id);
                is_hfpoint[id]=true;
				i++;
			}
			if(i<numOfHfpoint){
				numOfHfpoint=i;
				std::cout<<"real numOfHfpoint = "<<numOfHfpoint<<std::endl;
			}
		}		

        /*
		 *@description: load hfpoint
		 *@author: wanjingyi
		 *@date: 2021-03-01
		*/
		void load_hfpoint_and_qt(char* load_filename,int& numOfHfpoint,std::vector<int>& HFPoint,std::vector<bool>& is_hfpoint,std::vector<unsigned int>& queryTime,int& maxQueryTime,long long& totalQueryTime){//5%%
			std::cout<<"initial numOfHfpoint = "<<numOfHfpoint;
			std::ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int id;int t,i=0;
			char line[24];
			maxQueryTime=0;
			totalQueryTime=0;
			//read each line representing HFpoint to vector 
			while (in.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>id>>t;		
				queryTime[id]=t;
				totalQueryTime+=t;
				if(t>maxQueryTime) maxQueryTime=t;
				if(i>=numOfHfpoint)	continue;
				HFPoint.push_back(id);
                is_hfpoint[id]=true;
				i++;
			}
			if(i<numOfHfpoint){
				numOfHfpoint=i;
			}
			std::cout<<" real numOfHfpoint = "<<numOfHfpoint<<" maxQueryTime = "<<maxQueryTime<<" totalQueryTime = "<<totalQueryTime<<std::endl;
		}

        /*
		 *@description: load hfpoint
		 *@author: wanjingyi
		 *@date: 2021-01-16
		*/
		void load_hfpoint(char* load_filename,int& numOfHfpoint,std::vector<bool>& is_hfpoint,std::vector<int>& rank){//5%%
			std::cout<<"initial numOfHfpoint = "<<numOfHfpoint<<std::endl;
			std::ifstream in(load_filename);//input HFPoint file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<load_filename<<std::endl;}
			int id;int t,i=0;
			char line[24];
			//read each line representing HFpoint to vector 
			while (in.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>id>>t;		
				if(i>=numOfHfpoint)	break;
                is_hfpoint[id]=true;
                rank[id]=i;
				i++;
			}
			if(i<numOfHfpoint){
				numOfHfpoint=i;
				std::cout<<"real numOfHfpoint = "<<numOfHfpoint<<std::endl;
			}
		}

    	/**
		 * @Author: wanjingyi
		 * @description: load new node ids and original ids
		 * @param {*}
		 * @return {*}
		 */
		void load_id_map(char* newToOriginalFileName,std::vector<int>& newToOriginal,std::vector<int>& originalToNew,int& numOfOriginalVertices,int& numOfOverlayVertices){
            std::cout<<"*******************load id map start!*******************"<<std::endl;
			std::ifstream in(newToOriginalFileName);
			if(!in.is_open()) std::cerr<<newToOriginalFileName<<" cannot be opened!"<<std::endl;
            in>>numOfOriginalVertices>>numOfOverlayVertices;
			std::cout<<"numOfOriginalVertices="<<numOfOriginalVertices<<" numOfOverlayVertices="<<numOfOverlayVertices<<std::endl;
			originalToNew.resize(numOfOriginalVertices, -1);
			newToOriginal.resize(numOfOverlayVertices);
			int u,v;
            int i;
			//build the map relationship 
			for (i = 0; in >> u >> v;i++) {
				newToOriginal[u]=v;
				originalToNew[v]=u;
			}
			if(i!=numOfOverlayVertices) std::cerr<<"readlines error:i!=numOfOverlayVertices!"<<std::endl;
			in.close();
            std::cout<<"*******************load id map finished!*******************"<<std::endl;
		}		

		template<typename valueType>
		void load_workload_directed(char* queryPair_file, std::vector<valueType>& queryTime_s,vector<valueType>& queryTime_t,int n){
			//initilaize
			if(queryTime_s.empty()) queryTime_s.resize(n,0);
			if(queryTime_t.empty()) queryTime_t.resize(n,0);
			//read from file
			std::ifstream in(queryPair_file);//input HFPoint file to ifstream
			if(!in.is_open()) {std::cerr<<"Cannot open "<<queryPair_file<<std::endl;}
			int s,t;
			valueType qt;
			char line[24];
			while (in.getline(line,sizeof(line))){
				std::stringstream hp(line);
				hp>>s>>t>>qt;
				queryTime_s[s]+=qt;
				queryTime_t[t]+=qt;
			}
			in.close();
			std::cout<<"Load workload directed successfully!"<<std::endl;
			return;
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