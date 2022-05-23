/*
 * @Descripttion: h2h-index
 * @version: 
 * @Author: wanjingyi
 * @Date: 2021-02-20 22:11:22
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-10-16 11:21:44
 */
#ifndef _H2H_INDEX_H
#define _H2H_INDEX_H

#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <map>
#include <memory>
#include <algorithm>
#include <bitset>
#include <limits>
#include <climits>
#include <malloc.h>
#include <sys/time.h>
#include <xmmintrin.h> 
#include "./paras.h"
#include "./utils.h"
#include "./time_util.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define INF_WEIGHT SP_Constants::INF_WEIGHT

using namespace std;
using namespace time_util;
typedef vector<int>::iterator iterator_type;

struct h2h_index_t_p
{
	int* dis;
	int* pos;
	int height;
	int width;
	h2h_index_t_p(int h,int w):height(h),width(w){}
}__attribute__((aligned(64)));

struct h2h_bindex_t_p
{
	int* dis;
    int* r_dis;
	int* pos;
	int height;
	int width;
	h2h_bindex_t_p(int h,int w):height(h),width(w){}
}__attribute__((aligned(64)));

class Naive_rmq{
public:
    //2D array of precomputed
    //_arr[b][a] is the index of the minimum value in the range[a,a+b+1]
    vector<vector<int> > _arr;
    int _n;
    iterator_type _begin;
    iterator_type _end;
    Naive_rmq(iterator_type b,iterator_type e):_n(e-b),_begin(b),_end(e)
    {
        _arr.resize(_n);
        //fill(_arr.begin(),_arr.end(),vector<int>());
        fill_in();
    }

    int query(int u_index,int v_index){
        return _arr[v_index-u_index][u_index];
    }
protected:
    /**
    *Dynamic program to compute the answers to every possible query on the input
    **/
    void fill_in()
    {
        //The first level_arr[0] contains the answers to RMQ queries on interval of length 1,
        //which must just be the first element in the interval
        vector<int>& _arr_0=_arr[0];
        _arr_0.resize(_n);
        for(size_t i=0;i<_n;++i) _arr_0[i]=i;
        for(auto it=_arr.begin();it+1<_arr.end();++it){
            //Zips neighboring indexes of this level with a functor that chooses the index
            // producing the lesser value from each pair of indexes
            transform(it->begin(),it->end()-1,it->begin()+1,back_inserter(*(it+1)),
                      [this](const size_t& x,const size_t& y){
                        return _begin[x]<_begin[y]?x:y;
                      });
        }
    }
};

class Sparse_rmq{
public:
    //variables
    iterator_type _begin;
    iterator_type _end;
    //logn is the depth we need to precompute answers down to
    int _logn;
    //2D array of precompute answers,_arr[a][b] is the minimum value in the range [a,a+b]
    vector<vector<int> > _arr;
    int _n;
    //constructions
    Sparse_rmq(){
        _logn=0;
        _n=0;
    }

    Sparse_rmq(iterator_type b,iterator_type e):_n(e-b),_begin(b),_end(e){
        _logn=max(1,int(floor(log2(_n))));
        _arr.resize(_logn+1);
        fill_in();
    }

    void initialize(iterator_type b,iterator_type e){
        _n=e-b;
        _begin=b;
        _end=e;
        _logn=max(1,int(floor(log2(_n))));
        _arr.resize(_logn+1);
        fill_in();
    }

    int query(int u_index,int v_index) const
    {
        ++v_index;
        const int depth=max(0,int(floor(log2(v_index-u_index))));
        const int p_u=_arr[depth][u_index];
        const int p_v=_arr[depth][v_index-(1<<depth)];
        return _begin[p_u]<_begin[p_v]?p_u:p_v;
    }
protected:
    void fill_in()
    {
        //Each interval of zero length retuning the value at i
        vector<int>& _arr_0=_arr[0];
        _arr_0.resize(_n);
        for(int i=0;i<_n;++i) _arr_0[i]=i;
        //The depth goes up to logn-2^d
        for(int d=0;d<_logn;++d){
            //Interval is length 2^d,the up limit id (n-2^d)
            int width=1<<d;
            //From the next a array by zipping pairs of elements in the dth array
            //that are width apart,taking the lesser one
            
            transform(_arr[d].begin(),_arr[d].end()- width,
                      _arr[d].begin()+width,back_inserter(_arr[d+1]),
                      [this](const int& x, const int& y){
                        return _begin[x]<_begin[y]?x:y;
                      });
        }
    }
};

class LCA{
public:
    //lca variables
    vector<int> _euler;
    vector<int> _level;
    vector<int> _index;//get the index by node num
    int _n;
    //rmq variables
    int _logn;
    int _block_size;
    /**
    * Arrays of length 2n/lg(n) where the first array contains the minimum
    * element in the ith block of the input, and the second contains the
    * position of that element (as an offset from the beginning of the
    * original input, not from the beginning of the block).
    */
    vector<int> _super_array_vals;
    vector<int> _super_array_ids;
    typedef std::vector<int> block_identifier;//represent a normalized block
    /**
    *Tha map is responsible for each naive_rmq
    *and fro sure we only build once for each of sub block
    **/
    map<block_identifier,unique_ptr<Naive_rmq> > _sub_block_rmqs;
    /**
    *An array mapping sub blocks(by their indexes) to the naive rmq
    **/
    vector<Naive_rmq*> _sub_block_rmq_array;
    //The sparse RMQ implementation over _super_array_vals
    unique_ptr<Sparse_rmq> _super_rmq;
    //for (O(nlong),O(1)) lca query
    Sparse_rmq sparse_rmq;

    //***************constructions**************
    LCA(){
        _n=0;
        _logn=0;
        _block_size=0;
    }

    LCA(char* load_filename){
        //read variables from file
		load_lca_variables(load_filename);
        //construct_rmq_n();
        construct_rmq_nlogn();
    }

    void initialize(char* load_filename){
        //read variables from file
		load_lca_variables(load_filename);
        construct_rmq_nlogn();
    }

    void clear_nlogn(){
        vector<int> ().swap(_euler);
		_euler.clear();
        vector<int> ().swap(_level);
		_level.clear();
        vector<int> ().swap(_index);
		_index.clear();
        if(!sparse_rmq._arr.empty()){
            vector<vector<int> > ().swap(sparse_rmq._arr);
            sparse_rmq._arr.clear();
        }
    }

    #pragma g++ push_options
    #pragma optimize("O0")
    int lca_query(int u,int v){
        int u_index=_index[u];
        int v_index=_index[v];
        if(v_index<u_index){//swap
            int tmp=u_index;
            u_index=v_index;
            v_index=tmp;
        }
        //return _euler[rmq_query_n(u_index,v_index)];
        return _euler[rmq_query_nlogn(u_index,v_index)];
    }
    #pragma g++ pop_options

protected:
    void load_lca_variables(const char* load_filename){
        ifstream ifs(load_filename);
        int t;
        //load _euler
		int esize = 0;
        ifs.read((char*)&esize, sizeof(esize));
        //cout<<"seize="<<esize;//output test
        _euler.resize(esize);
        for(int i=0;i<esize;++i){
            ifs.read((char*)&t,sizeof(t));
            _euler[i]=t;
        }
        //load _level
        int lsize=0;
        ifs.read((char*)&lsize, sizeof(lsize));
        //cout<<"lsize="<<lsize;//output test
        _level.resize(lsize);
        _n=lsize;
        for(int i=0;i<lsize;++i){
            ifs.read((char*)&t,sizeof(t));
            _level[i]=t;
        }
        //load _index
        int isize=0;
        ifs.read((char*)&isize, sizeof(isize));
        //cout<<" isize="<<isize<<endl;//output test
        _index.resize(isize);
        for(int i=0;i<isize;++i){
            ifs.read((char*)&t,sizeof(t));
            _index[i]=t;
        }
        ifs.close();
        //ouput test
        // cout<<"_euler:"<<endl;//to be deleted
        // for(size_t i=0;i<_euler.size();++i) cout<<_euler[i]<<" ";//to be deleted
        // cout<<endl;//to be deleted
        // cout<<"_level:"<<endl;//to be deleted
        // for(size_t i=0;i<_level.size();++i) cout<<_level[i]<<" ";//to be deleted
        // cout<<endl;//to be deleted
        // cout<<"_index:"<<endl;//to be deleted
        // for(size_t i=0;i<_index.size();++i) cout<<_index[i]<<" ";//to be deleted
        // cout<<endl;//to be deleted
    }

    void test_sparse_rmq(){
        cout<<"LCA sparse rmq test:"<<endl;
        Sparse_rmq sparse_rmq(_level.begin(),_level.end());
        for(int i=0;i<numOfVertices;++i)
        {
            for(int j=i;j<numOfVertices;++j)
            {
                cout<<i<<"-"<<j<<":";
                int i_index=_index[i];
                int j_index=_index[j];
                if(j_index<i_index){
                    int tmp=i_index;
                    i_index=j_index;
                    j_index=tmp;
                }
                cout<<_euler[sparse_rmq.query(i_index,j_index)]<<endl;
            }
        }
    }

    void test_naive_rmq(){
        Naive_rmq naive_rmq(_level.begin(),_level.end());
        cout<<"LCA naive rmq test:"<<endl;
        for(int i=0;i<numOfVertices;++i)
        {
            for(int j=i;j<numOfVertices;++j)
            {
                int i_index=_index[i];
                int j_index=_index[j];
                if(j_index<i_index){
                    int tmp=i_index;
                    i_index=j_index;
                    j_index=tmp;
                }
                cout<<i<<"-"<<j<<":";
                cout<<_euler[naive_rmq.query(i_index,j_index)]<<endl;
            }
        }
    }

    int rmq_query_nlogn(int u_index,int v_index){
        return sparse_rmq.query(u_index,v_index);
    }

    int rmq_query_n(int u_index,int v_index)
    {

        int u_block_idx=int(u_index/_block_size);
        int u_offset=int(u_index%_block_size);
        int v_block_idx=int(v_index/_block_size);
        int v_offset=int(v_index%_block_size);
        //cout<<"u_index="<<u_index<<" v_index="<<v_index<<" u_block_idx="<<u_block_idx<<" u_offset="<<u_offset<<" v_block_idx="<<v_block_idx<<" v_offset="<<v_offset<<endl;//to be deleted
        Naive_rmq& u_naive=*_sub_block_rmq_array[u_block_idx];
        Naive_rmq& v_naive=*_sub_block_rmq_array[v_block_idx];
        int block_diff=v_block_idx-u_block_idx;
        if(block_diff==0){
            //cout<<"in the same block:"<<endl;//to be deleted
            //u and v are in the same block.One naive_rmq search suffices
            return(u_block_idx*_block_size)+u_naive.query(u_offset,v_offset);
        }else{
            const iterator_type u_block_end=min(_level.end(),_level.begin()+((u_block_idx+1)*_block_size));
            int u_min_idx=(u_block_idx*_block_size)+u_naive.query(u_offset,u_block_end-(_level.begin()+(_block_size*u_block_idx))-1);
            int v_min_idx=(v_block_idx*_block_size)+v_naive.query(0,v_offset);
            //cout<<"u_min_idx="<<u_min_idx<<" v_min_idx="<<v_min_idx<<endl;//to be deleted
            if(block_diff==1){
                //cout<<"in the adjacent block:"<<endl;//to be deleted
                //u and v are in adjacent block,not query super array
                //it doesn't handle zero-length intervals property
                return _level[u_min_idx]<_level[v_min_idx]?u_min_idx:v_min_idx;
            }else{
                //cout<<"in the interval block:"<<endl;//to be deleted
                //Full algorithm, using the sparse RMQ implementation
                //on the super array between u and v's blocks
                int super_idx=_super_rmq->query(u_block_idx+1,v_block_idx-1);
                const int& u_min_val=_level[u_min_idx];
                const int& v_min_val=_level[v_min_idx];
                if(u_min_val<v_min_val){
                    return u_min_val<_super_array_vals[super_idx]?u_min_idx:_super_array_ids[super_idx];
                }else{
                    return v_min_val<_super_array_vals[super_idx]?v_min_idx:_super_array_ids[super_idx];
                }
            }
        }
    }

    void construct_rmq_nlogn(){
        double qtime= GetCurrentTimeSec();
        sparse_rmq.initialize(_level.begin(),_level.end());
        qtime = GetCurrentTimeSec() - qtime;
        double avg_qtime = qtime/ (double)numOfVertices;
        cout << "rmq construct time:" << qtime * 1e6 <<  " microseconds" << endl;
        cout << "avg rmq construct time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
    }


    void construct_rmq_n(){
        double qtime= GetCurrentTimeSec();
        //block size
        _logn=max(1,int(floor(log2(_n))));
        _block_size=max(1,_logn/2);
        //cout<<"_logn="<<_logn<<" _block_size="<<_block_size<<endl;//to be deleted
        //initialize
        int block_length=ceil((float)_n/(float)_block_size);
        //cout<<" block_size="<<_block_size<<" block_length="<<block_length<<endl;//to be deleted
        _super_array_vals.reserve(block_length);
        _super_array_ids.reserve(block_length);
        //For each block,find the min element and normalize it and compute its naive rmq
        for(iterator_type block_begin=_level.begin();block_begin<_level.end();block_begin+=_block_size)
        {
            const iterator_type block_end=min(block_begin+_block_size,_level.end());
            //Find the min element of the block by brute force
            iterator_type block_min=min_element(block_begin,block_end);
            _super_array_vals.push_back(*block_min);
            _super_array_ids.push_back(block_min-_level.begin());
            //Compute the normalized block
            int init=*block_begin;
            vector<int> normalized_block(block_end-block_begin);
            transform(block_begin,block_end,normalized_block.begin(),[init](const int& val){return val-init;});
            //Find the normalized block int the map of RMQ structures,if not found, construct one
            unique_ptr<Naive_rmq>& naive_ptr=_sub_block_rmqs[normalized_block];
            if(!naive_ptr){
                naive_ptr.reset(new Naive_rmq(normalized_block.begin(),normalized_block.end()));
            }
            //Revord a pointer to the RMQ structure at this sub_block's index
            _sub_block_rmq_array.push_back(naive_ptr.get());
        }
        //cout<<"real block_length="<<_super_array_vals.size()<<" b-e="<<_super_array_vals.end()-_super_array_vals.begin()<<endl;//to be deleted
        //Construct the RMQ structure over the super array
        _super_rmq.reset(new Sparse_rmq(_super_array_vals.begin(),_super_array_vals.end()));
        qtime = GetCurrentTimeSec() - qtime;
        double avg_qtime = qtime/ (double)numOfVertices;
        cout << "rmq construct time:" << qtime * 1e6 <<  " microseconds" << endl;
        cout << "avg rmq construct time:" << avg_qtime * 1e6 <<  " microseconds" << endl;
    
    }

};

class H2H_Index{
public:
	h2h_index_t_p* h2h_index_p;
    h2h_bindex_t_p* h2h_bindex_p;//Add for directed
	int _numOfTreeNodes;
    LCA lca;//find the least common ancestor structure
    vector<double> _labelSize;
    vector<unsigned int> _queryTime;
    vector<bool> _HFPointFlag;
    vector<NodeID> _HFPoint;
    int _maxQueryTime{0};
    int _numOfHfpoint{0};
    long long _totalQueryTime{0};
	//**********constructions and deconstructions***********
	H2H_Index(){
		h2h_index_p=NULL;
        h2h_bindex_p=NULL;
		_numOfTreeNodes=0;
	}
    
	~H2H_Index(){
		// if(h2h_index_p!=NULL){
        //     for(int v=0;v<_numOfTreeNodes;++v){
        //         free(h2h_index_p[v].dis);
        //         free(h2h_index_p[v].pos);
        //     }
        //     free(h2h_index_p);
        // }
	}

    //**********functions****************
    // #pragma g++ push_options
    // #pragma optimize("O0")
    // int query_p(int s,int t){
    //     int lca_id=lca.lca_query(s,t);
    //     int distance = INF_WEIGHT;
    //     const h2h_index_t_p &idx_lca = h2h_index_p[lca_id];
    //     const h2h_index_t_p &idx_s = h2h_index_p[s];
    //     const h2h_index_t_p &idx_t = h2h_index_p[t];
    //     if(lca_id==s){
    //         return idx_t.dis[idx_lca.height];
    //     }else if(lca_id==t){
    //         return idx_s.dis[idx_lca.height];
    //     }else{
    //         int distance_t = INF_WEIGHT;
    //         _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
    //         _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
    //         _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
    //         for(int i=0;i<=idx_lca.width;++i){
    //             int pos_i=idx_lca.pos[i];
    //             distance_t=idx_s.dis[pos_i]+idx_t.dis[pos_i];
    //             if(distance_t<distance) distance=distance_t;
    //         }
    //         return distance;
    //     }
    //     return distance;
    // }
    // #pragma g++ pop_options

    int query_p(int s,int t){
        int lca_id=lca.lca_query(s,t);
        int distance = INF_WEIGHT;
        int distance_t = INF_WEIGHT;
        const h2h_index_t_p &idx_lca = h2h_index_p[lca_id];
        const h2h_index_t_p &idx_s = h2h_index_p[s];
        const h2h_index_t_p &idx_t = h2h_index_p[t];
        _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
        _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);

        for(int i=0;i<=idx_lca.width;++i){
            int pos_i=idx_lca.pos[i];
            distance_t=idx_s.dis[pos_i]+idx_t.dis[pos_i];
            if(distance_t<distance) distance=distance_t;
        }
       //cout<<s<<"-"<<t<<" lca="<<lca_id<<" dis="<<distance<<endl;//output test
        return distance;
    }

    int query_p_directed(int s,int t){
        int lca_id=lca.lca_query(s,t);
        int distance = INF_WEIGHT;
        int distance_t = INF_WEIGHT;
        const h2h_bindex_t_p &idx_lca = h2h_bindex_p[lca_id];
        const h2h_bindex_t_p &idx_s = h2h_bindex_p[s];
        const h2h_bindex_t_p &idx_t = h2h_bindex_p[t];
        _mm_prefetch(&idx_s.dis[0], _MM_HINT_T0);
        _mm_prefetch(&idx_t.dis[0], _MM_HINT_T0);
        _mm_prefetch(&idx_lca.pos[0], _MM_HINT_T0);
        for(int i=0;i<=idx_lca.width;++i){
            int pos_i=idx_lca.pos[i];
            distance_t=idx_s.dis[pos_i]+idx_t.r_dis[pos_i];
            if(distance_t<distance) distance=distance_t;
        }
       //cout<<s<<"-"<<t<<" lca="<<lca_id<<" dis="<<distance<<endl;//output test
        return distance;
    }

    void load_Index(char* load_index_filename,char* load_lca_filename){
        if(h2h_index_p){
            for(int v=0;v<_numOfTreeNodes;++v){
                free(h2h_index_p[v].dis);
                free(h2h_index_p[v].pos);
            }
            free(h2h_index_p);
        }
		h2h_index_p=NULL;
		
        ifstream ifs(load_index_filename);
        int isize=0,t=-1;
        ifs.read((char*)&isize, sizeof(isize));
        _numOfTreeNodes=isize;
        cout<<"_numOfTreeNodes="<<_numOfTreeNodes<<endl;
        h2h_index_p=(h2h_index_t_p*)memalign(64, _numOfTreeNodes * sizeof(h2h_index_t_p));
        for(int v=0;v<_numOfTreeNodes;++v)
        {
            h2h_index_t_p &idx=h2h_index_p[v];
            //pos
			ifs.read((char*)&isize, sizeof(isize));
            idx.width=isize;
			idx.pos = (int*)memalign(64, (isize+1) * sizeof(int));
            for(int i=0;i<=isize;++i){
                ifs.read((char*)&t, sizeof(t));
                idx.pos[i]=t;
            }
            //dis
            ifs.read((char*)&isize, sizeof(isize));
            idx.height=isize;
			idx.dis = (int*)memalign(64, (isize+1) * sizeof(int));
            for(int i=0;i<=isize;++i){
                ifs.read((char*)&t, sizeof(t));
                idx.dis[i]=t;
            }
        }
        //output test
        // for(int v=0;v<_numOfTreeNodes;++v){
        //     cout<<v<<":"<<endl;
        //     cout<<"pos:";
        //     for(int i=0;i<=h2h_index_p[v].width;++i) cout<<h2h_index_p[v].pos[i]<<" ";
        //     cout<<endl;
        //     cout<<"dis:";
        //     for(int i=0;i<=h2h_index_p[v].height;++i) cout<<h2h_index_p[v].dis[i]<<" ";
        //     cout<<endl;
        // }
        //build lca
        lca.initialize(load_lca_filename);
	}

    void load_Index(char* load_index_filename,char* load_lca_filename,double& time_rmq1){
        if(h2h_index_p){
            for(int v=0;v<_numOfTreeNodes;++v){
                free(h2h_index_p[v].dis);
                free(h2h_index_p[v].pos);
            }
            free(h2h_index_p);
        }
		h2h_index_p=NULL;
		
        ifstream ifs(load_index_filename);
        int isize=0,t=-1;
        ifs.read((char*)&isize, sizeof(isize));
        _numOfTreeNodes=isize;
        cout<<"_numOfTreeNodes="<<_numOfTreeNodes<<endl;
        h2h_index_p=(h2h_index_t_p*)memalign(64, _numOfTreeNodes * sizeof(h2h_index_t_p));
        for(int v=0;v<_numOfTreeNodes;++v)
        {
            h2h_index_t_p &idx=h2h_index_p[v];
            //pos
			ifs.read((char*)&isize, sizeof(isize));
            idx.width=isize;
			idx.pos = (int*)memalign(64, (isize+1) * sizeof(int));
            for(int i=0;i<=isize;++i){
                ifs.read((char*)&t, sizeof(t));
                idx.pos[i]=t;
            }
            //dis
            ifs.read((char*)&isize, sizeof(isize));
            idx.height=isize;
			idx.dis = (int*)memalign(64, (isize+1) * sizeof(int));
            for(int i=0;i<=isize;++i){
                ifs.read((char*)&t, sizeof(t));
                idx.dis[i]=t;
            }
        }
        //build lca
        time_rmq1=GetCurrentTimeSec();
        lca.initialize(load_lca_filename);
        time_rmq1=GetCurrentTimeSec()-time_rmq1;
	}

    void load_Index_directed(char* load_index_filename,char* load_lca_filename,double& time_rmq1){
        if(h2h_bindex_p){
            for(int v=0;v<_numOfTreeNodes;++v){
                free(h2h_bindex_p[v].dis);
                free(h2h_bindex_p[v].r_dis);
                free(h2h_bindex_p[v].pos);
            }
            free(h2h_bindex_p);
        }
		h2h_bindex_p=NULL;
		
        ifstream ifs(load_index_filename);
        int isize=0,t=-1;
        ifs.read((char*)&isize, sizeof(isize));
        _numOfTreeNodes=isize;
        cout<<"_numOfTreeNodes="<<_numOfTreeNodes<<endl;
        h2h_bindex_p=(h2h_bindex_t_p*)memalign(64, _numOfTreeNodes * sizeof(h2h_bindex_t_p));
        for(int v=0;v<_numOfTreeNodes;++v)
        {
            h2h_bindex_t_p &idx=h2h_bindex_p[v];
            //pos
			ifs.read((char*)&isize, sizeof(isize));
            idx.width=isize;
			idx.pos = (int*)memalign(64, (isize+1) * sizeof(int));
            for(int i=0;i<=isize;++i){
                ifs.read((char*)&t, sizeof(t));
                idx.pos[i]=t;
            }
            //rdis
            ifs.read((char*)&isize, sizeof(isize));
            idx.height=isize;
			idx.r_dis = (int*)memalign(64, (isize+1) * sizeof(int));
            for(int i=0;i<=isize;++i){
                ifs.read((char*)&t, sizeof(t));
                idx.r_dis[i]=t;
            }
            //dis
			idx.dis = (int*)memalign(64, (isize+1) * sizeof(int));
            for(int i=0;i<=isize;++i){
                ifs.read((char*)&t, sizeof(t));
                idx.dis[i]=t;
            }
        }
        //build lca
        time_rmq1=GetCurrentTimeSec();
        lca.initialize(load_lca_filename);
        time_rmq1=GetCurrentTimeSec()-time_rmq1;
	}

    void load_node_label_size(char* load_filename){
        _labelSize.resize(numOfVertices,0);
        ifstream in(load_filename);//input HFPoint file to ifstream
        if(!in.is_open()) {cerr<<"Cannot open "<<load_filename<<endl;}
        NodeID v; //each line representing the label size of each vertice
        double label_size;
        //read each line representing HFpoint to vector 
        int i;
        for (i = 0; in >> v >> label_size;i++) {
            //cout<<v<<":"<<label_size<<endl;//for debug
			_labelSize[v]=label_size;
		}
        if(i!=numOfVertices) cout<<"Load label size Error:i!=numOfVertices!"<<endl;
        in.close();
    }

	void analysis_and_append_result(char* analy_file_name,char* pointFreq_file_name,int hfRate=0){
        //load hfpoint
        _HFPoint.reserve(_numOfTreeNodes);
        _HFPointFlag.resize(_numOfTreeNodes,0);
        _queryTime.resize(_numOfTreeNodes,0);
        if(hfRate==0) _numOfHfpoint=_numOfTreeNodes;
        else _numOfHfpoint=static_cast<int>((double)_numOfTreeNodes*hfRate/(double)HF_DIVIDION);
		load_hfpoint_and_qt(pointFreq_file_name,_numOfHfpoint,_HFPoint,_HFPointFlag,_queryTime,_maxQueryTime,_totalQueryTime);
        //variables 
		double total_sum_size=0,total_ave_size=0;
		double hf_sum_size=0.0,hf_ave_size=0.0;
		double total_performance=0,total_performance_s=0,ratio;
		for(int v=0;v<_numOfTreeNodes;++v){
			double label_size=_labelSize[v];
			total_sum_size+=label_size;
			if(_HFPointFlag[v]) hf_sum_size+=label_size;
			ratio=(double)_queryTime[v]/((double)(_maxQueryTime*DIVISION_FACTOR));
			total_performance+=ratio*label_size;
		}
		total_ave_size= (double) total_sum_size/(double) _numOfTreeNodes;
		hf_ave_size= hf_sum_size/(double) _numOfHfpoint;
		total_performance_s=total_performance*((double)_maxQueryTime);
		//output result
		cout<<"_numOfTreeNodes = "<<_numOfTreeNodes<<" total_ave_size = "<<total_ave_size<<endl;
		cout<<"_numOfHFpoint = "<<_numOfHfpoint<<" hf_sum_size = "<<hf_sum_size<<" hf_ave_size = "<<hf_ave_size<<endl;
		cout<<"nomalization performance_result = "<<total_performance<<" standard performance_result = "<<total_performance_s<<endl;
		//append to file
		ofstream ofs(analy_file_name,ios::app|ios::out);//append way
		if(!ofs.is_open()) cout<<"Cannot open "<<analy_file_name<<endl;
		ofs.setf(ios::fixed);
		ofs.precision(4);
		ofs<<total_sum_size<<" "<<total_ave_size<<" "<<hf_sum_size<<" "<<hf_ave_size<<" "<<total_performance_s<<" ";
		ofs.close();
	}
};

#endif // !_H2H_INDEX_H