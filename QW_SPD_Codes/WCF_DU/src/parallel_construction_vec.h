/*
 * @Description: parallel pll vec model
 * @Author: wanjingyi
 * @Date: 2021-03-02 14:32:32
 * @LastEditTime: 2021-09-12 10:36:24
 */
#ifndef  PARALLEL_CONSTRUCTION_VEC_H
#define PARALLEL_CONSTRUCTION_VEC_H

#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <iostream>
#include <climits>
#include <xmmintrin.h>
#include <immintrin.h> 
#include <bitset>
#include <cmath>
#include "graph.h"
#include "omp.h"
#include "./time_util.h"
#include "./utils.h"

using std::vector;
using std::unordered_map;
using std::map;
using std::bitset;
using std::stable_sort;
using std::min;
using std::fill;
using namespace time_util;

namespace PADO_VEC{
    //**************************global typedefs*********************
    //typedef uint64_t idi; // unsinged long long
    typedef uint32_t idi; // unsigned int
    //typedef int weighti;
    typedef uint8_t weighti;
    //typedef int32_t weightiLarge;
    typedef int16_t weightiLarge;
    typedef uint8_t smalli;
    typedef uint32_t inti;
    const uint8_t WEIGHTI_MAX = UCHAR_MAX;
    const int16_t WEIGHTILARGE_MAX = SHRT_MAX;
    const uint8_t SMALLI_MAX = UCHAR_MAX;
    //**************************global typedefs*********************

    //**************************avx512_pado.h*********************
    //512位的AVX寄存器，如果整型Int是32位，那么一个寄存器能装16个整数，即同时操作16个整数
    const inti NUM_P_INT = 16;
    const __m512i INF_LARGE_v = _mm512_set1_epi32((int)WEIGHTILARGE_MAX);
    const __m512i LOW_TWO_BYTE_MASK = _mm512_set1_epi32((int)0xFFFF);
    const __m512i UNDEF_i32_v = _mm512_undefined_epi32();
    const __m256i INF_v_256i = _mm256_set1_epi16((int)-1);
    const __m512i INF_v = _mm512_set1_epi32((int)WEIGHTI_MAX);
    //**************************avx512_pado.h*********************

    // Structure for the type of label
    struct IndexType {
        vector<idi> vertices; // Vertices in the label, preresented as temperory ID
        vector<weightiLarge> distances;
        inti size() {
            return vertices.size();
        }
    }; //__attribute__((aligned(64)));

	// Structure for the type of temporary label
	struct ShortIndex {
		// Use a queue to store candidates
		vector<inti> candidates_que;
		inti end_candidates_que = 0;
		// Use a array to store distances of candidates; length of roots_size
		vector<weightiLarge> candidates_dists; // record the distances to candidates. If candidates_dists[c] = INF, then c is NOT a candidate
		// The candidates_dists is also used for distance query.

		// Use a queue to store temporary labels in this batch; use it so don't need to traverse vertices_dists.
		vector<inti> vertices_que; // Elements in vertices_que are roots_id, not real id
		idi end_vertices_que = 0;
		// Use an array to store distances to roots; length of roots_size
		vector<weightiLarge> vertices_dists; // labels_table

		// Usa a queue to record which roots have reached this vertex in this batch.
		// It is used for reset the dists_table
		vector<inti> reached_roots_que;
		idi end_reached_roots_que = 0;

		// Use a queue to store last inserted labels (IDs); distances are stored in vertices_dists.
		vector<inti> last_new_roots;
		idi end_last_new_roots = 0;

		ShortIndex() = default;
		explicit ShortIndex(idi num_roots) {
			candidates_que.resize(num_roots);
			candidates_dists.resize(num_roots, WEIGHTILARGE_MAX);
			last_new_roots.resize(num_roots);
			vertices_que.resize(num_roots);
			vertices_dists.resize(num_roots, WEIGHTILARGE_MAX);
			reached_roots_que.resize(num_roots);
		}
	}; //__attribute__((aligned(64)));

    // Structure of the public ordered index for distance queries.
    struct index_t {
        vector<NodeID> spt_v;
        vector<EdgeWeight> spt_d;
        NodeID size() {
            return spt_v.size();
        }
    };

    class WeightedVCPLLVEC{
    public:
        //***************variables*************
        vector<IndexType> L;
        vector<index_t> index_; // Ordered labels for original vertex ID
        inti BATCH_SIZE;
        //***************functions*************
        void reorder_labels(const vector<NodeID> &rank){
            int labels_count = 0;
            index_.resize(numOfVertices);
            for (NodeID v_id = 0; v_id < numOfVertices; ++v_id) {
                NodeID v_rank = rank[v_id];
                const IndexType &Lv = L[v_rank];
                NodeID size_labels = Lv.vertices.size();
                labels_count += size_labels;
                vector< std::pair<NodeID, EdgeWeight> > ordered_labels(size_labels);
                // Traverse v_id's all labels
                for (NodeID l_i = 0; l_i < size_labels; ++l_i) {
                    ordered_labels[l_i].first = Lv.vertices[l_i];
                    ordered_labels[l_i].second = Lv.distances[l_i];
                }
                // Sort
                sort(ordered_labels.begin(), ordered_labels.end());
                index_[v_id].spt_v.resize(size_labels+1);
                index_[v_id].spt_d.resize(size_labels+1);
                for (NodeID l_i = 0; l_i < size_labels; ++l_i) {
                    index_[v_id].spt_v[l_i]=ordered_labels[l_i].first;
                    index_[v_id].spt_d[l_i] = ordered_labels[l_i].second;
                }
                index_[v_id].spt_v[size_labels]=numOfVertices;
                index_[v_id].spt_d[size_labels]=INF_WEIGHT;
            }
            cout<<"Total label size="<<labels_count<<" Average label size="<<static_cast<double>(labels_count)/numOfVertices<<endl;
        }

        void save_labels(const char *save_filename)
        {
            string label_filename(save_filename);
            label_filename.append(".label");
            ofstream ofs(label_filename, ios::binary | ios::out);
            if(!ofs.is_open()){
                cout<<"Cannot open "<<label_filename<<endl;
                exit(EXIT_FAILURE);
            }
            // Store into file the number of vertices and the number of bit-parallel roots.
            ofs.write((const char*)&numOfVertices, sizeof(numOfVertices));
            for (NodeID v = 0; v < numOfVertices; ++v) {
                NodeID isize = index_[v].size();
                ofs.write((const char*)&isize, sizeof(isize));
                for (NodeID i = 0; i < index_[v].size(); ++i) {
                    ofs.write((const char*)&index_[v].spt_v[i], sizeof(index_[v].spt_v[i]));
                    ofs.write((const char*)&index_[v].spt_d[i], sizeof(index_[v].spt_d[i]));
                }
            }
            ofs.close();
        }

        void save_label_size(char* label_size_file,const vector<NodeID>& inv) {
            string labelSizefile(label_size_file);
            labelSizefile.append(".size");
            ofstream ofs(labelSizefile);
            if(!ofs.is_open()) {cerr<<"Cannot open "<<labelSizefile<<endl;}
            for (int i = 0; i < numOfVertices; ++i) {
                ofs << index_[i].size()-1 << endl;
            }
            ofs.close();
            //output the size by descending order,format:label_size
            string labelSizefile_prefix1(label_size_file);
            string labelSizefile1=labelSizefile_prefix1.append("_byOrder.size");
            ofstream out(labelSizefile1.c_str());
            if(!out.is_open()) {cerr<<"Cannot open "<<labelSizefile1<<endl;}
            for (int i = 0; i < numOfVertices; ++i) {
                out<<index_[inv[i]].size()-1 << endl;
            }		
            out.close();
        }

        void write_labels(char* write_filename,const vector<NodeID>& inv,bool isOrder=false)
        {
                string write_filename_prefix(write_filename);
                string write_filename1=write_filename_prefix.append(".list");
                ofstream ofs(write_filename1.c_str());
            for (NodeID v = 0; v < numOfVertices; ++v) 
            {
                NodeID isize = index_[v].size();
                ofs <<isize<<" "<<v;
                for (NodeID i = 0; i < index_[v].size(); ++i) {
                    ofs<<" "<<'('<<index_[v].spt_v[i]<<","<<index_[v].spt_d[i]<<")";
                }
                ofs<<endl;
            }
            ofs.close();

            //write labels with original NodeId in graph
            string write_filename2=write_filename_prefix.append("_original");
            ofstream out(write_filename2.c_str());
            for (NodeID v = 0; v < numOfVertices; ++v) 
            {
                NodeID isize = index_[v].size();
                out<<isize-1<<" "<<v;
                for (NodeID i = 0; i < index_[v].size()-1; ++i) {
                    out<<" "<<'('<<inv[index_[v].spt_v[i]]<<","<<index_[v].spt_d[i]<<")";
                }
                out<<endl;
            }			
            out.close();

            if(isOrder){
                string write_filename_prefix1(write_filename);
                string write_filename3=write_filename_prefix1.append(".list_order");
                ofstream out1(write_filename3.c_str());
                for (NodeID r= 0; r < numOfVertices; ++r) 
                {
                    NodeID v=inv[r];
                    NodeID isize = index_[v].size();
                    out1<<isize-1<<" "<<v;
                    for (NodeID i = 0; i < index_[v].size()-1; ++i) {
                        out1<<" "<<'('<<inv[index_[v].spt_v[i]]<<","<<index_[v].spt_d[i]<<")";
                    }
                    out1<<endl;
                }			
                out1.close();
                }
        }
    };

    class WeightedParaVertexCentricPLLVec :public WeightedVCPLLVEC{
    public:
        //****************variables**************
        int num_of_threads{5};
        //***************************constructions**************************
        WeightedParaVertexCentricPLLVec(){}
        
        WeightedParaVertexCentricPLLVec(WGraph& wgraph,int batchSize,int numThreads){
            num_of_threads=numThreads;
            BATCH_SIZE=batchSize;
            construct(wgraph);
        }

        ~WeightedParaVertexCentricPLLVec(){}
    
    protected:
        //*********************functions*********************
        // Function: return the position of the least significant bit that is set in a 32-bit integer
        // This magic code is from https://stackoverflow.com/a/757266/7252729
        inline inti get_position_lowest_set_1(uint32_t v) {
            static const inti MultiplyDeBruijnBitPosition[32] = 
            {
                0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
                31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
            };
            return MultiplyDeBruijnBitPosition[((uint32_t)((v & -v) * 0x077CB531U)) >> 27];
        }

        /**
         * @Author: wanjingyi
         * @description: 
         * @param {*}
         * @return {*}
         */
        template <typename T, typename Int>
        inline void TS_enqueue(
            vector<T> &queue,
            Int &end_queue,
            const T &e)
        {
            volatile Int old_i = end_queue;
            volatile Int new_i = old_i + 1;
            while (!CAS(&end_queue, old_i, new_i)) {
                old_i = end_queue;
                new_i = old_i + 1;
            }
            queue[old_i] = e;
        }

        /**
         * @Author: wanjingyi
         * @description: after finishing the label tables in the short_index, build the index according to it.And also reset the has_new_labels_queue
         * @param {*}
         * @return {*}
         */        
        inline void update_index(
            vector<IndexType> &L,
            vector<ShortIndex> &short_index,
            idi roots_start,
            vector<idi> &has_new_labels_queue,
            idi &end_has_new_labels_queue,
            vector<uint8_t> &has_new_labels)
        {
            // Traverse has_new_labels_queue for every vertex who got new labels in this batch
            #pragma omp parallel for
            for (idi i_q = 0; i_q < end_has_new_labels_queue; ++i_q) {
                idi v_id = has_new_labels_queue[i_q];
                has_new_labels[v_id] = 0; // Reset has_new_labels
                IndexType &Lv = L[v_id];
                ShortIndex &SI_v = short_index[v_id];
                inti bound_i_r = SI_v.end_vertices_que;
                for (inti i_r = 0; i_r < bound_i_r; ++i_r) {
                    idi r_root_id = SI_v.vertices_que[i_r];
                    weightiLarge dist = SI_v.vertices_dists[r_root_id];
                    if (WEIGHTILARGE_MAX == dist) {
                        continue;
                    }
                    SI_v.vertices_dists[r_root_id] = WEIGHTILARGE_MAX; // Reset v_id's vertices_dists
                    Lv.vertices.push_back(r_root_id + roots_start);
                    Lv.distances.push_back(dist);
                }
                SI_v.end_vertices_que = 0; // Clear v_id's vertices_que
            }
            end_has_new_labels_queue = 0; // Clear has_new_labels_queue
        }

        /**
         * @Author: wanjingyi
         * @description: reset distance table dists_table
         * @param {*}
         * @return {*}
         */        
        inline void reset_tables(
            vector<ShortIndex> &short_index,
            idi roots_start,
            inti roots_size,
            const vector<IndexType> &L,
            vector< vector<weightiLarge> > &dists_table,
            vector< vector<weightiLarge> > &roots_labels_buffer,
            vector<idi> &once_had_cand_queue,
            idi &end_once_had_cand_queue,
            vector<uint8_t> &once_had_cand)
        {
            // Reset roots_labels_buffer according to L (old labels)
            #pragma omp parallel for
            for (idi r_roots_id = 0; r_roots_id < roots_size; ++r_roots_id) {
                idi r_real_id = r_roots_id + roots_start;
                // Traverse labels of r
                const IndexType &Lr = L[r_real_id];
                idi bound_i_l = Lr.vertices.size();
                for (idi i_l = 0; i_l < bound_i_l; ++i_l) {
                    roots_labels_buffer[r_roots_id][Lr.vertices[i_l]] = WEIGHTILARGE_MAX;
                }
                roots_labels_buffer[r_roots_id][r_real_id] = WEIGHTILARGE_MAX;
            }

	        // Reset dists_table according to short_index[v].reached_roots_que
            #pragma omp parallel for
            for (idi i_q = 0; i_q < end_once_had_cand_queue; ++i_q) {
                idi v_id = once_had_cand_queue[i_q];
                once_had_cand[v_id] = false; // Reset once_had_cand
                ShortIndex &SI_v = short_index[v_id];
                // Traverse roots which have reached v_id
                inti bound_i_r = SI_v.end_reached_roots_que;
                for (inti i_r = 0; i_r < bound_i_r; ++i_r) {
                    dists_table[v_id][SI_v.reached_roots_que[i_r]] = WEIGHTILARGE_MAX;
                }
                SI_v.end_reached_roots_que = 0; // Clear v_id's reached_roots_que
            }
            end_once_had_cand_queue = 0; // Clear once_had_cand_queue.
        }

        /**
         * @Author: wanjingyi
         * @description: at the end of every batch, filter out wrong and redundant labels
         * @param {*}
         * @return {*}
         */        
        inline void filter_out_labels(
		const vector<idi> &has_new_labels_queue,
		idi end_has_new_labels_queue,
		vector<ShortIndex> &short_index,
		idi roots_start)
        {
            // Traverse has_new_labels_queue for every vertex who got new labels in this batch
            #pragma omp parallel for
            for (idi i_q = 0; i_q < end_has_new_labels_queue; ++i_q) {
                idi v_id = has_new_labels_queue[i_q];
                IndexType &Lv = L[v_id];
                ShortIndex &SI_v = short_index[v_id];
                inti bound_label_i = SI_v.end_vertices_que;
                // Traverse v_id's every new label
                for (inti i_c = 0; i_c < bound_label_i; ++i_c) {
                    idi c_root_id = SI_v.vertices_que[i_c];
                    weightiLarge dist_v_c = SI_v.vertices_dists[c_root_id];
                    ShortIndex &SI_c = short_index[c_root_id + roots_start];

                    // Vectorized Version
                    const __m512i c_root_id_v = _mm512_set1_epi32(c_root_id);
                    const __m512i dist_v_c_v = _mm512_set1_epi32(dist_v_c);
                    inti remainder_simd = bound_label_i % NUM_P_INT;
                    idi bound_i_r = bound_label_i - remainder_simd;
                    for (idi i_r = 0; i_r < bound_i_r; i_r += NUM_P_INT) {
                        // Other label r
                        __m512i r_root_id_v = _mm512_load_epi32(&SI_v.vertices_que[i_r]);
                        __mmask16 is_r_higher_ranked_m = _mm512_cmplt_epi32_mask(r_root_id_v, c_root_id_v);
                        if (!is_r_higher_ranked_m) {
                            continue;
                        }
                        // Distance v to r
                        __m512i dists_v_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_r_higher_ranked_m, r_root_id_v, &SI_v.vertices_dists[0], sizeof(weightiLarge));
                        dists_v_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_r_higher_ranked_m, dists_v_r_v, LOW_TWO_BYTE_MASK);
                        __mmask16 is_shorter_m = _mm512_cmple_epi32_mask(dists_v_r_v, dist_v_c_v);
                        if (!is_shorter_m) {
                            continue;
                        }
                        // Distance c to r
                        __m512i dists_c_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_shorter_m, r_root_id_v, &SI_c.vertices_dists[0], sizeof(weightiLarge));
                        dists_c_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_shorter_m, dists_c_r_v, LOW_TWO_BYTE_MASK);
                        is_shorter_m = _mm512_cmple_epi32_mask(dists_c_r_v, dist_v_c_v);
                        if (!is_shorter_m) {
                            continue;
                        }
                        // Distance compare
                        is_shorter_m = _mm512_mask_cmple_epi32_mask(is_shorter_m, _mm512_mask_add_epi32(INF_LARGE_v, is_shorter_m, dists_v_r_v, dists_c_r_v), dist_v_c_v);
                        if (is_shorter_m) {
                            int_mask_i32scatter_epiI(
                                    &SI_v.vertices_dists[0],
                                    is_shorter_m,
                                    c_root_id_v,
                                    WEIGHTILARGE_MAX);
                        }
                    }
                    if (remainder_simd) {
                        __mmask16 in_m = (__mmask16) ((uint16_t) 0xFFFF >> (NUM_P_INT - remainder_simd));
                        // Other label r
                        __m512i r_root_id_v = _mm512_mask_loadu_epi32(UNDEF_i32_v, in_m, &SI_v.vertices_que[bound_i_r]);
                        __mmask16 is_r_higher_ranked_m = _mm512_mask_cmplt_epi32_mask(in_m, r_root_id_v, c_root_id_v);
                        if (!is_r_higher_ranked_m) {
                            continue;
                        }
                        // Distance v to r
                        __m512i dists_v_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_r_higher_ranked_m, r_root_id_v, &SI_v.vertices_dists[0], sizeof(weightiLarge));
                        dists_v_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_r_higher_ranked_m, dists_v_r_v, LOW_TWO_BYTE_MASK);
                        __mmask16 is_shorter_m = _mm512_cmple_epi32_mask(dists_v_r_v, dist_v_c_v);
                        if (!is_shorter_m) {
                            continue;
                        }
                        // Distance c to r
                        __m512i dists_c_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_shorter_m, r_root_id_v, &SI_c.vertices_dists[0], sizeof(weightiLarge));
                        dists_c_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_shorter_m, dists_c_r_v, LOW_TWO_BYTE_MASK);
                        is_shorter_m = _mm512_cmple_epi32_mask(dists_c_r_v, dist_v_c_v);
                        if (!is_shorter_m) {
                            continue;
                        }
                        // Distance compare
                        is_shorter_m = _mm512_mask_cmple_epi32_mask(is_shorter_m, _mm512_mask_add_epi32(INF_v, is_shorter_m, dists_v_r_v, dists_c_r_v), dist_v_c_v);
                        if (is_shorter_m) {
                            int_mask_i32scatter_epiI(
                                    &SI_v.vertices_dists[0],
                                    is_shorter_m,
                                    c_root_id_v,
                                    WEIGHTILARGE_MAX);
                        }
                    }
                }
            }
        }

        /**
         * @Author: wanjingyi
         * @description: assign lanes of vindex in base_addr as int a.
         * @param {*}
         * @return {*}
         */        
        template <typename I> 
        inline void int_i32scatter_epiI(I *base_addr,__m512i vindex,I a)
        {
            int *tmp_vindex = (int *) &vindex;
            for (inti i = 0; i < 16; ++i) {
                base_addr[tmp_vindex[i]] = a;
            }
        }

        /**
         * @Author: wanjingyi
         * @description: according to mask k, assign lanes of vindex in base_addr as a
         * @param {*}
         * @return {*}
         */
        template <typename I>
        inline void int_mask_i32scatter_epiI(
                I *base_addr,
                __mmask16 k,
                __m512i vindex,
                I a)
        {
            int *tmp_vindex = (int *) &vindex;
            for (inti i = 0; i < 16; ++i) {
                uint16_t bit = (uint16_t) pow(2.0, i);
                if (k & bit) {
                    base_addr[tmp_vindex[i]] = a;
                }
            }
        }

        /**
         * @Author: wanjingyi
         * @description: function used for check
         * @param {*}
         * @return {*}
         */        
        inline weightiLarge distance_query(
            idi v_id,
            idi cand_root_id,
            const vector< vector<weightiLarge> > &roots_labels_buffer,
            const vector<ShortIndex> &short_index,
            const vector<IndexType> &L,
            idi roots_start,
            weightiLarge tmp_dist_v_c)
        {
            idi cand_real_id = cand_root_id + roots_start;
	        const IndexType &Lv = L[v_id];
            __m512i cand_real_id_v = _mm512_set1_epi32(cand_real_id);
            __m512i tmp_dist_v_c_v = _mm512_set1_epi32(tmp_dist_v_c);
            inti remainder_simd = Lv.vertices.size() % NUM_P_INT;
            idi bound_i_l = Lv.vertices.size() - remainder_simd;
            for (idi i_l = 0; i_l < bound_i_l; i_l += NUM_P_INT) {
                // Labels IDs
                __m512i r_v = _mm512_load_epi32(&Lv.vertices[i_l]);
                __mmask16 is_r_higher_ranked_m = _mm512_cmplt_epi32_mask(r_v, cand_real_id_v);
                if (!is_r_higher_ranked_m) {
                    continue;
                }
                // Labels dists
                __m512i dists_c_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_r_higher_ranked_m, r_v, &roots_labels_buffer[cand_root_id][0], sizeof(weightiLarge));
                dists_c_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_r_higher_ranked_m, dists_c_r_v, LOW_TWO_BYTE_MASK); // The least significant 2 byte is the weight.
                __mmask16 is_not_INF_m = _mm512_cmpneq_epi32_mask(dists_c_r_v, INF_LARGE_v);
                if (!is_not_INF_m) {
                    continue;
                }

                // Dist from v to c through label r
                __m256i tmp_dists_v_r_256i = _mm256_mask_loadu_epi16(INF_v_256i, is_not_INF_m, &Lv.distances[i_l]); // Weighti is 16-bit
                __m512i dists_v_r_v = _mm512_mask_cvtepi16_epi32(INF_LARGE_v, is_not_INF_m, tmp_dists_v_r_256i); // Convert to 32-bit
                __m512i label_dist_v_c_v = _mm512_mask_add_epi32(INF_LARGE_v, is_not_INF_m, dists_c_r_v, dists_v_r_v);
                __mmask16 is_label_dist_shorter_m = _mm512_mask_cmple_epi32_mask(is_not_INF_m, label_dist_v_c_v, tmp_dist_v_c_v);
                if (is_label_dist_shorter_m) {
                    // Need to return the shorter distance (might be equal)
                    inti index = get_position_lowest_set_1((uint32_t) is_label_dist_shorter_m);
                    return ((int *) &label_dist_v_c_v)[index]; // Return the distance
                }
            }
            if (remainder_simd) {
                __mmask16 in_m = (__mmask16) ((uint16_t) 0xFFFF >> (NUM_P_INT - remainder_simd));
                // Labels IDs
                __m512i r_v = _mm512_mask_loadu_epi32(UNDEF_i32_v, in_m, &Lv.vertices[bound_i_l]);
                __mmask16 is_r_higher_ranked_m = _mm512_mask_cmplt_epi32_mask(in_m, r_v, cand_real_id_v);
                if (is_r_higher_ranked_m) {
                    // Labels dists
                    __m512i dists_c_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_r_higher_ranked_m, r_v, &roots_labels_buffer[cand_root_id][0], sizeof(weightiLarge));
                    dists_c_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_r_higher_ranked_m, dists_c_r_v, LOW_TWO_BYTE_MASK); // The least significant 2 byte is the weight.
                    __mmask16 is_not_INF_m = _mm512_cmpneq_epi32_mask(dists_c_r_v, INF_LARGE_v);
                    if (is_not_INF_m) {
                        // Dist from v to c through label r
                        __m256i tmp_dists_v_r_256i = _mm256_mask_loadu_epi16(INF_v_256i, is_not_INF_m, &Lv.distances[bound_i_l]); // Weighti is 16-bit
                        __m512i dists_v_r_v = _mm512_mask_cvtepi16_epi32(INF_LARGE_v, is_not_INF_m, tmp_dists_v_r_256i); // Convert to 32-bit
                        __m512i label_dist_v_c_v = _mm512_mask_add_epi32(INF_LARGE_v, is_not_INF_m, dists_c_r_v, dists_v_r_v);
                        __mmask16 is_label_dist_shorter_m = _mm512_mask_cmple_epi32_mask(is_not_INF_m, label_dist_v_c_v, tmp_dist_v_c_v);
                        if (is_label_dist_shorter_m) {
                            // Need to return the shorter distance (might be equal)
                            inti index = get_position_lowest_set_1((uint32_t) is_label_dist_shorter_m);
                            return ((int *) &label_dist_v_c_v)[index]; // Return the distance
                        }
                    }
                }
            }

            // 2. Labels in short_index[v_id].vertices_que
            const ShortIndex &SI_v = short_index[v_id];
            const ShortIndex &SI_c = short_index[cand_root_id];

            const __m512i roots_start_v = _mm512_set1_epi32(roots_start);
            remainder_simd = SI_v.end_vertices_que % NUM_P_INT;
            idi bound_i_que = SI_v.end_vertices_que - remainder_simd;
            for (inti i_que = 0; i_que < bound_i_que; i_que += NUM_P_INT) {
                __m512i r_root_id_v = _mm512_load_epi32(&SI_v.vertices_que[i_que]);
                __m512i r_real_id_v = _mm512_add_epi32(r_root_id_v, roots_start_v);
                __mmask16 is_r_higher_ranked_m = _mm512_cmplt_epi32_mask(r_real_id_v, cand_real_id_v);
                if (!is_r_higher_ranked_m) {
                    continue;
                }
                // Dists from v to r
                __m512i dists_v_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_r_higher_ranked_m, r_root_id_v, &SI_v.vertices_dists[0], sizeof(weightiLarge));
                dists_v_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_r_higher_ranked_m, dists_v_r_v, LOW_TWO_BYTE_MASK);
                // Check r_root_id in cand_root_id's labels inserted in this batch
                {
                    // Dists from c to r
                    __m512i dists_c_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_r_higher_ranked_m, r_root_id_v, &SI_c.vertices_dists[0], sizeof(weightiLarge));
                    dists_c_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_r_higher_ranked_m, dists_c_r_v, LOW_TWO_BYTE_MASK);
                    __mmask16 is_not_INF_m = _mm512_mask_cmpneq_epi32_mask(is_r_higher_ranked_m, dists_c_r_v, INF_LARGE_v);
                    if (is_not_INF_m) {
                        // Label dist from v to c
                        __m512i label_dist_v_c_v = _mm512_mask_add_epi32(INF_LARGE_v, is_not_INF_m, dists_v_r_v, dists_c_r_v);
                        // Compare label dist to tmp dist
                        __mmask16 is_label_dist_shorter_m = _mm512_mask_cmple_epi32_mask(is_not_INF_m, label_dist_v_c_v, tmp_dist_v_c_v);
                        if (is_label_dist_shorter_m) {
                            // Need to return the shorter distance (might be equal)
                            inti index = get_position_lowest_set_1((uint32_t) is_label_dist_shorter_m);
                            return ((int *) &label_dist_v_c_v)[index]; // Return the distance
                        }
                    }
                }
            }
            if (remainder_simd) {
                __mmask16 in_m = (__mmask16) ((uint16_t) 0xFFFF >> (NUM_P_INT - remainder_simd));
                __m512i r_root_id_v = _mm512_mask_loadu_epi32(UNDEF_i32_v, in_m, &SI_v.vertices_que[bound_i_que]);
                __m512i r_real_id_v = _mm512_mask_add_epi32(UNDEF_i32_v, in_m, r_root_id_v, roots_start_v);
                __mmask16 is_r_higher_ranked_m = _mm512_mask_cmplt_epi32_mask(in_m, r_real_id_v, cand_real_id_v);
                if (is_r_higher_ranked_m) {
                    // Dists from v to r
                    __m512i dists_v_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_r_higher_ranked_m, r_root_id_v, &SI_v.vertices_dists[0], sizeof(weightiLarge));
                    dists_v_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_r_higher_ranked_m, dists_v_r_v, LOW_TWO_BYTE_MASK);
                    // Check r_root_id in cand_root_id's labels inserted in this batch
                    {
                        // Dists from c to r
                        __m512i dists_c_r_v = _mm512_mask_i32gather_epi32(INF_LARGE_v, is_r_higher_ranked_m, r_root_id_v, &SI_c.vertices_dists[0], sizeof(weightiLarge));
                        dists_c_r_v = _mm512_mask_and_epi32(INF_LARGE_v, is_r_higher_ranked_m, dists_c_r_v, LOW_TWO_BYTE_MASK);
                        __mmask16 is_not_INF_m = _mm512_mask_cmpneq_epi32_mask(is_r_higher_ranked_m, dists_c_r_v, INF_LARGE_v);
                        if (is_not_INF_m) {
                            // Label dist from v to c
                            __m512i label_dist_v_c_v = _mm512_mask_add_epi32(INF_LARGE_v, is_not_INF_m, dists_v_r_v, dists_c_r_v);
                            // Compare label dist to tmp dist
                            __mmask16 is_label_dist_shorter_m = _mm512_mask_cmple_epi32_mask(is_not_INF_m, label_dist_v_c_v, tmp_dist_v_c_v);
                            if (is_label_dist_shorter_m) {
                                // Need to return the shorter distance (might be equal)
                                inti index = get_position_lowest_set_1((uint32_t) is_label_dist_shorter_m);
                                return ((int *) &label_dist_v_c_v)[index]; // Return the distance
                            }
                        }
                    }
                }
            }
            return WEIGHTILARGE_MAX;
        }

        /**
         * @Author: wanjingyi
         * @description: Collect elements in the tmp_queue into the queue
         * @param {*}
         * @return {*}
         */        
        template <typename T>
        inline void collect_into_queue(
            vector<T> &tmp_queue,
            vector<idi> &offsets_tmp_queue, // the locations in tmp_queue for writing from tmp_queue
            vector<idi> &offsets_queue, // the locations in queue for writing into queue.
            idi num_elements, // total number of elements which need to be added from tmp_queue to queue
            vector<T> &queue,
            idi &end_queue)
        {
            if (0 == num_elements) {
                return;
            }
            idi i_bound = offsets_tmp_queue.size();
            #pragma omp parallel for
            for (idi i = 0; i < i_bound; ++i) {
                idi i_q_start = end_queue + offsets_queue[i];
                idi i_q_bound;
                if (i_bound - 1 != i) {
                    i_q_bound = end_queue + offsets_queue[i + 1];
                } else {
                    i_q_bound = end_queue + num_elements;
                }
                if (i_q_start == i_q_bound) {
                    // If the group has no elements to be added, then continue to the next group
                    continue;
                }
                idi end_tmp = offsets_tmp_queue[i];
                for (idi i_q = i_q_start; i_q < i_q_bound; ++i_q) {
                    queue[i_q] = tmp_queue[end_tmp++];
                }
            }
            end_queue += num_elements;
        }

        /**
         * @Author: wanjingyi
         * @description: Function to get the prefix sum of elements in offsets
         * @param {*}
         * @return {*}
         */        
        inline idi prefix_sum_for_offsets(vector<idi> &offsets)
        {
            idi size_offsets = offsets.size();
            if (1 == size_offsets) {
                idi tmp = offsets[0];
                offsets[0] = 0;
                return tmp;
            } else if (size_offsets < 2048) {
                idi offset_sum = 0;
                idi size = size_offsets;
                for (idi i = 0; i < size; ++i) {
                    idi tmp = offsets[i];
                    offsets[i] = offset_sum;
                    offset_sum += tmp;
                }
                return offset_sum;
            } else {
                // Parallel Prefix Sum, based on Guy E. Blelloch's Prefix Sums and Their Applications
                idi last_element = offsets[size_offsets - 1];
                idi size = 1 << ((idi) log2(size_offsets));
                idi tmp_element = offsets[size - 1];
                // Up-Sweep (Reduce) Phase
                idi log2size = log2(size);
                for (idi d = 0; d < log2size; ++d) {
                    idi by = 1 << (d + 1);
                    #pragma omp parallel for
                    for (idi k = 0; k < size; k += by) {
                        offsets[k + (1 << (d + 1)) - 1] += offsets[k + (1 << d) - 1];
                    }
                }
                // Down-Sweep Phase
                offsets[size - 1] = 0;
                for (idi d = log2(size) - 1; d != (idi) -1 ; --d) {
                    idi by = 1 << (d + 1);
                    #pragma omp parallel for
                    for (idi k = 0; k < size; k += by) {
                        idi t = offsets[k + (1 << d) - 1];
                        offsets[k + (1 << d) - 1] = offsets[k + (1 << (d + 1)) - 1];
                        offsets[k + (1 << (d + 1)) - 1] += t;
                    }
                }
                if (size != size_offsets) {
                    idi tmp_sum = offsets[size - 1] + tmp_element;
                    for (idi i = size; i < size_offsets; ++i) {
                        idi t = offsets[i];
                        offsets[i] = tmp_sum;
                        tmp_sum += t;
                    }
                }
                return offsets[size_offsets - 1] + last_element;
            }
        }

        /**
         * @Author: wanjingyi
         * @description: send messages to neighbors
         * @param {*}
         * @return {*}
         */        
        inline void send_messages(
            idi v_head,
            idi roots_start,
            WGraph& wgraph,
            vector< vector<weightiLarge> > &dists_table,
            vector<ShortIndex> &short_index,
            vector<idi> &tmp_has_cand_queue,
            idi &size_tmp_has_cand_queue,
            idi offset_tmp_queue,
            vector<uint8_t> &has_candidates,
            vector<idi> &tmp_once_had_cand_queue,
            idi &size_tmp_once_had_cand_queue,
            vector<uint8_t> &once_had_cand)
        {
            ShortIndex &SI_v_head = short_index[v_head];
            // Traverse v_head's every neighbor v_tail
            idi e_i_start = wgraph.vertices[v_head];
            idi e_i_bound = wgraph.vertices[v_head+1];
            for (idi e_i = e_i_start; e_i < e_i_bound; ++e_i) {
                idi v_tail = wgraph.edges[e_i].first;
                if (v_tail <= roots_start) { // v_tail has higher rank than any roots, then no roots can push new labels to it.
                    break;
                }
                weightiLarge weight_h_t = wgraph.edges[e_i].second;
                ShortIndex &SI_v_tail = short_index[v_tail];
                bool got_candidates = false; // A flag indicates if v_tail got new candidates.
                // Sequential Version
                idi bound_r_i = SI_v_head.end_last_new_roots;
                for (idi r_i = 0; r_i < bound_r_i; ++r_i) {
                    idi r_root_id = SI_v_head.last_new_roots[r_i]; // last inserted label r_root_id of v_head
                    if (v_tail <= r_root_id + roots_start) {
                        continue;
                    }
                    weightiLarge tmp_dist_r_t = SI_v_head.vertices_dists[r_root_id] + weight_h_t;
                    if (tmp_dist_r_t < dists_table[v_tail][r_root_id] && tmp_dist_r_t < SI_v_tail.candidates_dists[r_root_id]) {
                        // Mark r_root_id as a candidate of v_tail
                        if (CAS(&SI_v_tail.candidates_dists[r_root_id], WEIGHTILARGE_MAX, tmp_dist_r_t)) {
                            TS_enqueue(SI_v_tail.candidates_que, SI_v_tail.end_candidates_que, r_root_id);
                        }
                        volatile weightiLarge old_v = SI_v_tail.candidates_dists[r_root_id];
                        while (tmp_dist_r_t < old_v	&& !CAS(&SI_v_tail.candidates_dists[r_root_id], old_v, tmp_dist_r_t)) {
                            old_v = SI_v_tail.candidates_dists[r_root_id];
                        }
                        if (CAS(&dists_table[v_tail][r_root_id], WEIGHTILARGE_MAX, tmp_dist_r_t)) {
                            TS_enqueue(SI_v_tail.reached_roots_que, SI_v_tail.end_reached_roots_que, r_root_id);
                        }
                        old_v = dists_table[v_tail][r_root_id];
                        while (tmp_dist_r_t < old_v && !CAS(&dists_table[v_tail][r_root_id], old_v, tmp_dist_r_t)) {
                            old_v = dists_table[v_tail][r_root_id];
                        }
                        got_candidates = true;
                    }
                }
                    if (got_candidates) {
                    if (CAS(&has_candidates[v_tail], (uint8_t) 0, (uint8_t) 1)) {
                        tmp_has_cand_queue[offset_tmp_queue + size_tmp_has_cand_queue++] = v_tail;
                    }

                    if (CAS(&once_had_cand[v_tail], (uint8_t) 0, (uint8_t) 1)) {
                        tmp_once_had_cand_queue[offset_tmp_queue + size_tmp_once_had_cand_queue++] = v_tail;
                    }
                }
            }
            SI_v_head.end_last_new_roots = 0; // clear v_head's last_new_roots
        }

        /**
         * @Author: wanjingyi
         * @description: initialize all variables
         * @param {*}
         * @return {*}
         */        
        inline void initialize_tables(
			vector<ShortIndex> &short_index,
			vector< vector<weightiLarge> > &roots_labels_buffer,
			vector<idi> &active_queue,
			idi &end_active_queue,
			idi roots_start,
			inti roots_size,
			vector<IndexType> &L,
			vector<idi> &has_new_labels_queue,
			idi &end_has_new_labels_queue,
			vector<uint8_t> &has_new_labels)
        {
            idi roots_bound = roots_start + roots_size;
            // Distance Tables
            {
                // Traverse all roots
                #pragma omp parallel for
                for (idi r_root_id = 0; r_root_id < roots_size; ++r_root_id) {
                    idi r_real_id = r_root_id + roots_start;
                    const IndexType &Lr = L[r_real_id];
                    // Traverse r_real_id's all labels to initial roots_labels_buffer
                    // Sequential Version
                    idi l_i_bound = Lr.vertices.size();
                    _mm_prefetch(&Lr.vertices[0], _MM_HINT_T0);
                    _mm_prefetch(&Lr.distances[0], _MM_HINT_T0);
                    for (idi l_i = 0; l_i < l_i_bound; ++l_i) {
                        roots_labels_buffer[r_root_id][Lr.vertices[l_i]] = Lr.distances[l_i];
                    }

                    roots_labels_buffer[r_root_id][r_real_id] = 0;
                }
            }
            // Short_index
            {
                #pragma omp parallel for
                for (idi r_real_id = roots_start; r_real_id < roots_bound; ++r_real_id) {
                    ShortIndex &SI_r = short_index[r_real_id];
                    idi r_root_id = r_real_id - roots_start;
                    // Record itself as last inserted label distance
                    SI_r.vertices_que[SI_r.end_vertices_que++] = r_root_id;
                    SI_r.vertices_dists[r_root_id] = 0;
                    SI_r.last_new_roots[SI_r.end_last_new_roots++] = r_root_id;
                }
            }

            // Active queue
            {
                for (idi r_real_id = roots_start; r_real_id < roots_bound; ++r_real_id) {
                    active_queue[end_active_queue++] = r_real_id;
                }
            }

            // has_new_labels_queue: Put all roots into the has_new_labels_queue
            {
                for (idi r_real_id = roots_start; r_real_id < roots_bound; ++r_real_id) {
                    has_new_labels_queue[end_has_new_labels_queue++] = r_real_id;
                    has_new_labels[r_real_id] = 1;
                }
            }
        }

        /**
         * @Author: wanjingyi
         * @description: indexing in each batch
         * @param {*}
         * @return {*}
         */        
        inline void vertex_centric_labeling_in_batches(
			WGraph& wgraph,
			idi b_id,
			idi roots_start,
			inti roots_size,
			vector<IndexType> &L,
			vector<idi> &active_queue,
			idi end_active_queue,
			vector<uint8_t> &is_active,
			vector<idi> &has_cand_queue,
			idi end_has_cand_queue,
			vector<uint8_t> &has_candidates,
			vector< vector<weightiLarge> > &dists_table,
			vector< vector<weightiLarge> > &roots_labels_buffer,
			vector<ShortIndex> &short_index,
			vector<idi> &has_new_labels_queue,
			idi end_has_new_labels_queue,
			vector<uint8_t> &has_new_labels,
			vector<idi> &once_had_cand_queue,
			idi end_once_had_cand_queue,
			vector<uint8_t> &once_had_cand)
        {
            // At the beginning of a batch, initialize the labels L and distance buffer dist_matrix;
            initialize_tables(
                short_index,
                roots_labels_buffer,
                active_queue,
                end_active_queue,
                roots_start,
                roots_size,
                L,
                has_new_labels_queue,
                end_has_new_labels_queue,
                has_new_labels);

            while (0 != end_active_queue) {
                // First stage, sending distances.
                {
                    // Prepare for parallel processing the active_queue and adding to candidate_queue.
                    // Every vertex's offset location in tmp_candidate_queue
                    // It's used for every thread to write into tmp_candidate_queue and tmp_once_candidated_queue
                    vector<idi> offsets_tmp_queue(end_active_queue);
                    #pragma omp parallel for
                    for (idi i_queue = 0; i_queue < end_active_queue; ++i_queue) {
                        // Traverse all active vertices, get their out degrees.
                        offsets_tmp_queue[i_queue] = wgraph.vertices[active_queue[i_queue]+1]-wgraph.vertices[active_queue[i_queue]];
                        //offsets_tmp_queue[i_queue] = G.out_degrees[active_queue[i_queue]];
                    }
                    idi num_neighbors = prefix_sum_for_offsets(offsets_tmp_queue);
                    // every thread writes to tmp_candidate_queue at its offset location
                    vector<idi> tmp_has_cand_queue(num_neighbors);
                    // A vector to store the true number of pushed neighbors of every active vertex.
                    vector<idi> sizes_tmp_has_cand_queue(end_active_queue, 0);

                    vector<idi> tmp_once_had_cand_queue(num_neighbors);
                    vector<idi> sizes_tmp_once_had_cand_queue(end_active_queue, 0);


                    // Traverse the active queue, every active vertex sends distances to its neighbors
                    #pragma omp parallel for
                    for (idi i_queue = 0; i_queue < end_active_queue; ++i_queue) {
                        idi v_head = active_queue[i_queue];
                        is_active[v_head] = 0; // reset is_active

                        send_messages(
                                v_head,
                                roots_start,
                                wgraph,
                                dists_table,
                                short_index,
                                tmp_has_cand_queue,
                                sizes_tmp_has_cand_queue[i_queue],
                                offsets_tmp_queue[i_queue],
                                has_candidates,
                                tmp_once_had_cand_queue,
                                sizes_tmp_once_had_cand_queue[i_queue],
                                once_had_cand);
                    }

                    // According to sizes_tmp_candidate_queue, get the offset for inserting to the real queue
                    idi total_new = prefix_sum_for_offsets(sizes_tmp_has_cand_queue);
                    // Collect all candidate vertices from tmp_candidate_queue into candidate_queue.
                    collect_into_queue(
                                tmp_has_cand_queue,
                                offsets_tmp_queue, // the locations in tmp_queue for writing from tmp_queue
                                sizes_tmp_has_cand_queue, // the locations in queue for writing into queue.
                                total_new, // total number of elements which need to be added from tmp_queue to queue
                                has_cand_queue,
                                end_has_cand_queue);
                    total_new = prefix_sum_for_offsets(sizes_tmp_once_had_cand_queue);
                    collect_into_queue(
                            tmp_once_had_cand_queue,
                            offsets_tmp_queue,
                            sizes_tmp_once_had_cand_queue,
                            total_new,
                            once_had_cand_queue,
                            end_once_had_cand_queue);

                    end_active_queue = 0; // Set the active_queue empty
                }

                // Second stage, adding
                {
                    // Prepare for parallel processing the candidate_queue and adding to active_queue.
                    // Every vertex's offset location in tmp_active_queue is i_queue * roots_size
                    // It's used for every thread to write into tmp_candidate_queue and tmp_once_candidated_queue
                    vector<idi> offsets_tmp_queue(end_has_cand_queue);
                    #pragma omp parallel for
                    for (idi i_queue = 0; i_queue < end_has_cand_queue; ++i_queue) {
                        // Traverse all active vertices, get their out degrees.
                        // A ridiculous bug here. The v_id will, if any, only add itself to the active queue.
                        //offsets_tmp_queue[i_queue] = i_queue * roots_size;
                        offsets_tmp_queue[i_queue] = i_queue;
                    }
                    // every thread writes to tmp_active_queue at its offset location
                    vector<idi> tmp_active_queue(end_has_cand_queue);
                    // A vector to store the true number of pushed neighbors of every active vertex.
                    vector<idi> sizes_tmp_active_queue(end_has_cand_queue, 0);
                    // every thread writes to its offset location.
                    vector<idi> tmp_has_new_labels_queue(end_has_cand_queue);
                    // Store the size of every group recorded by every thread.
                    vector<idi> sizes_tmp_has_new_labels_queue(end_has_cand_queue, 0);

                    // Traverse vertices in the has_cand_queue to insert labels
                    #pragma omp parallel for
                    for (idi i_queue = 0; i_queue < end_has_cand_queue; ++i_queue) {
                        idi v_id = has_cand_queue[i_queue];
                        bool need_activate = 0;
                        ShortIndex &SI_v = short_index[v_id];
                        // Traverse v_id's all candidates
                        inti bound_cand_i = SI_v.end_candidates_que;
                        for (inti cand_i = 0; cand_i < bound_cand_i; ++cand_i) {
                            inti cand_root_id = SI_v.candidates_que[cand_i];
                            weightiLarge tmp_dist_v_c = SI_v.candidates_dists[cand_root_id];
                            // Distance check for pruning
                            weightiLarge query_dist_v_c;
                            if (WEIGHTILARGE_MAX ==
                                    (query_dist_v_c = distance_query(
                                            v_id,
                                            cand_root_id,
                                            //dists_table,
                                            roots_labels_buffer,
                                            short_index,
                                            L,
                                            roots_start,
                                            tmp_dist_v_c))) {
                                if (WEIGHTILARGE_MAX == SI_v.vertices_dists[cand_root_id]) {
                                    // Record cand_root_id as v_id's label
                                    SI_v.vertices_que[SI_v.end_vertices_que++] = cand_root_id;

                                }
                                // Record the new distance in the label table
                                SI_v.vertices_dists[cand_root_id] = tmp_dist_v_c;
                                SI_v.last_new_roots[SI_v.end_last_new_roots++] = cand_root_id;
                                need_activate = 1;
                            } else if (query_dist_v_c < tmp_dist_v_c){
                                // Update the dists_table
                                dists_table[v_id][cand_root_id] = query_dist_v_c;
                                // Need to send back the distance
                                // DEPRECATED! The correction should not be done here, because some shorter distance does not mean wrong label distances.
                                //					send_back(
                                //							v_id,
                                //							cand_root_id,
                                //							G,
                                //							dists_table,
                                //							short_index,
                                //							roots_start);
                            }
                        }
                        if (need_activate) {
                            if (!is_active[v_id]) {
                                is_active[v_id] = 1;
                                tmp_active_queue[offsets_tmp_queue[i_queue] + sizes_tmp_active_queue[i_queue]++] = v_id;
                            }
                            if (!has_new_labels[v_id]) {
                                has_new_labels[v_id] = 1;
                                tmp_has_new_labels_queue[offsets_tmp_queue[i_queue] + sizes_tmp_has_new_labels_queue[i_queue]++] = v_id;
                            }
                        }
                    }
                    // According to sizes_tmp_active_queue, get the offset for inserting to the real queue
                    idi total_new = prefix_sum_for_offsets(sizes_tmp_active_queue);
                    // Collect all candidate vertices from tmp_active_queue into active_queue.
                    collect_into_queue(
                                tmp_active_queue,
                                offsets_tmp_queue, // the locations in tmp_queue for writing from tmp_queue
                                sizes_tmp_active_queue, // the locations in queue for writing into queue.
                                total_new, // total number of elements which need to be added from tmp_queue to queue
                                active_queue,
                                end_active_queue);
                    // According to sizes_tmp_active_queue, get the offset for inserting to the real queue
                    total_new = prefix_sum_for_offsets(sizes_tmp_has_new_labels_queue);
                    // Collect all candidate vertices from tmp_candidate_queue into candidate_queue.
                    collect_into_queue(
                                tmp_has_new_labels_queue,
                                offsets_tmp_queue, // the locations in tmp_queue for writing from tmp_queue
                                sizes_tmp_has_new_labels_queue, // the locations in queue for writing into queue.
                                total_new, // total number of elements which need to be added from tmp_queue to queue
                                has_new_labels_queue,
                                end_has_new_labels_queue);
                    // Reset vertices' candidates_que and candidates_dists
                    #pragma omp parallel for
                    for (idi i_queue = 0; i_queue < end_has_cand_queue; ++i_queue) {
                        idi v_id = has_cand_queue[i_queue];
                        has_candidates[v_id] = 0; // reset has_candidates
                        ShortIndex &SI_v = short_index[v_id];

                        // Semi-Vectorized Version
                        inti remainder_simd = SI_v.end_candidates_que % NUM_P_INT;
                        inti bound_cand_i = SI_v.end_candidates_que - remainder_simd;
                        for (inti cand_i = 0; cand_i < bound_cand_i; cand_i += NUM_P_INT) {
                            //uint32_t inti >>> __m512i
                            __m512i cand_root_id_v = _mm512_load_epi32(&SI_v.candidates_que[cand_i]); // @suppress("Function cannot be resolved")
                            int_i32scatter_epiI(
                                    &SI_v.candidates_dists[0],
                                    cand_root_id_v,
                                    WEIGHTILARGE_MAX);
                        }
                        if (remainder_simd) {
                            __mmask16 in_m = (__mmask16) ((uint16_t) 0xFFFF >> (NUM_P_INT - remainder_simd));
                            __m512i cand_root_id_v = _mm512_mask_loadu_epi32(UNDEF_i32_v, in_m, &SI_v.candidates_que[bound_cand_i]);
                            int_mask_i32scatter_epiI(
                                    &SI_v.candidates_dists[0],
                                    in_m,
                                    cand_root_id_v,
                                    WEIGHTILARGE_MAX);
                        }

                        // Sequential Version
                        //			inti bound_cand_i = SI_v.end_candidates_que;
                        //			for (inti cand_i = 0; cand_i < bound_cand_i; ++cand_i) {
                        //				inti cand_root_id = SI_v.candidates_que[cand_i];
                        //				SI_v.candidates_dists[cand_root_id] = WEIGHTILARGE_MAX; // Reset candidates_dists after using in distance_query.
                        //			}
                        SI_v.end_candidates_que = 0; // Clear v_id's candidates_que
                    }
                    end_has_cand_queue = 0; // Set the has_cand_queue empty
                }
            }

            // Filter out wrong and redundant labels in this batch
            filter_out_labels(
                    has_new_labels_queue,
                    end_has_new_labels_queue,
                    short_index,
                    roots_start);

            // Reset dists_table and short_index
            reset_tables(
                    short_index,
                    roots_start,
                    roots_size,
                    L,
                    dists_table,
                    roots_labels_buffer,
                    once_had_cand_queue,
                    end_once_had_cand_queue,
                    once_had_cand);

            // Update the index according to labels_table
            update_index(
                    L,
                    short_index,
                    roots_start,
                    has_new_labels_queue,
                    end_has_new_labels_queue,
                    has_new_labels);
        }

        //************************************construct*******************************
        void construct(WGraph &wgraph)
        {
            L.resize(numOfVertices);
            idi remainer = numOfVertices % BATCH_SIZE;
            idi b_i_bound = numOfVertices / BATCH_SIZE;
            // Active queue
            vector<idi> active_queue(numOfVertices);
            idi end_active_queue = 0;
            vector<uint8_t> is_active(numOfVertices, 0);// is_active[v] is true means vertex v is in the active queue.
            // Queue of vertices having candidates
            vector<idi> has_cand_queue(numOfVertices);
            idi end_has_cand_queue = 0;
            vector<uint8_t> has_candidates(numOfVertices, 0); // has_candidates[v] is true means vertex v is in the queue has_cand_queue.
            // Distance table of shortest distance from roots to other vertices.
            vector< vector<weightiLarge> > dists_table(numOfVertices, vector<weightiLarge>(BATCH_SIZE, WEIGHTILARGE_MAX));
                // The distance table is roots_sizes by N. 
                // 1. record the shortest distance so far from every root to every vertex; (from v to root r)
                // 2. (DEPRECATED! Replaced by roots_labels_buffer) The distance buffer, recording old label distances of every root. It needs to be initialized every batch by labels of roots. (dists_table[r][l], r > l)
            // A tabel for roots' old label distances (which are stored in the dists_table in the last version)
            vector< vector<weightiLarge> > roots_labels_buffer(BATCH_SIZE, vector<weightiLarge>(numOfVertices, WEIGHTILARGE_MAX)); // r's old label distance from r to v
            // Every vertex has a ShortIndex object; the previous labels_table is now in ShortIndex structure
            vector<ShortIndex> short_index(numOfVertices, ShortIndex(BATCH_SIZE)); // Here the size of short_index actually is fixed because it is static.
                // Temporary distance table, recording in the current iteration the traversing distance from a vertex to a root.
                // The candidate table is replaced by the ShortIndex structure: every vertex has a queue and a distance array;
                // 1. the queue records last inserted labels.
                // 2. the distance array acts like a bitmap but restores distances.

            // A queue to store vertices which have got new labels in this batch. This queue is used for reset dists_table.
            vector<idi> has_new_labels_queue(numOfVertices);
            idi end_has_new_labels_queue = 0;
            vector<uint8_t> has_new_labels(numOfVertices, 0);

            // A queue to store vertices which have got candidates in this batch. This queue is used fore reset dists_table, accompanying with reached_roots_que in the ShortIndex.
            vector<idi> once_had_cand_queue(numOfVertices);
            idi end_once_had_cand_queue = 0;
            vector<uint8_t> once_had_cand(numOfVertices, 0);

            for (idi b_i = 0; b_i < b_i_bound; ++b_i) {
                vertex_centric_labeling_in_batches(
                        wgraph,
                        b_i,
                        b_i * BATCH_SIZE,
                        BATCH_SIZE,
                        L,
                        active_queue,
                        end_active_queue,
                        is_active,
                        has_cand_queue,
                        end_has_cand_queue,
                        has_candidates,
                        dists_table,
                        roots_labels_buffer,
                        short_index,
                        has_new_labels_queue,
                        end_has_new_labels_queue,
                        has_new_labels,
                        once_had_cand_queue,
                        end_once_had_cand_queue,
                        once_had_cand);
            }
            if (0 != remainer) {
                vertex_centric_labeling_in_batches(
                        wgraph,
                        b_i_bound,
                        b_i_bound * BATCH_SIZE,
                        remainer,
                        L,
                        active_queue,
                        end_active_queue,
                        is_active,
                        has_cand_queue,
                        end_has_cand_queue,
                        has_candidates,
                        dists_table,
                        roots_labels_buffer,
                        short_index,
                        has_new_labels_queue,
                        end_has_new_labels_queue,
                        has_new_labels,
                        once_had_cand_queue,
                        end_once_had_cand_queue,
                        once_had_cand);
            }
        }

    };
}

#endif // ! PARALLEL_CONSTRUCTION_VEC_H