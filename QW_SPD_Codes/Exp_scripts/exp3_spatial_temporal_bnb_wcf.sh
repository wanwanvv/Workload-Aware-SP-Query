#! /bin/bash
#用于批量构建混合方案的实验

#程序路径
EXE_INDEX="/home/wanwanvv/workspace/QW_SPD_Codes/WCF_DU/bin/wcf_du_run"
EXE_QUERY="/home/wanwanvv/workspace/QW_SPD_Codes/Query_processing/bin/query_run"
EXP_PREFIX="exp3_spatial_temporal"

#路经
DATADIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes/Datasets"
OUTPUTDIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes"
INDEX_PREFIX="WCF"

#数据集
DATASETS=("Manhattan" "NewYork" "Chengdu" "Haikou" "California" "Florida")

#划分参数
DAY_CNT_NUM=96
TIP_INDEX_PREFIX="BNB"

#构建索引参数
INDEXING_PARAS=(1) #(1 2)
DIRECTED_FLAG=(0) #(0 1)
DIRECTED_PREFIX=("Undirected") #("Undirected" "Directed")

#并行线程数
THREADS_NUM_PARAS=5

#查询参数
SHUFFLE="1"
NUMQUERY="1000000"

#程序计数器
CNT=0

#查询执行次数
QUERY_CNT=1

#建索引执行次数
INDEXING_CNT=1

for (( d=0;d<${#DIRECTED_FLAG[@]};d++))
do
	echo "==============================${DIRECTED_PREFIX[d]}=============================="
	output_path=${OUTPUTDIR_PREFIX}/${EXP_PREFIX}/${DIRECTED_PREFIX[d]}
	for(( i=0;i<${#INDEXING_PARAS[@]};i++))
	do
		echo "==============================${INDEX_PREFIX}_${INDEXING_PARAS[i]}=============================="
		output_prefix=${output_path}/${TIP_INDEX_PREFIX}_${INDEX_PREFIX}_${INDEXING_PARAS[i]}
		result_output_prefix=${output_prefix}/Result
		inter_file_dir=${output_prefix}/Inter_files
		label_dir=${output_prefix}/Labels
		label_size_dir=${output_prefix}/Label_size
		#重建结果文件夹
		if [ -d "${result_output_prefix}" ]; then  
			rm -rf ${result_output_prefix}
		fi
		mkdir ${result_output_prefix}
		result_file=${result_output_prefix}/results.txt
		echo "dataset time_interval_cnt total_labeling_time total_ordering_time total_compute_rmq_time query_cost_all_slices query_cost_avg_slice qtime_all_slices avg_qtime_avg_slice total_label_size" >> ${result_file}
		# ${#DATASETS[@]}
		for(( j=0;j<4;j++)) 
		do
			echo "==============================${DATASETS[j]}=============================="
			#重建中间结果文件夹
			if [ -d "${inter_file_dir}" ]; then  
				rm -rf ${inter_file_dir}
			fi
			mkdir ${inter_file_dir}
			#重建索引文件夹
			if [ -d "${label_dir}" ]; then  
				rm -rf ${label_dir}
			fi
			mkdir ${label_dir}
			#重建索引大小文件夹
			if [ -d "${label_size_dir}" ]; then  
				rm -rf ${label_size_dir}
			fi
			mkdir ${label_size_dir}
			dataset_prefix=${DATADIR_PREFIX}/${DATASETS[j]}
			graph_file=${dataset_prefix}/graph/${DATASETS[j]}.graph
			#查询频率文件夹
			query_freq_dir=${dataset_prefix}/${EXP_PREFIX}/freq_pred
			#查询workload文件夹
			query_pair_dir_ori=${dataset_prefix}/${EXP_PREFIX}/query
			if [ ! -d "${query_pair_dir_ori}" ]; then
				echo "${query_pair_dir_ori} doesn't exist!"  
				break
			fi
			query_pair_dir=${dataset_prefix}/${EXP_PREFIX}/query_bk
			if [ -d "${query_pair_dir}" ]; then  
				rm -rf ${query_pair_dir}
			fi
			#复制workload文件夹
		    	if [ -d "${query_pair_dir_ori}" ]; then  
				cp -r ${query_pair_dir_ori} ${query_pair_dir}
		    	fi
			echo -n "${DATASETS[j]} " >> ${result_file}
			tip_file=${dataset_prefix}/${EXP_PREFIX}/time_partitions/${TIP_INDEX_PREFIX}.txt
			paras_file=${dataset_prefix}/${EXP_PREFIX}/parameters/${INDEX_PREFIX}_${INDEXING_PARAS[i]}.paras
			if [ ! -d "${query_pair_dir}" ]; then  
				echo "${query_pair_dir} doesn't exist!"
				break
			fi
			for ((l=0; l<${INDEXING_CNT}; l++))
			do
				#执行批量构建索引
				${EXE_INDEX} -t -d ${DIRECTED_FLAG[d]} -w 1 -s 1 -o 0 -u 0 -m ${INDEXING_PARAS[i]} -n ${DAY_CNT_NUM} -c ${THREADS_NUM_PARAS} -g ${graph_file} -h ${query_freq_dir} -q ${query_pair_dir} -l ${label_dir} -b ${label_size_dir} -i ${inter_file_dir} -k ${tip_file} -p ${paras_file} -a ${result_file}
				let CNT++
			done
			for ((k=0; k<${QUERY_CNT}; k++))
			do
				#测试查询时间
				${EXE_QUERY} -l -d ${DIRECTED_FLAG[d]} -w 1 -p 1 -b 1 -o 0 -s 3 -x ${SHUFFLE} -n ${NUMQUERY} -m ${INDEXING_PARAS[i]} -i ${label_dir} -e ${label_size_dir} -v ${inter_file_dir} -y ${inter_file_dir} -a ${result_file} -q ${query_pair_dir}
				let CNT++
			done
			#输出size
			du -sh ${label_dir}|awk '{print $1}' >> ${result_file}
			#删除复制的workload文件夹
			rm -rf ${query_pair_dir}
		done
	done
done
echo Done.
echo ""
