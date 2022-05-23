#! /bin/bash
#用于生成WHP

DATADIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes/Datasets"
OUTPUTDIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes"
EXP_PREFIX="exp4_spatial_temporal_tune"

PLL=("WHP")
ORDERING_MODEL=2 
SHUFFLE="1"
NUMQUERY="10000"

EXE_INDEX="/home/wanwanvv/workspace/QW_SPD_Codes/PLL_WHP/bin/pll_whp_run"
EXE_QUERY="/home/wanwanvv/workspace/QW_SPD_Codes/Query_processing/bin/query_run"

#程序运行参数
DIRECTED_FLAG=(0)
DIRECTED_PREFIX=("Undirected")
BETWEENNESS_FILE=("original.betweenness")

#数据集
DATASETS=("Manhattan" "NewYork" "Chengdu" "Haikou" "Florida" "California")

#划分参数
DAY_CNT_NUM=96
TIP_INDEX_PREFIX="RL-TIP"

#计数程序执行次数
CNT=0

#查询执行次数
QUERY_CNT=1

#建索引执行次数
INDEXING_CNT=1

for (( d=0;d<${#DIRECTED_FLAG[@]};d++))
do
	echo "==============================${DIRECTED_PREFIX[d]}=============================="
	output_prefix=${OUTPUTDIR_PREFIX}/${EXP_PREFIX}/${DIRECTED_PREFIX[d]}/RL-TIP_${PLL}
	result_output_prefix=${output_prefix}/Result
	label_dir=${output_prefix}/Labels
	label_size_dir=${output_prefix}/Label_size
	if [ -d "${result_output_prefix}" ]; then  
		rm -rf ${result_output_prefix}
	fi
	mkdir ${result_output_prefix}
	query_result_file=${result_output_prefix}/query_results.txt
	label_result_file=${result_output_prefix}/label_results.txt
	echo "dataset time_interval_cnt total_ordering_timet total_hub_time total_labling_time total_label_size" >> ${label_result_file}
	echo "day_index cut_index query_cost average_query_time">> ${query_result_file}
	#for(( t=0;t<${#DATASETS[@]};t++))
	for(( j=2;j<3;j++))
	do
		echo "==============================${DATASETS[j]}=============================="
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
		betwenness_file=${dataset_prefix}/order_values/${BETWEENNESS_FILE[d]}
		query_freq_dir=${dataset_prefix}/${EXP_PREFIX}/freq_pred
		query_pair_dir=${dataset_prefix}/${EXP_PREFIX}/query_next
		tip_file=${dataset_prefix}/${EXP_PREFIX}/time_partitions/${TIP_INDEX_PREFIX}.txt
		paras_file=${dataset_prefix}/${EXP_PREFIX}/parameters/${PLL}.paras
		echo -n "${DATASETS[j]} " >> ${label_result_file}
		for(( k=0;k<${INDEXING_CNT};k++))
		do
			#执行批量建WHP程序(-q ${query_pair_dir} )
			${EXE_INDEX} -t -d ${DIRECTED_FLAG[d]} -w 1 -s 1 -u 0 -n ${DAY_CNT_NUM} -g ${graph_file} -h ${query_freq_dir} -l ${label_dir} -b ${label_size_dir} -p ${paras_file} -f ${betwenness_file} -k ${tip_file} -a ${label_result_file}
			let CNT++
		done
		#输出size
		du -sh ${label_dir}|awk '{print $1}' >> ${label_result_file}
		if [ ! -d "${query_pair_dir}" ]; then
			echo "${query_pair_dir} doesn't exist!"  
			break
		fi
		echo "${DATASETS[j]}" >> ${query_result_file}
		for(( k=0;k<${QUERY_CNT};k++))
		do	
			#测试查询时间
			${EXE_QUERY} -q -d ${DIRECTED_FLAG[d]} -w 1 -s 3 -b 1 -p 1 -o 1 -l ${label_dir} -t ${query_pair_dir} -a ${query_result_file} -z ${label_size_dir} -n ${NUMQUERY} -x ${SHUFFLE}
			let CNT++
		done
		#输出换行
		echo >> ${query_result_file}
	done
done
