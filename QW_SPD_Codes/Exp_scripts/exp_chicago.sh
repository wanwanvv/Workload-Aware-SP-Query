#! /bin/bash
#exp2_spatial

#路经
DATADIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes/Datasets"
OUTPUTDIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes"
EXP_PREFIX="exp2_spatial"

#程序运行参数
PLL=("BHP")
ORDERING_MODEL=(1)
SHUFFLE=1
NUMQUERY=1000000
DIRECTED_FLAG=(0)
DIRECTED_PREFIX=("Undirected")

EXE_INDEX="/home/wanwanvv/workspace/QW_SPD_Codes/PLL_WHP/bin/pll_whp_run"
EXE_QUERY="/home/wanwanvv/workspace/QW_SPD_Codes/Query_processing/bin/query_run"

#数据集
DATASETS=("Manhattan" "NewYork" "Chengdu" "Haikou" "Florida" "California")

#query数据范围
DAY_INDEX_START=20
DAY_INDEX_END=30

#索引构建实验次数
INDEXING_CNT=1

#程序计数
CNT=0

#查询次数
QUERY_CNT=1

for (( d=0;d<${#DIRECTED_FLAG[@]};d++))
do
	echo "==============================${DIRECTED_PREFIX[d]}=============================="
	for (( o=0;o<${#ORDERING_MODEL[@]};o++))
	do
		echo "==============================${PLL[o]}=============================="
		output_prefix=${OUTPUTDIR_PREFIX}/${EXP_PREFIX}/${DIRECTED_PREFIX[d]}/${PLL[o]}
		result_output_prefix=${output_prefix}/Result
		result_file=${result_output_prefix}/results.txt
		#if [ -f "${result_file}" ]; then  
			#rm -f ${result_file}
		#fi
		echo "dataset ordering_time labeling_time day20 day21 day22 day23 day24 day25 day26 day27 day28 day29 label_size" >> ${result_file}
		for(( t=0;t<4;t++))
		do
			echo "==============================${DATASETS[t]}=============================="
			dataset_prefix=${DATADIR_PREFIX}/${DATASETS[t]}
			label_file=${output_prefix}/Labels/${DATASETS[t]}.label
			label_size_file=${output_prefix}/Label_size/${DATASETS[t]}.size
			graph_file=${dataset_prefix}/graph/${DATASETS[t]}.graph
			query_freq_file=${dataset_prefix}/${EXP_PREFIX}/freq_pred/day${DAY_INDEX_START}-cut0.txt
			query_pair_dir=${dataset_prefix}/${EXP_PREFIX}/query
			#输出dataset skew_index
			echo -n "${DATASETS[t]} " >> ${result_file}
			for(( k=0;k<${INDEXING_CNT};k++))
			do
				#***********执行多次求平均值***********
				#执行建label程序
				${EXE_INDEX} -x -d ${DIRECTED_FLAG[d]} -w 1 -s 1 -m ${ORDERING_MODEL[o]} -g ${graph_file} -e ${label_file} -f ${result_file} -i ${label_size_file}
				#***********执行多次求平均值***********
				let CNT++
			done
			if [ ! -d "${query_pair_dir}" ]; then
				echo "${query_pair_dir} doesn't exist!"  
				break
			fi				
			for(( n=0;n<${QUERY_CNT};n++))
			do
				#测试查询时间
				${EXE_QUERY} -q -d ${DIRECTED_FLAG[d]} -w 1 -s 2 -b 1 -p 1 -l ${label_file} -t ${query_pair_dir} -a ${result_file} -z ${label_size_file} -n ${NUMQUERY} -x ${SHUFFLE}
				let CNT++
			done
			#输出size
			du -sh ${label_file}|awk '{print $1}' >> ${result_file}
		done
	done
	done
	echo Done.
echo ""
