#! /bin/bash
#用于生成WHP

DATADIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes/Datasets"
OUTPUTDIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes"
EXP_PREFIX="exp2_spatial"

PLL=("WHP")
ORDERING_MODEL=2 
SHUFFLE="1"
NUMQUERY="1000000"

EXE_INDEX="/home/wanwanvv/workspace/QW_SPD_Codes/PLL_WHP/bin/pll_whp_run"
EXE_QUERY="/home/wanwanvv/workspace/QW_SPD_Codes/Query_processing/bin/query_run"

#程序运行参数
DIRECTED_FLAG=(0)
DIRECTED_PREFIX=("Undirected")
BETWEENNESS_FILE=("original.betweenness")

#数据集
DATASETS=("Manhattan" "NewYork" "Chengdu" "Haikou" "Florida" "California")

#query数据范围
DAY_INDEX_START=20
DAY_INDEX_END=30

#权重参数(15)
FREQ_PARAS=(4 90 2 6)
TOTAL_PARA=100

#计数程序执行次数
CNT=0

#查询执行次数
QUERY_CNT=1

#建索引执行次数
INDEXING_CNT=1

for (( d=0;d<${#DIRECTED_FLAG[@]};d++))
do
	echo "==============================${DIRECTED_PREFIX[d]}=============================="
	output_prefix=${OUTPUTDIR_PREFIX}/${EXP_PREFIX}/${DIRECTED_PREFIX[d]}/${PLL}
	result_output_prefix=${output_prefix}/Result
	result_file=${result_output_prefix}/results.txt
	if [ -f "${result_file}" ]; then  
		rm -f ${result_file}
	fi
	echo "dataset ordering_time labeling_time day20 day21 day22 day23 day24 day25 day26 day27 day28 day29 label_size" >> ${result_file}
	for(( t=0;t<4;t++))
	do
		echo "==============================${DATASETS[t]}=============================="
		label_file=${output_prefix}/Labels/${DATASETS[t]}.label
		label_size_file=${output_prefix}/Label_size/${DATASETS[t]}.size
		dataset_prefix=${DATADIR_PREFIX}/${DATASETS[t]}
		graph_file=${dataset_prefix}/graph/${DATASETS[t]}.graph
		betwenness_file=${dataset_prefix}/order_values/${BETWEENNESS_FILE[d]}
		query_freq_file=${dataset_prefix}/${EXP_PREFIX}/freq_pred/day${DAY_INDEX_START}-cut0.txt
		query_pair_dir=${dataset_prefix}/${EXP_PREFIX}/query
		#参数
		bet_para=$[${TOTAL_PARA}-${FREQ_PARAS[t]}]
		#输出datset
		echo -n "${DATASETS[t]} " >> ${result_file}
		for(( k=0;k<${INDEXING_CNT};k++))
		do
			#执行建label程序
			${EXE_INDEX} -p -d ${DIRECTED_FLAG[d]} -w 1 -s 2 -g ${graph_file} -b ${betwenness_file} -e ${label_file} -f ${result_file} -q ${query_freq_file} -i ${label_size_file} -k ${FREQ_PARAS[t]} -l ${bet_para} -j 100
			let CNT++
		done
		if [ ! -d "${query_pair_dir}" ]; then
			echo "${query_pair_dir} doesn't exist!"  
			break
		fi
		for(( n=0;n<${QUERY_CNT};n++))
		do	
			#测试查询时间(-z ${label_size_file} )
			${EXE_QUERY} -q -d ${DIRECTED_FLAG[d]} -w 1 -s 2 -b 1 -p 1 -l ${label_file} -t ${query_pair_dir} -a ${result_file} -z ${label_size_file} -n ${NUMQUERY} -x ${SHUFFLE}
			let CNT++
		done
		#输出size
		du -sh ${label_file}|awk '{print $1}' >> ${result_file}
	done
done
