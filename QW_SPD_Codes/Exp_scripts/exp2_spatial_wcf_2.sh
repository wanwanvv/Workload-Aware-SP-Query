#! /bin/bash
#用于批量构建混合方案的实验

#程序路径
EXE_INDEX="/home/wanwanvv/workspace/QW_SPD_Codes/WCF_DU/bin/wcf_du_run"
EXE_QUERY="/home/wanwanvv/workspace/QW_SPD_Codes/Query_processing/bin/query_run"
EXP_PREFIX="exp2_spatial"

#路经
DATADIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes/Datasets"
OUTPUTDIR_PREFIX="/home/wanwanvv/workspace/QW_SPD_Codes"
INDEX_PREFIX="WCF"

#数据集
DATASETS=("Manhattan" "NewYork" "Chengdu" "Haikou" "California" "Florida")

#删除点参数
MAX_DEGREE=(10 65 35 8)

#overlay pll参数
FREQ_PARAS=(4 30 4 10)
TOTAL_PARA=100

#构建索引参数
INDEXING_PARAS=(2)
DIRECTED_FLAG=(0)
DIRECTED_PREFIX=("Undirected")

#并行线程数
THREADS_NUM_PARAS=5

#查询参数
SHUFFLE="1"
NUMQUERY="1000000"
DAY_INDEX_START=20
DAY_INDEX_END=30

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
		output_prefix=${output_path}/${INDEX_PREFIX}_${INDEXING_PARAS[i]}
		result_output_prefix=${output_prefix}/Result
		result_file=${result_output_prefix}/results.txt
		if [ -f "${result_file}" ]; then  
			rm -f ${result_file}
		fi
		echo "dataset labeling_time ordering_time compute_rmq_time day20 day21 day22 day23 day24 day25 day26 day27 day28 day29 label_size" >> ${result_file}
		for(( j=0;j<4;j++)) 
		do
			echo "==============================${DATASETS[j]}=============================="
			dataset_prefix=${DATADIR_PREFIX}/${DATASETS[j]}
			graph_file=${dataset_prefix}/graph/${DATASETS[j]}.graph
			#输出和中间文件一份
			inter_file_output_prefix=${output_prefix}/Inter_files/${DATASETS[j]}
			label_file=${output_prefix}/Labels/${DATASETS[j]}.label
			label_size_file=${output_prefix}/Label_size/${DATASETS[j]}.size
			lca_file=${inter_file_output_prefix}.lca
			isDeleted_file=${inter_file_output_prefix}.isDeleted
			#查询频率文件
			query_freq_file=${dataset_prefix}/${EXP_PREFIX}/freq_pred/day${DAY_INDEX_START}-cut0.txt
			#查询workload文件夹
			query_pair_dir=${dataset_prefix}/${EXP_PREFIX}/query
			#参数
			bet_para=$[${TOTAL_PARA}-${FREQ_PARAS[j]}]
			echo -n "${DATASETS[j]} " >> ${result_file}
			for ((l=0; l<${INDEXING_CNT}; l++))
			do
				#执行混合方案构建程序验证最大度对实验结果的影响(-b ${label_size_file})
				${EXE_INDEX} -z -d ${DIRECTED_FLAG[d]} -w 1 -s 2 -r 0 -o 0 -i ${INDEXING_PARAS[i]} -m ${MAX_DEGREE[j]} -k ${FREQ_PARAS[j]} -l ${bet_para} -j 100 -c ${THREADS_NUM_PARAS} -g ${graph_file} -q ${query_freq_file} -t ${inter_file_output_prefix} -a ${result_file} -e ${label_file} -b ${label_size_file}
				let CNT++
			done
			if [ ! -d "${query_pair_dir}" ]; then  
				echo "${query_pair_dir} doesn't exist!"
				break
			fi
			for ((k=0; k<${QUERY_CNT}; k++))
			do
				#测试查询时间
				${EXE_QUERY} -l -d ${DIRECTED_FLAG[d]} -w 1 -p 1 -b 1 -s 2 -x ${SHUFFLE} -n ${NUMQUERY} -m ${INDEXING_PARAS[i]} -i ${label_file} -v ${lca_file} -y ${isDeleted_file} -a ${result_file} -e ${label_size_file} -q ${query_pair_dir}
				let CNT++
			done
			#输出size
			du -sh ${label_file}|awk '{print $1}' >> ${result_file}
		done
	done
	echo Done.
	echo ""
done
