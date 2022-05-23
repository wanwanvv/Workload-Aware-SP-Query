#!/bin/bash

#make
make

updatedGraphFile="/data/d/linuxData/ManhattanProcess/Manhattan_update.txt"

outputIndexFile="/data/f/directedTest/Manhattan/Manhattan_H2H.index"
outputSTTable="/data/f/directedTest/Manhattan/Manhattan_ST_H2H.data"
outputBlockSTTable="/data/f/directedTest/Manhattan/Manhattan_ST_block_H2H.data"
outputQueryCost="/data/f/directedTest/Manhattan/queryCost_H2H.txt"
outputOriginalResult="/data/f/directedTest/Manhattan/originalResult-day0-cut7-H2H.txt"
outputResult="/data/f/directedTest/Manhattan/result-day0-cut7-H2H.txt"
numOfQueries=500000
queryPairFile="/data/d/linuxData/Datasets_exp/Manhattan/query/day0-cut7.txt"
HFVerticesFile="/data/d/linuxData/Datasets_exp/Manhattan/freq/day0-cut7.txt"

declare -i day dayy cut cutt gamma eta

# for day in {0..29}
#   do
#   for cut in {7..95}
#     do
    #H2H preprocessing
    ./bin/run -h -d -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult

#     #query test by H2H and O(nlogn) ST table
#     if [ $? == 0 ]; then     
#       echo "faild"
#     else
#       echo "success------------------------------------------------"

#       ./bin/run -q -d -n $numOfQueries -i $outputIndexFile -s $outputSTTable -q $queryPairFile -o $outputOriginalResult
#     fi

#     for gamma in {0..16}
#       do
#       for eta in {0..15}
#         do
#         #WH2H preprocessing
#         if [ $? == 0 ]; then     
#           echo "faild"
#         else
#           echo "success------------------------------------------------"
#           ./bin/run -w -d -r $gamma -e $eta -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult
#         fi
#         #query test by WH2H and O(nlogn) ST table
#         if [ $? == 0 ]; then     
#           echo "faild"
#         else
#           echo "success------------------------------------------------"
#           ./bin/run -q -d -n $numOfQueries -i $outputIndexFile -s $outputSTTable -q $queryPairFile -f $outputResult
#         fi
#       done
#     done
#     cutt=cut+1
#     if [ $cut == 95 ];then
#       cutt=7
#       echo $cutt
#     fi
#     queryPairFile=${queryPairFile//cut$cut/cut$cutt}
#     HFVerticesFile=${HFVerticesFile//cut$cut/cut$cutt}
#     outputOriginalResult=${outputOriginalResult//cut$cut/cut$cutt}
#     outputResult=${outputResult//cut$cut/cut$cutt}
#   done
#   dayy=day+1
#   queryPairFile=${queryPairFile//day$day/day$dayy}
#   HFVerticesFile=${HFVerticesFile//day$day/day$dayy}
#   outputOriginalResult=${outputOriginalResult//day$day/day$dayy}
#   outputResult=${outputResult//day$day/day$dayy}
#   echo $queryPairFile
#   echo $HFVerticesFile
# done

updatedGraphFile=${updatedGraphFile//Manhattan/Haikou}
outputIndexFile=${outputIndexFile//Manhattan/Haikou}
outputSTTable=${outputSTTable//Manhattan/Haikou}
outputBlockSTTable=${outputBlockSTTable//Manhattan/Haikou}
outputQueryCost=${outputQueryCost//Manhattan/Haikou}
outputOriginalResult=${outputOriginalResult//Manhattan/Haikou}
outputResult=${outputResult//Manhattan/Haikou}
queryPairFile=${queryPairFile//Manhattan/Haikou}
HFVerticesFile=${HFVerticesFile//Manhattan/Haikou}

queryPairFile=${queryPairFile//day30/day0}
HFVerticesFile=${HFVerticesFile//day30/day0}
outputOriginalResult=${outputOriginalResult//day30/day0}
outputResult=${outputResult//day30/day0}
queryPairFile=${queryPairFile//cut95/cut7}
HFVerticesFile=${HFVerticesFile//cut95/cut7}
outputOriginalResult=${outputOriginalResult//cut95/cut7}
outputResult=${outputResult//cut95/cut7}

# for day in {0..29}
#   do
#   for cut in {7..95}
#     do
    #H2H preprocessing
    ./bin/run -h -d -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult

#     #query test by H2H and O(nlogn) ST table
#     if [ $? == 0 ]; then     
#       echo "faild"
#     else
#       echo "success------------------------------------------------"

#       ./bin/run -q -d -n $numOfQueries -i $outputIndexFile -s $outputSTTable -q $queryPairFile -o $outputOriginalResult
#     fi

#     for gamma in {0..16}
#       do
#       for eta in {0..15}
#         do
#         #WH2H preprocessing
#         if [ $? == 0 ]; then     
#           echo "faild"
#         else
#           echo "success------------------------------------------------"
#           ./bin/run -w -d -r $gamma -e $eta -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult
#         fi
#         #query test by WH2H and O(nlogn) ST table
#         if [ $? == 0 ]; then     
#           echo "faild"
#         else
#           echo "success------------------------------------------------"
#           ./bin/run -q -d -n $numOfQueries -i $outputIndexFile -s $outputSTTable -q $queryPairFile -f $outputResult
#         fi
#       done
#     done
#     cutt=cut+1
#     if [ $cut == 95 ];then
#       cutt=7
#       echo $cutt
#     fi
#     queryPairFile=${queryPairFile//cut$cut/cut$cutt}
#     HFVerticesFile=${HFVerticesFile//cut$cut/cut$cutt}
#     outputOriginalResult=${outputOriginalResult//cut$cut/cut$cutt}
#     outputResult=${outputResult//cut$cut/cut$cutt}
#   done
#   dayy=day+1
#   queryPairFile=${queryPairFile//day$day/day$dayy}
#   HFVerticesFile=${HFVerticesFile//day$day/day$dayy}
#   outputOriginalResult=${outputOriginalResult//day$day/day$dayy}
#   outputResult=${outputResult//day$day/day$dayy}
#   echo $queryPairFile
#   echo $HFVerticesFile
# done

updatedGraphFile=${updatedGraphFile//Haikou/Chengdu}
outputIndexFile=${outputIndexFile//Haikou/Chengdu}
outputSTTable=${outputSTTable//Haikou/Chengdu}
outputBlockSTTable=${outputBlockSTTable//Haikou/Chengdu}
outputQueryCost=${outputQueryCost//Haikou/Chengdu}
outputOriginalResult=${outputOriginalResult//Haikou/Chengdu}
outputResult=${outputResult//Haikou/Chengdu}
queryPairFile=${queryPairFile//Haikou/Chengdu}
HFVerticesFile=${HFVerticesFile//Haikou/Chengdu}

queryPairFile=${queryPairFile//day30/day0}
HFVerticesFile=${HFVerticesFile//day30/day0}
outputOriginalResult=${outputOriginalResult//day30/day0}
outputResult=${outputResult//day30/day0}
queryPairFile=${queryPairFile//cut95/cut7}
HFVerticesFile=${HFVerticesFile//cut95/cut7}
outputOriginalResult=${outputOriginalResult//cut95/cut7}
outputResult=${outputResult//cut95/cut7}

# for day in {0..29}
#   do
#   for cut in {7..95}
#     do
    #H2H preprocessing
    ./bin/run -h -d -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult

#     #query test by H2H and O(nlogn) ST table
#     if [ $? == 0 ]; then     
#       echo "faild"
#     else
#       echo "success------------------------------------------------"

#       ./bin/run -q -d -n $numOfQueries -i $outputIndexFile -s $outputSTTable -q $queryPairFile -o $outputOriginalResult
#     fi

#     for gamma in {0..16}
#       do
#       for eta in {0..15}
#         do
#         #WH2H preprocessing
#         if [ $? == 0 ]; then     
#           echo "faild"
#         else
#           echo "success------------------------------------------------"
#           ./bin/run -w -d -r $gamma -e $eta -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult
#         fi
#         #query test by WH2H and O(nlogn) ST table
#         if [ $? == 0 ]; then     
#           echo "faild"
#         else
#           echo "success------------------------------------------------"
#           ./bin/run -q -d -n $numOfQueries -i $outputIndexFile -s $outputSTTable -q $queryPairFile -f $outputResult
#         fi
#       done
#     done
#     cutt=cut+1
#     if [ $cut == 95 ];then
#       cutt=7
#       echo $cutt
#     fi
#     queryPairFile=${queryPairFile//cut$cut/cut$cutt}
#     HFVerticesFile=${HFVerticesFile//cut$cut/cut$cutt}
#     outputOriginalResult=${outputOriginalResult//cut$cut/cut$cutt}
#     outputResult=${outputResult//cut$cut/cut$cutt}
#   done
#   dayy=day+1
#   queryPairFile=${queryPairFile//day$day/day$dayy}
#   HFVerticesFile=${HFVerticesFile//day$day/day$dayy}
#   outputOriginalResult=${outputOriginalResult//day$day/day$dayy}
#   outputResult=${outputResult//day$day/day$dayy}
#   echo $queryPairFile
#   echo $HFVerticesFile
# done

updatedGraphFile=${updatedGraphFile//Chengdu/NewYork}
outputIndexFile=${outputIndexFile//Chengdu/NewYork}
outputSTTable=${outputSTTable//Chengdu/NewYork}
outputBlockSTTable=${outputBlockSTTable//Chengdu/NewYork}
outputQueryCost=${outputQueryCost//Chengdu/NewYork}
outputOriginalResult=${outputOriginalResult//Chengdu/NewYork}
outputResult=${outputResult//Chengdu/NewYork}
queryPairFile=${queryPairFile//Chengdu/NewYork}
HFVerticesFile=${HFVerticesFile//Chengdu/NewYork}

queryPairFile=${queryPairFile//day30/day0}
HFVerticesFile=${HFVerticesFile//day30/day0}
outputOriginalResult=${outputOriginalResult//day30/day0}
outputResult=${outputResult//day30/day0}
queryPairFile=${queryPairFile//cut95/cut7}
HFVerticesFile=${HFVerticesFile//cut95/cut7}
outputOriginalResult=${outputOriginalResult//cut95/cut7}
outputResult=${outputResult//cut95/cut7}

# for day in {0..29}
#   do
#   for cut in {7..95}
#     do
    #H2H preprocessing
    ./bin/run -h -d -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult

#     #query test by H2H and O(nlogn) ST table
#     if [ $? == 0 ]; then     
#       echo "faild"
#     else
#       echo "success------------------------------------------------"

#       ./bin/run -q -d -n $numOfQueries -i $outputIndexFile -s $outputSTTable -q $queryPairFile -o $outputOriginalResult
#     fi

#     for gamma in {0..16}
#       do
#       for eta in {0..15}
#         do
#         #WH2H preprocessing
#         if [ $? == 0 ]; then     
#           echo "faild"
#         else
#           echo "success------------------------------------------------"
#           ./bin/run -w -d -r $gamma -e $eta -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult
#         fi
#         #query test by WH2H and O(nlogn) ST table
#         if [ $? == 0 ]; then     
#           echo "faild"
#         else
#           echo "success------------------------------------------------"
#           ./bin/run -q -d -n $numOfQueries -i $outputIndexFile -s $outputSTTable -q $queryPairFile -f $outputResult
#         fi
#       done
#     done
#     cutt=cut+1
#     if [ $cut == 95 ];then
#       cutt=7
#       echo $cutt
#     fi
#     queryPairFile=${queryPairFile//cut$cut/cut$cutt}
#     HFVerticesFile=${HFVerticesFile//cut$cut/cut$cutt}
#     outputOriginalResult=${outputOriginalResult//cut$cut/cut$cutt}
#     outputResult=${outputResult//cut$cut/cut$cutt}
#   done
#   dayy=day+1
#   queryPairFile=${queryPairFile//day$day/day$dayy}
#   HFVerticesFile=${HFVerticesFile//day$day/day$dayy}
#   outputOriginalResult=${outputOriginalResult//day$day/day$dayy}
#   outputResult=${outputResult//day$day/day$dayy}
#   echo $queryPairFile
#   echo $HFVerticesFile
# done


updatedGraphFile=${updatedGraphFile//NewYork/Florida}
outputIndexFile=${outputIndexFile//NewYork/Florida}
outputSTTable=${outputSTTable//NewYork/Florida}
outputBlockSTTable=${outputBlockSTTable//NewYork/Florida}
outputQueryCost=${outputQueryCost//NewYork/Florida}
outputOriginalResult=${outputOriginalResult//NewYork/Florida}
outputResult=${outputResult//NewYork/Florida}
queryPairFile=${queryPairFile//NewYork/Florida}
HFVerticesFile=${HFVerticesFile//NewYork/Florida}

    ./bin/run -h -d -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult

updatedGraphFile=${updatedGraphFile//Florida/California}
outputIndexFile=${outputIndexFile//Florida/California}
outputSTTable=${outputSTTable//Florida/California}
outputBlockSTTable=${outputBlockSTTable//Florida/California}
outputQueryCost=${outputQueryCost//Florida/California}
outputOriginalResult=${outputOriginalResult//Florida/California}
outputResult=${outputResult//Florida/California}
queryPairFile=${queryPairFile//Florida/California}
HFVerticesFile=${HFVerticesFile//Florida/California}

    ./bin/run -h -d -g $updatedGraphFile -v $HFVerticesFile -i $outputIndexFile -s $outputSTTable -b $outputBlockSTTable -q $outputQueryCost -o $outputOriginalResult -f $outputResult
