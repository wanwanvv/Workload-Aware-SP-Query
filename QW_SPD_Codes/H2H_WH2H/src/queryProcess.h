#pragma once
//#pragma GCC optimize("O0")
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <cassert>
#include <math.h>
#include <malloc.h>
#include <algorithm>

using namespace std;

//This class is used for query test. Read index, ST table and query pairs from files.
//And test the query performance.
class queryProcess
{
public:
  //Construct function. Read H2H or WH2H index.
  queryProcess(const char* inputH2HIndex, const bool& directed)
  {
    //Read index.
    ifstream ifsIndex(inputH2HIndex, ios::binary);

    //Read number of nodes.
    ifsIndex.read((char*)&numOfNodes, sizeof(numOfNodes));
    cout << numOfNodes << " vertices." << endl;

    //Read dis array.
    if (directed) {
      dis_in = (int_fast32_t**)malloc(numOfNodes * sizeof(int_fast32_t*));
      dis_out = (int_fast32_t**)malloc(numOfNodes * sizeof(int_fast32_t*));
      /*dis_in = (int_fast32_t**)memalign(32, numOfNodes * sizeof(int_fast32_t*));*/
      /*dis_out = (int_fast32_t**)memalign(32, numOfNodes * sizeof(int_fast32_t*));*/
      assert(dis_in && dis_out);
      for (int v = 0; v < numOfNodes; ++v)
      {
        int disSize = 0;
        ifsIndex.read((char*)&disSize, sizeof(disSize));
        dis_in[v] = (int_fast32_t*)malloc(disSize * sizeof(int_fast32_t));
        /*dis_in[v] = (int_fast32_t*)memalign(32, disSize * sizeof(int_fast32_t));*/
        assert(dis_in[v]);
        if (disSize > 0)
        {
          for (int i = 0; i < disSize; ++i)
          {
            int distance;
            ifsIndex.read((char*)&distance, sizeof(distance));
            dis_in[v][i] = distance;
          }
        }
        ifsIndex.read((char*)&disSize, sizeof(disSize));
        dis_out[v] = (int_fast32_t*)malloc(disSize * sizeof(int_fast32_t));
        /*dis_out[v] = (int_fast32_t*)memalign(32, disSize * sizeof(int_fast32_t));*/
        assert(dis_out[v]);
        if (disSize > 0)
        {
          for (int i = 0; i < disSize; ++i)
          {
            int distance;
            ifsIndex.read((char*)&distance, sizeof(distance));
            dis_out[v][i] = distance;
          }
        }
      }
    }
    else {
      dis = (int_fast32_t**)malloc(numOfNodes * sizeof(int_fast32_t*));
      /*dis = (int_fast32_t**)memalign(32, numOfNodes * sizeof(int_fast32_t*));*/
      assert(dis);
      for (int v = 0; v < numOfNodes; ++v)
      {
        int disSize = 0;
        ifsIndex.read((char*)&disSize, sizeof(disSize));
        dis[v] = (int_fast32_t*)malloc(disSize * sizeof(int_fast32_t));
        /*dis[v] = (int_fast32_t*)memalign(32, disSize * sizeof(int_fast32_t));*/
        assert(dis[v]);
        if (disSize > 0)
        {
          for (int i = 0; i < disSize; ++i)
          {
            int distance;
            ifsIndex.read((char*)&distance, sizeof(distance));
            dis[v][i] = distance;
          }
        }
      }
    }
    
    //Read pos array.
    pos = (int_fast32_t**)malloc(numOfNodes * sizeof(int_fast32_t*));
    /*pos = (int_fast32_t**)memalign(32, numOfNodes * sizeof(int_fast32_t*));*/
    assert(pos);
    for (int v = 0; v < numOfNodes; ++v)
    {
      int posSize;
      ifsIndex.read((char*)&posSize, sizeof(posSize));
      pos[v] = (int_fast32_t*)malloc((posSize + 2) * sizeof(int_fast32_t));
      /*pos[v] = (int_fast32_t*)memalign(32, (posSize + 2) * sizeof(int_fast32_t));*/
      assert(pos[v]);
      if (posSize > 0)
      {
        for (int i = 1; i <= posSize; ++i)
        {
          int position;
          ifsIndex.read((char*)&position, sizeof(position));
          pos[v][i] = position;
        }
      }
      pos[v][posSize + 1] = INT32_MAX;

      //The pos array in leaves is not useful.
      if (posSize == 0)
      {
        pos[v][0] = 0;
      }
      pos[v][0] = pos[v][posSize]; //We set the position of the node itself as the first element.
    }
    ifsIndex.close();
  }

  ~queryProcess()
  {
    // free(st_dep);
    // free(st_pos);
    // for (int i = 0; i < 2 * numOfNodes + 1; i++) {
    //   free(sparseTable[i]);
    // }
    // free(sparseTable);
    // for (int i = 0; i < numOfNodes; i++) {
    //   free(dis[i]);
    //   free(pos[i]);
    // }
    // free(pos);
    // free(dis);
  }

  //Read O(nlogn) ST table.
  void read_ST_logn(const char* inputST)
  {
    int sparseSize = 2 * numOfNodes + 1;
    sparseTable = (int_fast32_t**)malloc(sparseSize * sizeof(int_fast32_t*));
    /*sparseTable = (int_fast32_t**)memalign(32, sparseSize * sizeof(int_fast32_t*));*/
    assert(sparseTable);
    int length = log2(2 * numOfNodes + 1) + 1;
    for (int i = 0; i < (2 * numOfNodes + 1); i++)
    {
      sparseTable[i] = (int_fast32_t*)malloc(length * sizeof(int_fast32_t));
      /*sparseTable[i] = (int_fast32_t*)memalign(32, length * sizeof(int_fast32_t));*/
      assert(sparseTable[i]);
    }

    ifstream ifsST(inputST, ios::binary);
    assert(ifsST.is_open());
    for (int i = 0; i < sparseSize; i++)
    {
      for (int j = 0; j < length; j++)
      {
        int curSparseTable;
        ifsST.read((char*)&curSparseTable, sizeof(curSparseTable));
        sparseTable[i][j] = curSparseTable;
      }
    }

    //Read st_pos array.
    int posLength = numOfNodes;
    st_pos = (int_fast32_t*)malloc(posLength * sizeof(int_fast32_t));
    /*st_pos = (int_fast32_t*)memalign(32, posLength * sizeof(int_fast32_t));*/
    assert(st_pos);
    int posMAX = -1;
    for (int i = 0; i < posLength; i++)
    {
      int curST_Pos;
      ifsST.read((char*)&curST_Pos, sizeof(curST_Pos));
      if (posMAX < curST_Pos)
        posMAX = curST_Pos;
      st_pos[i] = curST_Pos;
    }

    //Read st_dep array.
    int depLength = numOfNodes;
    st_dep = (int_fast32_t*)malloc(depLength * sizeof(int_fast32_t));
    /*st_dep = (int_fast32_t*)memalign(32, depLength * sizeof(int_fast32_t));*/
    assert(st_dep);
    dep_MAX = 0; //Record position of the node with maximum depth.
    for (int i = 0; i < depLength; i++)
    {
      int curST_dep;
      ifsST.read((char*)&curST_dep, sizeof(curST_dep));
      if (curST_dep > st_dep[dep_MAX])
      {
        dep_MAX = i;
      }
      st_dep[i] = curST_dep;
    }
    ifsST.close();

    //We calculate all the log() value which will be used.
    logArray = (int_fast32_t*)malloc((posMAX + 2) * sizeof(int_fast32_t));
    /*logArray = (int_fast32_t*)memalign(32, (posMAX + 2) * sizeof(int_fast32_t));*/
    for (int i = 0; i < (posMAX + 2); i++)
    {
      logArray[i] = log2(i);
    }
  }

  //Read O(n) ST table.
  void read_ST_n(const char* inputST_block)
  {
    EulerLen = 2 * numOfNodes + 1;
    blockSize = log2(EulerLen) / 2;

    //Read globalST table.
    ifstream ifsST(inputST_block, ios::binary);
    assert(ifsST.is_open());
    globalST = (int_fast32_t**)malloc((EulerLen / blockSize + 1) * sizeof(int_fast32_t*));
    /*globalST = (int_fast32_t**)memalign(32, (EulerLen / blockSize + 1) * sizeof(int_fast32_t*));*/
    assert(globalST);
    int length = log2(EulerLen / blockSize + 1) + 1;
    for (int i = 0; i < (EulerLen / blockSize + 1); i++)
    {
      globalST[i] = (int_fast32_t*)malloc(length * sizeof(int_fast32_t));
      /*globalST[i] = (int_fast32_t*)memalign(32, length * sizeof(int_fast32_t));*/
      assert(globalST[i]);
    }
    for (int i = 0; i < (EulerLen / blockSize + 1); i++)
    {
      for (int j = 0; j < length; j++)
      {
        int curSparseTable;
        ifsST.read((char*)&curSparseTable, sizeof(curSparseTable));
        globalST[i][j] = curSparseTable;
      }
    }

    //Read insideST table.
    insideST = (int_fast32_t***)malloc((EulerLen / blockSize) * sizeof(int_fast32_t**));
    /*insideST = (int_fast32_t***)memalign(32, (EulerLen / blockSize) * sizeof(int_fast32_t**));*/
    int insideSize = log2(blockSize) + 1;
    for (int i = 0; i < EulerLen / blockSize - 1; i++)
    {
      insideST[i] = (int_fast32_t**)malloc(blockSize * sizeof(int_fast32_t*));
      /*insideST[i] = (int_fast32_t**)memalign(32, blockSize * sizeof(int_fast32_t*));*/
      for (int j = 0; j < blockSize; j++)
      {
        insideST[i][j] = (int_fast32_t*)malloc(insideSize * sizeof(int_fast32_t));
        /*insideST[i][j] = (int_fast32_t*)memalign(32, insideSize * sizeof(int_fast32_t));*/
      }
    }
    insideST[EulerLen / blockSize - 1] = (int_fast32_t**)malloc((blockSize + 1) * sizeof(int_fast32_t*));
    /*insideST[EulerLen / blockSize - 1] = (int_fast32_t**)memalign(32, (blockSize + 1) * sizeof(int_fast32_t*));*/
    insideSize = log2(blockSize + 1) + 1;
    for (int j = 0; j < (blockSize + 1); j++)
    {
      insideST[EulerLen / blockSize - 1][j] = (int_fast32_t*)malloc(insideSize * sizeof(int_fast32_t));
      /*insideST[EulerLen / blockSize - 1][j] = (int_fast32_t*)memalign(32, (insideSize + 1) * sizeof(int_fast32_t));*/
    }
    insideSize = log2(blockSize) + 1;
    for (int i = 0; i < EulerLen / blockSize - 1; i++)
    {
      for (int j = 0; j < blockSize; j++)
      {
        for (int k = 0; k < insideSize; k++)
        {
          int curInsideST;
          ifsST.read((char*)&curInsideST, sizeof(curInsideST));
          insideST[i][j][k] = curInsideST;
        }
      }
    }
    insideSize = log2((blockSize + 1)) + 1;
    for (int j = 0; j < (blockSize + 1); j++)
    {
      for (int k = 0; k < insideSize; k++)
      {
        int curInsideST;
        ifsST.read((char*)&curInsideST, sizeof(curInsideST));
        insideST[EulerLen / blockSize - 1][j][k] = curInsideST;
      }
    }

    //Read st_pos array.
    int posLength = numOfNodes;
    st_pos = (int_fast32_t*)malloc(posLength * sizeof(int_fast32_t));
    /*st_pos = (int_fast32_t*)memalign(32, posLength * sizeof(int_fast32_t));*/
    assert(st_pos);
    int posMAX = -1;
    for (int i = 0; i < posLength; i++)
    {
      int curST_Pos;
      ifsST.read((char*)&curST_Pos, sizeof(curST_Pos));
      if (posMAX < curST_Pos)
        posMAX = curST_Pos;
      st_pos[i] = curST_Pos;
    }

    //Read st_dep array.
    int depLength = numOfNodes;
    st_dep = (int_fast32_t*)malloc(depLength * sizeof(int_fast32_t));
    /*st_dep = (int_fast32_t*)memalign(32, depLength * sizeof(int_fast32_t));*/
    assert(st_dep);
    dep_MAX = 0; //Record position of the node with maximum depth.
    for (int i = 0; i < numOfNodes + 1; i++)
    {
      int curST_dep;
      ifsST.read((char*)&curST_dep, sizeof(curST_dep));
      if (curST_dep > st_dep[dep_MAX])
      {
        dep_MAX = i;
      }
      st_dep[i] = curST_dep;
    }
    ifsST.close();

    //We calculate all the log() value which will be used.
    logArray = (int_fast32_t*)malloc((posMAX + 2) * sizeof(int_fast32_t));
    /*logArray = (int_fast32_t*)memalign(32, (posMAX + 2) * sizeof(int_fast32_t));*/
    for (int i = 0; i < (posMAX + 2); i++)
    {
      logArray[i] = log2(i);
    }
  }

  //Calculate Log2() quickly.
  unsigned int log2_fast(unsigned int x)
  {
    /*unsigned int ret;
    __asm__ __volatile__("bsrl %1, %%eax"
      : "=a"(ret)
      : "m"(x));
    return ret;*/
    return log2(x);
  }

  //Read query pairs from file. Randomly shuffle queries. And increase the number of queries to num.
  void queriesGenarator(const int& num, const char* queryPair_file)
  {
    TESTNUM = num;
    queries.reserve(TESTNUM * 1.1);
    int m, n, fre;
    ifstream ifsQuery(queryPair_file);
    for (; ifsQuery >> m >> n >> fre;)
    {
      if (m != n)
      {
        for (int i = 0; i < fre; i++)
        {
          queries.push_back(make_pair(m, n));
        }
      }
    }
    int queriesPairSize = queries.size();
    if (queriesPairSize == 0)
    {
      cout << "The queries are wrong!" << endl;
      return;
    }

    //Increase queries.
    while (queries.size() < TESTNUM)
    {
      for (int i = 0; i < queriesPairSize; i++)
      {
        queries.push_back(queries[i]);
      }
    }
    
    //Randomly shuffle.
    std::random_shuffle(queries.begin(), queries.end());

    TESTNUM = queries.size();
    warmup = queriesPairSize;
  }

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Test the performance by O(nlogn) ST table on indirected graph.
  void test_ST_logn()
  {
    int testSum = queries.size() - warmup;
    int totalNum = queries.size();
    double result = 0;
    for (int j = 0; j < 3; j++)
    {
      //warmup
      for (int i = 0; i < warmup; ++i)
      {
        int s = queries[i].first;
        int t = queries[i].second;
        query_logArray_logn(s, t);
      }
      //test
      double qTime = clock();
      for (int i = warmup; i < totalNum; ++i)
      {
        int s = queries[i].first;
        int t = queries[i].second;
        query_logArray_logn(s, t);
      }
      double endtime = clock();
      result += (endtime - qTime) / (double)(testSum);
    }
    result /= 3.0;
    cout << "The average response time is : " << result << " microseconds" << endl;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Test the performance by O(nlogn) ST table on directed graph.
  void test_ST_logn_directed()
  {
    int testSum = queries.size() - warmup;
    int totalNum = queries.size();
    double result = 0;
    for (int j = 0; j < 3; j++)
    {
      //warmup
      for (int i = 0; i < warmup; ++i)
      {
        int s = queries[i].first;
        int t = queries[i].second;
        query_logArray_logn_directed(s, t);
      }
      //test
      double qTime = clock();
      for (int i = warmup; i < totalNum; ++i)
      {
        int s = queries[i].first;
        int t = queries[i].second;
        query_logArray_logn_directed(s, t);
      }
      double endtime = clock();
      result += (endtime - qTime) / (double)(testSum);
    }
    result /= 3.0;
    cout << "The average response time is : " << result << " microseconds" << endl;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Test the performance by O(n) ST table.
  void test_ST_n()
  {
    int testSum = queries.size() - warmup;
    int totalNum = queries.size();
    double result = 0;
    for (int j = 0; j < 3; j++)
    {
      //warmup
      for (int i = 0; i < warmup; ++i)
      {
        int s = queries[i].first;
        int t = queries[i].second;
        query_ST_On(s, t);
      }
      //test
      double qTime = clock();
      for (int i = warmup; i < totalNum; ++i)
      {
        int s = queries[i].first;
        int t = queries[i].second;
        query_ST_On(s, t);
      }
      double endtime = clock();
      result += (endtime - qTime) / (double)(testSum);
    }
    result /= 3.0;
    cout << "The average response time is : " << result << " microseconds" << endl;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Test the performance by O(n) ST table.
  void test_ST_n_directed()
  {
    int testSum = queries.size() - warmup;
    int totalNum = queries.size();
    double result = 0;
    for (int j = 0; j < 3; j++)
    {
      //warmup
      for (int i = 0; i < warmup; ++i)
      {
        int s = queries[i].first;
        int t = queries[i].second;
        query_ST_On_directed(s, t);
      }
      //test
      double qTime = clock();
      for (int i = warmup; i < totalNum; ++i)
      {
        int s = queries[i].first;
        int t = queries[i].second;
        query_ST_On_directed(s, t);
      }
      double endtime = clock();
      result += (endtime - qTime) / (double)(testSum);
    }
    result /= 3.0;
    cout << "The average response time is : " << result << " microseconds" << endl;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Search LCA by O(nlogn) ST table using log2_fast().
  int search_LCA_ST(int x, int y)
  {
    int ans = 0;
    if (st_pos[x] > st_pos[y])
      swap(x, y);
    int l1 = st_pos[x], l2 = st_pos[y];
    int len = log2_fast(l2 - l1 + 1);
    int calLen = 1 << len;
    int_fast32_t st_l1_len = sparseTable[l1][len];
    int_fast32_t st_l2_callen = sparseTable[l2 - calLen + 1][len];
    ans = st_l1_len;
    if (st_dep[st_l1_len] > st_dep[st_l2_callen])
      ans = st_l2_callen;
    return ans;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Search LCA by O(nlogn) ST table using logArray.
  int search_LCA_ST_logArray(int x, int y)
  {
    int ans = 0;
    if (st_pos[x] > st_pos[y])
      swap(x, y);
    int l1 = st_pos[x], l2 = st_pos[y];
    int len = logArray[l2 - l1 + 1];
    int calLen = 1 << len;
    int_fast32_t st_l1_len = sparseTable[l1][len];
    int_fast32_t st_l2_callen = sparseTable[l2 - calLen + 1][len];
    ans = st_l1_len;
    if (st_dep[st_l1_len] > st_dep[st_l2_callen])
      ans = st_l2_callen;
    return ans;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Search LCA by O(n) ST table.
  int search_LCA_ST_On(int x, int y)
  {
    int ans = 0;
    if (st_pos[x] > st_pos[y])
      swap(x, y);
    int l1 = st_pos[x], l2 = st_pos[y];
    int xBlock = l1 / blockSize;
    int yBlock = l2 / blockSize;
    int l1_Location = l1; //The location of l1 inside block.
    if (xBlock != 0)
    {
      l1_Location = l1 % (xBlock * blockSize);
    }
    int l2_Location = l2; //The location of l2 inside block.
    if (yBlock != 0)
    {
      l2_Location = l2 % (yBlock * blockSize);
    }
    if (xBlock == yBlock)
    { //If the two block are the same, return directly.
      return search_insideST_Block(xBlock, l1_Location, l2_Location);
    }
    int ans1 = search_globalST_Block(xBlock + 1, yBlock - 1);             //Search LCA in the block between two blocks.
    int ans2 = search_insideST_Block(xBlock, l1_Location, blockSize - 1); //Search LCA in block xBlock.
    int ans3 = search_insideST_Block(yBlock, 0, l2_Location);             //Search LCA in block yBlock.
    ans = (st_dep[ans1] < st_dep[ans2]) ? ans1 : ans2;
    ans = (st_dep[ans] < st_dep[ans3]) ? ans : ans3;
    return ans;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Search LCA in globalST where l1 < l2.
  int search_globalST_Block(const int& l1, const int& l2)
  {
    //Process search between blocks where l1 = l2 + 1.
    if (l1 > l2)
    {
      return dep_MAX;
    }
    if (l1 == l2)
    {
      return globalST[l1][0];
    }
    int ans = 0;
    int len = logArray[l2 - l1 + 1];
    if (st_dep[globalST[l1][len]] < st_dep[globalST[l2 - (1 << len) + 1][len]])
      ans = globalST[l1][len];
    else
      ans = globalST[l2 - (1 << len) + 1][len];
    return ans;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Search LCA in insideST table where l1 < l2.
  int search_insideST_Block(const int& block, const int& l1, const int& l2)
  {
    //Process search between blocks where l1 = l2 + 1.
    if (l1 > l2)
    {
      return dep_MAX;
    }
    if (l1 == l2)
    {
      return insideST[block][l1][0];
    }
    int ans = 0;
    int len = logArray[l2 - l1 + 1];
    if (st_dep[insideST[block][l1][len]] < st_dep[insideST[block][l2 - (1 << len) + 1][len]])
      ans = insideST[block][l1][len];
    else
      ans = insideST[block][l2 - (1 << len) + 1][len];
    return ans;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Query by O(nlogn) ST table using log2_fast() in indirected graph.
  int query(int& m, int& n)
  {
    if (m == n) {
      return 0;
    }
    int ans = search_LCA_ST(m, n);
    if (ans == m)
    { //If the ancestor is m, return distance directly.
      return dis[n][pos[ans][0]];
    }
    if (ans == n)
    { //If the ancestor is n, return distance directly.
      return dis[m][pos[ans][0]];
    }
    //return 0;
    int distance{ INT32_MAX };
    int_fast32_t* curPos = pos[ans];
    int_fast32_t* dis_n = dis[n];
    int_fast32_t* dis_m = dis[m];
    for (int i = 1;; i++)
    {
      if (curPos[i] == INT32_MAX)
      { //curPos结束位置是INT32_MAX标志
        break;
      }
      int curDis = dis_m[curPos[i]] + dis_n[curPos[i]];
      if (distance > curDis)
      {
        distance = curDis;
      }
    }
    return distance;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Query by O(nlogn) ST table using log2_fast() in directed graph.
  int query_directed(int& m, int& n)
  {
    if (m == n) {
      return 0;
    }
    int ans = search_LCA_ST(m, n);
    if (ans == m)
    { //If the ancestor is m, return distance directly.
      return dis_in[n][pos[ans][0]];
    }
    if (ans == n)
    { //If the ancestor is n, return distance directly.
      return dis_out[m][pos[ans][0]];
    }
    //return 0;
    int distance{ INT32_MAX };
    int_fast32_t* curPos = pos[ans];
    int_fast32_t* dis_n = dis_in[n];
    int_fast32_t* dis_m = dis_out[m];
    for (int i = 1;; i++)
    {
      if (curPos[i] == INT32_MAX)
      { //curPos结束位置是INT32_MAX标志
        break;
      }
      int curDis = dis_m[curPos[i]] + dis_n[curPos[i]];
      if (distance > curDis)
      {
        distance = curDis;
      }
    }
    return distance;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Query by O(nlogn) ST table using logArray in indirected graph.
  int query_logArray_logn(int& m, int& n)
  {
    if (m == n) {
      return 0;
    }
    int ans = search_LCA_ST_logArray(m, n);

    //Return distance in O(1) time if there is an ancestor-child relationship.
    if (ans == m)
    {
      return dis[n][pos[ans][0]];
    }
    if (ans == n)
    {
      return dis[m][pos[ans][0]];
    }
    int distance{ INT32_MAX };
    int_fast32_t* curPos = pos[ans];
    int_fast32_t* dis_n = dis[n];
    int_fast32_t* dis_m = dis[m];
    for (int i = 1;; i++)
    {
      if (curPos[i] == INT32_MAX)
      {
        break;
      }
      int curDis = dis_m[curPos[i]] + dis_n[curPos[i]];
      if (distance > curDis)
      {
        distance = curDis;
      }
    }
    return distance;
  };
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Query by O(nlogn) ST table using logArray in directed graph.
  int query_logArray_logn_directed(int& m, int& n)
  {
    if (m == n) {
      return 0;
    }
    int ans = search_LCA_ST_logArray(m, n);

    //Return distance in O(1) time if there is an ancestor-child relationship.
    if (ans == m)
    {
      return dis_in[n][pos[ans][0]];
    }
    if (ans == n)
    {
      return dis_out[m][pos[ans][0]];
    }
    int distance{ INT32_MAX };
    int_fast32_t* curPos = pos[ans];
    int_fast32_t* dis_n = dis_in[n];
    int_fast32_t* dis_m = dis_out[m];
    for (int i = 1;; i++)
    {
      if (curPos[i] == INT32_MAX)
      {
        break;
      }
      int curDis = dis_m[curPos[i]] + dis_n[curPos[i]];
      if (distance > curDis)
      {
        distance = curDis;
      }
    }
    return distance;
  };
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Query by O(n) ST table in indirected graph.
  int query_ST_On(int& m, int& n)
  {
    if (m == n) {
      return 0;
    }
    int ans = search_LCA_ST_On(m, n);

    //Return distance in O(1) time if there is an ancestor-child relationship.
    if (ans == m)
    {
      return dis[n][pos[ans][0]];
    }
    if (ans == n)
    {
      return dis[m][pos[ans][0]];
    }
    int distance{ INT32_MAX };
    int_fast32_t* curPos = pos[ans];
    int_fast32_t* dis_n = dis[n];
    int_fast32_t* dis_m = dis[m];
    for (int i = 1;; i++)
    {
      if (curPos[i] == INT32_MAX)
      {
        break;
      }
      int curDis = dis_m[curPos[i]] + dis_n[curPos[i]];
      if (distance > curDis)
      {
        distance = curDis;
      }
    }
    return distance;
  }
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize("O0")
  //Query by O(n) ST table in directed graph.
  int query_ST_On_directed(int& m, int& n)
  {
    if (m == n) {
      return 0;
    }
    int ans = search_LCA_ST_On(m, n);

    //Return distance in O(1) time if there is an ancestor-child relationship.
    if (ans == m)
    {
      return dis_in[n][pos[ans][0]];
    }
    if (ans == n)
    {
      return dis_out[m][pos[ans][0]];
    }
    int distance{ INT32_MAX };
    int_fast32_t* curPos = pos[ans];
    int_fast32_t* dis_n = dis_in[n];
    int_fast32_t* dis_m = dis_out[m];
    for (int i = 1;; i++)
    {
      if (curPos[i] == INT32_MAX)
      {
        break;
      }
      int curDis = dis_m[curPos[i]] + dis_n[curPos[i]];
      if (distance > curDis)
      {
        distance = curDis;
      }
    }
    return distance;
  }
#pragma GCC pop_options

private:
  int numOfNodes = 0;
  int TESTNUM;                    //number of queries
  int warmup;                     //number of queries for warmup
  vector<pair<int, int>> queries; //query pairs
  int_fast32_t** dis;             //the dis array in Index for indirected graph
  int_fast32_t** dis_in;          //the dis_in array in Index for directed graph
  int_fast32_t** dis_out;         //the dis_out array in Index for directed graph
  int_fast32_t** pos;             //the pos array in Index
  int_fast32_t** sparseTable;     //ST table
  int_fast32_t* st_pos;           //Save the position of each node in Euler sequence.
  int_fast32_t* st_dep;           //Save the height of each node in tree decomposition.
  int_fast32_t* logArray;         //log2() array

  int_fast32_t** globalST;  //ST table between blocks
  int_fast32_t*** insideST; //ST table inside block
  int EulerLen;             //length of Euler sequence
  int blockSize;            //The size of blocks in O(n) ST table.
  int dep_MAX;
};