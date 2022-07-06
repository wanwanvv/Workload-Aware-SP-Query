# Worlolad-Aware-SP-Query
Source codes of **ICDE-2022** research paper "Workload-Aware Shortest Path Distance Querying in Road Networks".

# Introduction
论文《负载感知的路网最短路径查询》中用C++实现的主要算法和流程框架<br>
<p align="center">
<img src=".\img\framework.png" height = "200" alt="" align=center />
<br><br>
<b>Figure 1.</b> The architecture of WCF.
</p>

# Problem
最短路径距离查询是指给定一个网络上的起点和终点，要求返回这两点间的最短路径距离。其作为一项基本操作，是许多基于位置服务的应用中的**基础构建模块**，例如GPS导航、打车服务、物流运输、自动驾驶等。随着移动互联网的发展，大量的并发查询请求在不同时间段内组成不同的查询负载，给服务器资源带来了巨大的压力。**如何在尽可能小的资源占用的情况下，快速高效地处理最短路径查询负载，面临巨大挑战。**

# Motivation
查询负载具有的时空特性：
* 空间偏斜性：二八定律，少部分点被大量查询(hot vertices)；
* 时间局部性查询具有时间连续性

<p align="center">
<img src=".\img\motivation.png" height = "200" alt="" align=center />
<br><br>
<b>Figure 2.</b> Motivation.
</p>

# Environment
- Ubuntu 18.04 AMD-64G<br>
- gcc-7.5 g++-7.5<br>

# PLL_WHP
### 算法介绍
* 作为对比的State-of-the-art算法PLL
  + 节点排序扑(结构特征例如度大小、路径重要性、中介中心度等)
  + Dijkstra搜索
  + 添加标签索引
  + 两步剪枝
      
* 论文中基于此改进的wPLL算法
  + 节点排序(路网结构特征和查询负载特征)
  + Dijkstra搜索
  + 添加标签索引(高频点更早被处理，减小索引大小)
  + 两步剪枝
      
### Usage
    详情请见.\PLL_WHP\README

## H2H_WH2H
### 算法介绍
    * 作为对比的State-of-the-art算法H2H
      基于树结构的最短距离查询方法，按照节点度什序进行节点删除，并通过权值更新和填边保证删除点后剩余图的正确性，构建完树结构后按照从上到下的方式计算每个节点的标签索引，然后利用标签索引进行最短路径查询
      + 树分解/节点删除
      + 连接形成树结构
      + 索引构建
      + 计算LCA(运用ST表获得O(1)时间复杂度)
      + 最短距离查询
      
    * 论文中基于此改进的WH2H算法
      节点重要性综合考虑*查询频率*和*节点结构重要性*，使得形成的树结构中高频点更靠近根结点，减小高频点查询开销。
      + 节点排序
      + 树分解
      + 分块技术
      
<p align="center">
<img src=".\img\tree.png" height = "300" alt="" align=center />
<br><br>
<b>Figure 3.</b> 树分解.
</p>

### Usage
    详情请见.\H2H_WH2H\README

## WCF_DU
### 实现算法
    * 基于层级结构的混合方法WCF
      综合利用查询频率和路网结构提取出层次结构，在高层结构上利用基于节点排序的方法来减小高频查询点的索引大小，从而减小总的查询开销，然后在包含大部分低频查询点的低层次结构上利用基于树分解的方法来减小索引构建开销。WCF包含以下三个步骤：
      + 构建层级结构(core-forest，core为距离保留图，forest-边缘为多棵树结构)
      + 索引构建(由于森林结构互斥可用多线程并行加速)
      + 最短路径距离查询

<p align="center">
<img src=".\img\wcf.png" height = "200" alt="" align=center />
<br><br>
<b>Figure 4.</b> Core-forest结构.
</p>

### Usage
    详情请见.\WCF_DU\README

## 数据结构
#### 图结构存储
定义在graph.h中，可以读取有权/无权和有向/无向图，提供去除自环和重边功能，为了优化空间和访问速度，将所有边(u,v,w)按照起始点u排序依次存储在一维数组edges中，用二维邻接表adj只存储每个节点最后一条邻接边在edges中的索引位置

<p align="center">
<img src=".\img\graph.png" height = "200" alt="" align=center />
<br><br>
<b>Figure 5.</b> Graph example.
</p>

#### 标签索引
定义在labels.h中，通过struct来实现，支持各种类型的图，构建阶段主要利用vector，查询阶段采用指针，并进行字节对齐，和__mm_prefetch等优化手段。为了提高查询速度，提供缓存功能。

#### 自定义堆
定义在heap.h中，模板参数log_k-几叉树实现, k_t-键值key类型, id_t-值value类型，提供update()增加或者更新key对应的value、top()返回堆顶元素、extract_min()获得堆顶元素并删除该元素。<br>
数据结构：用一个一维vector数组elements来存储所有的元素，然后用一个一维vector数组positions来存储键值为key的元素在elements中的位置。<br>
用于实现节点排序功能。<br>

### 最短路径查询基础算法
定义在graph_search.h中，支持各种类型图，提供基础的BFS、Dijkstra、双向BFS以及双向Dijkstra，并提供多个点的并行计算。

## 文件存储
主要使用二进制文件来存储索引减小磁盘开销

## Exp_scripts
  复现论文中实验结果的批量Linus shell运行脚本，数据集放在Datasets文件夹下。
