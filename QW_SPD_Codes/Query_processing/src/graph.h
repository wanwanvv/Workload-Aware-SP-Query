/*
 * @Descripttion: define the graph class and util functions of graph
 * @version: 1.0
 * @Author: wanjingyi
 * @Date: 2021-02-09 15:37:10
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-10-15 09:16:58
 */

/*
This file referenced:
An Experimental Study on Hub Labeling based Shortest Path Algorithms [Experiments and Analyses]

Authors: Ye Li, Leong Hou U, Man Lung Yiu, Ngai Meng Kou
Contact: yb47438@umac.mo
Affiliation: University of Macau

The MIT License (MIT)

Copyright (c) 2016 University of Macau

All rights reserved. 
*/

#pragma once
#ifndef GRAPH_H
#define GRAPH_H

#include<vector>
#include<fstream>
#include <algorithm>
#include<cstdint> 
#include "paras.h"
#include<iostream>
#include<assert.h>
#include<map>


using namespace std;

//typedef unsigned int NodeID;
typedef int NodeID;
typedef long long EdgeID;
//typedef unsigned int EdgeID;
typedef pair<NodeID, NodeID> Edge;
//typedef double EdgeWeight;
typedef  int EdgeWeight;
typedef pair<NodeID, EdgeWeight> NodeEdgeWeightPair;



#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG 

#endif