/*
 * @Description: variables used statement
 * @Author: wanjingyi
 * @Date: 2021-01-26 11:37:47
 * @LastEditTime: 2021-10-14 10:11:42
 */

#pragma once
#ifndef PARAS_H
#define PARAS_H
#include<cstdint> 
#include <limits>

#define DIVISION_FACTOR 1
#define HF_DIVIDION 10000
#define AMPLIFY_FACTOR 1

namespace SP_Constants {
	bool DIRECTED_FLAG = false;
	bool WEIGHTED_FLAG = false;
	int numOfVertices = 0;
	int numOfEdges = 0;
	const extern int INF_WEIGHT = std::numeric_limits<int>::max() / 3;
}

#endif