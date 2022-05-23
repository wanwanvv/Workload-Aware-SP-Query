/*
 * @Description: variables used statement
 * @Author: wanjingyi
 * @Date: 2021-01-26 11:37:47
 * @LastEditTime: 2021-09-18 11:55:45
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

namespace Processing{
        template<typename T>
	    struct  calcCoefficient//struct to store the relative coefficients of the importance of a node
        {
            T deg_mult=0;//coefficient of node degree
            T freq_mult=0;//coefficient of node query frequency
            T cov_mult=0;//coefficient of coverage
            T dep_mult=0;//coefficient of the distance from current node to last choose node denote the uniformity
            T bet_mult=0;//coefficient of the betwenness of each node
            bool is_deg_mult=false;
            bool is_freq_mult=false;
            bool is_cov_mult=false;
            bool is_dep_mult=false;
            bool is_bet_mult=false;
            calcCoefficient():deg_mult(0),freq_mult(0),cov_mult(0),dep_mult(0),is_deg_mult(false),is_freq_mult(false),is_cov_mult(false),is_dep_mult(false),is_bet_mult(false){}
            calcCoefficient(T de,T frq,T cov,T dep,T bet,bool is_deg, bool is_freq, bool is_cov, bool is_dep, bool is_bet):deg_mult(de),freq_mult(frq),cov_mult(cov),dep_mult(dep),bet_mult(bet),is_deg_mult(is_deg),is_freq_mult(is_freq),is_cov_mult(is_cov),is_dep_mult(is_dep),is_bet_mult(is_bet){}
        };
}

#endif