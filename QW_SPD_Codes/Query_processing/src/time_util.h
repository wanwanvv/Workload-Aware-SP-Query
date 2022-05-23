/*
 * @Description:  compute time util functions
 * @Author: wanjingyi
 * @Date: 2021-01-26 11:37:47
 * @LastEditTime: 2021-10-14 10:15:16
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

#ifndef TIME_UTIL_H
#define TIME_UTIL_H

#include <sys/time.h>
#include<map>

namespace time_util {

	//double GetCurrentTimeSec() { return 0; }

	
	double GetCurrentTimeSec() {
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec + tv.tv_usec * 1e-6;
	}
	
}
#endif