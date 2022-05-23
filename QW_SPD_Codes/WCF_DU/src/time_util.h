/*
 * @Descripttion: 
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-09-13 11:03:52
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-09-13 11:04:18
 */
/*
 * @Description:  compute time util functions
 * @Author: wanjingyi
 * @Date: 2021-01-26 11:37:47
 * @LastEditTime: 2021-09-09 16:49:09
 */
#pragma once
#ifndef TIME_UTIL_H
#define TIME_UTIL_H

#include <sys/time.h>
#include <map>

namespace time_util {

	//double GetCurrentTimeSec() { return 0; }
	double GetCurrentTimeSec() {
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec + tv.tv_usec * 1e-6;
	}
}

#endif