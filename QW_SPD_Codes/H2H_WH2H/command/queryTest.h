//#pragma GCC optimize("O0")
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <cassert>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "../src/command.h"
#include "../src/queryProcess.h"

/*
Read H2H or WH2H index, O(nlogn) or O(n) ST table, and query pairs from files.
Then test the query performance of index.
The output is the average query time.
*/

using namespace std;

namespace command {
	class queryTest :public Command {
		int main(int argc, char* argv[]) {
			int TESTNUM = 100000;
			char inputH2HIndex[255] = "";
			char inputST[255] = "";
			char inputST_Block[255] = "";
			char queryPair_file[255] = "";
			bool ST = true;
			bool directed{ false }; //directed graph or indirected graph
			if (argc < 9) {
				exit_with_help();
			}
			for (int i = 2; i < argc; i++) {
				if (argv[i][0] != '-') break;
				if (argv[i][1] == 'd') {
					i++;
					directed = true;
				}
				if (++i >= argc)
					exit_with_help();
				switch (argv[i - 1][1]) {
				case 'n':
					TESTNUM = atoi(argv[i]);
					break;
				case 'i':
					strcpy(inputH2HIndex, argv[i]);
					break;
				case 's':
					strcpy(inputST, argv[i]);
					break;
				case 'b':
					strcpy(inputST_Block, argv[i]);
					ST = false;
					break;
				case 'q':
					strcpy(queryPair_file, argv[i]);
					break;
				default:
					exit_with_help();
				}
			}

			queryProcess process(inputH2HIndex, directed);

			//Read O(nlogn) ST table or O(n) ST table.
			if (ST)
				process.read_ST_logn(inputST);
			else
				process.read_ST_n(inputST_Block);

			// for (int i = 0; i < 25; i++) {
			// 	for (int j = 0; j < 25; j++) {
			// 		int dis2 = process.query_logArray_logn_directed(i, j);
			// 		int dis3 = process.query_ST_On_directed(i, j);
			// 		if (dis2 != dis3 || dis2 == INT32_MAX || dis3 == INT32_MAX) {
			// 			cout << "wrong!" << endl;
			// 			cout << i << " " << j << " " << dis2 << " " << dis3 << endl;
			// 			//break;
			// 		}
			// 	}
			// }

			//Read query pairs from files, shuffle queries, and increase queries to TESTNUM.
 			process.queriesGenarator(TESTNUM, queryPair_file);

			//performance test
			if(directed){
				if (ST)
					process.test_ST_logn_directed();
				else
					process.test_ST_n_directed();
			}
			else{
				if (ST)
				process.test_ST_logn();
				else
					process.test_ST_n();
			}
			return 1;
		}
	};
}