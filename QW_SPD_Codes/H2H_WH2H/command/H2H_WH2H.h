//#pragma GCC optimize(3)

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <time.h>
#include "../src/command.h"
#include "../src/node.h"
#include "../src/tree.h"
#include "../src/edge.h"
#include "../src/graph.h"

using namespace std;

/*
Graph is preprocessed by H2H or WH2H in this command.
The output is the index, ST table, and Query Cost of H2H or WH2H.
*/

#define MAXLISTSIZE 10000 //the maximum length of inverted list, if the inverted list out of flow during runtime, you may enlarge the number.
#define MAXDEGREE 100		//the maximum degree of the graph

namespace command {
	class H2H_WH2H :public Command {
		int main(int argc, char* argv[]) {

			bool WH2H{ false };
			bool directed{ false };
			float gamma{ 0.1 };
			int BLOCKSIZE{ 5 };

			char graph_file[255] = "";
			char HFVertex_file[255] = "";
			char outputH2HIndex[255] = "";
			char outputST[255] = "";
			char outputST_block[255] = "";
			char outputQueryCost[255] = "";

			if (argc < 12) {
				exit_with_help();
			}
			if (argv[1][1] == 'w')
				WH2H = true;
			for (int i = 2; i < argc; i++) {
				if (argv[i][0] != '-') break;
				if (argv[i][1] == 'd') {
					i++;
					directed = true;
				}
				if (++i >= argc)
					exit_with_help();
				switch (argv[i - 1][1]) {
				case 'r':
					gamma = atof(argv[i]);
					break;
				case 'e':
					BLOCKSIZE = atoi(argv[i]);
					break;
				case 'g':
					strcpy(graph_file, argv[i]);
					break;
				case 'v':
					strcpy(HFVertex_file, argv[i]);
					break;
				case 'i':
					strcpy(outputH2HIndex, argv[i]);
					break;
				case 's':
					strcpy(outputST, argv[i]);
					break;
				case 'b':
					strcpy(outputST_block, argv[i]);
					break;
				case 'q':
					strcpy(outputQueryCost, argv[i]);
					break;
				default:
					exit_with_help();
				}
			}
			assert(gamma > 0 && gamma < 1);
			graph vertexProcess(graph_file, MAXLISTSIZE, gamma, BLOCKSIZE, MAXDEGREE, directed);

			if (WH2H) {
				vertexProcess.readHFVertexData(HFVertex_file);//Read high frequency point data when using WH2H.
			}

			//remove vertices by H2H or WH2H
			time_t t_start = clock();
			if (WH2H) {
				if (directed) {
					vertexProcess.removeVertices_WH2H_directed();
				}
				else {
					vertexProcess.removeVertices_WH2H();
				}
			}
			else {
				if (directed) { //directed graph
					vertexProcess.removeVertices_H2H_directed();
				}
				else { //undirected graph
					vertexProcess.removeVertices_H2H();
				}
			}
			time_t t_end = clock();
			cout << "The time of removing vertices is: " << (t_end - t_start) / 1000 << " ms." << endl;

			//connect all nodes by child-parent relationship
			t_start = clock();
			tree H2HTree(vertexProcess.nodes, vertexProcess.level, directed);
			t_end = clock();
			cout << "The time of creating tree is: " << (t_end - t_start) / 1000 << " ms." << endl;
			

			//Create index by tree structure.
			t_start = clock();
			H2HTree.createLabel(vertexProcess.nodes, NULL, directed, vertexProcess.dis_dijkstra);
			t_end = clock();
			cout << "The time of creating index is: " << (t_end - t_start) / 1000 << " ms." << endl;


			//Calculate tree height, tree-width, average tree-height, and average tree-width.
			cout << "============================================" << endl;
			cout << "The Tree Height is: " << H2HTree.getHeight(H2HTree.getRoot()) << endl;
			int avgHeight{ 0 };
			int maxWidth{ 0 };
			int avgWidth{ 0 };
			for (int i = 0; i < vertexProcess.numOfVertices; i++) {
				avgHeight += vertexProcess.nodes[i].height;
				if (maxWidth < vertexProcess.nodes[i].pos.size()) {
					maxWidth = vertexProcess.nodes[i].pos.size();
				}
				if (vertexProcess.nodes[i].child.size() != 0) {
					avgWidth += vertexProcess.nodes[i].pos.size();
				}
			}

			cout << "The tree-width is: " << maxWidth << endl;
			cout << "The avg tree-height is: " << (float)avgHeight / vertexProcess.numOfVertices << endl;
			cout << "The avg tree-width is: " << (float)avgWidth / vertexProcess.numOfVertices << endl;

			//Generate ST table for the search of LCA in O(1) time. 
			t_start = clock();
			// We generate two ST table, the space complexity is O(nlogn) for the first one,
			// and O(n) for the second one. We finally use the first one for queries because 
			// it is a little faster than the second one.
			H2HTree.createST(vertexProcess.nodes);
			t_end = clock();
			cout << "The time of creating ST table is: " << (t_end - t_start) / 1000 << " ms." << endl;
			//H2HTree.search_LCA_ST_On(0, 2);
			//H2HTree.search_LCA_ST_On(0, 22);

			//for (int i = 0; i < vertexProcess.numOfVertices; i++) {
			//	for (int j = 0; j < vertexProcess.numOfVertices; j++) {
			//		int dis1 = vertexProcess.dis_dijkstra[i][j];
			//		int dis2 = H2HTree.query_ST_directed(vertexProcess.nodes, i, j);
			//		int dis3= H2HTree.query_ST_On_directed(vertexProcess.nodes, i, j);
			//		if (dis1 != dis3 || dis2 != dis1) {
			//			cout << "wrong!" << endl;
			//			cout << i << " " << j << " " << dis2 << " " << dis3 << endl;
			//			//break;
			//		}
			//		/*if (dis1 != dis3) {
			//			break;
			//		}*/
			//	}
			//}


			//Output index, ST table, and query cost.
			// H2HTree.calQueryCost(vertexProcess.nodes);
			// H2HTree.outputQueryCost(outputQueryCost, vertexProcess.nodes);
			H2HTree.save_labels(outputH2HIndex, vertexProcess.nodes, directed);
			H2HTree.save_st_binary(outputST, outputST_block);

			return 1;
		}
	};
}