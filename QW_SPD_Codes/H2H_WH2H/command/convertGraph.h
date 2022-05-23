#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <cstring>
#include <string>
#include "../src/command.h"
#include "../src/edge.h"

using namespace std;

/*
Any graph is converted to an undirected graph in this command.
The output is an undirected graph used as the input of H2H or WH2H.
*/

namespace command {
	class convertGraph :public Command {
	public:

		int main(int argc, char* argv[]) {

			char graph_file[255] = "";
			char output_file[255] = "";

			if (argc < 6) {
				exit_with_help();
			}
			for (int i = 2; i < argc; i++) {
				if (argv[i][0] != '-') break;
				if (++i >= argc)
					exit_with_help();
				switch (argv[i - 1][1]) {
				case 'g':
					strcpy(graph_file, argv[i]);
					break;
				case 'o':
					strcpy(output_file, argv[i]);
					break;
				default:
					exit_with_help();
				}
			}

			edge temp;
			ifstream ifs(graph_file);
			int numOfNodes{ 0 };
			int numOfEdges{ 0 };

			//get number of nodes from the graph
			for (; ifs >> temp.u >> temp.v >> temp.weight;) {
				numOfEdges++;
				if (temp.u > numOfNodes) {
					numOfNodes = temp.u;
				}
				if (temp.v > numOfNodes) {
					numOfNodes = temp.v;
				}
			}
			numOfNodes++;

			//save edges by the order of the first node
			vector<list<edge>> es(numOfNodes);
			ifs.close();
			ifs.open(graph_file);

			//read edges from the graph
			for (; ifs >> temp.u >> temp.v >> temp.weight;) {
				temp.u;
				temp.v;
				es[temp.u].push_back(temp);
			}

			//make all the edges unique
			for (int i = 0; i < numOfNodes; i++) {
				es[i].sort([](edge p1, edge p2) {return p1.v < p2.v; });
				es[i].unique([](edge p1, edge p2) {return p1.v == p2.v; });
			}

			//convert a directed graph to an undirected graph
			for (int i = 0; i < numOfNodes; i++) {
				for (list<edge>::iterator it = es[i].begin(); it != es[i].end(); it++) {
					list<edge>::iterator itv = es[it->v].begin();
					for (; itv != es[it->v].end(); itv++) {
						if (itv->v == it->u)
							break;
					}
					if (itv == es[it->v].end()) {
						temp.u = it->v;
						temp.v = it->u;
						temp.weight = it->weight;
						es[it->v].push_back(temp);
					}
				}
			}

			//test whether there is an edge (u,v) where u==v
			for (int i = 0; i < numOfNodes; i++) {
				for (list<edge>::iterator it = es[i].begin(); it != es[i].end(); it++) {
					if (it->u == it->v) {
						es[i].erase(it);
					}
				}
			}

			//output graph
			ofstream ofs(output_file);
			ofs << numOfNodes << " " << numOfEdges << endl;
			for (int i = 0; i < numOfNodes; i++) {
				if (es[i].size() == 0) {
					i = i;
				}
				es[i].sort([](const edge e1, const edge e2) {return e1.v < e2.v; });
				for (list<edge>::iterator it = es[i].begin(); it != es[i].end(); it++) {
					ofs << it->u << " " << it->v << " " << it->weight << endl;
				}
			}
			return 1;
		}
	};
}