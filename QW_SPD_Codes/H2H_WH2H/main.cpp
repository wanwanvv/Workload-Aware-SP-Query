#include <iostream>
#include <fstream>
using namespace std;

#include "command/convertGraph.h"
#include "command/H2H_WH2H.h"
#include "command/queryTest.h"

void exit_with_help()
{
	printf("Command Usage: (the first parameter specifies the procedure be executed)\n\n");
	printf("-------------------------------------------------------------------\n");
	printf("Convert graph into undirected graph and test the correctness:\n");
	printf("ShortestPathTest.exe -c -g [graphFileName] -o [outputGraphFileName]\n\n");
	printf("-------------------------------------------------------------------\n");
	printf("H2H preprocessing:\n\n");
	printf("ShortestPathTest.exe -h -g [convertedGraphFileName] -i [outputIndexFileName] \n\
			-s [outputSTFileName] -b [outputSTBlockFileName] -q [outputQueryCostFileName]\n\n");

	printf("WH2H preprocessing:\n\n");
	printf("ShortestPathTest.exe -w -r [gammaValue] -e [etaValue] -g [convertedGraphFileName] \n\
			-v [HighFrequencyVerticesFileName] -i [outputIndexFileName] \n\
			-s [outputSTFileName] -b [outputSTBlockFileName] -q [outputQueryCostFileName]\n\n");

	printf("-------------------------------------------------------------------\n");
	printf("Query processing: (make sure preprocessing is finished)\n\n");
	printf("ShortestPathTest.exe -q -n [numberOfTestPairs]-i [inputIndexFileName] \n\
			-s [inputSTFileName] / -b [inputSTBlockFileName] -q [inputQueryPairFileName]\n\
			Note that just use \"-s\" or \"-b\" as ST table\n\n");
	exit(1);
}

/** The main program. */

int main(int argc, char* argv[])
{
	// The program is controlled by command-line arguments. The order of those
	// arguments is important. The first argument specifies the Command-
	// class that is used.

	int opt = 'm';
	if (argc > 1) {
		switch (argv[1][1]) {
		case 'c':
			opt = 'c';
			break;
		case 'h':
			opt = 'h';
			break;
		case 'w':
			opt = 'h';
			break;
		case  'q':
			opt = 'q';
			break;

		default:
			exit_with_help();
			break;
		}
	}

	Command* m = NULL;
	int result = 0;
	switch (opt) {

		/****************Graph data transform***************************************************/
	case 'c':
		m = new command::convertGraph();
		break;

		/****************Preprocessing**********************************************************/

		/* ************* *
		* H2H or WH2H*
		* ************* */
	case 'h':
		m = new command::H2H_WH2H();
		break;
		/***********************Query Testing************************************************/

		/* ******************** *
		* Query Testing        *
		* ******************** */

	case 'q':
		m = new command::queryTest();
		break;

	default:
		break;

	}
	if (m != NULL) {
		result = m->main(argc, argv);
		delete m;
	}
	else {
		exit_with_help();
	}
	return result;
}