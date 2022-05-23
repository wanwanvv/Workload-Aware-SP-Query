/*
 * @Description: 
 * @Author: wanjingyi
 * @Date: 2021-01-26 14:31:17
 * @LastEditTime: 2022-02-15 11:45:53
 */
#include <iostream>
#include <stdlib.h> 
#include <fstream>

#include "command/IndexProcessing.h"
#include "command/constructionLabeling.h"
#include "command/queryDistance.h"
#include "command/tipIndexing.h"
#include "command/generateData.h"

using namespace std;

void exit_with_help(){
    printf("Command Usage: (the first parameter specifies the procedure be executed)\n");
    printf("-------------------------------------------------------------------\n");
    printf("(1) -x: Indexing by based ordering schemes:\n");
    printf("\tpll_whp_run -x -d [directedGraphFlag] -w [weightedGraphFlag] -s [experimentFlag] -m [orderingSchemes] [-a [s_beta]] -g [graphFileName] -e [outputLabelFileName] -i [outputLabelSizeFileName] [-f [outputResultDirName]] [-q [queryFreqFileName]] \n");
    printf("-------------------------------------------------------------------\n");
    printf("(2) -p: Indexing by custom ordering scheme(WHP):\n");
    printf("\tpll_whp_run -p -d [directedGraphFlag] -w [weightedGraphFlag] -s [experimentFlag] -g [graphFileName] -e [exportLabelFileName] [-q [queryFreqFileName]] [-i [labelSizeFileName]] [-f [outputResultDirName]] [-b [betwennessFileName]] [-k [freqParameter]] [-l [betParameter]] \n");
    printf("-------------------------------------------------------------------\n");
    printf("(3) -q: Query shortest path distance using 2-hop labels:\n");
    printf("\tpll_whp_run -q -d [directedGraphFlag] -w [weightedGraphFlag] -s [experimentFlag] -p [dataType] -l [inputLabelFileName] [-t [queryPairFileName]] [-i [queryFreqFileName]] [-a [outputResultDirName]] [-f [outputDistanceFileName]] [-z [inputLabelSizeFileName]] [-u [multiQuery]] [-n [numOfQuery]] [-x [shuffleFlag]] \n");
    printf("-------------------------------------------------------------------\n");
    exit(1);
}

/** The main program. */
int main(int argc,char* argv[])
{
    int opt = 'm';
    if (argc > 1){
        switch (argv[1][1]){
            case 'x':
                opt = 'x';
                break;
             case 'q':
                opt = 'q';
                break;
            case 'p':
                opt = 'p';
                break;
            case 't':
                opt = 't';
                break;
            case 'r':
                opt = 'r';
                break;
            default:
                exit_with_help();
                break;
        }
    }
    if (argc == 1) exit_with_help();

    Command* m = NULL;
    int result = 0;
    switch (opt) {
        case 'x':
            m = new command::IndexProcessing();
            break;
        case 'p':
            m = new command::ConstructLabel();
            break;
        case 'q':
            m = new command::QueryProcessing();
            break;
        case 't':
            m=new command::TIPIndexing();
            break;
        case 'r':
            m=new command::GenerateData();
            break;
        default:
            break;
    }

    if (m != NULL) {
        result = m->main(argc, argv);
        delete m;
    } else {
        exit_with_help();
    }
    return result;

    return 0;
}//main
