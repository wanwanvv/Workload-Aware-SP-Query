/*
 * @Description: 
 * @Author: wanjingyi
 * @Date: 2021-01-26 14:31:17
 * @LastEditTime: 2021-11-03 15:33:12
 */
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "./command/deleteNode.h"
#include "./command/integratedIndexing.h"
#include "./command/parallel_indexingLabeling.h"
#include "./command/integratedQuery.h"
#include "./command/tipIndexing.h"

using namespace std;

void exit_with_help(){
    printf("Command Usage: (the first parameter specifies the procedure be executed)\n");
    printf("-------------------------------------------------------------------\n");
    printf("wcf_du_run -e -d [directedGraphFlag] -w [weightedGraphFlag] -s [experimentFlag] -g [graphFileName] \n -h [queryFreqFileName] -o [outputDirName] -m [removeMaxDegree] [-i [deleteModel]] [-r [hfRate]] \n [-t [deleteNodeRate]] [-b [numOfBfsHop]] [-l [initialDegreeOfLocalSearch]] \n");
    printf("wcf_du_run -z -d [directedGraphFlag] -w [weightedGraphFlag] -o [coreOrderingModel] -s [experimentFlag] \n -i [indexingModel] -m [removeMaxDegree] -k [freqParameter] -l [betParameter] [-n [bestRemoveDegree]] \n-c [numOfThreads] -g [graphFileName] -q [queryFreqFileName] -t [interFileOutputPrefix] -a [outputResultPrefix] -e [outputLabelFileName] -b [outputLabelSizeFileName] [-f [inputOrderFileName]]\n");
    printf("\twcf_du_run -l -d [directedGraphFlag] -w [weightedGraphFlag] -p [dataType] -b [batchQueryFlag] -s [experimentFlag] \n [-x [shuffleFlag]] [-n [numOfQuery]] [-u [multiQuery]] -m [indexingModel] \n -i [outputLabelFileName] -v [lcaFileName] -y [isDeletedFileName] [-a [outputResultPrefix]] [-q [queryDataFileName]] [-f [outputDistanceFileName]] [-j [matchPrefix]]\n");
    printf("bin/wcf_du_run -d [directedGraphFlag] -w [weightedGraphFlag] -g [graphFileName] [-o [orderingModel]] [-b [batchSize]] -a [outputDirName] [-n [numOfThreads]] \n");
    printf("wcf_du_run -u -d [directedGraphFlag] -w [weightedGraphFlag] -s [updateModel] -o [coreOrderingModel] [-x [shuffleFlag]]\n [-c [numOfThreads]] -e [multiQuery] -m [indexingModel] -y [dayIndex] -n [numOfTimeSlice] -t [greedyUpdateThreshold] -g [graphFileName] [-p [parametersFileName]] \n-h [queryFreqDirName] -q [queryDataDirName] -l [outputLabelDirName] -b [outputLabelSizeDirName] -i [interFileOutputDirName] [-a [outputResultPrefix]] [-k [updateIndexFileName]]\n");
    exit(1);
}

/** The main program. */
int main(int argc,char* argv[])
{
    int opt = 'm';
    if (argc > 1){
        switch (argv[1][1]){
            case 'e':
                opt = 'e';
                break;
            case 'z':
                opt = 'z';
                break;
            case 'l':
                opt = 'l';
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
        case 'e':
            m = new command::DeleteNode();
            break;
        case 'z':
            m=new command::IntegrationIndexing();
            break;
        case 'l':
            m=new command::IntegrationQuery();
            break;
        case 't':
            m=new command::TIPIndexing();
            break;
        case 'r':
            m=new command::ParallelIndexingLabeling();
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
