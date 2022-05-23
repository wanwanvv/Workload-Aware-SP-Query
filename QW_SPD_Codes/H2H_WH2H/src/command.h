
#ifndef _COMMAND_H
#define _COMMAND_H

#include <math.h>

class Command {
public:
  virtual ~Command() {};
  virtual int main(int argc, char* argv[]) = 0;
  void exit_with_help()
  {
    printf("Command Usage: (the first parameter specifies the procedure be executed)\n\n");
    printf("-------------------------------------------------------------------\n");
    printf("Convert graph into undirected graph and test the correctness:\n");
    printf("ShortestPathTest.exe -c -g [graphFileName] -o [outputGraphFileName]\n");
    printf("-------------------------------------------------------------------\n");
    printf("H2H preprocessing:\n\n");
    printf("ShortestPathTest.exe -h (-d) -g [convertedGraphFileName] -i [outputIndexFileName] \n\
                    -s [outputSTFileName] -b [outputSTBlockFileName] -q [outputQueryCostFileName]\n\n");

    printf("WH2H preprocessing:\n\n");
    printf("ShortestPathTest.exe -w (-d) -r [gammaValue] -e [etaValue] -g [convertedGraphFileName] -v [HighFrequencyVerticesFileName] -i [outputIndexFileName] \n\
                    -s [outputSTFileName] -b [outputSTBlockFileName] -q [outputQueryCostFileName]\n\n");

    printf("-------------------------------------------------------------------\n");
    printf("Query processing: (make sure preprocessing is finished)\n\n");
    printf("ShortestPathTest.exe -q (-d) -n [numberOfTestPairs]-i [iutputIndexFileName] \n\
                    -s [iutputSTFileName] / -b [iutputSTBlockFileName] -q [iutputQueryPairFileName]\n\
                    Note that just use \"-s\" or \"-b\" as ST table");
    exit(1);
  }
};

#endif // _COMMAND_H
