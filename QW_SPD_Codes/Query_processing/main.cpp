/*
 * @Descripttion: 
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-10-14 08:34:37
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-10-16 17:27:11
 */
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "./command/queryDistance.h"
#include "./command/integratedQuery.h"

using namespace std;

void exit_with_help(){
    printf("Command Usage: (the first parameter specifies the procedure be executed)\n");
    printf("-------------------------------------------------------------------\n");
    printf("Parameter explanation:\n");
    printf("-------------------------------------------------------------------\n");
    printf("Examples:\n");
    printf("-------------------------------------------------------------------\n");
    exit(1);
}

/** The main program. */
int main(int argc,char* argv[])
{
    int opt = 'm';
    if (argc > 1){
        switch (argv[1][1]){
             case 'q':
                opt = 'q';
                break;
             case 'l':
                opt = 'l';
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
        case 'q':
            m = new command::QueryProcessing();
            break;
        case 'l':
            m=new command::IntegrationQuery();
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
}
