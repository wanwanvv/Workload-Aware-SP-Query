/*
 * @Descripttion: 
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-11-03 15:49:04
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-11-05 09:11:44
 */

#ifndef _COMMAND_TIP_INDEXING
#define _COMMAND_TIP_INDEXING
#include <cstring>
#include <iostream>
#include "../command.h"
#include "../src/time_util.h"
#include "../src/paras.h"
#include "../src/timeInterval.h"

using namespace std;

namespace command{
    class TIPIndexing:public Command{
        public:
            void exit_with_help(){
                printf("Usage:\n");
                printf("\tsspexp_run -u -d [directedGraphFlag] -w [weightedGraphFlag] -s [specialFlag] -m [indexingSchemes] -o [OrderingSchemes] -g [graphFileName] \n -e [exportLabelFileName] [-q [query_freq_file]] [-h [HFPoint_file]] [-p label_list_file] [-i [label_size_file]] [-f size_analysis_file\n [-b betweenness_filename] [-c coverage_filename] [-r high_frequency rate default-5%%] [-j -k -l -u -v coeffient of params(degree,frequency,betwenness,coverage,depth;0~10)]\n");
                printf("-------------------------------------------------------------------\n");
            }

            int main(int argc, char* argv[]){
                char graphFileName[255] = "";
                char labelDirName[255] = "";
                char labelSizeDirName[255]="";
                char pointFreqDirName[255] = "";
                char analyDirName[255] = "";
                char parasFileName[255] = "";
                char updateIndexFileName[255] = "";
                char queryDataDirName[255] = "";
                char orderWeightFileName[255]=""; //output prder weight filename
                char betwennessFileName[255]=""; //a betwenness file computed before
                char coverageFileName[255]="";//a coverageness file computed before
                int day_slice_cnt=0;
                int t_directed_flag = 0;
                int t_weighted_flag = 0;
                int t_special_flag=0;
                int graphType=0; //0-full connected 1-u is start and t is destination
                double updateThreshold=0;

                for(int i = 2; i < argc; i++){
                    if(argv[i][0] != '-') break;
                    if(++i >= argc){
                        exit_with_help();
                        std::cout<<"error:++i >= argc"<<std::endl;
                    }
                    switch (argv[i-1][1]){
                        case 'd':
                            t_directed_flag = atoi(argv[i]);
                            break;
                        case 'w':
                            t_weighted_flag = atoi(argv[i]);
                            break;
                        case 's':
                            t_special_flag = atoi(argv[i]);
                            break;
                        case 'n':
                            day_slice_cnt = atoi(argv[i]);
                            break;
                        case 'u':
                            graphType = atoi(argv[i]);
                            break;
                        case 'g':
                            strcpy(graphFileName,argv[i]);
                            std::cout<<"graphFileName = "<<graphFileName<<std::endl;
                            break;
                        case 'l':
                            strcpy(labelDirName, argv[i]);
                            std::cout<<"labelDirName = "<<labelDirName<<std::endl;
                            break;
                        case 'a':
                            strcpy(analyDirName, argv[i]);
                            if(*analyDirName!='\0'){
                                std::cout<<"analyDirName = "<<analyDirName<<std::endl;
                            }
                            break;
                        case 'p':
                            strcpy(parasFileName,argv[i]);
                            std::cout<<"parasFileName = "<<parasFileName<<std::endl;
                            break;
                        case 'h':
                            strcpy(pointFreqDirName, argv[i]);
                            std::cout<<"pointFreqDirName = "<<pointFreqDirName<<std::endl;
                            break;
                        case 'q':
                            strcpy(queryDataDirName, argv[i]);
                            std::cout<<"queryDataDirName = "<<queryDataDirName<<std::endl;
                            break;
                        case 'b':
                            strcpy(labelSizeDirName, argv[i]);
                            std::cout<<"labelSizeDIrName = "<<labelSizeDirName<<std::endl;
                            break;
                        case 'k':
                            strcpy(updateIndexFileName,argv[i]);
                            std::cout<<"updateIndexFileName = "<<updateIndexFileName<<std::endl;
                            break;
                        case 'o':
                            strcpy(orderWeightFileName, argv[i]);
                            std::cout<<"orderWeightFileName = "<<orderWeightFileName<<endl;
                            break;
                        case 'f':
                            strcpy(betwennessFileName, argv[i]);
                            std::cout<<"betwennessFileName = "<<betwennessFileName<<endl;
                            break;
                        case 'c':
                            strcpy(coverageFileName, argv[i]);
                            std::cout<<"coverageFileName = "<<coverageFileName<<endl;
                            break;
                        case 'e':
                            updateThreshold=atof(argv[i]);
                            std::cout<<"updateThreshold = "<<updateThreshold<<std::endl;
                            break;
                        default:
                            exit_with_help();
                            std::cout<<"default!"<<std::endl;
                            break;
                    }
                }

                if (t_directed_flag == 1)
                    DIRECTED_FLAG = true;
                if (t_weighted_flag == 1)
                    WEIGHTED_FLAG = true;

                if(t_special_flag==0){//GTP
                    std::cout<<"t_special_flag=0 using GTP......"<<endl;
                    TimeInterval(graphFileName,pointFreqDirName,parasFileName,
                                                analyDirName,labelDirName,labelSizeDirName,
                                                betwennessFileName,coverageFileName,orderWeightFileName,
                                                graphType,updateThreshold,queryDataDirName);
                }else if(t_special_flag==1){//RL-TIP
                    std::cout<<"t_special_flag=1 using RL-TIP...... "<<endl;
                    TimeInterval(graphFileName,pointFreqDirName,updateIndexFileName,
                                                parasFileName,analyDirName,labelDirName,
                                                labelSizeDirName,betwennessFileName,coverageFileName,
                                                orderWeightFileName,day_slice_cnt,graphType,
                                                queryDataDirName);
                }

                return 0;    
            }//main ends
    };
}

#endif