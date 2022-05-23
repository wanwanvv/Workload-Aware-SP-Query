/*
 * @Descripttion: 
 * @version: 
 * @Author: Wan Jingyi
 * @Date: 2021-10-29 14:31:09
 * @LastEditors: Wan Jingyi
 * @LastEditTime: 2021-11-03 16:15:23
 */
#ifndef _COMMAND_TIP_INDEXING
#define _COMMAND_TIP_INDEXING
#include <cstring>
#include <iostream>
#include "../command.h"
#include "../src/paras.h"
#include "../src/timeInterval.h"

#define numOfVertices SP_Constants::numOfVertices
#define numOfEdges SP_Constants::numOfEdges
#define WEIGHTED_FLAG SP_Constants::WEIGHTED_FLAG  
#define DIRECTED_FLAG SP_Constants::DIRECTED_FLAG

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
                char outputDirName[255] = "";
                char parasFileName[255] = "";
                char updateIndexFileName[255] = "";
                char givenOrderFileName[255]="";
                char queryDataDirName[255] = "";
                int day_slice_cnt=0;
                int t_directed_flag = 0;
                int t_weighted_flag = 0;
                int t_special_flag=0;
                int t_index_model=0;
                int t_ordering_flag=0;
                int num_threads=0;
                int index_type=0;//0-WHP, 1-WCF
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
                        case 'm':
                            t_index_model = atoi(argv[i]);
                            break;
                        case 'n':
                            day_slice_cnt = atoi(argv[i]);
                            break;
                        case 'o':
                            t_ordering_flag = atoi(argv[i]);
                            break;
                        case 'u':
                            graphType = atoi(argv[i]);
                            break;
                        case 'c':
                           num_threads = atoi(argv[i]);
                           std::cout<<"num_threads="<<num_threads<<std::endl;
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
                        case 'i':
                            strcpy(outputDirName,argv[i]);
                            std::cout<<"outputDirName="<<outputDirName<<std::endl;
                            break;
                        case 'b':
                            strcpy(labelSizeDirName, argv[i]);
                            std::cout<<"labelSizeDIrName = "<<labelSizeDirName<<std::endl;
                            break;
                        case 'k':
                            strcpy(updateIndexFileName,argv[i]);
                            std::cout<<"updateIndexFileName = "<<updateIndexFileName<<std::endl;
                            break;
                        case 'f':
                            strcpy(givenOrderFileName,argv[i]);
                            std::cout<<"givenOrderFileName="<<givenOrderFileName<<std::endl;
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
                                                outputDirName,givenOrderFileName,t_index_model,
                                                t_ordering_flag,num_threads,graphType,
                                                updateThreshold,queryDataDirName);
                }else if(t_special_flag==1){//RL-TIP
                    std::cout<<"t_special_flag=1 using RL-TIP...... "<<endl;
                    TimeInterval(graphFileName,pointFreqDirName,updateIndexFileName,
                                                parasFileName,analyDirName,labelDirName,
                                                labelSizeDirName,outputDirName,givenOrderFileName,
                                                t_index_model,t_ordering_flag,num_threads,
                                                day_slice_cnt,graphType,queryDataDirName);
                }
                
                return 0;
            }
            //main ends
    };
}


#endif // !_COMMAND_TIP_INDEXING
