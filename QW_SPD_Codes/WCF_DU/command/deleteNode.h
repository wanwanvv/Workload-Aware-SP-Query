#ifndef _COMMAND_DELETE_NODE
#define _COMMAND_DELETE_NODE

#include <cstring>
#include <iostream>
#include "../command.h"
#include "../src/paras.h"
#include "../src/time_util.h"
#include "../src/overlay_graph.h"

using namespace time_util; 

namespace command{
    class DeleteNode:public Command
    {
    public:
        void exit_with_help(){
            printf("Usage:\n");
            printf("wcf_du_run -e -d [directedGraphFlag] -w [weightedGraphFlag] -s [experimentFlag] -g [graphFileName] \n -h [queryFreqFileName] -o [outputDirName] -m [removeMaxDegree] [-i [deleteModel]] [-r [hfRate]] \n [-t [deleteNodeRate]] [-b [numOfBfsHop]] [-l [initialDegreeOfLocalSearch]] \n");
            printf("-------------------------------------------------------------------\n");
            printf("Parameter explanation:\n");
            printf("[deleteNodeRate]: user-defined vale of deleted vertices divided by 1000 [0,1000], default-0 means reserve all high frequent vertices\n");
            printf("[deleteModel] = 0: elete nodes by degree \
 			\t 1: delete nodes by query freq and degree \
 			\t 2: delete nodes by query freq and degree optimized with local search \
			\t 3: delete nodes by query freq and degree optimized with local search (parameters for\n");
            printf("[removeMaxDegree]: the max degree of nodes to be removed (>=1)\n");
            printf("numOfBfsHop]: the bfs search hop limitation for local search (optimization)\n");
            printf("[initialDegreeOfLocalSearch]: the degree to start local search (optimization)\n");
            printf("[hfRate]: user-defined vale of high frequent vertices divided by 10000 [0,10000], default-0 means all read vertices\n");
            printf("-------------------------------------------------------------------\n");
            exit(1);
        }

        int main(int argc, char *argv[]){
            int t_directed_flag = 0;
            int t_weighted_flag = 0;
            int t_experiment_flag = 0;
            char graph_file[255];//original graph
            char HFPoint_file[255];//hfPoint filename
            char output_dir_name[255] ;//delete lf point and add shortcut filename
            char analysisFileName[255];//experiment modelmoutput result
            int deleteModel = 0;//delete node way 
            int hfRate=0; //high frequency point rate default-0 means all read vertices
            int deleteNodeRate=0;//deleteNode ratio default-0 means delete non-hfpoint
            int maxDegree=5;//max node degree to be deleted default-5
            int bfsHop_num=1;//bfs search threshold default=-1
            int localSearchDegree=INF_WEIGHT;//when >=degree then start  local search
            int graphType=0; //0-full connected 1-u is start and t is destination

            // if(argc < 14)  exit_with_help();
            for(int i = 2; i < argc; i++){
                if(argv[i][0] != '-') break;
                if(++i >= argc) exit_with_help();
                switch (argv[i-1][1]){
                    case 'g':
                        strcpy(graph_file,argv[i]);
                        std::cout<<"graph_file = "<<graph_file<<std::endl;
                        break;
                    case 'h':
                        strcpy(HFPoint_file,argv[i]);
                        std::cout<<"HFPoint_file = "<<HFPoint_file<<std::endl;
                        break;
                    case 'o':
                        strcpy(output_dir_name,argv[i]);
                        std::cout<<"output_dir_name ="<<output_dir_name<<std::endl;
                        break;                 
                    case 'a':
                        strcpy(analysisFileName,argv[i]);
                        std::cout<<"analysisFileName="<<analysisFileName<<std::endl;
                        break;    
                    case 'd':
                        t_directed_flag = atoi(argv[i]);
                        break;
                    case 'w':
                        t_weighted_flag = atoi(argv[i]);
                        break;
                    case 's':
                        t_experiment_flag = atoi(argv[i]);
                        break;
                    case 'm':
                        maxDegree =atoi(argv[i]);
                        std::cout<<"maxDegree = "<<maxDegree<<std::endl;
                        break;
                    case 'i':
                        deleteModel =atoi(argv[i]);
                        break;
                    case 'c':
                        graphType = atoi(argv[i]);
                        std::cout<<"graphType = "<<graphType<<std::endl;
                        break;
                    case 'r':
                        hfRate = atoi(argv[i]);
                        std::cout<<"hfRate = "<<hfRate<<std::endl;
                        break;
                    case 'b':
                        bfsHop_num = atoi(argv[i]);
                        std::cout<<"bfsHop_num = "<<bfsHop_num<<std::endl;
                        break;
                    case 'l':
                        localSearchDegree = atoi(argv[i]);
                        std::cout<<"localSearchDegree = "<<localSearchDegree<<std::endl;
                        break;
                    case 't':
                        deleteNodeRate = atoi(argv[i]);
                        std::cout<<"deleteNodeRate = "<<deleteNodeRate<<std::endl;
                        break;
                    default:
                        exit_with_help();
                        break;
                }
            }

            if (t_directed_flag == 1) DIRECTED_FLAG = true;
            if (t_weighted_flag == 1) WEIGHTED_FLAG = true;

            //check
            //if(isDeleteByFreq&&*HFPoint_file=='\0') exit_with_help();

            if (DIRECTED_FLAG == true){
                if( WEIGHTED_FLAG == true){
                        std::cout<<"directed weighted graph..."<<endl;
                        Node_deletion node_deletion;
                        node_deletion.load_and_copy_original_graph(graph_file,graphType);
                        if (numOfVertices == 0 || numOfEdges == 0){
                            cout << "Corruptted graph file";
                            return 0;
                        }
                        if(t_experiment_flag==0){
                            cout << "experiment_flag=0!"<<endl;
                            node_deletion.deleteNode_directed(deleteModel,maxDegree,deleteNodeRate,hfRate,HFPoint_file);
                            node_deletion.generate_overlay_graph_directed(output_dir_name,graphType);
                        }else if(t_experiment_flag==1){
                            cout << "experiment_flag=1!"<<endl;
                        }
                        return 0;
                }else{}
            }else{//Undirected
                if( WEIGHTED_FLAG == true){//Weighted
                    std::cout<<"undirected weighted graph..."<<endl;
                    Node_deletion node_deletion;
                    node_deletion.load_and_copy_original_graph(graph_file,graphType);
                    if (numOfVertices == 0 || numOfEdges == 0){
                        cout << "Corruptted graph file" << endl;
                        return 0;
                    }
                    if(t_experiment_flag==0){
                        node_deletion.deleteNode(deleteModel ,maxDegree,deleteNodeRate,hfRate,HFPoint_file,bfsHop_num,localSearchDegree);
                        node_deletion.generate_overlay_graph(output_dir_name);
                    }else if(t_experiment_flag==1){
                        node_deletion.deleteNode(deleteModel ,maxDegree,deleteNodeRate,hfRate,HFPoint_file,bfsHop_num,localSearchDegree,analysisFileName);
                        node_deletion.generate_overlay_graph(output_dir_name,analysisFileName);
                    }

                    return 0;
                }else{}
            }
            return 0;
        }//main end


    };//DeleteNode class end
    
}

#endif