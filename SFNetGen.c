#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/*

This code generates a scale-free network of N nodes and whose power-law
distribution has an exponent of beta. These parameters are introduced below.
Usage:
	./SFNetGen [filename] [TFSteps]
The program receives 2 parameters:
filename -- the name of the output file containing the adjacency matrix of the network
TFSteps -- the code performs a triad formation step to maximize the clustering coefficient
	of the network.

The program outputs several files:
[filename] contains the adjacency matrix of the created network
[filename]_CCs shows the clustering coefficient of every node in the network
[filename]_inDD contains the in-degree distribution of the network. (This should follow
	a Poisson distribution roughly)
[filename]_inDegs shows the in-degree of every node in the network
[filename]_outDD contains the out-degree distribution of the network. This follows a power-
	law distribution when plotted.
[filename]_outDegs shows the out-degree of every node in the network
[filename]_T is the adjacency matrix of the transpose of the network created. Here, the
	direction of edges get reversed. Thus, out-degrees become in-degrees, and vice versa.
 
*/

struct nodeA{
    short value;
    struct nodeA *right;
};

struct nodeB{
    short value;
    struct nodeB *right;
    struct nodeB *left;
    struct nodeB *up;
    struct nodeB *down;
};

struct nodeC{
    float value;
    struct nodeC *right;
};


int main(int argc, char *argv[]){
	int N = 1024;
	float C = 200;
	float beta = -1.0;
	int x0 = 1;
    int TFSteps = atoi(argv[2]);
    int node_1ClusteringBoostStep = 0;
    
	int outDD;
	int nodeCount;
	int remNodes = N;
	int i,j,k;
    int ren, col;
	int nodes2connect, nodes2create;
    
    float pInOrOutN = 0.5;
    float r;
    
    int numClusteringSteps;
    int flagListCreated;
    int availableInN, availableOutN;
    
    int flagRepeated;
    float meanCC;
    int numNeighbors;
    int numOnes;
    
    struct nodeA *nodeList;
    struct nodeB *Net;
    struct nodeA *inDD;
    struct nodeA *outDeg;
    struct nodeA *inDeg;
    struct nodeA *inNeighbors, *outNeighbors, *inNLast, *outNLast, *Neighbors;
    struct nodeC *CCs;
    
    struct nodeA *auxA1, *auxA2, *auxA3, *auxA4, *auxA5, *auxA6;
    struct nodeB *auxB1, *auxB2, *auxB3, *auxB4;
    struct nodeC *auxC1, *auxC2;
    
    
    int triad[3];
    int triadIndex;
    int flagTriad;
    int successful1NodeTriads;
    
    char filename[30];
    //filename[0] = 0;
    
    strcpy(filename,argv[1]);
	FILE *outDegDist;
	if((outDegDist = fopen(strcat(filename,"_outDD"),"wt")) == NULL) printf("file error\n");
    
    strcpy(filename,argv[1]);
	FILE *inDegDist;
	if((inDegDist = fopen(strcat(filename,"_inDD"),"wt")) == NULL) printf("file error\n");
    
	FILE *SFNet;
	if((SFNet = fopen(argv[1],"wt")) == NULL) printf("file error\n");
    
    strcpy(filename, argv[1]);
	FILE *SFNetT; //transpose of SFNet
	if((SFNetT = fopen(strcat(filename,"_T"),"wt")) == NULL) printf("file error\n");
    
    strcpy(filename,argv[1]);
	FILE *inDegs;
	if((inDegs = fopen(strcat(filename,"_inDegs"),"wt")) == NULL) printf("file error\n");
    
    strcpy(filename,argv[1]);
	FILE *outDegs;
	if((outDegs = fopen(strcat(filename,"_outDegs"),"wt")) == NULL) printf("file error\n");
    
    strcpy(filename,argv[1]);
    FILE *NetCCs;
	if((NetCCs = fopen(strcat(filename,"_CCs"),"wt")) == NULL) printf("file error\n");    
    
	srand48(time(NULL));
    
    //initialize
    //create networks:
    // 1) nodeList
    nodes2create = N;
    nodeList = (struct nodeA*) malloc (sizeof(struct nodeA));
    if(nodeList == 0){
        printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
        return 1;
    }
    nodeList->value = 0;
    nodeList->right = NULL;
    auxA1 = nodeList;
    nodes2create--;
    do{
        auxA2 = (struct nodeA*) malloc (sizeof(struct nodeA));
        if(auxA2 == 0){
            printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
            return 1;
        }
        auxA2->value = 0;
        auxA2->right = NULL;
        auxA1->right = auxA2;
        auxA1 = auxA2;
        nodes2create--;        
    } while(nodes2create > 0);
    printf("nodeList created!\n\n");
    
    // 2) inDD
    nodes2create = N+1;
    inDD = (struct nodeA*) malloc (sizeof(struct nodeA));
    if(inDD == 0){
        printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
        return 1;
    }
    inDD->value = 0;
    inDD->right = NULL;
    auxA1 = inDD;
    nodes2create--;
    do{
        auxA2 = (struct nodeA*) malloc (sizeof(struct nodeA));
        if(auxA2 == 0){
            printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
            return 1;
        }
        auxA2->value = 0;
        auxA2->right = NULL;
        auxA1->right = auxA2;
        auxA1 = auxA2;
        nodes2create--;
    } while(nodes2create > 0);
    printf("inDD created!\n\n");
    
    // 3) outDeg
    nodes2create = N;
    outDeg = (struct nodeA*) malloc (sizeof(struct nodeA));
    if(outDeg == 0){
        printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
        return 1;
    }
    outDeg->value = 0;
    outDeg->right = NULL;
    auxA1 = outDeg;
    nodes2create--;
    do{
        auxA2 = (struct nodeA*) malloc (sizeof(struct nodeA));
        if(auxA2 == 0){
            printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
            return 1;
        }
        auxA2->value = 0;
        auxA2->right = NULL;
        auxA1->right = auxA2;
        auxA1 = auxA2;
        nodes2create--;        
    } while(nodes2create > 0); 
    printf("outDeg created!\n\n");
    
    // 4) inDeg
    nodes2create = N;
    inDeg = (struct nodeA*) malloc (sizeof(struct nodeA));
    if(inDeg == 0){
        printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
        return 1;
    }
    inDeg->value = 0;
    inDeg->right = NULL;
    auxA1 = inDeg;
    nodes2create--;
    do{
        auxA2 = (struct nodeA*) malloc (sizeof(struct nodeA));
        if(auxA2 == 0){
            printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
            return 1;
        }
        auxA2->value = 0;
        auxA2->right = NULL;
        auxA1->right = auxA2;
        auxA1 = auxA2;
        nodes2create--;
    } while(nodes2create > 0);    
    printf("inDeg created!\n\n");
    
    // 5) Net
    ren = 1;
    col = 1;
    nodes2create = N*N;
    Net = (struct nodeB*) malloc (sizeof(struct nodeB));
    if(Net == 0){
        printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
        return 1;
    }
    Net->value = 0;
    Net->right = NULL;
    Net->left = NULL;
    Net->up = NULL;
    Net->down = NULL;
    auxB1 = Net;
    auxB3 = Net;
    nodes2create--;
    col++;
    do{
        auxB4 = (struct nodeB*) malloc (sizeof(struct nodeB));
        if(auxB4 == 0){
            printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
            return 1;
        }
        auxB4->value = 0;
        auxB4->right = NULL;
        auxB4->left = NULL;
        auxB4->up = NULL;
        auxB4->down = NULL;
        if(ren == 1){
            auxB1->right = auxB4;
            auxB4->left = auxB1;
        }
        else{
            if(col == 1){
                auxB2 = auxB3;
                auxB3 = auxB4;
            }
            else{
                auxB1->right = auxB4;
                auxB4->left = auxB1;
            }
            auxB2->down = auxB4;
            auxB4->up = auxB2;
            auxB2 = auxB2->right;
        }
        auxB1 = auxB4;
        col++;
        if(col > N){
            col = 1;
            ren++;
        }
        nodes2create--;
    } while(nodes2create > 0);
    printf("Net created!\n\n");

    // 6) CCs
    nodes2create = N;
    CCs = (struct nodeC*) malloc (sizeof(struct nodeC));
    if(CCs == 0){
        printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
        return 1;
    }
    CCs->value = 0.0;
    CCs->right = NULL;
    auxC1 = CCs;
    nodes2create--;
    do{
        auxC2 = (struct nodeC*) malloc (sizeof(struct nodeC));
        if(auxC2 == 0){
            printf("ERROR: Out of memory! Remaining nodes: %d\n",nodes2create);
            return 1;
        }
        auxC2->value = 0.0;
        auxC2->right = NULL;
        auxC1->right = auxC2;
        auxC1 = auxC2;
        nodes2create--;        
    } while(nodes2create > 0);
    printf("CCs list created!\n\n");

    
//////////////////////////////////////////////////////////////////////////////    

	//main routine: compute degree distribution and connect nodes
	nodeCount = 0;
    triadIndex = 0;
    flagTriad = 0;
    successful1NodeTriads = 0;
	do{
		//routine to get degree:
		//outDD = (int)((float)(N*C) * pow((float)x0,beta));
        outDD = (int)(C * pow((float)x0,beta));
		nodeCount += outDD;
		//fprintf(outDegDist, "%d %d %d\n", x0, outDD, nodeCount);
		fprintf(outDegDist, "%d %d\n", x0, outDD);
        
		//routine to connect according to degree:
		do{
			//choose a node at random from available nodes as node-origin:
			printf("Choosing node-origin");

            do{
                printf("*");
                if(x0 == 1 && node_1ClusteringBoostStep == 1 && triadIndex > 0 && flagTriad == 1){
                    i = triad[triadIndex];
                    flagTriad = 0;
                }
                else{
                    i = (int)(drand48()*(double) N);
                    if(x0 == 1 && node_1ClusteringBoostStep == 1 && triadIndex > 0){
                        triadIndex = 0;
                    }
                }
                //move to position i+1 in the nodeList
                col = 0;
                auxA1 = nodeList;
                auxA2 = outDeg;
                while(col < i){
                    auxA1 = auxA1->right;
                    auxA2 = auxA2->right;
                    col++;
                }
            } while((auxA1->value) == 1);            
            if(x0 == 1 && node_1ClusteringBoostStep == 1){
                switch(triadIndex){
                    case 0:
                        triad[0] = i;
                        break;
                    case 2:
                        flagTriad = 1;
                        break;
                }
            }

            printf("\n");
			//connect node-origin to another x0 nodes:
			nodes2connect = x0;
			do{
				//choose node-destination:
				printf("Choosing node-destination");
                
                do{
                    printf("+");
                    if(x0 == 1 && node_1ClusteringBoostStep == 1 && triadIndex == 2 && flagTriad == 1){
                        j = triad[0];
                        flagTriad = 0;
                    }
                    else
                        j = (int)(drand48()*(double) N);
                    ren = 0;
                    col = 0;
                    auxB1 = Net;
                    auxA3 = inDeg;
                    auxA4 = outDeg;
                    //move pointer to row of node origin:
                    while(ren < i){
                        auxB1 = auxB1->down;
                        ren++;
                    }
                    //move pointer to column of node destination
                    //and the other pointer to inDeg and outDeg of node destination:
                    while(col < j){
                        auxB1 = auxB1->right;
                        auxA3 = auxA3->right;
                        auxA4 = auxA4->right;
                        col++;
                    }
                } while(i == j || (auxB1->value) == 1);
                if(x0 == 1 && node_1ClusteringBoostStep == 1){
                    if(triadIndex == 2){
                        triadIndex = 0;
                        successful1NodeTriads++;
                    }
                    else{
                        triadIndex++;
                        triad[triadIndex] = j;
                        flagTriad = 1;
                    }
                }
				printf("\n");
				//connect both:
				auxB1->value = 1;
                //increment outdegree of node origin:
                auxA2->value++;
                //increment indegree of node destination:
                auxA3->value++;
                
				nodes2connect--;
                
                //here goes some code to implement tunable clustering!
                if(nodes2connect > 0 && TFSteps > 0){
                    numClusteringSteps = TFSteps;
                    availableInN = 0;
                    availableOutN = 0;
                    
                    if(auxA3->value > 0 || auxA4->value >0){                                                
                        if(auxA3->value > 0){
                            ren = 0;
                            col = 0;                            
                            //j has inNeighbors, populate list with inNeighbors of j:
                            auxB1 = Net;
                            flagListCreated = 0;
                            while(col < j){
                                auxB1 = auxB1->right;
                                col++;
                            }
                            for(ren = 0; ren < N; ren++){
                                if(auxB1->value == 1 || auxB1->value == -1){
                                    if(flagListCreated == 0){
                                        inNeighbors = (struct nodeA*) malloc (sizeof(struct nodeA));
                                        if(inNeighbors == 0){
                                            printf("ERROR: Out of memory to create inNeighbors list!\n");
                                            return 1;
                                        }
                                        inNeighbors->value = ren;
                                        inNeighbors->right = NULL;
                                        inNLast = inNeighbors;
                                        flagListCreated = 1;
                                    }
                                    else{
                                        auxA6 = (struct nodeA*) malloc (sizeof(struct nodeA));
                                        if(auxA6 == 0){
                                            printf("ERROR: Out of memory to create inNeighbors list!\n");
                                            return 1;
                                        }
                                        auxA6->value = ren;
                                        auxA6->right = NULL;
                                        inNLast->right = auxA6;
                                        inNLast = auxA6;
                                    }
                                }
                                auxB1 = auxB1->down;
                            }
                            availableInN = 1;
                        }
                        
                        if(auxA4->value > 0){
                            ren = 0;
                            col = 0;                            
                            //j has outNeighbors, populate list with outNeighbors of j:
                            auxB1 = Net;
                            flagListCreated = 0;
                            while(ren < j){
                                auxB1 = auxB1->down;
                                ren++;
                            }
                            for(col = 0; col < N; col++){
                                if(auxB1->value == 1 || auxB1->value == -1){
                                    if(flagListCreated == 0){
                                        outNeighbors = (struct nodeA*) malloc (sizeof(struct nodeA));
                                        if(outNeighbors == 0){
                                            printf("ERROR: Out of memory to create outNeighbors list!\n");
                                            return 1;
                                        }
                                        outNeighbors->value = col;
                                        outNeighbors->right = NULL;
                                        outNLast = outNeighbors;
                                        flagListCreated = 1;
                                    }
                                    else{
                                        auxA6 = (struct nodeA*) malloc (sizeof(struct nodeA));
                                        if(auxA6 == 0){
                                            printf("ERROR: Out of memory to create outNeighbors list!\n");
                                            return 1;
                                        }
                                        auxA6->value = col;
                                        auxA6->right = NULL;
                                        outNLast->right = auxA6;
                                        outNLast = auxA6;
                                    }
                                }
                                auxB1 = auxB1->right;
                            }
                            availableOutN = 1;
                        }
                        //at this point we have populated neighbor lists for j
                        while(nodes2connect > 0 && numClusteringSteps > 0 && (availableInN == 1 || availableOutN == 1)){
                            if(availableInN == 1 && availableOutN == 1){
                                //flip a (biased?) coin:
                                r = drand48();
                                if(r <= pInOrOutN){
                                    //choose an outNeighbor!
                                    k = outNeighbors->value;
                                    if(outNeighbors == outNLast)
                                        availableOutN = 0;
                                    else
                                        outNeighbors = outNeighbors->right;
                                }
                                else{
                                    //choose an inNeighbor!
                                    k = inNeighbors->value;
                                    if(inNeighbors == inNLast)
                                        availableInN = 0;
                                    else
                                        inNeighbors = inNeighbors->right;
                                }
                            }
                            else{
                                if(availableOutN == 1){
                                    //choose an outNeighbor!
                                    k = outNeighbors->value;
                                    if(outNeighbors == outNLast)
                                        availableOutN = 0;
                                    else
                                        outNeighbors = outNeighbors->right;                                    
                                }
                                if(availableInN == 1){
                                    //choose an inNeighbor!
                                    k = inNeighbors->value;
                                    if(inNeighbors == inNLast)
                                        availableInN = 0;
                                    else
                                        inNeighbors = inNeighbors->right;                                    
                                }
                            }
                            ren = 0;
                            col = 0;
                            auxB1 = Net;
                            auxA3 = inDeg;//at this moment we lose track of inDeg of j and focus on k
                            auxA4 = outDeg;//at this moment we lose track of outDeg of j and focus on k
                            //move pointer to row of node origin:
                            while(ren < i){
                                auxB1 = auxB1->down;
                                ren++;
                            }
                            //move pointer to column of node destination
                            //and the other pointer to inDeg and outDeg of node destination:
                            while(col < k){
                                auxB1 = auxB1->right;
                                auxA3 = auxA3->right;
                                auxA4 = auxA4->right;
                                col++;
                            }
                            if(i != k && (auxB1->value) == 0){
                                //connect i to k:
                                printf(".");
                                auxB1->value = 1;
                                auxA2->value++;//increments outDeg of i
                                auxA3->value++;//increments inDeg of k
                                nodes2connect--;
                                numClusteringSteps--;
                            }
                        }
                        printf("\n");
                        //at the end, free memory!
                    }//else, do nothing, continue with the rest of the routine
                }//here ends code for tunable clustering
			} while(nodes2connect > 0);
			printf("Node %d is now connected to another %d nodes. Remaining nodes: %d\n", i, x0, --remNodes);
            auxA1->value = 1;
			outDD--;
		} while(outDD > 0 && nodeCount < N);
		//ends routine to connect!
        
		x0++;		
	} while(x0 <= N && nodeCount < N);
	//end of main routine
    //UP TO HERE THE NETWORK IS CREATED!
    
///////////////////////////////////////////////////////////////////////////////
    
    //Code to compute clustering coeffs per node!
    auxC1 = CCs;//CCs
    auxB1 = Net;//scan the net vertically (inDeg)
    auxB2 = Net;//scan the net horizontally (outDeg)
    auxA1 = outDeg;//outDegree per node
    auxA2 = inDeg; //inDegree per node
    meanCC = 0.0;
    
    printf("\n");

    for(i = 0; i < N; i++){
        if((auxA1->value + auxA2->value) > 1){
            numNeighbors = 0;
            numOnes = 0;
            flagListCreated = 0;
            //gather outNeighbors of i:
            if(auxA1->value > 0){
                col = 0;
                auxB3 = auxB2;
                for(col = 0; col < N; col++){
                    if(auxB3->value == 1 || auxB3->value == -1){
                        if(flagListCreated == 0){
                            Neighbors = (struct nodeA*) malloc (sizeof(struct nodeA));
                            if(Neighbors == 0){
                                printf("ERROR: Out of memory when creating neighbor list for node %d!\n", i);
                                return 1;
                            }
                            Neighbors->value = col;
                            Neighbors->right = NULL;
                            auxA3 = Neighbors;
                            flagListCreated = 1;
                            numNeighbors++;
                        }
                        else{
                            auxA4 = (struct nodeA*) malloc (sizeof(struct nodeA));
                            if(auxA4 == 0){
                                printf("ERROR: Out of memory when creating neighbor list for node %d!\n", i);
                                return 1;                            
                            }
                            auxA4->value = col;
                            auxA4->right = NULL;
                            auxA3->right = auxA4;
                            auxA3 = auxA4;
                            numNeighbors++;
                        }
                    }
                    auxB3 = auxB3->right;
                }            
            }
            //gather inNeighbors of i:
            if(auxA2->value > 0){
                ren = 0;
                auxB3 = auxB1;
                for(ren = 0; ren < N; ren++){
                    if(auxB3->value == 1 || auxB3->value == -1){
                        if(flagListCreated == 0){
                            Neighbors = (struct nodeA*) malloc (sizeof(struct nodeA));
                            if(Neighbors == 0){
                                printf("ERROR: Out of memory when creating neighbor list for node %d!\n", i);
                                return 1;                            
                            }
                            Neighbors->value = ren;
                            Neighbors->right = NULL;
                            auxA3 = Neighbors;
                            flagListCreated = 1;
                            numNeighbors++;
                        }
                        else{
                            //first check for repeated neighbors, ie. those who are counted as outNeighbors of i
                            auxA4 = Neighbors;
                            flagRepeated = 0;
                            do{
                                if(auxA4->value == ren){
                                    flagRepeated = 1;
                                    break;
                                }
                                auxA4 = auxA4->right;
                            }while(auxA4 != NULL);
                            if(flagRepeated == 0){
                                auxA4 = (struct nodeA*) malloc (sizeof(struct nodeA));
                                if(auxA4 == 0){
                                    printf("ERROR: Out of memory when creating neighbor list for node %d!\n", i);
                                    return 1;                                
                                }
                                auxA4->value = ren;
                                auxA4->right = NULL;
                                auxA3->right = auxA4;
                                auxA3 = auxA4;
                                numNeighbors++;
                            }
                        }
                    }
                    auxB3 = auxB3->down;
                }
            }//up to this point we have a list with all neighbors of i
            //now, look for -existing- connections among neighbors of i:
            auxA3 = Neighbors;
            while(auxA3 != NULL){
                auxA4 = Neighbors;
                while(auxA4 != NULL){
                    if(auxA3 != auxA4){
                        j = auxA3->value;
                        k = auxA4->value;
                        //now, check if (j,k) = 1 in the adjacency matrix!
                        auxB3 = Net;
                        ren = 0;
                        col = 0;
                        while(ren < j){
                            auxB3 = auxB3->down;
                            ren++;
                        }
                        while(col < k){
                            auxB3 = auxB3->right;
                            col++;
                        }
                        if(auxB3->value == -1 || auxB3->value == 1)
                            numOnes++;
                    }
                    auxA4 = auxA4->right;
                }
                auxA3 = auxA3->right;
            }
            auxC1->value = (float) numOnes / (numNeighbors*(numNeighbors - 1));
        }
        else
            auxC1->value = 0.0;
        printf("CC(%d) = %g\n",i, auxC1->value);
        meanCC += auxC1->value;
        
        auxB1 = auxB1->right;
        auxB2 = auxB2->down;
        
        auxC1 = auxC1->right;
        auxA1 = auxA1->right;//outDeg
        auxA2 = auxA2->right;//inDeg
    }
    meanCC = (float) meanCC / N;
    printf("\nMean CC: %g\n", meanCC);
    
///////////////////////////////////////////////////////////////////////////////
    
    printf("\nSuccesful 1-node triads: %d\n", successful1NodeTriads);
    
	printf("\nClustering Steps: %d\n", TFSteps);
    
    //compute in-degree distribution:        
    auxA1 = inDeg;
    for(col = 0; col < N; col++){
        i = 0;
        auxA2 = inDD;
        while(i < (auxA1->value)){
            auxA2 = auxA2->right;
            i++;
        }
        auxA2->value++;
        auxA1 = auxA1->right;        
    }
    
	//print final network, CCs and degrees    
    auxB2 = Net;
    auxA1 = inDeg;
    auxA2 = outDeg;
    auxA3 = inDD;
    auxC1 = CCs;
    
    for(ren = 0; ren < N; ren++){
        auxB1 = auxB2;
        for(col = 0; col < N; col++){
            fprintf(SFNet, "%d", auxB1->value);
            if((col+1) < N){
                fprintf(SFNet, " ");
            }
            auxB1 = auxB1->right;
        }
        
        fprintf(inDegs, "d(%d) = %d\n", ren, auxA1->value);
        auxA1 = auxA1->right;
        fprintf(outDegs, "d(%d) = %d\n", ren, auxA2->value);
        auxA2 = auxA2->right;
        fprintf(inDegDist, "%d %d\n", ren, auxA3->value);
        auxA3 = auxA3->right;
        fprintf(NetCCs, "CC(%d) = %g\n", ren, auxC1->value);
        auxC1 = auxC1->right;
        
        auxB2 = auxB2->down;
        fprintf(SFNet, "\n");
    }
    fprintf(inDegDist, "%d %d\n", ren, auxA3->value);

    auxB2 = Net;
    for(col = 0; col < N; col++){
        auxB1 = auxB2;
        for(ren = 0; ren < N; ren++){
            fprintf(SFNetT, "%d", auxB1->value);
            if((ren+1) < N){
                fprintf(SFNetT, " ");
            }
            auxB1 = auxB1->down;
        }
        auxB2 = auxB2->right;
        fprintf(SFNetT, "\n");
    }
    
    //routine to check if theres a node with indegree or outdegree zero:
    auxA3 = inDeg;
    auxA4 = outDeg;
    for(i = 0; i < N; i++){
        if(auxA3->value == 0 || auxA4->value == 0)
            printf("Node %d is not reachable or transmitter!\n", i);
        auxA3 = auxA3->right;
        auxA4 = auxA4->right;
    }
    
	fclose(outDegDist);
	fclose(inDegDist);
	fclose(SFNet);
	fclose(SFNetT);
	fclose(inDegs);
	fclose(outDegs);
    fclose(NetCCs);

    
	printf("\n");
    
    //free(nodeList);
    
    return 0;
}