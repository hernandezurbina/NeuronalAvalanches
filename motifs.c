#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*

There are 13 different motif types in a directed network.
This code counts the presence of each motif in a given directed network.
./motifs [file] [numNodes]

[file] -- name of file containing the adjacency matrix of the network
[numNodes] -- number of nodes in the network

output:
[file]_motifsDistribution -- file describing the motif count for each type

*/

int main(int argc, char *argv[]) {
	int N = atoi(argv[2]);
	int Net[N][N];
	int i,j,k;
	int m,n;
	int triad[3][3];
	char ch[1];
	int outDeg[3], inDeg[3];
	int inDD[3], outDD[3];
	int motifDist[14];
	int node;
    char filename[50];

	FILE * net;
	if ((net=fopen(argv[1],"r"))==NULL) printf("file error\n");

    strcpy(filename, argv[1]);
	FILE *motifsDistribution;
	if((motifsDistribution = fopen(strcat(filename,"_motifsDistribution"),"wt")) == NULL) printf("file error\n");
    

	//read from the net file specifying topology
	i = 0;
	j = 0;

	ch[0] = getc(net);
	while (ch[0] != EOF){
		if(ch[0] != ' '){
			if(ch[0] == '\n'){
				i += 1;
				j = 0;
			}
			else{
				Net[i][j] = atoi(ch);
				j += 1;
			}
		}
	 	ch[0] = getc(net);
	}
	fclose(net);

	//initialize distribution of motifs:
	for(i = 0; i < 14; i++){
		motifDist[i] = 0;
	}

	/*printf("\n");
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			printf("%d ", Net[i][j]);
		}
		printf("\n");
	}
	printf("\n");*/

	//starts main routine
	for(i = 0; i < (N-2); i++){
		for(j = i+1; j < (N-1); j++){
			for(k = j+1; k < N; k++){
				//initialize histograms for this run:
				for(m = 0; m < 3; m++){
					inDD[m] = 0;
					outDD[m] = 0;
				}

				//printf("motif: %d %d %d\n", i, j, k);
				triad[0][0] = Net[i][i];
				triad[0][1] = Net[i][j];
				triad[0][2] = Net[i][k];
				triad[1][0] = Net[j][i];
				triad[1][1] = Net[j][j];
				triad[1][2] = Net[j][k];
				triad[2][0] = Net[k][i];
				triad[2][1] = Net[k][j];
				triad[2][2] = Net[k][k];
				for(m = 0; m < 3; m++){
					inDeg[m] = 0;
					outDeg[m] = 0;
					for(n = 0; n < 3; n++){
						//printf("%d ", triad[m][n]);
						inDeg[m] += triad[n][m];
						outDeg[m] += triad[m][n];
					}
					inDD[inDeg[m]]++; //fill histogram
					outDD[outDeg[m]]++; //fill histogram
					//printf("\n");
				}
				//printf("outDD: ");
				/*for(m = 0; m < 3; m++){
					printf("[%d]: %d ", m, outDD[m]);
				}
				printf("\ninDD: ");
				for(m = 0; m < 3; m++){
					printf("[%d]: %d ", m, inDD[m]);
				}
				printf("\n\n");*/
				//fill histogram for motif distribution:
				if(outDD[0] == 2 && outDD[1] == 0 && outDD[2] == 1 && inDD[0] == 1 && inDD[1] == 2 && inDD[2] == 0){
					motifDist[1]++;
					continue;
				}
				if(outDD[0] == 1 && outDD[1] == 2 && outDD[2] == 0 && inDD[0] == 1 && inDD[1] == 2 && inDD[2] == 0){
					motifDist[2]++;
					continue;
				}
				if(outDD[0] == 1 && outDD[1] == 1 && outDD[2] == 1 && inDD[0] == 0 && inDD[1] == 3 && inDD[2] == 0){
					motifDist[3]++;
					continue;
				}
				if(outDD[0] == 1 && outDD[1] == 2 && outDD[2] == 0 && inDD[0] == 2 && inDD[1] == 0 && inDD[2] == 1){
					motifDist[4]++;
					continue;
				}
				if(outDD[0] == 1 && outDD[1] == 1 && outDD[2] == 1 && inDD[0] == 1 && inDD[1] == 1 && inDD[2] == 1){
					motifDist[5]++;
					continue;
				}
				if(outDD[0] == 1 && outDD[1] == 0 && outDD[2] == 2 && inDD[0] == 0 && inDD[1] == 2 && inDD[2] == 1){
					motifDist[6]++;
					continue;
				}
				if(outDD[0] == 0 && outDD[1] == 3 && outDD[2] == 0 && inDD[0] == 1 && inDD[1] == 1 && inDD[2] == 1){
					motifDist[7]++;
					continue;
				}
				if(outDD[0] == 0 && outDD[1] == 2 && outDD[2] == 1 && inDD[0] == 0 && inDD[1] == 2 && inDD[2] == 1){
					//which node has outdegree 2:
					for(m = 0; m < 3; m++){
						if(outDeg[m] == 2){
							node = m;
							break;
						}
					}					
					if(inDeg[node] == 2){
						motifDist[8]++;
					}
					else{
						motifDist[10]++;
					}
					continue;
				}
				if(outDD[0] == 0 && outDD[1] == 3 && outDD[2] == 0 && inDD[0] == 0 && inDD[1] == 3 && inDD[2] == 0){
					motifDist[9]++;
					continue;
				}
				if(outDD[0] == 0 && outDD[1] == 2 && outDD[2] == 1 && inDD[0] == 1 && inDD[1] == 0 && inDD[2] == 2){
					motifDist[11]++;
					continue;
				}
				if(outDD[0] == 0 && outDD[1] == 1 && outDD[2] == 2 && inDD[0] == 0 && inDD[1] == 1 && inDD[2] == 2){
					motifDist[12]++;
					continue;
				}
				if(outDD[0] == 0 && outDD[1] == 0 && outDD[2] == 3 && inDD[0] == 0 && inDD[1] == 0 && inDD[2] == 3){
					motifDist[13]++;
					continue;
				}
			}
		}
	}//end main routine

	//output distribution
	for(i = 1; i < 14; i++){
		fprintf(motifsDistribution, "%d %d\n", i, motifDist[i]);
	}

	fclose(motifsDistribution);
	return(1);
}
