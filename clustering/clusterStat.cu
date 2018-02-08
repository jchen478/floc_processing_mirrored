#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "clusterStat.h"

using namespace std;

__global__ void clusterStat(float **var, int **intVar){

	int nfib = *intVar[6]; 
	int *cluster = intVar[35]; 
	int *clusterAccess = intVar[36]; 
	int *clusterCount = intVar[37]; 

	int m = threadIdx.x + blockIdx.x*blockDim.x; 
	int mm; 

	clusterCount[m] = 0; 
	if (clusterAccess[m] == 1){
		for (mm = 0; mm < nfib; mm++){
			clusterCount[m] += cluster[m*nfib + mm]; 
		}
		if(clusterCount[m] > 1)
			printf("m %6d count %6d\n", m, clusterCount[m]); 
	}
}
