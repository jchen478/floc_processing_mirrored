#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include "clusteringZero.h"

using namespace std;

__global__ void clusteringZero(float **var, int **intVar){

	int nfib = *intVar[6]; 
	int *cluster = intVar[35]; 
	int *clusterAccess = intVar[36]; 
	
	int m = threadIdx.x+blockIdx.x*blockDim.x;
	
	int ind;

	clusterAccess[m] = 1; 
	for(ind = 0; ind<nfib; ind++){
		cluster[m*nfib+ind] = 0; 
		if(ind == m) cluster[m*nfib+ind] = 1;
	}
}
