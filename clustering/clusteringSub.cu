#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "clusteringSub.h"

using namespace std;

__global__ void clusteringSub(float **var, int **intVar){

	int nfib = *intVar[6]; 
	int *cluster = intVar[35]; 
	int *clusterAccess = intVar[36]; 

	int m, n, ind, N; 
	int m_start, limit, dum; 
	bool combine; 

	int tid = threadIdx.x + blockIdx.x*blockDim.x; 

	dum = nfib / (blockDim.x*gridDim.x); 

	m_start = tid*dum; 
	limit = (tid + 1)*dum; 

	combine = false; 

	for (m = m_start; m < limit; m++){
		if(clusterAccess[m] == 1){
			N = 0; 
			for (ind = 0; ind < nfib; ind++){
				N += cluster[m*nfib + ind]; 
			}
			if (N == 1){
				clusterAccess[m] = 0; 
			}
		}
	}


	for (m = m_start; m < limit; m++){
		if (clusterAccess[m] == 1){
			for (n = m_start; n < limit; n++){
				if ((clusterAccess[n] == 1) && (m != n)){
					for (ind = 0; ind < nfib; ind++){
						if ((cluster[n*nfib + ind] == 1) && (cluster[m*nfib + ind] == cluster[n*nfib + ind])){
							combine = true;
							clusterAccess[n] = 0;
							break;
						}
					}
					if (combine){
						cluster[m*nfib + n] = 1;
						for (ind = m_start; ind < nfib; ind++){
							if (cluster[n*nfib + ind] == 1){
								cluster[m*nfib + ind] = 1;
							}
						}
						n = m_start;
						combine = false;
					}
				}
			}
		}
	}
}
