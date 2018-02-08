#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "clustering.h"

using namespace std;

__global__ void clustering(float **var, int **intVar){

	int nfib = *intVar[6]; 
	int *cluster = intVar[35]; 
	int *clusterAccess = intVar[36]; 

	int m, n, ind; 
	bool combine; 

	combine = false; 
	
	for (m = 0; m < nfib; m++){
		if (clusterAccess[m] == 1){
			for (n = 0; n < nfib; n++){
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
						for (ind = 0; ind < nfib; ind++){						
							if (cluster[n*nfib + ind] == 1){
								cluster[m*nfib + ind] = 1;
							}
						}
						n = 0;
						combine = false; 
					}
				}
			}
		}
	}
}
