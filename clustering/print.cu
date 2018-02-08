#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include "print.h"
#include <math.h>



__global__ void print(float **var, int **intVar){

	int nfib = *intVar[6]; 
	int *cluster = intVar[35]; 
	int m, n; 

	for (m = 120; m < 160; m++){
		for (n = 120; n < 160 ; n++){
			printf("%2d", cluster[m*nfib + n]); 
		}
		printf("\n"); 
	}
	

}
