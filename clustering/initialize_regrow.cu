#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include "initialize_regrow.h"

using namespace std;

__global__ void initialize_regrow(float **var, int **intVar){

	int nseg = *intVar[7];
	float rp = *var[114];

	float *rcmx = var[0];
	float *rcmy = var[1];
	float *rcmz = var[2];
	float *rx = var[3];
	float *ry = var[4];
	float *rz = var[5];
	float *px = var[18];
	float *py = var[19];
	float *pz = var[20];

	float pjx, pjy, pjz, pkx, pky, pkz; 
	int mi, m, k, j; 

	m = threadIdx.x + blockDim.x * blockIdx.x;

	/*if (m == 0){
		printf("regrow nseg %d rp %f\n", nseg, rp); 
	}*/

	mi = m*nseg;

	// Regrow the fibers		
	pjx = 0.0; pjy = 0.0; pjz = 0.0;
	for (j = 1; j<nseg; j++){
		pkx = 0.0; pky = 0.0; pkz = 0.0;
		for (k = 1; k<(j); k++){
			pkx = pkx + *(px + m*nseg + k);
			pky = pky + *(py + m*nseg + k);
			pkz = pkz + *(pz + m*nseg + k);
		}
		pjx = pjx + (*(px + m*nseg) + *(px + m*nseg + j) + 2 * pkx);
		pjy = pjy + (*(py + m*nseg) + *(py + m*nseg + j) + 2 * pky);
		pjz = pjz + (*(pz + m*nseg) + *(pz + m*nseg + j) + 2 * pkz);
	}
	*(rx + mi) = *(rcmx + m) - (rp / (float)nseg)*pjx;
	*(ry + mi) = *(rcmy + m) - (rp / (float)nseg)*pjy;
	*(rz + mi) = *(rcmz + m) - (rp / (float)nseg)*pjz;
}
