#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include "initialize_regrow2.h"

using namespace std;

__global__ void initialize_regrow2(float **var, int **intVar){

	int    nseg = *intVar[7];
	float rp = *var[114];

	float *rx = var[3];
	float *ry = var[4];
	float *rz = var[5];
	float *px = var[18];
	float *py = var[19];
	float *pz = var[20];

	float pjx, pjy, pjz;
	int mi, m, i, mth;

	int tid = threadIdx.x + blockDim.x*blockIdx.x;

	m = tid / (nseg - 1);
	i = tid - m * (nseg - 1) + 1;
	mi = m*nseg + i;

	// Regrow the fibers
	pjx = 0.0; pjy = 0.0; pjz = 0.0;
	for (mth = 1; mth<(i); mth++){
		pjx = pjx + *(px + m*nseg + mth);
		pjy = pjy + *(py + m*nseg + mth);
		pjz = pjz + *(pz + m*nseg + mth);
	}
	*(rx + mi) = *(rx + m*nseg) + rp**(px + m*nseg) + 2.0*rp*pjx + rp**(px + mi);
	*(ry + mi) = *(ry + m*nseg) + rp**(py + m*nseg) + 2.0*rp*pjy + rp**(py + mi);
	*(rz + mi) = *(rz + m*nseg) + rp**(pz + m*nseg) + 2.0*rp*pjz + rp**(pz + mi);

}
