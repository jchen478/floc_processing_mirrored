#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "delta_twist_zero.h"

using namespace std;

__global__ void delta_twist_zero(float **var, int **intVar){

	int *ncnt = intVar[1];
	int *ncpf = intVar[2];
	int  nfib = *intVar[6];
	int  nseg = *intVar[7];
	int *bin = intVar[12];
	int nxbin = *intVar[15];
	int nybin = *intVar[16];
	int nzbin = *intVar[17];
	int *status = intVar[19];
	int *nc = intVar[21];

	float *fx = var[60];
	float *fy = var[61];
	float *fz = var[62];
	float *tx = var[63];
	float *ty = var[64];
	float *tz = var[65];
	float *fcx = var[66];
	float *fcy = var[67];
	float *fcz = var[68];
	float *tcx = var[69];
	float *tcy = var[70];
	float *tcz = var[71];

	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	int tid2 = threadIdx.x + blockDim.x * blockIdx.x;
	int tid5 = threadIdx.x + blockDim.x * blockIdx.x;
	int tid8 = threadIdx.x + blockDim.x * blockIdx.x;

	status[tid5] = 0;

	fx[tid5] = 0.0;
	fy[tid5] = 0.0;
	fz[tid5] = 0.0;
	tx[tid5] = 0.0;
	ty[tid5] = 0.0;
	tz[tid5] = 0.0;
	fcx[tid5] = 0.0;
	fcy[tid5] = 0.0;
	fcz[tid5] = 0.0;
	tcx[tid5] = 0.0;
	tcy[tid5] = 0.0;
	tcz[tid5] = 0.0;
	nc[tid5] = 0;
	ncpf[tid5] = 0;
	ncnt[tid5] = 0;

	while (tid8 < nxbin*nybin*nzbin){
		bin[tid8] = 0;
		tid8 += blockDim.x*gridDim.x;
	}

	if (tid == 0){
		*intVar[3] = 0; // num_groups
		*intVar[4] = 0; // overs
		*intVar[5] = 0; // total_contacts
	}

}
