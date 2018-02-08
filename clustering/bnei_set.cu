#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "bnei_set.h"

using namespace std;

__global__ void bnei_set(float **var, int **intVar){

	int *bnei = intVar[14];
	int nxbin = *intVar[15];
	int nybin = *intVar[16];
	int nzbin = *intVar[17];
	int bdimx = *intVar[22];
	int bdimy = *intVar[23];
	int bdimz = *intVar[24];

	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	int xbin, ybin, zbin, xcen, ycen, zcen;

	int xpos, ypos, zpos;
	int xind, yind, zind;
	int xdiff, ydiff, zdiff;
	int ind;

	int tid2;

	xcen = (bdimx - 1) / 2;
	ycen = (bdimy - 1) / 2;
	zcen = (bdimz - 1) / 2;

	zbin = tid / (nxbin*nybin);
	ybin = (tid - zbin*nxbin*nybin) / nxbin;
	xbin = tid - ybin*nxbin - zbin*nxbin*nybin;

	for (tid2 = 0; tid2 < bdimx*bdimy*bdimz; tid2++){

		zpos = tid2 / (bdimx*bdimy);
		ypos = (tid2 - zpos*bdimx*bdimy) / bdimx;
		xpos = tid2 - ypos*bdimx - zpos*bdimx*bdimy;

		xdiff = xpos - xcen;
		xind = xbin + xdiff;
		if (xind < 0){
			xind += nxbin;
		}
		if (xind >= nxbin){
			xind -= nxbin;
		}
		ydiff = ypos - ycen;
		yind = ybin + ydiff;
		if (yind < 0){
			yind += nybin;
		}
		if (yind >= nybin){
			yind -= nybin;
		}
		zdiff = zpos - zcen;
		zind = zbin + zdiff;
		if (zind < 0){
			zind += nzbin;
		}
		if (zind >= nzbin){
			zind -= nzbin;
		}

		ind = xind + yind*nxbin + zind*nxbin*nybin;
		if ((xbin + ybin*nxbin + zbin*nxbin*nybin + tid2 * nxbin*nybin*nzbin) >= 0 && (xbin + ybin*nxbin + zbin*nxbin*nybin + tid2 * nxbin*nybin*nzbin) < nxbin*nybin*nzbin*bdimx*bdimy*bdimz){
			bnei[xbin + ybin*nxbin + zbin*nxbin*nybin + tid2 * nxbin*nybin*nzbin] = ind;
		}
		else{
			printf("tid2 %4d\n", tid2);
		}
	}
}
