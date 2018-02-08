#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "cell.h"

using namespace std;

__global__ void cell(float **var, int **intVar){
		
	int *bin = intVar[12]; 
	int *list = intVar[13]; 
	int nxbin = *intVar[15];
	int nybin = *intVar[16];
	int nzbin = *intVar[17];
	int maxBin = *intVar[28]; 
	int *bnum = intVar[31]; 

	float *rx = var[3]; 
	float *ry = var[4]; 
	float *rz = var[5]; 
	float sidex = *var[119]; 
	float sidey = *var[120]; 
	float sidez = *var[121]; 
	float delta_rx = *var[138];
	float dx = *var[108]; 
	float dy = *var[109]; 
	float dz = *var[110]; 

	int xbin, ybin, zbin, old, corx, cory, corz;
	float rxmi, rymi, rzmi; 


	int mi = threadIdx.x + blockDim.x*blockIdx.x; 

	cory = roundf(ry[mi] / sidey);
	corz = roundf(rz[mi] / sidez);
	rxmi = rx[mi] - corz*delta_rx; 
	corx = roundf(rxmi / sidex); 
	rymi = ry[mi] - cory*sidey; 
	rzmi = rz[mi] - corz*sidez; 
	rxmi = rxmi - corx*sidex; 

	xbin = int(floorf(rxmi / dx)) + nxbin / 2; 
	ybin = int(floorf(rymi / dy)) + nybin / 2;
	zbin = int(floorf(rzmi / dz)) + nzbin / 2;

	if (xbin == nxbin){
		xbin--; 
	}
	if (ybin == nybin){
		ybin--;
	}
	if (zbin == nzbin){
		zbin--;
	}
	if (xbin == -1){
		xbin++; 
	}
	if (ybin == -1){
		ybin++; 
	}
	if (zbin == -1){
		zbin++; 
	}
	if (xbin < 0 || ybin < 0 || zbin < 0 || xbin >= nxbin || ybin >= nybin || zbin >= nzbin){
		printf("error in cell: bin index out of bound\n"); 
		printf("mi = %4d (%3d %3d %3d) rx ry rz %10f %10f %10f \n", mi, xbin, ybin, zbin, rx[mi], ry[mi], rz[mi]);
	}

	old = atomicAdd(bin + xbin + ybin*nxbin + zbin*nxbin*nybin, 1);
	list[xbin + ybin*nxbin + zbin*nxbin*nybin + old*nxbin*nybin*nzbin] = mi;
	bnum[mi] = xbin + ybin*nxbin + zbin*nxbin*nybin;
	if (old > maxBin){
		printf("error in cell: number of segments per bin %4d exceeds maxBin\n", old);
	}
}