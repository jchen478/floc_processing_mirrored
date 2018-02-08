#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "link.h"

using namespace std;

__global__ void link(float **var, int **intVar){

	const int npcn = 2000;

	int nseg = *intVar[7];
	int *bin = intVar[12];
	int *list = intVar[13];
	int *bnei = intVar[14];
	int nxbin = *intVar[15];
	int nybin = *intVar[16];
	int nzbin = *intVar[17];
	int *bnum = intVar[31];
	int *potCon = intVar[32];
	int *potConSize = intVar[33];


	int mi = blockIdx.x;
	int tid = threadIdx.x;
	int miBin = bnum[mi];
	int nextBin = bnei[miBin + tid*nxbin*nybin*nzbin];
	int nextTot = bin[nextBin];

	int nj, i, pos;

	if (tid == 0){
		potConSize[mi] = 0;
	}
	__syncthreads();

	for (i = 0; i < nextTot; i++){
		nj = list[nextBin + i*nxbin*nybin*nzbin];
		// skip when segments are connected
		if (mi >= nj){
			continue;
		}
		if ((mi - (mi / nseg)*nseg) != 0 && nj == mi - 1){
			continue;
		}
		if ((nj - (nj / nseg)*nseg) != 0 && mi == nj - 1){
			continue;
		}
		if ((mi - (mi / nseg)*nseg) != nseg - 1 && nj == mi + 1){
			continue;
		}
		if ((nj - (nj / nseg)*nseg) != nseg - 1 && mi == nj + 1){
			continue;
		}
		//printf("link mi nj %4d %4d\n", mi, nj); 
		pos = atomicAdd(potConSize + mi, 1);
		if (pos >= npcn - 1) printf("allocate more space for potCon pos %4d\n", pos);
		potCon[mi * npcn + pos] = nj;
	}
}
