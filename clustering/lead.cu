#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "lead.h"

using namespace std;


__global__ void lead(float **var, int **intVar){

	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	int *ncpf = intVar[2];
	// ncpf[tid] - number of segments in direct contact with segment tid
	int *clist = intVar[18];
	// clist[tid*maxCon+i] - the i th fiber in contact with segment tid
	int *status = intVar[19];
	// status[tid] - whether the list of fiber tid has been accessed
	int *lead_clist = intVar[20];
	// lead_clist[tid*maxGr+k] - the k th fiber in group led by tid
	int *nc = intVar[21];
	// nc[tid] - number of fiber in group led by tid
	int maxCon = *intVar[25];
	// maxCon - maximum number of fibers contacting one fiber
	int maxGr = *intVar[27];
	// maxGr - maximum number of fibers in a group

	int pos, fiber, locfiber, i, j, k; 
	bool add;

	pos = 0;
	add = false;

	for (i = 0; i < ncpf[tid]; i++){
		fiber = clist[tid*maxCon + i];
		if (fiber < tid){
			status[tid] = 1;
			return;
		}
		lead_clist[maxGr*tid + nc[tid]] = clist[tid*maxCon + i];
		nc[tid]++;
	}
	while (pos != nc[tid]){
		fiber = lead_clist[maxGr*tid + pos];
		for (j = 0; j < ncpf[fiber]; j++){
			locfiber = clist[fiber*maxCon + j];
			if (locfiber < tid){
				status[tid] = 1;
				return;
			}
			if (locfiber == tid){
				continue;
			}
			add = true;
			for (k = 0; k < nc[tid]; k++){
				if (locfiber == lead_clist[maxGr*tid + k]){
					add = false;
					break;
				}
			}
			if (add){
				lead_clist[maxGr*tid + nc[tid]] = locfiber;
				nc[tid]++;
				if (nc[tid] >= maxGr){
					printf("error in lead: increase maxGr\n"); 
				}
			}
		}
		pos++;
	}

}