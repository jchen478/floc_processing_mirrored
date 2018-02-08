#include <cuda_runtime_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

__global__ void shift(float *rx, float *ry, float *rz,
	float *sidex, float *sidey, float *sidez);
__global__ void print(int *dim1, int *dim2, int *sliceMax);
__global__ void zeroBin(int *bin); 
__global__ void binning(int *bin, float *rx, float *ry, float *rz,
	float *d_dx, float *d_dy, float *d_dz,
	int *d_nx, int *d_ny, int *d_nz);
__global__ void flocSizex(int *bin, int *binFill,
	int *d_nx, int *d_ny, int *d_nz, int *d_maxGap, int *sliceMax);
__global__ void flocSizey(int *bin, int *binFill,
	int *d_nx, int *d_ny, int *d_nz, int *d_maxGap, int *sliceMax);
__global__ void flocSizez(int *bin, int *binFill,
	int *d_nx, int *d_ny, int *d_nz, int *d_maxGap, int *sliceMax);
__global__ void maxX(int *sliceMax, int *d_nx, int *d_ny, int *d_nz, int *flocDim); 
__global__ void maxY(int *sliceMax, int *d_nx, int *d_ny, int *d_nz, int *flocDim);
__global__ void maxZ(int *sliceMax, int *d_nx, int *d_ny, int *d_nz, int *flocDim);

int main(){

	// Relevant files
	FILE *Parameters;
	FILE *rxfile, *ryfile, *rzfile;
	FILE *ClusterInfo, *Cluster_results, *Cluster;
	FILE *Cluster_size;

	Parameters = fopen("Parameters.in", "r");
	rxfile = fopen("../rx.txt", "rb");
	ryfile = fopen("../ry.txt", "rb");
	rzfile = fopen("../rz.txt", "rb");
	ClusterInfo = fopen("ClusterInfo.in", "r");
	Cluster_results = fopen("Cluster_results.txt", "r");
	Cluster = fopen("Cluster.txt", "r");
	Cluster_size = fopen("Cluster_size.txt", "w");

	// Variables in input files
	int nfib, nseg, config_write, idum;
	int nFloc, floc_cutoff, maxGap;
	float dt, strain, mustat, mukin, dum;
	float rps, sidex, sidey, sidez;

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &rps);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &mustat);
	fscanf(Parameters, "%*[^\n]%f", &mukin);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dt);
	fscanf(Parameters, "%*[^\n]%f", &strain);
	fscanf(Parameters, "%*[^\n]%f", &sidex);
	fscanf(Parameters, " %f", &sidey);
	fscanf(Parameters, " %f", &sidez);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &config_write);

	// Read in cluster info
	fscanf(ClusterInfo, "%d", &idum);
	fscanf(ClusterInfo, "%*[^\n]%d", &nFloc);
	fscanf(ClusterInfo, "%*[^\n]%d", &floc_cutoff);
	fscanf(ClusterInfo, "%*[^\n]%d", &maxGap);

	// Read in number of fibers in each floc
	int *nfibFloc, *fiberId;
	nfibFloc = (int*)malloc(nFloc*sizeof(int)); 
	fiberId = (int*)malloc(nfib*sizeof(int));
	for (int i = 0; i < nFloc; i++){
		fscanf(Cluster_results, "%d", &idum);
		fscanf(Cluster_results, " %d", nfibFloc + i);
	}

	// Divide domain into bins
	float dx, dy, dz;
	int nx, ny, nz, *d_maxGap;
	int *d_bin, *d_binFill; 
	int *d_sliceMaxX, *d_sliceMaxY, *d_sliceMaxZ;
	int *d_flocDim, *flocDim; 
 
	dx = rps; 	
	dy = rps;	
	dz = rps;
	nx = int(ceil(sidex / dx)); 
	ny = int(ceil(sidey / dy));
	nz = int(ceil(sidez / dz));
	dx = sidex / float(nx); 
	dy = sidey / float(ny); 
	dz = sidez / float(nz);
	flocDim = (int*)malloc(3*sizeof(int)); 
	cudaMalloc((void**)&d_maxGap, sizeof(int));
	cudaMalloc((void**)&d_flocDim, 3*sizeof(int));
	cudaMalloc((void**)&d_sliceMaxX, ny*nz*sizeof(int));
	cudaMalloc((void**)&d_sliceMaxY, nx*nz*sizeof(int));
	cudaMalloc((void**)&d_sliceMaxZ, nx*ny*sizeof(int));
	cudaMalloc((void**)&d_bin, nx*ny*nz*sizeof(int)); 
	cudaMalloc((void**)&d_binFill, nx*ny*nz*sizeof(int));
	cudaMemcpy(d_maxGap, &maxGap, sizeof(int), cudaMemcpyHostToDevice); 


	printf("binning parameters: %4d %4d %4d %10.6f %10.6f %10.6f\n",
		nx, ny, nz, dx, dy, dz); 
	

	// Read in configuration at the last frame //
	float *rx, *ry, *rz;
	float *rxC, *ryC, *rzC; 
	int step, nConfig;
	rx = (float*)malloc(nfib*nseg*sizeof(float));
	ry = (float*)malloc(nfib*nseg*sizeof(float));
	rz = (float*)malloc(nfib*nseg*sizeof(float));
	rxC = (float*)malloc(nfib*nseg*sizeof(float));
	ryC = (float*)malloc(nfib*nseg*sizeof(float));
	rzC = (float*)malloc(nfib*nseg*sizeof(float));
	nConfig = int(strain / dt / float(config_write)) + 1;
	
	for (step = 0; step < nConfig; step++){
		fread(&dum, sizeof(float), 1, rxfile);
		fread(rx, sizeof(float), nfib*nseg, rxfile);
		fread(&dum, sizeof(float), 1, ryfile);
		fread(ry, sizeof(float), nfib*nseg, ryfile);
		fread(&dum, sizeof(float), 1, rzfile);
		fread(rz, sizeof(float), nfib*nseg, rzfile);
	}

	// Copy variables to device
	float *d_sidex, *d_sidey, *d_sidez;
	float *d_rxC, *d_ryC, *d_rzC;
	float *d_dx, *d_dy, *d_dz; 
	int *d_nx, *d_ny, *d_nz;
	cudaMalloc((void**)&d_sidex, sizeof(float));
	cudaMalloc((void**)&d_sidey, sizeof(float));
	cudaMalloc((void**)&d_sidez, sizeof(float));
	cudaMalloc((void**)&d_rxC, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_ryC, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_rzC, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_dx, sizeof(float));
	cudaMalloc((void**)&d_dy, sizeof(float));
	cudaMalloc((void**)&d_dz, sizeof(float));
	cudaMalloc((void**)&d_nx, sizeof(int));
	cudaMalloc((void**)&d_ny, sizeof(int));
	cudaMalloc((void**)&d_nz, sizeof(int));
	cudaMemcpy(d_sidex, &sidex, sizeof(float), cudaMemcpyHostToDevice); 
	cudaMemcpy(d_sidey, &sidey, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sidez, &sidez, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dx, &dx, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dy, &dy, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dz, &dz, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_nx, &nx, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ny, &ny, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_nz, &nz, sizeof(int), cudaMemcpyHostToDevice);

	// temporary variables
	int f, m, n, mi, ni, i, flocCount; 

	// loop through the number of flocs
	for (f = 0; f < nFloc; f++){

		// read the fiber indexes belonging to group f
		fscanf(Cluster, " %d", &idum);
		fscanf(Cluster, " %d", &idum);
		for (m = 0; m < nfibFloc[f]; m++){
			fscanf(Cluster, " %d", &fiberId[m]);
		}
		// only calculate size for flocs larger than cutoff
		if (nfibFloc[f] >= floc_cutoff){

			flocCount = nfibFloc[f]; 

			// 1. obtain coordinates of fibers in floc f
			for (m = 0; m < flocCount; m++){
				n = fiberId[m]; 
				for (i = 0; i < nseg; i++){
					mi = m*nseg + i; 
					ni = n*nseg + i; 
					rxC[mi] = rx[ni];
					ryC[mi] = ry[ni];
					rzC[mi] = rz[ni];
				}
			}

			// 2. copy coordinates to device
			cudaMemcpy(d_rxC, rxC, flocCount*nseg*sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(d_ryC, ryC, flocCount*nseg*sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(d_rzC, rzC, flocCount*nseg*sizeof(float), cudaMemcpyHostToDevice);

			// 3. shift coordinates so that all are in the box and positive 
			shift << <flocCount,nseg >> >(d_rxC, d_ryC, d_rzC, d_sidex, d_sidey, d_sidez);

			// 4. zero overall bin 
			zeroBin << <nx*ny, nz >> >(d_bin); 

			// 5. place segments into bin
			binning << < flocCount,nseg >> >(d_bin, d_rxC, d_ryC, d_rzC,
				d_dx, d_dy, d_dz, d_nx, d_ny, d_nz); 

			// 6x. zero bins to contain info for filled bin
			zeroBin << <nx*ny, nz >> >(d_binFill);

			// 7x. fill in holes in the floc and obtain a maximum matrix
			flocSizex << <ny, nz >> >(d_bin, d_binFill, 
				d_nx, d_ny, d_nz, d_maxGap, d_sliceMaxX); 

			// 8x. find maximum of that matrix 
			maxX << <1, nz >> >(d_sliceMaxX, d_nx, d_ny, d_nz, d_flocDim); 
			
			// 6y. zero bins to contain info for filled bin
			zeroBin << <nx*ny, nz >> >(d_binFill);

			// 7y. fill in holes in the floc and obtain a maximum matrix
			flocSizey << <nx, nz >> >(d_bin, d_binFill,
				d_nx, d_ny, d_nz, d_maxGap, d_sliceMaxY);

			//print << <1, 1 >> >(d_nx, d_nz, d_sliceMaxY);
			// x - dim1 = y; dim2 = z
			// y - dim1 = x; dim2 = z
			// z - dim1 = x; dim2 = y 

			//cudaDeviceSynchronize();
			//printf("after maximizing \n"); 

			// 8y. find maximum of that matrix 
			maxY << <1, nz >> >(d_sliceMaxY, d_nx, d_ny, d_nz, d_flocDim);

			//print << <1, 1 >> >(d_nx, d_nz, d_sliceMaxY);
			// x - dim1 = y; dim2 = z
			// y - dim1 = x; dim2 = z
			// z - dim1 = x; dim2 = y 

			// 6z. zero bins to contain info for filled bin
			zeroBin << <nx*ny, nz >> >(d_binFill);

			// 7z. fill in holes in the floc and obtain a maximum matrix
			flocSizez << <nx, ny >> >(d_bin, d_binFill,
				d_nx, d_ny, d_nz, d_maxGap, d_sliceMaxZ);

			// 8y. find maximum of that matrix 
			maxZ << <1, ny >> >(d_sliceMaxZ, d_nx, d_ny, d_nz, d_flocDim);
			
			// 9. copy maximum info back to host
			cudaMemcpy(flocDim, d_flocDim, 3 * sizeof(int), cudaMemcpyDeviceToHost); 
			fprintf(Cluster_size, "%6d %4d %6.2f %6.2f %6.2f %8.2f %10.2f %10.2f %10.2f %6.2f %6.2f %6.2f %5d %5d %5d %4d %8.2f %8.2f %8.2f\n",
				nfib, nseg, mustat, mukin, rps, strain, sidex, sidey, sidez, dx, dy, dz, f, floc_cutoff, flocCount, maxGap,
				dx*float(flocDim[0]), dy*float(flocDim[1]), dz*float(flocDim[2]));
			
		}
	}


	// close files
	fclose(Parameters); 
	fclose(Cluster_size);
	fclose(ClusterInfo); fclose(Cluster); fclose(Cluster_results);
	fclose(rxfile); fclose(ryfile); fclose(rzfile); 
	
	// free memory
	free(rx); free(ry); free(rz);
	free(rxC); free(ryC); free(rzC);
	free(nfibFloc); free(fiberId); 
	cudaFree(d_bin); cudaFree(d_binFill); cudaFree(d_maxGap); 
	cudaFree(d_sliceMaxX); cudaFree(d_sliceMaxY); cudaFree(d_sliceMaxZ);
	cudaFree(d_flocDim);
	cudaFree(d_dx); cudaFree(d_dy); cudaFree(d_dz);
	cudaFree(d_nx); cudaFree(d_ny); cudaFree(d_nz);
	cudaFree(d_rxC); cudaFree(d_ryC); cudaFree(d_rzC);
	cudaFree(d_sidex); cudaFree(d_sidey); cudaFree(d_sidez);

	return 0;
}

__global__ void shift(float *rx, float *ry, float *rz,
	float *d_sidex, float *d_sidey, float *d_sidez){

	int mi = threadIdx.x + blockDim.x*blockIdx.x;
	float sidex, sidey, sidez;

	sidex = *d_sidex;
	sidey = *d_sidey;
	sidez = *d_sidez;

	// shift coordinates inside box
	if (rx[mi] < -sidex / 2.0 || rx[mi] > sidex / 2){
		rx[mi] -= copysignf(sidex,rx[mi]);
	}
	if (ry[mi] < -sidey / 2.0 || ry[mi] > sidey / 2){
		ry[mi] -= copysignf(sidey, ry[mi]);
	}
	if (rz[mi] < -sidez / 2.0 || rz[mi] > sidez / 2){
		rz[mi] -= copysignf(sidez, rz[mi]);
	}
	// shift coordinates so all are positive
	rx[mi] += sidex / 2.0;
	ry[mi] += sidey / 2.0;
	rz[mi] += sidez / 2.0;
}
__global__ void zeroBin(int *bin){

	int b = threadIdx.x + blockDim.x*blockIdx.x;
	bin[b] = 0;
}
__global__ void print(int *dim1, int *dim2, int *sliceMax){
	
	// x - dim1 = y; dim2 = z
	// y - dim1 = x; dim2 = z
	// z - dim1 = x; dim2 = y 

	int i, j;

	for (j = 0; j < *dim2; j++){
		for (i = 0; i < *dim1; i++){
			printf("%2d ", sliceMax[i + j * *dim1]);
		}
		printf("\n");
	}

}
__global__ void binning(int *bin, float *rx, float *ry, float *rz,
	float *d_dx, float *d_dy, float *d_dz, 
	int *d_nx, int *d_ny, int *d_nz){

	int mi = threadIdx.x + blockDim.x*blockIdx.x;

	float dx, dy, dz;
	int nx, ny, nz; 
	int xloc, yloc, zloc;
	dx = *d_dx; 
	dy = *d_dy; 
	dz = *d_dz; 
	nx = *d_nx; 
	ny = *d_ny; 
	nz = *d_nz; 

	xloc = floorf(rx[mi] / dx); 
	yloc = floorf(ry[mi] / dy);
	zloc = floorf(rz[mi] / dz);

	atomicAdd(bin + xloc + yloc*nx + zloc*nx*ny, 1); 
	
	// check for out of bound access
	if (xloc < 0 || xloc >= nx){
		printf("x index out of bound %4d\n", xloc); 
	}
	if (yloc < 0 || yloc >= ny){
		printf("y index out of bound %4d\n", yloc);
	}
	if (zloc < 0 || zloc >= nz){
		printf("x index out of bound %4d\n", zloc);
	}
}
__global__ void flocSizex(int *bin, int *binFill,
	int *d_nx, int *d_ny, int *d_nz, int *d_maxGap, int *sliceMax){

	int nx, ny, nz, maxGap;
	nx = *d_nx;
	ny = *d_ny;
	nz = *d_nz;
	maxGap = *d_maxGap;

	// new variables
	int xloc, yloc, zloc;
	int sum, oldBin, newBin, ii;
	bool access; 
	yloc = blockIdx.x;
	zloc = threadIdx.x;

	// zero maximum 
	sliceMax[yloc + zloc*ny] = 0; 

	// find sum of x array at yloc and zloc
	sum = 0;
	for (xloc = 0; xloc < nx; xloc++){
		sum += bin[xloc + yloc*nx + zloc*nx*ny]; 
	}

	access = false; 
	// only perform calculations if there are fibers in array
	if (sum != 0){
		oldBin = 0;
		newBin = 0; 
		// find the position of the previous occupied bin
		for (xloc = nx - 1; xloc >= nx - maxGap; xloc--){
			if (bin[xloc + yloc*nx + zloc*nx*ny] != 0){
				binFill[xloc + yloc*nx + zloc*nx*ny] = 1;
				oldBin = xloc; 
				access = true; 
				break;
			}
		}
		// fill in the holes from oldBin to newBin if within maxGap
		for (xloc = 0; xloc <= nx; xloc++){
			if (bin[xloc + yloc*nx + zloc*nx*ny] != 0){
				binFill[xloc + yloc*nx + zloc*nx*ny] = 1; 
				newBin = xloc; 
				if (access){
					if ((nx - 1) - oldBin + newBin <= maxGap){
						for (ii = oldBin; ii < nx; ii++){
							binFill[ii + yloc*nx + zloc*nx*ny] = 1;
						}
						for (ii = 0; ii < newBin; ii++){
							binFill[ii + yloc*nx + zloc*nx*ny] = 1;
						}
					}
					access = false; 
				}
				else if (bin[oldBin + yloc*nx + zloc*nx*ny] != 0){
					if (newBin - oldBin - 1 <= maxGap){
						for (ii = oldBin + 1; ii < newBin; ii++){
							binFill[ii + yloc*nx + zloc*nx*ny] = 1;
						}
					}
				}
				oldBin = newBin;
			}
		}
		// update maximum from this perspective
		sum = 0;
		for (xloc = 0; xloc < nx; xloc++){
			sum += binFill[xloc + yloc*nx + zloc*nx*ny];
		}
		sliceMax[yloc + zloc*ny] = sum;
		if (sum >= nx){
			printf("yloc zloc %4d %4d\n", yloc, zloc); 
		}
	}
}
__global__ void maxX(int *sliceMax, int *d_nx, int *d_ny, int *d_nz, int *flocDim){
	
	int nx, ny, nz;
	nx = *d_nx;
	ny = *d_ny;
	nz = *d_nz;

	int zloc = threadIdx.x; 
	int yloc, z; 

	for (yloc = 1; yloc < ny; yloc++){
		if (sliceMax[yloc + zloc*ny] > sliceMax[zloc*ny]){
			sliceMax[zloc*ny] = sliceMax[yloc + zloc*ny]; 
		}
	}

	__syncthreads(); 

	if (zloc == 0){
		for (z = 1; z < nz; z++){
			if (sliceMax[0] < sliceMax[z*ny]){
				sliceMax[0] = sliceMax[z*ny];
			}
		}
		flocDim[0] = sliceMax[0]; 
	}

}

__global__ void flocSizey(int *bin, int *binFill,
	int *d_nx, int *d_ny, int *d_nz, int *d_maxGap, int *sliceMax){

	int nx, ny, nz, maxGap;
	nx = *d_nx;
	ny = *d_ny;
	nz = *d_nz;
	maxGap = *d_maxGap;

	// new variables
	int xloc, yloc, zloc;
	int sum, oldBin, newBin, ii;
	bool access;
	xloc = blockIdx.x;
	zloc = threadIdx.x;

	// zero maximum 
	sliceMax[xloc + zloc*nx] = 0;

	// find sum of y array at xloc and zloc
	sum = 0;
	for (yloc = 0; yloc < ny; yloc++){
		sum += bin[xloc + yloc*nx + zloc*nx*ny];		
	}

	access = false;
	// only perform calculations if there are fibers in array
	if (sum != 0){
		oldBin = 0;
		newBin = 0;
		// find the position of the previous occupied bin
		for (yloc = ny - 1; yloc >= ny - maxGap; yloc--){
			if (bin[xloc + yloc*nx + zloc*nx*ny] != 0){
				binFill[xloc + yloc*nx + zloc*nx*ny] = 1;
				oldBin = yloc;
				access = true;
				break;
			}
		}
		// fill in the holes from oldBin to newBin if within maxGap
		for (yloc = 0; yloc < ny; yloc++){
			if (bin[xloc + yloc*nx + zloc*nx*ny] != 0){
				binFill[xloc + yloc*nx + zloc*nx*ny] = 1;
				newBin = yloc;
				if (access){
					if ((ny - 1) - oldBin + newBin <= maxGap){
						for (ii = oldBin; ii < ny; ii++){
							binFill[xloc + ii*nx + zloc*nx*ny] = 1;
						}
						for (ii = 0; ii < newBin; ii++){
							binFill[xloc + ii*nx + zloc*nx*ny] = 1;
						}
					}
					access = false;
				}
				else if(bin[xloc + oldBin*nx+zloc*nx*ny] != 0){
					if (newBin - oldBin - 1 <= maxGap){
						for (ii = oldBin + 1; ii < newBin; ii++){
							binFill[xloc + ii*nx + zloc*nx*ny] = 1;
						}
					}
				}
				oldBin = newBin;
			}
		}

		// update maximum from this perspective
		sum = 0;
		for (yloc = 0; yloc < ny; yloc++){
			sum += binFill[xloc + yloc*nx + zloc*nx*ny];
		}
		sliceMax[xloc + zloc*nx] = sum;
		if (sum >= ny){
			printf("xloc zloc %4d %4d\n", xloc, zloc);
		}
	}
}
__global__ void maxY(int *sliceMax, int *d_nx, int *d_ny, int *d_nz, int *flocDim){

	int nx, ny, nz;
	nx = *d_nx;
	ny = *d_ny;
	nz = *d_nz;

	int zloc = threadIdx.x;
	int xloc, z;

	for (xloc = 1; xloc < nx; xloc++){
		if (sliceMax[xloc + zloc*nx] > sliceMax[zloc*nx]){
			sliceMax[zloc*nx] = sliceMax[xloc + zloc*nx];
		}
	}

	__syncthreads();

	if (zloc == 0){
		for (z = 1; z < nz; z++){
			if (sliceMax[0] < sliceMax[z*nx]){
				sliceMax[0] = sliceMax[z*nx];
			}
		}
		flocDim[1] = sliceMax[0];
	}

}
__global__ void flocSizez(int *bin, int *binFill,
	int *d_nx, int *d_ny, int *d_nz, int *d_maxGap, int *sliceMax){

	int nx, ny, nz, maxGap;
	nx = *d_nx;
	ny = *d_ny;
	nz = *d_nz;
	maxGap = *d_maxGap;

	// new variables
	int xloc, yloc, zloc;
	int sum, oldBin, newBin, ii;
	bool access;
	xloc = blockIdx.x;
	yloc = threadIdx.x;


	// zero maximum 
	sliceMax[xloc + yloc*nx] = 0;

	// find sum of y array at xloc and yloc
	sum = 0;
	for (zloc = 0; zloc < nz; zloc++){
		sum += bin[xloc + yloc*nx + zloc*nx*ny];
	}

	access = false;
	// only perform calculations if there are fibers in array
	if (sum != 0){
		oldBin = 0;
		newBin = 0;
		// find the position of the previous occupied bin
		for (zloc = nz - 1; zloc >= nz - maxGap; zloc--){
			if (bin[xloc + yloc*nx + zloc*nx*ny] != 0){
				binFill[xloc + yloc*nx + zloc*nx*ny] = 1;
				oldBin = zloc;
				access = true;
				break;
			}
		}
		// fill in the holes from oldBin to newBin if within maxGap
		for (zloc = 0; zloc <= nz; zloc++){
			if (bin[xloc + yloc*nx + zloc*nx*ny] != 0){
				binFill[xloc + yloc*nx + zloc*nx*ny] = 1;
				newBin = zloc;
				if (access){
					if ((nz - 1) - oldBin + newBin <= maxGap){
						for (ii = oldBin; ii < nz; ii++){
							binFill[xloc + yloc*nx + ii*nx*ny] = 1;
						}
						for (ii = 0; ii < newBin; ii++){
							binFill[xloc + yloc*nx + ii*nx*ny] = 1;
						}
					}
					access = false;
				}
				else if (bin[xloc + yloc*nx + oldBin*nx*ny] != 0){
					if (newBin - oldBin - 1 <= maxGap){
						for (ii = oldBin + 1; ii < newBin; ii++){
							binFill[xloc + yloc*nx + ii*nx*ny] = 1;
						}
					}
				}
				oldBin = newBin;
			}
		}
		// update maximum from this perspective
		sum = 0;
		for (zloc = 0; zloc < nz; zloc++){
			sum += binFill[xloc + yloc*nx + zloc*nx*ny];
		}
		sliceMax[xloc + yloc*nx] = sum;
		if (sum >= nz){
			printf("xloc yloc %4d %4d\n", xloc, yloc);
		}
	}
}
__global__ void maxZ(int *sliceMax, int *d_nx, int *d_ny, int *d_nz, int *flocDim){

	int nx, ny, nz;
	nx = *d_nx;
	ny = *d_ny;
	nz = *d_nz;

	int yloc = threadIdx.x;
	int xloc, y;

	for (xloc = 1; xloc < nx; xloc++){
		if (sliceMax[xloc + yloc*nx] > sliceMax[yloc*nx]){
			sliceMax[yloc*nx] = sliceMax[xloc + yloc*nx];
		}
	}

	__syncthreads();

	if (yloc == 0){
		for (y = 1; y < ny; y++){
			if (sliceMax[0] < sliceMax[y*nx]){
				sliceMax[0] = sliceMax[y*nx];
			}
		}
		flocDim[2] = sliceMax[0];
	}
}
