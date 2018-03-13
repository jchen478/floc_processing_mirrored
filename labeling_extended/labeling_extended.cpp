
#include "/usr/local/intel/2017/compilers_and_libraries_2017.3.191/linux/mkl/include/mkl.h"
#include "/usr/local/intel/2017/compilers_and_libraries_2017.3.191/linux/mkl/include/mkl_lapacke.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std; 

void binning(int nfib, int nseg, float *rcmx, float *rcmy, float *rcmz, 
	float dx, float dy, float dz, int nx, int ny, int nz, int maxBin,
	int *bin, int *binList, int *binFib);
bool contact(float rxmi, float rymi, float rzmi,
	float rxnj, float rynj, float rznj,
	float pxmi, float pymi, float pzmi,
	float pxnj, float pynj, float pznj,
	float sidex, float sidey, float sidez,
	float contact_cutoff, float rep_cutoff, float delta_rx, float rp, int m, int i, int n, int j);
void checkBins(int m, int nseg, int binInd, int* contactList, int* contactLabel, int *nCon,
	float *rx, float *ry, float *rz, float *px, float *py, float *pz,
	float sidex, float sidey, float sidez,
	float contact_cutoff, float rep_cutoff, float delta_rx, float rp, int *bin, int *binList,
	int *fibLabel, int maxBin); 
void parallel_sort(float sx, float sy, float sz, float pxmi, float pymi, float pzmi,
	float pxnj, float pynj, float pznj, float pdotp, float rp,
	float *xmin, float *ymin);
void gyration(float *flocRcmx, float *flocRcmy, float *flocRcmz, float *Rxx, float *Rxy, float *Rxz, float *Ryy, float *Ryz, float *Rzz,
	float *rcmx, float *rcmy, float *rcmz, 
	float *rx, float *ry, float *rz, int ind, int *nfibFloc, int *flocList, 
	int maxNfib, int nfib, int nseg);
void MLE(float *d1, float *d2, float *d3, float *V1, float *V2, float *V3,   
	float *rx, float *ry, float *rz, int ind, int *nfibFloc, int *flocList, 
	int maxNfib, int nfib, int nseg);


int main(){

	// Relevant files
	FILE *Parameters, *center_mass;
	FILE *rxfile, *ryfile, *rzfile;
	FILE *pxfile, *pyfile, *pzfile;
	FILE *Cluster_results, *Cluster;
	FILE *Cluster_fibID;
	FILE *Cluster_results_max, *Cluster_max;
	FILE *Cluster_properties;

	Parameters = fopen("Parameters.in", "r");
	center_mass = fopen("center_mass.txt", "rb");
	rxfile = fopen("rx.txt", "rb");
	ryfile = fopen("ry.txt", "rb");
	rzfile = fopen("rz.txt", "rb");
	pxfile = fopen("px.txt", "rb");
	pyfile = fopen("py.txt", "rb");
	pzfile = fopen("pz.txt", "rb");
	Cluster_results = fopen("Cluster_results.txt", "w");
	Cluster = fopen("Cluster.txt", "w");
	Cluster_results_max = fopen("Cluster_results_max.txt", "w");
	Cluster_max = fopen("Cluster_max.txt", "w");
	Cluster_fibID = fopen("Cluster_fibID.txt", "w");
	Cluster_properties = fopen("Cluster_properties.txt", "w");

	// Variables in input files
	int nfib, nseg;
	float rps, sidex, sidey, sidez, dum;
	float contact_cutoff, rep_cutoff;
	float delta_rx = 0.0; 
	float sidexo, sideyo, sidezo;

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &rps);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &contact_cutoff);
	fscanf(Parameters, "%*[^\n]%f", &rep_cutoff);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &sidex);
	fscanf(Parameters, " %f", &sidey);
	fscanf(Parameters, " %f", &sidez);
	fscanf(Parameters, "%*[^\n]%f", &dum);

	nfib *= 27;
	sidexo = sidex; 
	sideyo = sidey; 
	sidezo = sidez; 
	sidex *= 3.0;
	sidey *= 3.0;
	sidez *= 3.0;
	sidex += float(nseg)*rps; 
	sidey += float(nseg)*rps;
	sidez += float(nseg)*rps;

	// array allocation parameters
	int maxFloc = 20000; // max number of flocs
	int maxNfib = nfib; // maximum fibers in a floc
	int maxCon = 300;	// maximum number of contacts with one fiber
	int maxBin = 5500;	// maximum number of whole fibers in a bin

	// binning
	float dx, dy, dz;
	int nx, ny, nz, *bin, *binList, *binFib;
	dx = 2.0*float(nseg)*rps + 2.0 + 2.0*contact_cutoff;
	dy = 2.0*float(nseg)*rps + 2.0 + 2.0*contact_cutoff;
	dz = 2.0*float(nseg)*rps + 2.0 + 2.0*contact_cutoff;
	contact_cutoff = (contact_cutoff + 2.0)*(contact_cutoff + 2.0); 
	rep_cutoff = (rep_cutoff + 2.0)*(rep_cutoff + 2.0);
	nx = int(floorf(sidex/dx));
	ny = int(floorf(sidey/dy));
	nz = int(floorf(sidez/dz));
	dx = sidex/float(nx);	
	dy = sidey/float(ny);	
	dz = sidez/float(nz);	

	binFib = (int*)malloc(nfib*nseg*sizeof(int));
	bin = (int*)malloc(nx*ny*nz*sizeof(int));
	binList = (int*)malloc(nx*ny*nz*maxBin*sizeof(int));

	// clustering array allocation
	int *nfibFloc, *flocList, *fibLabel;
	fibLabel = (int*)malloc(nfib*sizeof(int));
	nfibFloc = (int*)malloc(maxFloc*sizeof(int));
	flocList = (int*)malloc(maxFloc*maxNfib*sizeof(int));

	int binInd, xloc, yloc, zloc, ii, xloc2, yloc2, zloc2;
	int nCon, minLabel, labelInd, oldLabel, iter; 
	int *contactList, *contactLabel;
	bool newLabel;
	contactList = (int*)malloc(maxCon*sizeof(int));
	contactLabel = (int*)malloc(maxCon*sizeof(int));

	// Read in configuration at the last frame //
	float *rx, *ry, *rz;
	float *rcmx, *rcmy, *rcmz;
	float *px, *py, *pz;
	int m, i, n, j;
	rcmx = (float*)malloc(nfib*sizeof(float));
	rcmy = (float*)malloc(nfib*sizeof(float));
	rcmz = (float*)malloc(nfib*sizeof(float));
	rx = (float*)malloc(nfib*nseg*sizeof(float));
	ry = (float*)malloc(nfib*nseg*sizeof(float));
	rz = (float*)malloc(nfib*nseg*sizeof(float));
	px = (float*)malloc(nfib*nseg*sizeof(float));
	py = (float*)malloc(nfib*nseg*sizeof(float));
	pz = (float*)malloc(nfib*nseg*sizeof(float));
	
	fread(&dum, sizeof(float), 1, center_mass); 
	fread(rcmx, sizeof(float), nfib, center_mass); 
	fread(rcmy, sizeof(float), nfib, center_mass);
	fread(rcmz, sizeof(float), nfib, center_mass);
	fread(&dum, sizeof(float), 1, rxfile);
	fread(rx, sizeof(float), nfib*nseg, rxfile);
	fread(&dum, sizeof(float), 1, ryfile);
	fread(ry, sizeof(float), nfib*nseg, ryfile);
	fread(&dum, sizeof(float), 1, rzfile);
	fread(rz, sizeof(float), nfib*nseg, rzfile);
	fread(&dum, sizeof(float), 1, pxfile);
	fread(px, sizeof(float), nfib*nseg, pxfile);
	fread(&dum, sizeof(float), 1, pyfile);
	fread(py, sizeof(float), nfib*nseg, pyfile);
	fread(&dum, sizeof(float), 1, pzfile);
	fread(pz, sizeof(float), nfib*nseg, pzfile);
	
	// initialize binning arrays
	for (i=0; i < nx*ny*nz; i++){
		bin[i] = 0; 
		for (j = 0; j < maxBin; j++){
			binList[i*maxBin+j] = 0;
		}
	}
	for (i = 0; i < nfib; i++){
		binFib[i] = 0; 
		fibLabel[i] = -1; 
	}
	for (i = 0; i < maxFloc; i++){
		nfibFloc[i] = 0;
		for (j = 0; j<maxNfib; j++){
			flocList[i*maxNfib+j] = -1;
		}
	}

	binning(nfib,nseg,rcmx,rcmy,rcmz,dx,dy,dz,nx,ny,nz,maxBin,bin,binList,binFib);	 

	labelInd = 0; 
	// loop through all fibers
	for (m = 0; m < nfib; m++){

		if (rcmx[m] == -1)
			continue; 

		// find bin index of central bin
		binInd = binFib[m];
		zloc = binInd / (nx*ny);
		yloc = (binInd - nx*ny*zloc) / nx;
		xloc = binInd - yloc*nx - zloc*nx*ny; 

		nCon = 0;
		// check central bin for contact
		checkBins (m,nseg,binInd,contactList,contactLabel,
			&nCon,rx,ry,rz,px,py,pz,sidex,sidey,sidez,
			contact_cutoff,rep_cutoff,delta_rx, rps, bin, binList, fibLabel, maxBin);

		// check neighboring bin for contact
		// case 1
		xloc2 = xloc + 1;
		yloc2 = yloc;
		zloc2 = zloc; 
		if (xloc2 != nx){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny; 
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}
		
		// case 2
		xloc2 = xloc;
		yloc2 = yloc + 1;
		zloc2 = zloc;
		if (yloc2 != ny){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}
		
		// case 3
		xloc2 = xloc;
		yloc2 = yloc;
		zloc2 = zloc + 1;
		if (zloc2 != nz){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 4
		xloc2 = xloc - 1;
		yloc2 = yloc;
		zloc2 = zloc;
		if (xloc2 >= 0){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}
		
		// case 5
		xloc2 = xloc;
		yloc2 = yloc - 1;
		zloc2 = zloc;
		if (yloc2 >= 0){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 6
		xloc2 = xloc;
		yloc2 = yloc;
		zloc2 = zloc - 1;
		if (zloc2 >= 0){ 
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 7
		xloc2 = xloc + 1;
		yloc2 = yloc + 1;
		zloc2 = zloc;
		if ((xloc2 != nx) && (yloc2 != ny)){ 
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 8
		xloc2 = xloc + 1;
		yloc2 = yloc - 1;
		zloc2 = zloc;
		if ((xloc2 != nx) && (yloc2 >= 0)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 9
		xloc2 = xloc - 1;
		yloc2 = yloc + 1;
		zloc2 = zloc;
		if ((xloc2 >= 0) && (yloc2 != ny)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}
		
		// case 10
		xloc2 = xloc - 1;
		yloc2 = yloc - 1;
		zloc2 = zloc;
		if ((xloc2 >= 0) && (yloc2 >= 0)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}
		
		// case 11
		xloc2 = xloc + 1;
		yloc2 = yloc;
		zloc2 = zloc + 1;
		if ((xloc2 != nx) && (zloc2 != nz)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}
		
		// case 12
		xloc2 = xloc - 1;
		yloc2 = yloc;
		zloc2 = zloc + 1;
		if ((xloc2 >= 0) && (zloc2 != nz)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 13
		xloc2 = xloc + 1;
		yloc2 = yloc;
		zloc2 = zloc - 1;
		if ((xloc2 != nx) && (zloc2 >= 0)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 14
		xloc2 = xloc - 1;
		yloc2 = yloc;
		zloc2 = zloc - 1;
		if ((xloc2 >= 0) && (zloc2 >= 0)){ 
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 15
		xloc2 = xloc;
		yloc2 = yloc + 1; 
		zloc2 = zloc + 1;
		if ((yloc2 != ny) && (zloc2 != nz)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 16
		xloc2 = xloc;
		yloc2 = yloc + 1;
		zloc2 = zloc - 1;
		if ((yloc2 != ny) && (zloc2 >= 0)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 18
		xloc2 = xloc;
		yloc2 = yloc - 1;
		zloc2 = zloc + 1;
		if ((yloc2 >= 0) && (zloc2 != nz)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 18
		xloc2 = xloc;
		yloc2 = yloc - 1;
		zloc2 = zloc - 1;
		if ((yloc2 >= 0) && (zloc2 >= 0)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 19
		xloc2 = xloc + 1;
		yloc2 = yloc + 1;
		zloc2 = zloc + 1;
		if ((xloc2 != nx) && (yloc2 != ny) && (zloc2 != nz)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 20
		xloc2 = xloc + 1;
		yloc2 = yloc + 1;
		zloc2 = zloc - 1;
		if ((xloc2 != nx) && (yloc2 != ny) && (zloc2 >= 0)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 21
		xloc2 = xloc + 1;
		yloc2 = yloc - 1;
		zloc2 = zloc + 1;
		if ((xloc2 != nx) && (yloc2 >= 0) && (zloc2 != nz)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 22
		xloc2 = xloc + 1;
		yloc2 = yloc - 1;
		zloc2 = zloc - 1;
		if ((xloc2 != nx) && (yloc2 >= 0) && (zloc2 >= 0)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 23
		xloc2 = xloc - 1;
		yloc2 = yloc + 1;
		zloc2 = zloc + 1;
		if ((xloc2 >= 0) && (yloc2 != ny) && (zloc2 != nz)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 24
		xloc2 = xloc - 1;
		yloc2 = yloc + 1;
		zloc2 = zloc - 1;
		if ((xloc2 >= 0) && (yloc2 != ny) && (zloc2 >= 0)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 25
		xloc2 = xloc - 1;
		yloc2 = yloc - 1;
		zloc2 = zloc + 1;
		if ((xloc2 >= 0) && (yloc2 >= 0) && (zloc2 != nz)){
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		// case 26
		xloc2 = xloc - 1;
		yloc2 = yloc - 1;
		zloc2 = zloc - 1;
		if ((xloc2 >= 0) && (yloc2 >= 0) && (zloc2 >= 0)){ 
			binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
			checkBins(m, nseg, binInd, contactList, contactLabel,
				&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
				contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);
		}

		//////////////// Labeling //////////////
		minLabel = labelInd; 
		newLabel = true;
		// when the fiber has contacts
		if (nCon != 0){
			// find the smallest label within the group of contacting fibers
			for (i = 0; i < nCon; i++){
				// don't need new label
				if (contactLabel[i] > -1){
					if (contactLabel[i] < minLabel){
						minLabel = contactLabel[i];
					}
					newLabel = false;
				}
			}
			if (newLabel){
				minLabel = labelInd;
				labelInd++;
				if (minLabel >= maxFloc) printf("allocate more space for storing flocs, i.e. maxFloc\n"); 
			}
			// label all fibers in group
			// add fibers to list of fibers in a floc
			oldLabel = fibLabel[m]; 
			if (oldLabel == -1){
				flocList[minLabel*maxNfib + nfibFloc[minLabel]] = m;
				nfibFloc[minLabel]++;
			}
			// combine to new group
			else if (oldLabel != minLabel && nfibFloc[oldLabel] > 0){
				iter = nfibFloc[oldLabel];
				for (ii = 0; ii < iter ; ii++){
					flocList[minLabel*maxNfib + nfibFloc[minLabel]] =
						flocList[oldLabel*maxNfib + ii];
					fibLabel[flocList[oldLabel*maxNfib + ii]] = minLabel;
					nfibFloc[minLabel]++;
					nfibFloc[oldLabel] = -1; 
				}
			}
			fibLabel[m] = minLabel;
			for (i = 0; i < nCon; i++){
				n = contactList[i];
				oldLabel = fibLabel[n];
				// add a new fiber to group
				if (oldLabel == -1){
					flocList[minLabel*maxNfib + nfibFloc[minLabel]] = n;
					nfibFloc[minLabel]++;
				}
				// combine to new group
				else if (oldLabel != minLabel && nfibFloc[oldLabel] > 0){
					iter = nfibFloc[oldLabel];	
					for (ii = 0; ii < iter; ii++){
						flocList[minLabel*maxNfib + nfibFloc[minLabel]] =
							flocList[oldLabel*maxNfib + ii];
						fibLabel[flocList[oldLabel*maxNfib + ii]] = minLabel;
						nfibFloc[minLabel]++;
						nfibFloc[oldLabel] = -1;
					}
				}
				fibLabel[n] = minLabel;
			}		
		}		
	}

	// the largest floc
	int maxfloc, ind;
	maxfloc = 0;

	// printf floc list
	int nCluster = 0; 
	for (i = 0; i < maxFloc; i++){
		if (nfibFloc[i] >= 2){
			//printf("%4d %4d: ", i, nfibFloc[i]);
			fprintf(Cluster, "%6d %6d ", i, nfibFloc[i]);
			fprintf(Cluster_results, "%6d %6d\n", nCluster, nfibFloc[i]);
			nCluster++; 
			for (j = 0; j < nfibFloc[i]; j++){
				//printf("%4d ", flocList[i*maxNfib + j]); 
				fprintf(Cluster, "%8d ", flocList[i*maxNfib + j]);
			}
			//printf("\n"); 
			fprintf(Cluster, "\n");
		}
		if (nfibFloc[i] > maxfloc){
			ind = i; 
			maxfloc = nfibFloc[i]; 
		}
		
	}
	fprintf(Cluster_results, "Number of clusters: \n");
	fprintf(Cluster_results, "%6d\n", nCluster);

	for (m = 0; m < nfib; m++){
		fprintf(Cluster_fibID, "%8d %4d\n", m + 1, fibLabel[m]);
	}

	for (j = 0; j < nfibFloc[ind]; j++){
		fprintf(Cluster_max, "%8d\n", flocList[ind*maxNfib + j]);
	}
	fprintf(Cluster_results_max, "%6d %6d\n", ind, nfibFloc[ind]);

	float Rxx, Rxy, Rxz; 
	float Ryy, Ryz, Rzz; 
	float flocRcmx, flocRcmy, flocRcmz; 
	Rxx = 0.0; Rxy = 0.0; Rxz = 0.0; 
	Ryy = 0.0; Ryz = 0.0; Rzz = 0.0; 
	gyration(&flocRcmx, &flocRcmy, &flocRcmz, &Rxx, &Rxy, &Rxz, &Ryy, &Ryz, &Rzz,
		rcmx, rcmy, rcmz, rx, ry, rz, ind, nfibFloc, flocList, maxNfib, nfib, nseg); 

	float R[9], W[3];
	int info;
	R[0] = Rxx;
	R[1] = Rxy;
	R[2] = Rxz;
	R[3] = 0.0; 
	R[4] = Ryy; 
	R[5] = Ryz;
	R[6] = 0.0;
	R[7] = 0.0; 
	R[8] = Rzz; 
	
	info = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'V','U',3,R,3,W);
	if (info != 0){
		fprintf(Cluster_properties, "Floc index\n %4d\n",ind);
		fprintf(Cluster_properties, "Number of fibers in floc\n %6d\n", nfibFloc[ind]);
		fprintf(Cluster_properties, "Floc center of mass\n");
		fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n", flocRcmx, flocRcmy, flocRcmz);
		fprintf(Cluster_properties,"matrix is not solvable: info %4d\n", info);
		fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n%15.6f %15.6f %15.6f\n%15.6f %15.6f %15.6f\n",
			Rxx, Rxy, Rxz, Rxy, Ryy, Ryz, Rxz, Ryz, Rzz);
		fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n%15.6f %15.6f %15.6f\n%15.6f %15.6f %15.6f\n",
			R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8]);
		return 1;
	}
	
	float d1, d2, d3; 
	float V1[3], V2[3], V3[3]; 
	V1[0] = R[0]; V2[0] = R[1]; V3[0] = R[2]; 
	V1[1] = R[3]; V2[1] = R[4]; V3[1] = R[5];
	V1[2] = R[6]; V2[2] = R[7]; V3[2] = R[8];
	
	// find the maximum linear extent in each principal direction
	MLE(&d1, &d2, &d3, V1, V2, V3, rx, ry, rz, ind, nfibFloc, flocList, maxNfib, nfib, nseg); 

	int per;
	
	per = 0; 	
	if (d1 >= sidexo || d2 >= sideyo || d3 >= sidezo)
		per = 1; 

	fprintf(Cluster_properties, "Binning properties\n");
	fprintf(Cluster_properties, "=======================================\n");
	fprintf(Cluster_properties, "Original Sidex Sidey Sidez\n%15.6f %15.6f %15.6f\n", sidexo, sideyo, sidezo); 
	fprintf(Cluster_properties, "Sidex Sidey Sidez\n%15.6f %15.6f %15.6f\n", sidex, sidey, sidez); 
	fprintf(Cluster_properties, "dx dy dz\n%15.6f %15.6f %15.6f\n", dx, dy, dz);
	fprintf(Cluster_properties, "nx ny nz\n%4d %4d %4d\n", nx, ny, nz);
	fprintf(Cluster_properties, "\nCluster properties\n");
	fprintf(Cluster_properties, "=======================================\n");
	fprintf(Cluster_properties, "Floc index\n %4d\n",ind);
	fprintf(Cluster_properties, "Number of fibers in floc\n %6d\n", nfibFloc[ind]);
	fprintf(Cluster_properties, "Floc center of mass\n");
	fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n", flocRcmx, flocRcmy, flocRcmz);
	fprintf(Cluster_properties, "Gyration tensor\n");
	fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n%15.6f %15.6f %15.6f\n%15.6f %15.6f %15.6f\n",
		Rxx, Rxy, Rxz, Rxy, Ryy, Ryz, Rxz, Ryz, Rzz);
	fprintf(Cluster_properties, "Principle moments\n");
	fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n",W[0],W[1],W[2]);
	fprintf(Cluster_properties, "Principal directions\n");
	fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n%15.6f %15.6f %15.6f\n%15.6f %15.6f %15.6f\n",
		R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7],R[8]);
	fprintf(Cluster_properties, "Principle moments square root\n");
	fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n",sqrt(W[0]),sqrt(W[1]),sqrt(W[2]));
	fprintf(Cluster_properties, "Radius of gyration\n");
	fprintf(Cluster_properties, "%15.6f\n",sqrt(W[0]+W[1]+W[2]));
	fprintf(Cluster_properties, "Maximum linear extent in principal directions\n");
	fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n",d1,d2,d3);
	fprintf(Cluster_properties, "Percolation (1: yes, 0: no)\n");
	fprintf(Cluster_properties, "%4d\n",per);
	fprintf(Cluster_properties, "Asphericity\n");
	fprintf(Cluster_properties, "%15.6f\n",W[2]-0.5*(W[0]+W[1]));
	fprintf(Cluster_properties, "Acylindricity\n");
	fprintf(Cluster_properties, "%15.6f\n",W[1]-W[0]);
	fprintf(Cluster_properties, "Relative shape anisotropy k2\n");
	fprintf(Cluster_properties, "%15.6f\n",1.5*(W[1]*W[1]+W[2]*W[2]+W[0]*W[0])/(W[0]+W[1]+W[2])-0.5);
	fprintf(Cluster_properties, "Invariants\n");
	fprintf(Cluster_properties, "%15.6f %15.6f %15.6f\n",
		(W[0]+W[1]+W[2])/W[2],(W[0]*W[1]+W[1]*W[2]+W[0]*W[2])/(W[2]*W[2]),W[0]*W[1]*W[2]/(W[2]*W[2]*W[2]));

	
	free(rcmx); free(rcmy); free(rcmz);
	free(rx); free(ry); free(rz);
	free(px); free(py); free(pz);
	free(bin); free(binList); free(binFib);  
	free(nfibFloc); free(flocList); 
	free(fibLabel); 	
	free(contactList); 
	free(contactLabel); 

	fclose(Parameters); 
	fclose(center_mass); 
	fclose(Cluster); 
	fclose(Cluster_results);
	fclose(Cluster_max);
	fclose(Cluster_properties); 
	fclose(Cluster_results_max);
	fclose(Cluster_fibID);
	fclose(rxfile);
	fclose(rzfile);
	fclose(ryfile);
	fclose(pxfile);
	fclose(pzfile);
	fclose(pyfile);

	return 0;
}

void MLE(float *d1, float *d2, float *d3, float *V1, float *V2, float *V3, 
	float *rx, float *ry, float *rz, int ind, int *nfibFloc, int *flocList, 
	int maxNfib, int nfib, int nseg){

	int m, n, i, j, k, c, mi, nj, flocCount; 
	float dd1, dd2, dd3, dx, dy, dz;

	flocCount = nfibFloc[ind]; 
	*d1 = 0.0; 
	*d2 = 0.0; 
	*d3 = 0.0; 

	for (k = 0; k < flocCount-1; k++){
		m = flocList[ind*maxNfib + k];
		for (c = k+1; c < flocCount; c++){
			n = flocList[ind*maxNfib + c];  
			for (i = 0; i < nseg; i++){
				mi = m*nseg+i;
				for (j = 0; j < nseg; j++){
					nj = n*nseg + j; 
					dx = rx[mi] - rx[nj]; 
					dy = ry[mi] - ry[nj];
					dz = rz[mi] - rz[nj]; 
					// find dot product
					dd1 = abs(dx*V1[0]+dy*V1[1]+dz*V1[2]);
					dd2 = abs(dx*V2[0]+dy*V2[1]+dz*V2[2]);
					dd3 = abs(dx*V3[0]+dy*V3[1]+dz*V3[2]);
					if (dd1 > *d1) *d1 = dd1;
					if (dd2 > *d2) *d2 = dd2;
					if (dd3 > *d3) *d3 = dd3;
				}
			}
		}
	}
}

void gyration(float *flocRcmx, float *flocRcmy, float *flocRcmz, 
	float *Rxx, float *Rxy, float *Rxz, float *Ryy, float *Ryz, float *Rzz,
	float *rcmx, float *rcmy, float *rcmz,
	float *rx, float *ry, float *rz, int ind, int *nfibFloc, int *flocList,
	int maxNfib, int nfib, int nseg){

	float dx, dy, dz; 
	int m, i, k, mi, flocCount; 

	flocCount = nfibFloc[ind]; 
	*flocRcmx = 0.0; 
	*flocRcmy = 0.0; 
	*flocRcmz = 0.0; 

	// find floc center of mass
	for (k = 0; k < flocCount; k++){
		m = flocList[ind*maxNfib + k];
		*flocRcmx += rcmx[m]; 
		*flocRcmy += rcmy[m];
		*flocRcmz += rcmz[m];
	}
	*flocRcmx /= float(flocCount); 
	*flocRcmy /= float(flocCount);
	*flocRcmz /= float(flocCount);

	for (k = 0; k < flocCount; k++){
		m = flocList[ind*maxNfib + k];
		for (i = 0; i < nseg; i++){
			mi = m*nseg + i; 
			dx = (rx[mi] - *flocRcmx); 
			dy = (ry[mi] - *flocRcmy);
			dz = (rz[mi] - *flocRcmz);
			*Rxx += dx*dx; 
			*Ryy += dy*dy; 
			*Rzz += dz*dz; 
			*Rxy += dx*dy; 
			*Rxz += dx*dz; 
			*Ryz += dy*dz; 
		}
	}
	*Rxx /= float(flocCount*nseg); 
	*Ryy /= float(flocCount*nseg);
	*Rzz /= float(flocCount*nseg);
	*Rxy /= float(flocCount*nseg);
	*Rxz /= float(flocCount*nseg);
	*Ryz /= float(flocCount*nseg);
}

void checkBins(int m, int nseg, int binInd, int* contactList, int* contactLabel, int *nCon,
	float *rx, float *ry, float *rz, float *px, float *py, float *pz,
	float sidex, float sidey, float sidez,
	float contact_cutoff, float rep_cutoff, float delta_rx,
	float rp, int *bin, int *binList,
	int *fibLabel,int maxBin){

	int nof, b, n, mi, nj, i, j; 
	bool isContact, repeat;

	// determine contacts wihtin central bin
	nof = bin[binInd];
	for (b = 0; b < nof; b++){
		n = binList[binInd*maxBin + b];
		isContact = false;
		repeat = false; 
		for (i = 0; i < *nCon; i++){
			if (n == contactList[i]){
				repeat = true; 
			}
		}
		// excluding comparison with itself
		if (!repeat && n != m){
			// loop through segments
			for (i = 0; i < nseg; i++){
				mi = m*nseg + i;
				for (j = 0; j < nseg; j++){
					nj = n*nseg + j;					
						isContact =
							contact(rx[mi], ry[mi], rz[mi], rx[nj], ry[nj], rz[nj],
							px[mi], py[mi], pz[mi], px[nj], py[nj], pz[nj],
							sidex, sidey, sidez, contact_cutoff, rep_cutoff, 
							delta_rx, rp, m, i, n, j);
					if (isContact){
						break;
					}
				}
				if (isContact) break;
			}
			if (isContact){
				// save to list of contacting fibers with m
				contactList[*nCon] = n;
				contactLabel[*nCon] = fibLabel[n];
				*nCon = *nCon + 1;
			}
		}
	}
}

bool contact(float rxmi, float rymi, float rzmi,
	float rxnj, float rynj, float rznj,
	float pxmi, float pymi, float pzmi,
	float pxnj, float pynj, float pznj,
	float sidex, float sidey, float sidez,
	float contact_cutoff, float rep_cutoff, float delta_rx, 
	float rp, int m, int i, int n, int j){

	float sxx, syy, szz, corx, cory, corz;
	float pdotp, xmin, ymin, dx, dy, dz, sep;
	float xi[9], yj[9], sep_tmp, gij;
	int ipos, ith; 

	sxx = rxnj - rxmi;
	syy = rynj - rymi;
	szz = rznj - rzmi;
	pdotp = pxmi*pxnj + pymi*pynj + pzmi*pznj;
	xmin = (-(pxnj * sxx + pynj * syy + pznj * szz)* pdotp
		+ (pxmi * sxx + pymi * syy + pzmi * szz))
		/ (1.0 - pdotp*pdotp);
	ymin = ((pxmi * sxx + pymi * syy + pzmi * szz)* pdotp
		- (pxnj * sxx + pynj * syy + pznj * szz))
		/ (1.0 - pdotp*pdotp);

	dx = rxnj + ymin*pxnj - rxmi - xmin*pxmi;
	dy = rynj + ymin*pynj - rymi - xmin*pymi;
	dz = rznj + ymin*pznj - rzmi - xmin*pzmi;
	sep = dx*dx + dy*dy + dz*dz;

	ipos = 8;
	yj[0] = rp;
	xi[0] = pxmi*sxx + pymi*syy + pzmi*szz + yj[0] * pdotp;
	yj[1] = -rp;
	xi[1] = pxmi*sxx + pymi*syy + pzmi*szz + yj[1] * pdotp;
	xi[2] = rp;
	yj[2] = -(pxnj*sxx + pynj*syy + pznj*szz) + xi[2] * pdotp;
	xi[3] = -rp;
	yj[3] = -(pxnj*sxx + pynj*syy + pznj*szz) + xi[3] * pdotp;
	xi[4] = rp;    yj[4] = rp;
	xi[5] = rp;    yj[5] = -rp;
	xi[6] = -rp;   yj[6] = rp;
	xi[7] = -rp;   yj[7] = -rp;
	xi[8] = xmin;  yj[8] = ymin;

		// Check if segments are parallel
	if (fabsf(pdotp*pdotp - 1.0) <= 1.0e-6) {
		parallel_sort(sxx, syy, szz, pxmi, pymi, pzmi,
			pxnj, pynj, pznj, pdotp, rp, &xmin, &ymin);
		sep = (sxx + ymin*pxnj - xmin*pxmi)*(sxx + ymin*pxnj - xmin*pxmi) +
			(syy + ymin*pynj - xmin*pymi)*(syy + ymin*pynj - xmin*pymi) +
			(szz + ymin*pznj - xmin*pzmi)*(szz + ymin*pznj - xmin*pzmi);
	}
	else if (sep < rep_cutoff && (fabsf(xmin) >= rp || fabsf(ymin) >= rp)){
		sep = 1000.0;
		// check which end-side or end-end separation
		// is the smallest
		for (ith = 0; ith < 8; ith++){
			sep_tmp = (sxx + yj[ith] * pxnj - xi[ith] * pxmi)*(sxx + yj[ith] * pxnj - xi[ith] * pxmi) +
				(syy + yj[ith] * pynj - xi[ith] * pymi)*(syy + yj[ith] * pynj - xi[ith] * pymi) +
				(szz + yj[ith] * pznj - xi[ith] * pzmi)*(szz + yj[ith] * pznj - xi[ith] * pzmi);
			if (sep_tmp < sep && fabsf(xi[ith]) <= rp && fabsf(yj[ith]) <= rp){
				sep = sep_tmp;
				ipos = ith;
			}
		}
		xmin = xi[ipos];
		ymin = yj[ipos];
	}
	gij = sqrtf(sep);
	if (sep < contact_cutoff){
		return true; 
	}
	return false; 
}
void binning(int nfib, int nseg, float *rcmx, float *rcmy, float *rcmz, 
	float dx, float dy, float dz, int nx, int ny, int nz, int maxBin,
	int *bin, int *binList, int *binFib){

	int mi; 
	int xloc, yloc, zloc;
	int ind;

	for (mi = 0; mi < nfib; mi++){

		if (rcmx[mi] == -1){
			continue; 
		}
		
		xloc = floorf(rcmx[mi] / dx); 
		yloc = floorf(rcmy[mi] / dy);
		zloc = floorf(rcmz[mi] / dz);
		// check for out of bound access
		if (xloc < 0 || xloc >= nx){
			if (xloc == nx){
				 xloc--;
				 break;
			}
			printf("x index out of bound: xloc %4d dx %f nx %4d rcmx %f \n", xloc,dx,nx,rcmx[mi]); 
		}
		if (yloc < 0 || yloc >= ny){
			if (yloc == ny){
				 yloc--;
				 break;
			}
			printf("y index out of bound: yloc %4d dy %f ny %4d rcmy %f \n", yloc,dy,ny,rcmy[mi]); 
		}
		if (zloc < 0 || zloc >= nz){
			if (zloc == nz){
				 zloc--;
				 break;
			}
			printf("z index out of bound: zloc %4d dz %f nz %4d rcmz %f \n", zloc,dz,nz,rcmz[mi]); 
		}
		ind = xloc+yloc*nx+zloc*ny*nz;	
		binList[ind*maxBin+bin[ind]] = mi;
		bin[ind] = bin[ind] + 1;	
		binFib[mi] = ind;
		if(bin[ind] >= maxBin) printf("larger binlist!\n");
	}
}

void parallel_sort(float sx, float sy, float sz, float pxmi, float pymi, float pzmi,
	float pxnj, float pynj, float pznj, float pdotp, float rp,
	float *xmin, float *ymin){

	float posneg, pn2, dist, sijp, sijm, sjip, sjim;

	// The different end point to fiber contact points
	sijp = pxmi*sx + pymi*sy + pzmi*sz + rp*pdotp;
	sijm = pxmi*sx + pymi*sy + pzmi*sz - rp*pdotp;
	sjip = -(pxnj*sx + pynj*sy + pznj*sz) + rp*pdotp;
	sjim = -(pxnj*sx + pynj*sy + pznj*sz) - rp*pdotp;

	// for fiber i
	if (fabsf(sijp) < fabsf(sijm)){
		*xmin = sijp;
		posneg = 1.0;
	}
	else if (fabsf(sijp) > fabsf(sijm)){
		*xmin = sijm;
		posneg = -1.0;
	}
	else{
		*xmin = 0.0;
		posneg = 0.0;
	}
	if (*xmin >= rp){
		*xmin = rp;
	}
	if (*xmin <= -rp){
		*xmin = -rp;
	}
	// for fiber j
	if (fabsf(sjip) < fabsf(sjim)){
		*ymin = sjip;
	}
	else if (fabsf(sjip) > fabsf(sjim)){
		*ymin = sjim;
	}
	else{
		*ymin = 0.0;
		posneg = 0.0;
	}
	if (*ymin >= rp){
		*ymin = rp;
	}
	if (*ymin <= -rp){
		*ymin = -rp;
	}	
	if (fabsf(*xmin) < rp && fabsf(*ymin) < rp){
		if (pdotp > 0.0){
			pn2 = 1.0;
		}
		else{
			pn2 = -1.0;
		}
		dist = (rp + posneg**xmin) / 2.0;
		*xmin = *xmin - posneg*dist;
		*ymin = *ymin + posneg*pn2*dist;		
	}
}
