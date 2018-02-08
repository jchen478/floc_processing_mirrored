#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void shift(int nfib, int nseg, float *rx, float *ry, float *rz, 
	float *rcmx, float *rcmy, float *rcmz,
	float sidex, float sidey, float sidez);
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

int main(){



	// Relevant files
	FILE *Parameters, *center_mass;
	FILE *rxfile, *ryfile, *rzfile;
	FILE *pxfile, *pyfile, *pzfile;
	FILE *Cluster_results, *Cluster;
	FILE *Cluster_fibID;

	Parameters = fopen("Parameters.in", "r");
	center_mass = fopen("../../center_mass.txt", "rb");
	rxfile = fopen("../../rx.txt", "rb");
	ryfile = fopen("../../ry.txt", "rb");
	rzfile = fopen("../../rz.txt", "rb");
	pxfile = fopen("../../px.txt", "rb");
	pyfile = fopen("../../py.txt", "rb");
	pzfile = fopen("../../pz.txt", "rb");
	Cluster_results = fopen("Cluster_results.txt", "w");
	Cluster = fopen("Cluster.txt", "w");
	Cluster_fibID = fopen("Cluster_fibID.txt", "w");

	// Variables in input files
	int nfib, nseg, config_write;
	float dt, strain, mustat, mukin, dum;
	float rps, sidex, sidey, sidez;
	float contact_cutoff, rep_cutoff;
	float delta_rx; 

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &rps);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &mustat);
	fscanf(Parameters, "%*[^\n]%f", &mukin);
	fscanf(Parameters, "%*[^\n]%f", &contact_cutoff);
	fscanf(Parameters, "%*[^\n]%f", &rep_cutoff);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dt);
	fscanf(Parameters, "%*[^\n]%f", &strain);
	fscanf(Parameters, "%*[^\n]%f", &sidex);
	fscanf(Parameters, " %f", &sidey);
	fscanf(Parameters, " %f", &sidez);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &config_write);

	delta_rx = strain*sidex; 
	delta_rx -= lroundf(delta_rx / sidex)*sidex;

	// array allocation parameters
	int maxFloc = 5000; // max number of flocs
	int maxNfib = nfib; // maximum fibers in a floc
	int maxCon = 500;	// maximum number of contacts with one fiber
	int maxBin = 1000;	// maximum number of whole fibers in a bin

	// binning
	float dx, dy, dz;
	int nx, ny, nz, *bin, *binList, *binFib;
	dx = 1.5*float(nseg)*rps + 2.0*contact_cutoff;
	dy = 1.5*float(nseg)*rps + 2.0*contact_cutoff;
	dz = 1.5*float(nseg)*rps + 2.0*contact_cutoff;
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
	int step, nConfig;
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
	nConfig = int(strain / dt / float(config_write)) + 1;
	
	for (step = 0; step < nConfig; step++){
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
	}
	
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

	shift(nfib,nseg,rx,ry,rz,rcmx,rcmy,rcmz,sidex,sidey,sidez);	
	binning(nfib,nseg,rcmx,rcmy,rcmz,dx,dy,dz,nx,ny,nz,maxBin,bin,binList,binFib);	 

	labelInd = 0; 
	// loop through all fibers
	for (m = 0; m < nfib; m++){

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
		if (xloc2 == nx) xloc2 -= nx; 
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny; 
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 2
		xloc2 = xloc;
		yloc2 = yloc + 1;
		zloc2 = zloc;
		if (yloc2 == ny) yloc2 -= ny;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 3
		xloc2 = xloc;
		yloc2 = yloc;
		zloc2 = zloc + 1;
		if (zloc2 == nz) zloc2 -= nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 4
		xloc2 = xloc - 1;
		yloc2 = yloc;
		zloc2 = zloc;
		if (xloc2 < 0) xloc2 += nx;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 5
		xloc2 = xloc;
		yloc2 = yloc - 1;
		zloc2 = zloc;
		if (yloc2 < 0) yloc2 += ny;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 6
		xloc2 = xloc;
		yloc2 = yloc;
		zloc2 = zloc - 1;
		if (zloc2 < 0) zloc2 += nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 7
		xloc2 = xloc + 1;
		yloc2 = yloc + 1;
		zloc2 = zloc;
		if (xloc2 == nx) xloc2 -= nx;
		if (yloc2 == ny) yloc2 -= ny;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 8
		xloc2 = xloc + 1;
		yloc2 = yloc - 1;
		zloc2 = zloc;
		if (xloc2 == nx) xloc2 -= nx;
		if (yloc2 < 0) yloc2 += ny;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 9
		xloc2 = xloc - 1;
		yloc2 = yloc + 1;
		zloc2 = zloc;
		if (xloc2 < 0) xloc2 += nx;
		if (yloc2 == ny) yloc2 -= ny;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 10
		xloc2 = xloc - 1;
		yloc2 = yloc - 1;
		zloc2 = zloc;
		if (xloc2 < 0) xloc2 += nx;
		if (yloc2 < 0) yloc2 += ny;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 11
		xloc2 = xloc + 1;
		yloc2 = yloc;
		zloc2 = zloc + 1;
		if (xloc2 == nx) xloc2 -= nx;
		if (zloc2 == nz) zloc2 -= nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 12
		xloc2 = xloc - 1;
		yloc2 = yloc;
		zloc2 = zloc + 1;
		if (xloc2 < 0) xloc2 += nx;
		if (zloc2 == nz) zloc2 -= nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 13
		xloc2 = xloc + 1;
		yloc2 = yloc;
		zloc2 = zloc - 1;
		if (xloc2 == nx) xloc2 -= nx;
		if (zloc2 < 0) zloc2 += nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 14
		xloc2 = xloc - 1;
		yloc2 = yloc;
		zloc2 = zloc - 1;
		if (xloc2 < 0) xloc2 += nx;
		if (zloc2 < 0) zloc2 += nz; 
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 15
		xloc2 = xloc;
		yloc2 = yloc + 1; 
		zloc2 = zloc + 1;
		if (yloc2 == ny) yloc2 -= ny;
		if (zloc2 == nz) zloc2 -= nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 16
		xloc2 = xloc;
		yloc2 = yloc + 1;
		zloc2 = zloc - 1;
		if (yloc2 == ny) yloc2 -= ny;
		if (zloc2 < 0) zloc2 += nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 18
		xloc2 = xloc;
		yloc2 = yloc - 1;
		zloc2 = zloc + 1;
		if (yloc2 < 0) yloc2 += ny;
		if (zloc2 == nz) zloc2 -= nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 18
		xloc2 = xloc;
		yloc2 = yloc - 1;
		zloc2 = zloc - 1;
		if (yloc2 < 0) yloc2 += ny;
		if (zloc2 < 0) zloc2 += nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 19
		xloc2 = xloc + 1;
		yloc2 = yloc + 1;
		zloc2 = zloc + 1;
		if (xloc2 == nx) xloc2 -= nx;
		if (yloc2 == ny) yloc2 -= ny;
		if (zloc2 == nz) zloc2 -= nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 20
		xloc2 = xloc + 1;
		yloc2 = yloc + 1;
		zloc2 = zloc - 1;
		if (xloc2 == nx) xloc2 -= nx;
		if (yloc2 == ny) yloc2 -= ny;
		if (zloc2 < 0) zloc2 += nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 21
		xloc2 = xloc + 1;
		yloc2 = yloc - 1;
		zloc2 = zloc + 1;
		if (xloc2 == nx) xloc2 -= nx;
		if (yloc2 < 0) yloc2 += ny;
		if (zloc2 == nz) zloc2 -= nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 22
		xloc2 = xloc + 1;
		yloc2 = yloc - 1;
		zloc2 = zloc - 1;
		if (xloc2 == nx) xloc2 -= nx;
		if (yloc2 < 0) yloc2 += ny;
		if (zloc2 < 0) zloc2 += nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 23
		xloc2 = xloc - 1;
		yloc2 = yloc + 1;
		zloc2 = zloc + 1;
		if (xloc2 < 0) xloc2 += nx;
		if (yloc2 == ny) yloc2 -= ny;
		if (zloc2 == nz) zloc2 -= nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 24
		xloc2 = xloc - 1;
		yloc2 = yloc + 1;
		zloc2 = zloc - 1;
		if (xloc2 < 0) xloc2 += nx;
		if (yloc2 == ny) yloc2 -= ny;
		if (zloc2 < 0) zloc2 += nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 25
		xloc2 = xloc - 1;
		yloc2 = yloc - 1;
		zloc2 = zloc + 1;
		if (xloc2 < 0) xloc2 += nx;
		if (yloc2 < 0) yloc2 += ny;
		if (zloc2 == nz) zloc2 -= nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

		// case 26
		xloc2 = xloc - 1;
		yloc2 = yloc - 1;
		zloc2 = zloc - 1;
		if (xloc2 < 0) xloc2 += nx;
		if (yloc2 < 0) yloc2 += ny;
		if (zloc2 < 0) zloc2 += nz;
		binInd = xloc2 + yloc2*nx + zloc2*nx*ny;
		checkBins(m, nseg, binInd, contactList, contactLabel,
			&nCon, rx, ry, rz, px, py, pz, sidex, sidey, sidez,
			contact_cutoff, rep_cutoff, delta_rx, rps, bin, binList, fibLabel, maxBin);

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

	// printf floc list
	int nCluster = 0; 
	for (i = 0; i < maxFloc; i++){
		if (nfibFloc[i] > 0){
			fprintf(Cluster, "%6d %6d ", i, nfibFloc[i]);
			fprintf(Cluster_results, "%6d %6d\n", i, nfibFloc[i]);
			nCluster++; 
			for (j = 0; j < nfibFloc[i]; j++){
				fprintf(Cluster, "%8d ", flocList[i*maxNfib + j]);
			}
			fprintf(Cluster, "\n");
		}
		
	}
	fprintf(Cluster_results, "Number of clusters: \n");
	fprintf(Cluster_results, "%6d\n", nCluster);

	for (m = 0; m < nfib; m++){
		fprintf(Cluster_fibID, "%8d %4d\n", m + 1, fibLabel[m]);
	}

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
	fclose(Cluster_fibID);
	fclose(rxfile);
	fclose(rzfile);
	fclose(ryfile);
	fclose(pxfile);
	fclose(pzfile);
	fclose(pyfile);


	return 0;
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
					if (isContact) break;
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
	float rxmi_shift, rymi_shift, rzmi_shift;
	float xi[9], yj[9], sep_tmp, gij;
	int ipos, ith; 

	// find minimum image (for shear flow system)
	sxx = rxnj - rxmi;
	syy = rynj - rymi;
	szz = rznj - rzmi;
	cory = roundf(syy / sidey);
	corz = roundf(szz / sidez);
	sxx = sxx - corz*delta_rx;
	corx = roundf(sxx / sidex);
	sxx = sxx - corx*sidex;
	syy = syy - cory*sidey;
	szz = szz - corz*sidez;
	rxmi_shift = rxnj - sxx;
	rymi_shift = rynj - syy;
	rzmi_shift = rznj - szz;
	pdotp = pxmi*pxnj + pymi*pynj + pzmi*pznj;
	xmin = (-(pxnj * sxx + pynj * syy + pznj * szz)* pdotp
		+ (pxmi * sxx + pymi * syy + pzmi * szz))
		/ (1.0 - pdotp*pdotp);
	ymin = ((pxmi * sxx + pymi * syy + pzmi * szz)* pdotp
		- (pxnj * sxx + pynj * syy + pznj * szz))
		/ (1.0 - pdotp*pdotp);

	dx = rxnj + ymin*pxnj - rxmi_shift - xmin*pxmi;
	dy = rynj + ymin*pynj - rymi_shift - xmin*pymi;
	dz = rznj + ymin*pznj - rzmi_shift - xmin*pzmi;
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

void shift(int nfib, int nseg, float *rx, float *ry, float *rz, 
	float *rcmx, float *rcmy, float *rcmz,
	float sidex, float sidey, float sidez){

	int mi;
	
	for (mi = 0; mi < nfib; mi++){
		// shift coordinates inside box
		if (rcmx[mi] < -sidex / 2.0 || rcmx[mi] > sidex / 2){
			rcmx[mi] -= copysignf(sidex, rcmx[mi]);
		}
		if (rcmy[mi] < -sidey / 2.0 || rcmy[mi] > sidey / 2){
			rcmy[mi] -= copysignf(sidey, rcmy[mi]);
		}
		if (rcmz[mi] < -sidez / 2.0 || rcmz[mi] > sidez / 2){
			rcmz[mi] -= copysignf(sidez, rcmz[mi]);
		}
		// shift coordinates so all are positive
		rcmx[mi] += sidex / 2.0;
		rcmy[mi] += sidey / 2.0;
		rcmz[mi] += sidez / 2.0;
	}
}
void binning(int nfib, int nseg, float *rcmx, float *rcmy, float *rcmz, 
	float dx, float dy, float dz, int nx, int ny, int nz, int maxBin,
	int *bin, int *binList, int *binFib){

	int mi; 
	int xloc, yloc, zloc;
	int ind;

	for (mi = 0; mi < nfib; mi++){
		
		xloc = floorf(rcmx[mi] / dx); 
		yloc = floorf(rcmy[mi] / dy);
		zloc = floorf(rcmz[mi] / dz);
		if (xloc == nx) xloc--; 
		if (yloc == ny) yloc--; 
		if (zloc == nz) zloc--; 
		// check for out of bound access
		if (xloc < 0 || xloc >= nx){
			printf("x index out of bound: mi %6d xloc %4d dx %f nx %4d rcmx %f \n", mi, xloc,dx,nx,rcmx[mi]); 
		}
		if (yloc < 0 || yloc >= ny){
			printf("y index out of bound: mi %6d yloc %4d dy %f ny %4d rcmy %f \n", mi, yloc,dy,ny,rcmy[mi]); 
		}
		if (zloc < 0 || zloc >= nz){
			printf("z index out of bound: mi %6d zloc %4d dz %f nz %4d rcmz %f \n", mi, zloc,dz,nz,rcmz[mi]); 
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
