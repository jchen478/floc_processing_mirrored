#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <math.h>

void regrow(float *rx, float *ry, float *rz, float *px, float *py, float *pz,
	    float *rcmx, float *rcmy, float *rcmz, float rps, int nfib, int nseg);
void shift(float *rx, float *ry, float *rz, float sidex, float sidey, float sidez,
           int nfib, int nseg);
void extend(float *rx, float *ry, float *rz, float *rxE, float *ryE, float *rzE,
            float sidex, float sidey, float sidez, float cellx, float celly,float cellz,
            int pos, int nfib, int nseg); 
float find_min(float *rx, int nfib, int nseg);
void shift_origin(float *rx, float offset, int nfib, int nseg); 

int main(){

	// Relevant files
	FILE *Parameters;
	FILE *Centers_of_Mass, *Euler_Parameters;
	FILE *rxfile, *ryfile, *rzfile;
	FILE *pxfile, *pyfile, *pzfile;
	FILE *center_mass;

	Parameters = fopen("Parameters.in", "r");
	Centers_of_Mass = fopen("Centers_of_Mass.in", "r");
	Euler_Parameters = fopen("Euler_Parameters.in", "r");
	center_mass = fopen("center_mass.txt", "wb");
	rxfile = fopen("rx.txt", "wb");
	ryfile = fopen("ry.txt", "wb");
	rzfile = fopen("rz.txt", "wb");
	pxfile = fopen("px.txt", "wb");
	pyfile = fopen("py.txt", "wb");
	pzfile = fopen("pz.txt", "wb");

	// Variables in input files
	int nfib, nseg;
	float dum, rps, sidex, sidey, sidez;
	float contact_cutoff;

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &rps);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &contact_cutoff);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &sidex);
	fscanf(Parameters, " %f", &sidey);
	fscanf(Parameters, " %f", &sidez);


	// Read in configuration at the last frame //
	float *rcmx, *rcmy, *rcmz, *rx, *ry, *rz;
	float *q0, *q1, *q2, *q3, *px, *py, *pz;
	float *rcmxE, *rcmyE, *rcmzE;
	float *rxE, *ryE, *rzE, *pxE, *pyE, *pzE; 
	int m, i, n, j, mi, pos, idum1, idum2;
	float minVal; 

	rcmx = (float*)malloc(nfib*sizeof(float));
	rcmy = (float*)malloc(nfib*sizeof(float));
	rcmz = (float*)malloc(nfib*sizeof(float));
	rcmxE = (float*)malloc(27*nfib*sizeof(float));
	rcmyE = (float*)malloc(27*nfib*sizeof(float));
	rcmzE = (float*)malloc(27*nfib*sizeof(float));
	rx = (float*)malloc(nfib*nseg*sizeof(float));
	ry = (float*)malloc(nfib*nseg*sizeof(float));
	rz = (float*)malloc(nfib*nseg*sizeof(float));
	px = (float*)malloc(nfib*nseg*sizeof(float));
	py = (float*)malloc(nfib*nseg*sizeof(float));
	pz = (float*)malloc(nfib*nseg*sizeof(float));
	q0 = (float*)malloc(nfib*nseg*sizeof(float));
	q1 = (float*)malloc(nfib*nseg*sizeof(float));
	q2 = (float*)malloc(nfib*nseg*sizeof(float));
	q3 = (float*)malloc(nfib*nseg*sizeof(float));
	rxE = (float*)malloc(27*nfib*nseg*sizeof(float));
	ryE = (float*)malloc(27*nfib*nseg*sizeof(float));
	rzE = (float*)malloc(27*nfib*nseg*sizeof(float));
	pxE = (float*)malloc(27*nfib*nseg*sizeof(float));
	pyE = (float*)malloc(27*nfib*nseg*sizeof(float));
	pzE = (float*)malloc(27*nfib*nseg*sizeof(float));
	
	// Read in Centers of Mass, Euler Parameters
	for (m = 0; m < nfib; m++){
		fscanf(Centers_of_Mass, "%d %f %f %f ",
			&idum1, rcmx + m, rcmy + m, rcmz + m);
		for (i = 0; i < nseg; i++){
			mi = m*nseg + i;
			fscanf(Euler_Parameters, "%d %d %f %f %f %f",
				&idum1, &idum2, q0 + mi, q1 + mi, q2 + mi, q3 + mi);
			dum = sqrtf(q0[mi]*q0[mi]+q1[mi]*q1[mi]+q2[mi]*q2[mi]+q3[mi]*q3[mi]);
			q0[mi] /= dum; 
			q1[mi] /= dum; 
			q2[mi] /= dum; 
			q3[mi] /= dum; 
			px[mi] = 2.0*(q1[mi]*q3[mi] + q0[mi]*q2[mi]);
			py[mi] = 2.0*(q3[mi]*q2[mi] - q0[mi]*q1[mi]);
			pz[mi] = 2.0*(q0[mi]*q0[mi] + q3[mi]*q3[mi]) - 1.0;
			dum = sqrtf(px[mi]*px[mi] + py[mi]*py[mi] + pz[mi]*pz[mi]);
			px[mi] /= dum; 
			py[mi] /= dum; 
			pz[mi] /= dum; 
		}
	}

	// shift centers of mass to positive numbers
	shift(rcmx, rcmy, rcmz, sidex, sidey, sidez, nfib, 1);
	// regrow fiber segments from centers of mass
	regrow(rx,ry,rz,px,py,pz,rcmx,rcmy,rcmz,rps,nfib,nseg);
	
	// find the minimum position in each direction
	minVal = find_min(rx, nfib, nseg); 
	// shift by minimum position so that 
	shift_origin(rx, -1.0*minVal, nfib, nseg); 
	
	minVal = find_min(ry, nfib, nseg);
	shift_origin(ry, -1.0*minVal, nfib, nseg);

	minVal = find_min(rz, nfib, nseg);
	shift_origin(rz, -1.0*minVal, nfib, nseg);

	// extend segment position and orientation
	pos = 0; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,0.0,0.0,0.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,1.0,0.0,0.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,2.0,0.0,0.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,0.0,1.0,0.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,0.0,2.0,0.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,0.0,0.0,1.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,0.0,0.0,2.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,1.0,1.0,0.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,2.0,1.0,0.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,1.0,2.0,0.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,2.0,2.0,0.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,1.0,0.0,1.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,2.0,0.0,1.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,1.0,0.0,2.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,2.0,0.0,2.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,0.0,1.0,1.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,0.0,2.0,1.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,0.0,1.0,2.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,0.0,2.0,2.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,1.0,1.0,1.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,1.0,1.0,2.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,1.0,2.0,1.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,1.0,2.0,2.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,2.0,1.0,1.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,2.0,1.0,2.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,2.0,2.0,1.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));
	pos += nfib*nseg; 
	extend(rx,ry,rz,rxE,ryE,rzE,sidex,sidey,sidez,2.0,2.0,2.0,pos,nfib,nseg); 
	memcpy(pxE + pos, px, nfib*nseg*sizeof(float));
	memcpy(pyE + pos, py, nfib*nseg*sizeof(float));
	memcpy(pzE + pos, pz, nfib*nseg*sizeof(float));

	pos = 0; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,0.0,0.0,0.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,1.0,0.0,0.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,2.0,0.0,0.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,0.0,1.0,0.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,0.0,2.0,0.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,0.0,0.0,1.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,0.0,0.0,2.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,1.0,1.0,0.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,2.0,1.0,0.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,1.0,2.0,0.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,2.0,2.0,0.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,1.0,0.0,1.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,2.0,0.0,1.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,1.0,0.0,2.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,2.0,0.0,2.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,0.0,1.0,1.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,0.0,2.0,1.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,0.0,1.0,2.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,0.0,2.0,2.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,1.0,1.0,1.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,1.0,1.0,2.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,1.0,2.0,1.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,1.0,2.0,2.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,2.0,1.0,1.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,2.0,1.0,2.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,2.0,2.0,1.0,pos,nfib,1); 
	pos += nfib*1; 
	extend(rcmx,rcmy,rcmz,rcmxE,rcmyE,rcmzE,sidex,sidey,sidez,2.0,2.0,2.0,pos,nfib,1); 

	// write to file
	float tprint = 0.0;
	fwrite(&tprint, sizeof(float), 1, center_mass);
	fwrite(rcmxE, sizeof(float), 27*nfib, center_mass);
	fwrite(rcmyE, sizeof(float), 27*nfib, center_mass);
	fwrite(rcmzE, sizeof(float), 27*nfib, center_mass);

	fwrite(&tprint, sizeof(float), 1, pxfile);
	fwrite(&tprint, sizeof(float), 1, pyfile);
	fwrite(&tprint, sizeof(float), 1, pzfile);
	fwrite(&tprint, sizeof(float), 1, rxfile);
	fwrite(&tprint, sizeof(float), 1, ryfile);
	fwrite(&tprint, sizeof(float), 1, rzfile);

	fwrite(pxE, sizeof(float), 27*nfib*nseg, pxfile);
	fwrite(pyE, sizeof(float), 27*nfib*nseg, pyfile);
	fwrite(pzE, sizeof(float), 27*nfib*nseg, pzfile);
	fwrite(rxE, sizeof(float), 27*nfib*nseg, rxfile);
	fwrite(ryE, sizeof(float), 27*nfib*nseg, ryfile);
	fwrite(rzE, sizeof(float), 27*nfib*nseg, rzfile);
	
	free(rx); free(ry); free(rz);
	free(rcmx); free(rcmy); free(rcmz);
	free(px); free(py); free(pz);
	free(q0); free(q1); free(q2); free(q3);
	free(pxE); free(pyE); free(pzE);
	free(rxE); free(ryE); free(rzE);
	free(rcmxE); free(rcmyE); free(rcmzE);

	fclose(Parameters); 
	fclose(Euler_Parameters);
	fclose(Centers_of_Mass);
	fclose(center_mass);
	fclose(rxfile);
	fclose(rzfile);
	fclose(ryfile);
	fclose(pxfile);
	fclose(pzfile);
	fclose(pyfile);

	return 0;
}


void shift_origin(float *rx, float offset, int nfib, int nseg){

	int m, i;

	for (m = 0; m < nfib; m++){
		for (i = 0; i < nseg; i++){
			rx[m*nseg + i] += offset; 
		}
	}
}

float find_min(float *rx, int nfib, int nseg){

	float minVal = 0.0; 
	int m, i; 

	for (m = 0; m < nfib; m++){
		for (i = 0; i < nseg; i++){
			if (rx[m*nseg + i] < minVal){
				minVal = rx[m*nseg + i]; 
			}
		}
	}
	return minVal; 
}

void extend(float *rx, float *ry, float *rz, float *rxE, float *ryE, float *rzE,
            float sidex, float sidey, float sidez, float cellx, float celly,float cellz,
            int pos, int nfib, int nseg){

	int m, i, mi, ind; 
	
	for (m = 0; m < nfib; m++){
		for (i = 0; i < nseg; i++){
			mi = m*nseg+i; 
			ind = pos + mi;
			rxE[ind] = rx[mi] + cellx*sidex;
			ryE[ind] = ry[mi] + celly*sidey;
			rzE[ind] = rz[mi] + cellz*sidez;
		}
	}
}

void shift(float *rx, float *ry, float *rz, float sidex, float sidey, float sidez,
           int nfib, int nseg){

	int m, i, mi; 
	
	for (m = 0; m < nfib; m++){
		for (i = 0; i < nseg; i++){
			mi = m*nseg+i; 
			// shift segments outside of box into box
			if (rx[mi] < -sidex/2.0 || rx[mi] > sidex/2.0)
				rx[mi] -= copysign(sidex,rx[mi]);
			if (ry[mi] < -sidey/2.0 || ry[mi] > sidey/2.0)
				ry[mi] -= copysign(sidey,ry[mi]);
			if (rz[mi] < -sidez/2.0 || rz[mi] > sidez/2.0)
				rz[mi] -= copysign(sidez,rz[mi]);
			// shift segments so that all coord. are positive
			rx[mi] += sidex/2.0;
			ry[mi] += sidey/2.0;
			rz[mi] += sidez/2.0;
		}
	}
}
void regrow(float *rx, float *ry, float *rz, float *px, float *py, float *pz,
	    float *rcmx, float *rcmy, float *rcmz, float rps, int nfib, int nseg){

	int m, i, j, k, mi, mth;
	float pjx, pjy, pjz; 
	float pkx, pky, pkz; 
	// regrow fiber segments
	for (m = 0; m < nfib; m++){
		for (i = 0; i < nseg; i++){
			mi = m*nseg + i; 
			pjx = 0.0; pjy = 0.0; pjz = 0.0;
			if (i == 0){
				for (j = 1; j < nseg; j++){
					pkx = 0.0; pky = 0.0; pkz = 0.0;
					for (k = 1; k < j; k++){
						pkx += px[m*nseg + k];
						pky += py[m*nseg + k];
						pkz += pz[m*nseg + k];
					}
					pjx += px[m*nseg] + px[m*nseg + j] + 2 * pkx;
					pjy += py[m*nseg] + py[m*nseg + j] + 2 * pky;
					pjz += pz[m*nseg] + pz[m*nseg + j] + 2 * pkz;
				}
				rx[mi] = rcmx[m] - rps / float(nseg) * pjx; 
				ry[mi] = rcmy[m] - rps / float(nseg) * pjy; 
				rz[mi] = rcmz[m] - rps / float(nseg) * pjz; 
			}
			else{
				for (mth = 1; mth < i; mth++){
					pjx += px[m*nseg + mth];
					pjy += py[m*nseg + mth];
					pjz += pz[m*nseg + mth];
				}
				rx[mi] = rx[m*nseg] + rps*px[m*nseg] + 2.0*rps*pjx + rps*px[mi];
				ry[mi] = ry[m*nseg] + rps*py[m*nseg] + 2.0*rps*pjy + rps*py[mi];
				rz[mi] = rz[m*nseg] + rps*pz[m*nseg] + 2.0*rps*pjz + rps*pz[mi];
			}
		}
	} 
}
