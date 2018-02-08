//
//  main.cpp
//  flexfric
//
//  Created by Jing-Yao Chen on 7/12/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h> /* malloc, calloc, free, exit */
#include <cmath>
#include <time.h>
#include <cuda_runtime_api.h>

// function linking
#include "setVarArray.h"
#include "initialize.h"
#include "initialize_regrow.h"
#include "initialize_regrow2.h"
#include "bnei_set.h"
#include "cell.h"
#include "link.h"
#include "delta_twist_zero.h"
#include "total_contacts_zero.h"
#include "contact.h"
#include "clustering.h"
#include "clusterStat.h"
#include "clusteringSub.h"
#include "clusteringZero.h"
#include "print.h"

using namespace std;


#define gpuErrchk(ans) {gpuAssert((ans),__FILE__,__LINE__);}
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true){

	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) {
			getchar();
			exit(code);
		}
	}
}

int main(void) {

	const int npcn = 2000; 

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////
	FILE *Parameters, *Centers_of_Mass, *Euler_Parameters, *Equilibrium_Angles;
	FILE *README; 
	FILE *Cluster_results, *Cluster, *Cluster_fibID; 

	// open files //
	Parameters = fopen("Parameters.in", "r");
	Centers_of_Mass = fopen("Centers_of_Mass.in", "r");
	Euler_Parameters = fopen("Euler_Parameters.in", "r");
	Equilibrium_Angles = fopen("Equilibrium_Angles.in", "r");
	README = fopen("README.txt", "w");
	Cluster_results = fopen("Cluster_results.txt", "w");
	Cluster = fopen("Cluster.txt", "w");
	Cluster_fibID = fopen("Cluster_fibID.txt", "w");

	int nfib, nseg, config_write, maxCon, maxGr, maxBin;
	int fac, bdimx, bdimy, bdimz, fiberPerBlock, contact_write;
	float rp, kb, mu_stat, mu_kin;
	float contact_cutoff, rep_cutoff, over_cut;
	float dt, strain, sidex, sidey, sidez;
	float fstar, fact, Astar, decatt, elf, delta_rx;
	float duxdx, duydx, duzdx, duxdy, duydy;
	float duzdy, duxdz, duydz, duzdz;
	float fraction_rp, dx, dy, dz;
	int nfibGrid, nfibBlock;
	int blasGrid, blasBlock;
	int stress_write;
	float dum; 


	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &rp);
	fscanf(Parameters, "%*[^\n]%f", &kb);
	fscanf(Parameters, "%*[^\n]%f", &mu_stat);
	fscanf(Parameters, "%*[^\n]%f", &mu_kin);
	fscanf(Parameters, "%*[^\n]%f", &contact_cutoff);
	fscanf(Parameters, "%*[^\n]%f", &rep_cutoff);
	fscanf(Parameters, "%*[^\n]%f", &over_cut);
	fscanf(Parameters, "%*[^\n]%f", &dt);
	fscanf(Parameters, "%*[^\n]%f", &strain);
	fscanf(Parameters, "%*[^\n]%f", &sidex);
	fscanf(Parameters, " %f", &sidey);
	fscanf(Parameters, " %f", &sidez);
	fscanf(Parameters, "%*[^\n]%f", &fraction_rp);
	fscanf(Parameters, "%*[^\n]%d", &config_write);
	fscanf(Parameters, "%*[^\n]%d", &contact_write);
	fscanf(Parameters, "%*[^\n]%f", &fstar);
	fscanf(Parameters, "%*[^\n]%f", &fact);
	fscanf(Parameters, "%*[^\n]%f", &Astar);
	fscanf(Parameters, "%*[^\n]%f", &decatt);
	fscanf(Parameters, "%*[^\n]%f", &delta_rx);
	fscanf(Parameters, "%*[^\n]%f", &duxdx);
	fscanf(Parameters, "%*[^\n]%f", &duydx);
	fscanf(Parameters, "%*[^\n]%f", &duzdx);
	fscanf(Parameters, "%*[^\n]%f", &duxdy);
	fscanf(Parameters, "%*[^\n]%f", &duydy);
	fscanf(Parameters, "%*[^\n]%f", &duzdy);
	fscanf(Parameters, "%*[^\n]%f", &duxdz);
	fscanf(Parameters, "%*[^\n]%f", &duydz);
	fscanf(Parameters, "%*[^\n]%f", &duzdz);
	fscanf(Parameters, "%*[^\n]%d", &fac);
	fscanf(Parameters, "%*[^\n]%f", &elf);
	fscanf(Parameters, "%*[^\n]%f", &dx);
	fscanf(Parameters, "%*[^\n]%f", &dy);
	fscanf(Parameters, "%*[^\n]%f", &dz);
	fscanf(Parameters, "%*[^\n]%d", &bdimx);
	fscanf(Parameters, "%*[^\n]%d", &bdimy);
	fscanf(Parameters, "%*[^\n]%d", &bdimz);
	fscanf(Parameters, "%*[^\n]%d", &fiberPerBlock);
	fscanf(Parameters, "%*[^\n]%d", &maxCon);
	fscanf(Parameters, "%*[^\n]%d", &maxGr);
	fscanf(Parameters, "%*[^\n]%d", &maxBin);
	fscanf(Parameters, "%*[^\n]%d", &nfibGrid);
	fscanf(Parameters, "%*[^\n]%d", &nfibBlock);
	fscanf(Parameters, "%*[^\n]%d", &blasGrid);
	fscanf(Parameters, "%*[^\n]%d", &blasBlock);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &stress_write);

	////////////////////////////////////////
	//             allocation             //
	////////////////////////////////////////

	int *d_nfib, *d_nseg, *d_maxCon, *d_maxGr, *d_fac;
	int *d_bdimx, *d_bdimy, *d_bdimz, *d_maxBin;
	float *d_rp, *d_kb, *d_mu_stat, *d_mu_kin;
	float *d_contact_cutoff, *d_rep_cutoff, *d_over_cut;
	float *d_dt, *d_sidex, *d_sidey, *d_sidez;
	float *d_fstar, *d_fact, *d_Astar, *d_decatt, *d_elf;
	float *d_delta_rx, *d_duxdx, *d_duxdy, *d_duxdz, *d_duydx;
	float *d_duydy, *d_duydz, *d_duzdx, *d_duzdy, *d_duzdz;
	float *d_dx, *d_dy, *d_dz;

	cudaMalloc((void**)&d_nfib, sizeof(int));
	cudaMalloc((void**)&d_nseg, sizeof(int));
	cudaMalloc((void**)&d_maxCon, sizeof(int));
	cudaMalloc((void**)&d_maxGr, sizeof(int));
	cudaMalloc((void**)&d_maxBin, sizeof(int));
	cudaMalloc((void**)&d_fac, sizeof(int));
	cudaMalloc((void**)&d_bdimx, sizeof(int));
	cudaMalloc((void**)&d_bdimy, sizeof(int));
	cudaMalloc((void**)&d_bdimz, sizeof(int));
	cudaMalloc((void**)&d_rp, sizeof(float));
	cudaMalloc((void**)&d_kb, sizeof(float));
	cudaMalloc((void**)&d_mu_stat, sizeof(float));
	cudaMalloc((void**)&d_mu_kin, sizeof(float));
	cudaMalloc((void**)&d_contact_cutoff, sizeof(float));
	cudaMalloc((void**)&d_rep_cutoff, sizeof(float));
	cudaMalloc((void**)&d_over_cut, sizeof(float));
	cudaMalloc((void**)&d_dt, sizeof(float));
	cudaMalloc((void**)&d_sidex, sizeof(float));
	cudaMalloc((void**)&d_sidey, sizeof(float));
	cudaMalloc((void**)&d_sidez, sizeof(float));
	cudaMalloc((void**)&d_fstar, sizeof(float));
	cudaMalloc((void**)&d_fact, sizeof(float));
	cudaMalloc((void**)&d_Astar, sizeof(float));
	cudaMalloc((void**)&d_decatt, sizeof(float));
	cudaMalloc((void**)&d_elf, sizeof(float));
	cudaMalloc((void**)&d_delta_rx, sizeof(float));
	cudaMalloc((void**)&d_duxdx, sizeof(float));
	cudaMalloc((void**)&d_duxdy, sizeof(float));
	cudaMalloc((void**)&d_duxdz, sizeof(float));
	cudaMalloc((void**)&d_duydx, sizeof(float));
	cudaMalloc((void**)&d_duydy, sizeof(float));
	cudaMalloc((void**)&d_duydz, sizeof(float));
	cudaMalloc((void**)&d_duzdx, sizeof(float));
	cudaMalloc((void**)&d_duzdy, sizeof(float));
	cudaMalloc((void**)&d_duzdz, sizeof(float));
	cudaMalloc((void**)&d_dx, sizeof(float));
	cudaMalloc((void**)&d_dy, sizeof(float));
	cudaMalloc((void**)&d_dz, sizeof(float));

	// variables on host
	float pi = 3.14159265;
	//float *rcmx = (float*)malloc(nfib* sizeof(float));
	float *rcmx, *rcmy, *rcmz, *rx, *ry, *rz;
	float *q0, *q1, *q2, *q3, *px, *py, *pz;
	float *ux, *uy, *uz, *wx, *wy, *wz;
	int *num_groups, *total_contacts;

	cudaMallocHost((void**)&rcmx, nfib*sizeof(float));
	cudaMallocHost((void**)&rcmy, nfib*sizeof(float));
	cudaMallocHost((void**)&rcmz, nfib*sizeof(float));
	cudaMallocHost((void**)&rx, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&ry, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&rz, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&q0, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&q1, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&q2, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&q3, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&px, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&py, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&pz, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&ux, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&uy, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&uz, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&wx, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&wy, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&wz, nfib*nseg*sizeof(float));
	cudaMallocHost((void**)&num_groups, sizeof(int));
	cudaMallocHost((void**)&total_contacts, sizeof(int));

	float C0, C1, C2, C3, C4;
	float Omega_x, Omega_y, Omega_z, kt;
	float E11, E12, E13, E22, E23, E33;
	float Xa, Ya, Xc, Yc, Yh;
	float *thetaeq = (float*)malloc(nfib*nseg*sizeof(float));
	float *phieq = (float*)malloc(nfib*nseg*sizeof(float));
	float re, ecc;
	float  nL3, volfrac, consistency;
	// re - effective segment aspect ratio
	// ecc - eccentricity
	// time_steps - number of time steps
	// num_groups - number of groups
	// overs - number of overlapping contacts	
	// nL3 - dimensionless concentration ("n-L-cubed")
	// volfrac - volume fraction
	// consist - consistency = volfrac/2.6

	// floating point variables on device
	float *d_rcmx, *d_rcmy, *d_rcmz, *d_rx, *d_ry, *d_rz;
	float *d_q0, *d_q1, *d_q2, *d_q3;
	float *d_q0dot, *d_q1dot, *d_q2dot, *d_q3dot;
	float *d_qe0, *d_qe1, *d_qe2, *d_qe3;
	float *d_px, *d_py, *d_pz;
	float *d_ucmx, *d_ucmy, *d_ucmz;
	float *d_ucox, *d_ucoy, *d_ucoz;
	float *d_ux, *d_uy, *d_uz;
	float *d_uxfl, *d_uyfl, *d_uzfl;
	float *d_wx, *d_wy, *d_wz;
	// rcmx,rcmy,rcmz - fiber center of mass
	// rx,ry,rz - segment centers
	// q0,q1,q2,q3 - segment euler parameters
	// q0dot... - time derivatives of Euler Parameters
	// qe1,... - time derivatives of Euler Parameters from previous time
	// px,py,pz - segment orientational vectors
	// ucmx,ucmy,ucmz - fiber center of mass velocity
	// ucox,ucoy,ucoz -  fiber center of mass velocity from previous time
	// ux,uy,uz - segment velocity
	// uxfl,uyfl,uzfl - ambient flow field velocity
	// wx,wy,wz - segment angular velocity	 

	cudaMalloc((void**)&d_rcmx, nfib*sizeof(float));
	cudaMalloc((void**)&d_rcmy, nfib*sizeof(float));
	cudaMalloc((void**)&d_rcmz, nfib*sizeof(float));
	cudaMalloc((void**)&d_rx, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_ry, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_rz, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_q0, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_q1, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_q2, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_q3, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_q0dot, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_q1dot, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_q2dot, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_q3dot, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_qe0, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_qe1, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_qe2, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_qe3, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_px, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_py, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_pz, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_ucmx, nfib*sizeof(float));
	cudaMalloc((void**)&d_ucmy, nfib*sizeof(float));
	cudaMalloc((void**)&d_ucmz, nfib*sizeof(float));
	cudaMalloc((void**)&d_ucox, nfib*sizeof(float));
	cudaMalloc((void**)&d_ucoy, nfib*sizeof(float));
	cudaMalloc((void**)&d_ucoz, nfib*sizeof(float));
	cudaMalloc((void**)&d_ux, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_uy, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_uz, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_uxfl, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_uyfl, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_uzfl, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_wx, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_wy, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_wz, nfib*nseg*sizeof(float));

	float *d_R11, *d_R12, *d_R13, *d_R21, *d_R22, *d_R23;
	float *d_R11eq, *d_R12eq, *d_R13eq, *d_R21eq, *d_R22eq, *d_R23eq;
	float *d_R31eq, *d_R32eq, *d_R33eq;
	float *d_fx, *d_fy, *d_fz, *d_tx, *d_ty, *d_tz;
	float *d_fcx, *d_fcy, *d_fcz, *d_tcx, *d_tcy, *d_tcz;
	float *d_fbx, *d_fby, *d_fbz, *d_tbx, *d_tby, *d_tbz;
	// R11,R12,... - segment rotation matrix
	// R11eq,... - equilibrium rotation matrix
	// Xx,Xy,Xz - intrafiber constraint forces
	// Yx,Yy,Yz - restoring bending/twisting torques
	// fx,fy,fz - frictional forces on segment
	// tx,ty,tz - frictional torques on segment
	// fbx,fby,fbz - segment body force (i.e. gravity,etc.)
	// tbx,tby,tbz - segment body torque
	// fcx,fcy,fcz - segment colloidal repulsive force
	// tcx,tcy,tcz - segment colloidal repulsive torque

	cudaMalloc((void**)&d_R11, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R12, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R13, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R21, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R22, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R23, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R11eq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R12eq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R13eq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R21eq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R22eq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R23eq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R31eq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R32eq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_R33eq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_fx, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_fy, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_fz, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_tx, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_ty, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_tz, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_fcx, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_fcy, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_fcz, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_tcx, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_tcy, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_tcz, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_fbx, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_fby, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_fbz, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_tbx, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_tby, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_tbz, nfib*nseg*sizeof(float));

	float *d_D1, *d_D2, *d_D3;
	float *d_A11, *d_A12, *d_A13, *d_A23, *d_A22, *d_A33;
	float *d_C11, *d_C12, *d_C13, *d_C23, *d_C22, *d_C33;
	float *d_C0, *d_C1, *d_C2, *d_C3, *d_C4;
	float *d_g, *d_Omega_x, *d_Omega_y, *d_Omega_z;
	// D1,D2,D3 - "Jeffery" term in hydrodynamic torque
	// A11,... - A inverse term in notes
	// C11,... - C inverse term in notes
	// C0,C1... - constant groupings
	// g - gap between fiber centers at contact point
	// Omega_... - ambient flow field angular velocity

	cudaMalloc((void**)&d_D1, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_D2, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_D3, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_A11, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_A12, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_A13, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_A23, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_A22, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_A33, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_C11, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_C12, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_C13, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_C23, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_C22, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_C33, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_C0, sizeof(float));
	cudaMalloc((void**)&d_C1, sizeof(float));
	cudaMalloc((void**)&d_C2, sizeof(float));
	cudaMalloc((void**)&d_C3, sizeof(float));
	cudaMalloc((void**)&d_C4, sizeof(float));
	cudaMalloc((void**)&d_g, 1 * sizeof(float));
	cudaMalloc((void**)&d_Omega_x, sizeof(float));
	cudaMalloc((void**)&d_Omega_y, sizeof(float));
	cudaMalloc((void**)&d_Omega_z, sizeof(float));

	float *d_kt, *d_E11, *d_E12, *d_E13, *d_E22, *d_E23, *d_E33;
	float *d_neighb_cutoff, *d_Ya, *d_Yc, *d_Yh;
	float *d_thetaeq, *d_phieq;
	// kt - twisting constant	
	// E11,E12,... - rate of strain tensor
	// Xa,Ya,Xc,Yc,Yh - scalar resistance functions	
	// A_fric, P_fric - matrix and RHS for friction contact
	// thetaeq/phieq - angles associated with the equil. configuration

	cudaMalloc((void**)&d_kt, sizeof(float));
	cudaMalloc((void**)&d_E11, sizeof(float));
	cudaMalloc((void**)&d_E12, sizeof(float));
	cudaMalloc((void**)&d_E13, sizeof(float));
	cudaMalloc((void**)&d_E22, sizeof(float));
	cudaMalloc((void**)&d_E23, sizeof(float));
	cudaMalloc((void**)&d_E33, sizeof(float));
	cudaMalloc((void**)&d_neighb_cutoff, sizeof(float));
	cudaMalloc((void**)&d_Ya, sizeof(float));
	cudaMalloc((void**)&d_Yc, sizeof(float));
	cudaMalloc((void**)&d_Yh, sizeof(float));
	cudaMalloc((void**)&d_thetaeq, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_phieq, nfib*nseg*sizeof(float));

	float *nxV, *nyV, *nzV, *gV;
	float *GijxV, *GijyV, *GijzV;
	float *GjixV, *GjiyV, *GjizV;

	cudaMalloc((void**)&nxV, nfib*nseg*maxCon*sizeof(float));
	cudaMalloc((void**)&nyV, nfib*nseg*maxCon*sizeof(float));
	cudaMalloc((void**)&nzV, nfib*nseg*maxCon*sizeof(float));
	cudaMalloc((void**)&gV, nfib*nseg*maxCon*sizeof(float));
	cudaMalloc((void**)&GijxV, nfib*nseg*maxCon*sizeof(float));
	cudaMalloc((void**)&GijyV, nfib*nseg*maxCon*sizeof(float));
	cudaMalloc((void**)&GijzV, nfib*nseg*maxCon*sizeof(float));
	cudaMalloc((void**)&GjixV, nfib*nseg*maxCon*sizeof(float));
	cudaMalloc((void**)&GjiyV, nfib*nseg*maxCon*sizeof(float));
	cudaMalloc((void**)&GjizV, nfib*nseg*maxCon*sizeof(float));

	int *cluster, *clusterAccess, *d_clusterCount, *clusterCount, *singleCluster; 
	cudaMalloc((void**)&cluster, nfib*nfib*sizeof(float));
	cudaMalloc((void**)&clusterAccess, nfib*sizeof(float));
	cudaMalloc((void**)&d_clusterCount, nfib*sizeof(float));
	clusterCount = (int*)malloc(nfib*sizeof(int));
	singleCluster = (int*)malloc(nfib*sizeof(int));

	// calculate bin parameters
	int nxbin, nybin, nzbin;
	nxbin = int(floorf(sidex / dx));
	nybin = int(floorf(sidey / dy));
	nzbin = int(floorf(sidez / dz));
	printf("dx dy dz %10f %10f %10f\n", dx, dy, dz);
	printf("nxbin nybin nzbin %4d %4d %4d\n", nxbin, nybin, nzbin);
	if (nxbin % 2 != 0){
		nxbin--;
	}
	if (nybin % 2 != 0){
		nybin--;
	}
	if (nzbin % 2 != 0){
		nzbin--;
	}
	dx = sidex / float(nxbin);
	dy = sidey / float(nybin);
	dz = sidez / float(nzbin);
	printf("dx dy dz %10f %10f %10f\n", dx, dy, dz);
	printf("nxbin nybin nzbin %4d %4d %4d\n", nxbin, nybin, nzbin);

	// integer variables on device
	int *d_ifiber, *d_ncnt, *d_ncpf, *d_num_groups, *d_overs;
	int *d_total_contacts;
	int *bin, *list, *bnei, *d_nxbin, *d_nybin, *d_nzbin;
	int *clist, *status, *lead_clist, *nc, *clist_pos;
	int *bnum, *potCon, *potConSize, *groupId;
	// ifiber - index of fibers in contact
	// ncnt - number of contacts in groups
	// ncpf - number of fibers contacting fiber 
	// num_groups - number of groups in contact
	// overs - number of overlapping contacts
	// total_contacts - number of contacts in system
	// step - time step
	// csrRow - indices to first nonzero value in ith block row
	// csrCol - col indices of corresponding value in bsrVal
	// bin - number of fibers in bin[x]
	// list - list of fibers in each bin
	// bnei - neighbor of each bin
	// nxbin - number of bins in x direction
	// clist - list of fibers in contact
	// status - marks whether a fiber is in charge of solving
	//          the system of contacts
	// lead_clist - list of all fibers in contacting group
	// nc - number of contacts in groups 
	// clist_pos - position of contact in clist

	cudaMalloc((void**)&d_ifiber, nfib*nseg * 2 * maxGr*sizeof(int));
	cudaMalloc((void**)&d_ncnt, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&d_ncpf, (nfib*nseg)*sizeof(int));
	cudaMalloc((void**)&d_num_groups, sizeof(int));
	cudaMalloc((void**)&d_overs, sizeof(int));
	cudaMalloc((void**)&d_total_contacts, sizeof(int));
	cudaMalloc((void**)&bin, nxbin*nybin*nzbin*sizeof(int));
	cudaMalloc((void**)&list, maxBin * nxbin*nybin*nzbin*sizeof(int));
	cudaMalloc((void**)&bnei, bdimx*bdimy*bdimz*nxbin*nybin*nzbin*sizeof(int));
	cudaMalloc((void**)&d_nxbin, sizeof(int));
	cudaMalloc((void**)&d_nybin, sizeof(int));
	cudaMalloc((void**)&d_nzbin, sizeof(int));
	cudaMalloc((void**)&clist, nfib*nseg*maxCon*sizeof(int));
	cudaMalloc((void**)&status, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&lead_clist, nfib*nseg*maxGr*sizeof(int));
	cudaMalloc((void**)&nc, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&clist_pos, nfib*nseg*maxGr*sizeof(int));
	cudaMalloc((void**)&bnum, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&potCon, nfib*nseg * npcn * sizeof(int));
	cudaMalloc((void**)&potConSize, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&groupId, nfib*nseg*sizeof(int));

	if (cudaSuccess != cudaGetLastError()){ printf("after int variables !\n"); getchar(); }

	// arrays of device pointers 
	const int numVar = 166;
	const int numIntVar = 38;
	// numVar - number of float variables
	// numIntVar - number of integer variables
	float **var;
	int **intVar;
	cudaMalloc((void**)&var, numVar*sizeof(float*));
	cudaMalloc((void**)&intVar, numIntVar*sizeof(int*));

	// kernel launch parameters
	if (nfib < fiberPerBlock) fiberPerBlock = nfib;
	int numThreads = fiberPerBlock * nseg;
	int numBlocks = (nfib * nseg) / numThreads;
	// numThreads - number of threads per block
	// numBlocks  - number of blocks 

	// temporary variables
	int idum1, idum2, m, i, mi, mm;

	int nCluster; 

	int *clusterString = (int*)malloc(nfib*sizeof(int)); 

	////////////////////////////////////////
	//           more input files         //
	////////////////////////////////////////

	// Read in Centers of Mass, Euler Parameters, Equilibrium Angles //
	for (m = 0; m < nfib; m++){

		fscanf(Centers_of_Mass, "%d %f %f %f ",
			&idum1, rcmx + m, rcmy + m, rcmz + m);
		for (i = 0; i < nseg; i++){
			mi = m*nseg + i;
			fscanf(Euler_Parameters, "%d %d %f %f %f %f",
				&idum1, &idum2, q0 + mi, q1 + mi, q2 + mi, q3 + mi);
			if (i > 0){
				fscanf(Equilibrium_Angles, "%d %d %f %f", &idum1, &idum2, thetaeq + mi, phieq + mi);
			}
		}
	}
	fclose(Parameters);
	fclose(Centers_of_Mass);
	fclose(Euler_Parameters);
	fclose(Equilibrium_Angles);

	if (cudaSuccess != cudaGetLastError()){ printf("after reading input files!\n"); getchar(); }


	////////////////////////////////////////
	//       constant initialization      //
	////////////////////////////////////////

	elf = 0.0;
	// Rate of strain tensor for the flow field
	E11 = 0.5*(duxdx + duxdx);
	E12 = 0.5*(duydx + duxdy);
	E13 = 0.5*(duzdx + duxdz);
	E22 = 0.5*(duydy + duydy);
	E23 = 0.5*(duzdy + duydz);
	E33 = 0.5*(duzdz + duzdz);
	// Fluid vorticity vector
	Omega_x = 0.5*(duzdy - duydz);
	Omega_y = 0.5*(duxdz - duzdx);
	Omega_z = 0.5*(duydx - duxdy);
	// Resistance functions
	re = fraction_rp*rp;
	ecc = sqrtf(re*re - 1.0) / re;
	Xa = (8.0*powf(ecc, 3.0) / 3.0) / (-2.0*ecc + (1.0 + ecc*ecc)*logf((1.0 + ecc) / (1.0 - ecc)));
	Ya = (16.0*powf(ecc, 3.0) / 3.0) / (2.0*ecc + (3.0*ecc*ecc - 1.0)*logf((1.0 + ecc) / (1.0 - ecc)));
	Xc = (4.0*powf(ecc, 3.0) / 3.0)*(1.0 - ecc*ecc) / (2.0*ecc - (1.0 - ecc*ecc)*logf((1.0 + ecc) / (1.0 - ecc)));
	Yc = (4.0*powf(ecc, 3.0) / 3.0)*(2.0 - ecc*ecc) / (-2.0*ecc + (1.0 + ecc*ecc)*logf((1.0 + ecc) / (1.0 - ecc)));
	Yh = (4.0*powf(ecc, 5.0) / 3.0) / (-2.0*ecc + (1.0 + ecc*ecc)*logf((1.0 + ecc) / (1.0 - ecc)));
	// Other constants //
	kt = 0.67*kb; // twisting constant
	C0 = 3.0 / (4.0*rp);
	C1 = 3.0 / (4.0*rp*rp);
	C2 = 3.0 / (4.0*Yc);
	C3 = 1.0 / Xa - 1.0 / Ya;
	C4 = 1.0 / Xc - 1.0 / Yc;
	contact_cutoff = powf((contact_cutoff + 2.0), 2.0);
	rep_cutoff = powf((rep_cutoff + 2.0), 2.0);
	nL3 = float(nfib)*powf((float(2.0)*float(nseg)*rp), 3.0) / (sidex*sidey*sidez);
	volfrac = float(nfib)*float(nseg)*pi*(float(2.0)*rp) / (sidex*sidey*sidez);
	consistency = volfrac / float(2.6000);

	printf("nL3 %10f\nvolfrac %10f\nconsistency %10f\n", nL3, volfrac, consistency);

	////////////////////////////////////////
	//             README output          //
	////////////////////////////////////////

	// Make the README file to expflain which run this is
	fprintf(README, "Fibers with kinetic friction\n\n");
	fprintf(README, "Number of Fibers: %d\n", nfib);
	fprintf(README, "Number of Segments: %d\n", nseg);
	fprintf(README, "Aspect Ratio of a Segment: %.15f\n", rp);
	fprintf(README, "Aspect Ratio of a Fiber: %.15f\n", rp*nseg);
	fprintf(README, "Time Step: %.15E\n", dt);
	fprintf(README, "Total Strain: %.15E\n", strain);
	fprintf(README, "Box Side X: %.15f\n", sidex);
	fprintf(README, "Box Side Y: %.15f\n", sidey);
	fprintf(README, "Box Side Z: %.15f\n", sidez);
	fprintf(README, "Coefficients of Friction: %.15f %.15f\n", mu_stat, mu_kin);
	fprintf(README, "Bending Constant: %.15f\n", kb);
	fprintf(README, "Twisting Constant: %.15f\n", kt);
	fprintf(README, "concentration, nL3: %.15f\n", nL3);
	fprintf(README, "Volume fraction: %.15E\n", volfrac);
	fprintf(README, "Consistency: %.15E\n", consistency);
	fclose(README);

	////////////////////////////////////////
	//            Memory to Device        //
	////////////////////////////////////////

	cudaMemcpy(d_dx, &dx, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dy, &dy, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dz, &dz, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_nxbin, &nxbin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_nybin, &nybin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_nzbin, &nzbin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_bdimx, &bdimx, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_bdimy, &bdimy, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_bdimz, &bdimz, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_maxCon, &maxCon, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_maxGr, &maxGr, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_maxBin, &maxBin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_nfib, &nfib, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_nseg, &nseg, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_rp, &rp, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_kb, &kb, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_mu_stat, &mu_stat, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_mu_kin, &mu_kin, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_over_cut, &over_cut, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dt, &dt, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sidex, &sidex, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sidey, &sidey, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sidez, &sidez, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_fstar, &fstar, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_fact, &fact, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Astar, &Astar, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_decatt, &decatt, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_delta_rx, &delta_rx, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_duxdx, &duxdx, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_duxdy, &duxdy, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_duxdz, &duxdz, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_duydx, &duydx, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_duydy, &duydy, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_duydz, &duydz, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_duzdx, &duzdx, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_duzdy, &duzdy, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_duzdz, &duzdz, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_fac, &fac, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_elf, &elf, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_E11, &E11, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_E12, &E12, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_E13, &E13, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_E22, &E22, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_E23, &E23, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_E33, &E33, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Omega_x, &Omega_x, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Omega_y, &Omega_y, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Omega_z, &Omega_z, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Ya, &Ya, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Yc, &Yc, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Yh, &Yh, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_kt, &kt, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_C0, &C0, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_C1, &C1, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_C2, &C2, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_C3, &C3, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_C4, &C4, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_contact_cutoff, &contact_cutoff, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_rep_cutoff, &rep_cutoff, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_rcmx, rcmx, nfib*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_rcmy, rcmy, nfib*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_rcmz, rcmz, nfib*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_q0, q0, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_q1, q1, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_q2, q2, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_q3, q3, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_thetaeq, thetaeq, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_phieq, phieq, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);

	////////////////////////////////////////
	//           Initialization           //
	////////////////////////////////////////

	// stores device pointers into double pointer arrays
	setVarArray << <1, 1 >> >(var, intVar, d_ifiber, d_ncnt, d_ncpf,
		d_num_groups, d_overs, d_total_contacts, d_nfib, d_nseg, d_fac, 
		bin, list, bnei, d_nxbin, d_nybin, d_nzbin,
		clist, status, lead_clist, nc, d_bdimx, d_bdimy, d_bdimz, d_maxCon, clist_pos, d_maxGr, d_maxBin,
		bnum, potCon, potConSize, groupId,
		cluster, clusterAccess, d_clusterCount,
		d_rcmx, d_rcmy, d_rcmz, d_rx, d_ry, d_rz, d_q0, d_q1, d_q2, d_q3, d_q0dot, d_q1dot,
		d_q2dot, d_q3dot, d_qe0, d_qe1, d_qe2, d_qe3, d_px, d_py, d_pz, d_ucmx, d_ucmy,
		d_ucmz, d_ucox, d_ucoy, d_ucoz, d_ux, d_uy, d_uz, d_uxfl, d_uyfl, d_uzfl, d_wx, d_wy,
		d_wz, d_R11, d_R12, d_R13, d_R21, d_R22, d_R23, d_R11eq, d_R12eq, d_R13eq, d_R21eq,
		d_R22eq, d_R23eq, d_R31eq, d_R32eq, d_R33eq, 
		d_fx, d_fy, d_fz, d_tx, d_ty, d_tz, d_fcx, d_fcy, d_fcz,
		d_tcx, d_tcy, d_tcz, d_fbx, d_fby, d_fbz, d_tbx, d_tby, d_tbz, 
		d_C0, d_C1, d_C2, d_C3, d_C4, d_g, d_duxdx, d_duxdy, d_duxdz, d_duydx, d_duydy,
		d_duydz, d_duzdx, d_duzdy, d_duzdz, d_dx, d_dy, d_dz,
		d_Omega_x, d_Omega_y, d_Omega_z, d_rp, d_kb,
		d_kt, d_dt, d_over_cut, d_sidex, d_sidey, d_sidez, d_E11, d_E12, d_E13, d_E22, d_E23,
		d_E33, d_contact_cutoff, d_rep_cutoff, d_neighb_cutoff, d_Ya, d_Yc, d_Yh, d_thetaeq, d_phieq, d_delta_rx, d_mu_stat, d_mu_kin,
		d_fstar, d_fact, d_Astar, d_decatt, d_elf, GijxV, GijyV, GijzV, GjixV,
		GjiyV, GjizV, nxV, nyV, nzV, gV);
	

	// initialize variables
	initialize << <numBlocks, numThreads >> >(var, intVar);	
	if (cudaSuccess != cudaGetLastError()){ printf("after initialize!\n"); getchar(); exit(0); }

	initialize_regrow << <nfibGrid, nfibBlock >> >(var, intVar);

	// if (nseg > 1)
	initialize_regrow2 << <nfib / fiberPerBlock, (nseg - 1) * fiberPerBlock >> >(var, intVar);
	if (cudaSuccess != cudaGetLastError()){ printf("after initialize regrow 2!\n"); getchar();  exit(0); }

	gpuErrchk(cudaDeviceSynchronize()); 

	// set list of neighbor for each bin
	bnei_set << < nzbin*nybin, nxbin >> > (var, intVar);

	gpuErrchk(cudaDeviceSynchronize());

	// 1. zero sorting, friction, X forces variables
	delta_twist_zero << <numBlocks, numThreads >> >(var, intVar);

	gpuErrchk(cudaDeviceSynchronize());

	
	// 5. Linked cell sort
	//// a. put fibers into bins
	clusteringZero<<<nfibGrid,nfibBlock>>>(var,intVar); 

	gpuErrchk(cudaDeviceSynchronize());

	cell << < numBlocks, numThreads >> > (var, intVar);
	//// b. link with fibers in neighbor bins
	link << < nfib*nseg, bdimx*bdimy*bdimz >> > (var, intVar);
	contact << <numBlocks, numThreads >> >(var, intVar);

	//clusteringSub << <1, 32 >> >(var, intVar);

	//clusteringSub << <1, 16 >> >(var, intVar);
	
	clustering << <1,1 >> >(var, intVar); 

	clusterStat << <nfibGrid, nfibBlock >> >(var, intVar);

	cudaDeviceSynchronize(); 

	cudaMemcpy(clusterCount, d_clusterCount, nfib*sizeof(int), cudaMemcpyDeviceToHost);


	for (m = 0; m < nfib; m++){
		clusterString[m] = 0; 
	}

	fprintf(Cluster_results, "Cluster ID, Number of fibers in cluster: \n");
	
	nCluster = 0; 
	for (m = 0; m < nfib; m++){
		if (clusterCount[m] > 1){
			nCluster++; 
			fprintf(Cluster_results, "%6d %6d\n", nCluster, clusterCount[m]);
			fprintf(Cluster, "%6d %6d ", nCluster, clusterCount[m]);
			//printf("%6d %6d \n", nCluster, clusterCount[m]);
			cudaMemcpy(singleCluster, cluster + m*nfib, nfib*sizeof(int), cudaMemcpyDeviceToHost); 
			for (mm = 0; mm < nfib; mm++){
				if (singleCluster[mm] != 0){
					clusterString[mm] = nCluster;
					fprintf(Cluster, "%8d ", mm);
					//printf("%8d ", mm);
				}				
			}
			fprintf(Cluster, "\n");
			//printf("\n");
		}
	}

	fprintf(Cluster_results, "Number of clusters: \n");
	fprintf(Cluster_results, "%6d\n", nCluster);

	for (m = 0; m < nfib; m++){
		fprintf(Cluster_fibID, "%8d %4d\n", m + 1, clusterString[m]); 
	}

	fclose(Cluster); fclose(Cluster_results); fclose(Cluster_fibID); 

	// free host memory
	//free(rcmx); 
	cudaFreeHost(rcmx); cudaFreeHost(rcmy); cudaFreeHost(rcmz);
	cudaFreeHost(rx); cudaFreeHost(ry); cudaFreeHost(rz);
	cudaFreeHost(px); cudaFreeHost(py); cudaFreeHost(pz);
	cudaFreeHost(ux); cudaFreeHost(uy); cudaFreeHost(uz);
	cudaFreeHost(wx); cudaFreeHost(wy); cudaFreeHost(wz);
	cudaFreeHost(q0); cudaFreeHost(q1); cudaFreeHost(q2); cudaFreeHost(q3);
	cudaFreeHost(num_groups); cudaFreeHost(total_contacts);
	free(phieq); free(thetaeq); 
	free(clusterString); free(clusterCount); free(singleCluster);

	// free device memory
	cudaFree(var);       cudaFree(intVar);
	cudaFree(d_phieq);   cudaFree(d_thetaeq);
	cudaFree(d_rcmx);    cudaFree(d_rcmy);    cudaFree(d_rcmz);
	cudaFree(d_rx);	     cudaFree(d_ry);      cudaFree(d_rz);
	cudaFree(d_q0);      cudaFree(d_q1);      cudaFree(d_q2);	 cudaFree(d_q3);
	cudaFree(d_q0dot);   cudaFree(d_q1dot);   cudaFree(d_q2dot); cudaFree(d_q3dot);
	cudaFree(d_qe0);     cudaFree(d_qe1);     cudaFree(d_qe2);	 cudaFree(d_qe3);
	cudaFree(d_px);      cudaFree(d_py);      cudaFree(d_pz);
	cudaFree(d_ucmx);    cudaFree(d_ucmy);    cudaFree(d_ucmz);
	cudaFree(d_ucox);    cudaFree(d_ucoy);    cudaFree(d_ucoz);
	cudaFree(d_ux);      cudaFree(d_uy);      cudaFree(d_uz);
	cudaFree(d_uxfl);    cudaFree(d_uyfl);    cudaFree(d_uzfl);
	cudaFree(d_wx);      cudaFree(d_wy);      cudaFree(d_wz);
	cudaFree(d_R11);     cudaFree(d_R12);     cudaFree(d_R13);
	cudaFree(d_R21);     cudaFree(d_R22);     cudaFree(d_R23);
	cudaFree(d_R11eq);   cudaFree(d_R12eq);   cudaFree(d_R13eq);
	cudaFree(d_R21eq);   cudaFree(d_R22eq);   cudaFree(d_R23eq);
	cudaFree(d_R31eq);   cudaFree(d_R32eq);   cudaFree(d_R33eq);
	cudaFree(d_fx);      cudaFree(d_fy);      cudaFree(d_fz);
	cudaFree(d_tx);      cudaFree(d_ty);      cudaFree(d_tz);
	cudaFree(d_fcx);     cudaFree(d_fcy);     cudaFree(d_fcz);
	cudaFree(d_tcx);     cudaFree(d_tcy);     cudaFree(d_tcz);
	cudaFree(d_fbx);     cudaFree(d_fby);     cudaFree(d_fbz);
	cudaFree(d_tbx);     cudaFree(d_tby);     cudaFree(d_tbz);
	cudaFree(d_ifiber);  cudaFree(d_ncnt);    cudaFree(d_ncpf);
	cudaFree(d_Omega_x); cudaFree(d_Omega_y); cudaFree(d_Omega_z);
	cudaFree(d_rp);      cudaFree(d_kb);      cudaFree(d_kt);
	cudaFree(d_sidex);   cudaFree(d_sidey);   cudaFree(d_sidez);
	cudaFree(d_E11);     cudaFree(d_E12);     cudaFree(d_E13);
	cudaFree(d_E22);     cudaFree(d_E23);     cudaFree(d_E33);
	cudaFree(d_Ya);      cudaFree(d_Yc);	  cudaFree(d_Yh);
	cudaFree(d_Astar);   cudaFree(d_decatt);  cudaFree(d_delta_rx);
	cudaFree(d_mu_stat); cudaFree(d_mu_kin);  cudaFree(d_contact_cutoff);
	cudaFree(d_fstar);   cudaFree(d_fact);    cudaFree(d_rep_cutoff);
	cudaFree(d_elf);	 cudaFree(d_g);       cudaFree(d_neighb_cutoff);
	cudaFree(d_overs);   cudaFree(d_dt);      cudaFree(d_num_groups);
	cudaFree(d_nfib);    cudaFree(d_nseg);    cudaFree(d_total_contacts);
	cudaFree(d_fac);     cudaFree(d_maxCon);  cudaFree(d_over_cut);
	cudaFree(d_duxdx);   cudaFree(d_duxdy);   cudaFree(d_duxdz);
	cudaFree(d_duydx);   cudaFree(d_duydy);   cudaFree(d_duydz);
	cudaFree(d_duzdx);   cudaFree(d_duzdy);   cudaFree(d_duzdz);
	cudaFree(d_C0);		 cudaFree(d_C1);	  cudaFree(d_C2);
	cudaFree(d_C3);		 cudaFree(d_C4); 
	cudaFree(nc);  
	cudaFree(bin);       cudaFree(list);      cudaFree(bnei);
	cudaFree(d_nxbin);   cudaFree(d_nybin);   cudaFree(d_nzbin);
	cudaFree(clist);     cudaFree(status);    cudaFree(lead_clist);
	cudaFree(d_dx);      cudaFree(d_dy);      cudaFree(d_dz);
	cudaFree(d_bdimx);   cudaFree(d_bdimy);   cudaFree(d_bdimz);
	cudaFree(clist_pos); cudaFree(nxV);       cudaFree(nyV);
	cudaFree(nzV);       cudaFree(gV);        cudaFree(GijxV);
	cudaFree(GijyV);     cudaFree(GijzV);     cudaFree(GjixV);
	cudaFree(GjiyV);     cudaFree(GjizV);     cudaFree(d_maxGr);
	cudaFree(d_maxBin); 
	cudaFree(groupId);  
	cudaFree(bnum);      cudaFree(potCon);	  cudaFree(potConSize);
	cudaFree(cluster); cudaFree(clusterAccess); cudaFree(d_clusterCount); 
	cudaDeviceReset();

	return 0;
}
