#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include "initialize.h"

using namespace std;

__global__ void initialize(float **var, int **intVar){

	int  nfib = *intVar[6];
	int  nseg = *intVar[7];

	float *q0 = var[6];
	float *q1 = var[7];
	float *q2 = var[8];
	float *q3 = var[9];
	float *q0dot = var[10];
	float *q1dot = var[11];
	float *q2dot = var[12];
	float *q3dot = var[13];
	float *qe0 = var[14];
	float *qe1 = var[15];
	float *qe2 = var[16];
	float *qe3 = var[17];
	float *px = var[18];
	float *py = var[19];
	float *pz = var[20];
	float *ux = var[27];
	float *uy = var[28];
	float *uz = var[29];
	float *uxfl = var[30];
	float *uyfl = var[31];
	float *uzfl = var[32];
	float *wx = var[33];
	float *wy = var[34];
	float *wz = var[35];
	float *R11 = var[36];
	float *R12 = var[37];
	float *R13 = var[38];
	float *R21 = var[39];
	float *R22 = var[40];
	float *R23 = var[41];
	float *R11eq = var[42];
	float *R12eq = var[43];
	float *R13eq = var[44];
	float *R21eq = var[45];
	float *R22eq = var[46];
	float *R23eq = var[47];
	float *R31eq = var[48];
	float *R32eq = var[49];
	float *R33eq = var[50];
	float *tx = var[63];
	float *ty = var[64];
	float *tz = var[65];
	float *fcx = var[66];
	float *fcy = var[67];
	float *fcz = var[68];
	float *tcx = var[69];
	float *tcy = var[70];
	float *tcz = var[71];
	float *fbx = var[72];
	float *fby = var[73];
	float *fbz = var[74];
	float *tbx = var[75];
	float *tby = var[76];
	float *tbz = var[77];
	float elf = *var[148];
	float *thetaeq = var[136];
	float *phieq = var[137];

	float q0mi, q1mi, q2mi, q3mi;
	float pxmi, pymi, pzmi, thetaeqmi, phieqmi;

	int mi = threadIdx.x + blockIdx.x*blockDim.x;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	int m, i;
	float dum;
	// temporary variables

	m = mi / nseg;
	i = mi - m*nseg;
	//printf("accessing initialize\n"); 
	// zero varaibles
	ux[mi] = 0.0; uy[mi] = 0.0; uz[mi] = 0.0;
	wx[mi] = 0.0; wy[mi] = 0.0; wz[mi] = 0.0;
	q0dot[mi] = 0.0; q1dot[mi] = 0.0; q2dot[mi] = 0.0; q3dot[mi] = 0.0;
	qe0[mi] = 0.0; qe1[mi] = 0.0; qe2[mi] = 0.0; qe3[mi] = 0.0;
	uxfl[mi] = 0.0; uyfl[mi] = 0.0; uzfl[mi] = 0.0;
	fcx[mi] = 0.0; fcy[mi] = 0.0; fcz[mi] = 0.0;
	tx[mi] = 0.0; ty[mi] = 0.0; tz[mi] = 0.0;
	tcx[mi] = 0.0; tcy[mi] = 0.0; tcz[mi] = 0.0;

	q0mi = q0[mi]; q1mi = q1[mi]; q2mi = q2[mi]; q3mi = q3[mi];

	// Normalize Euler Parameters
	dum = sqrtf(q0mi*q0mi + q1mi*q1mi + q2mi*q2mi + q3mi*q3mi);

	q0mi /= dum;
	q1mi /= dum;
	q2mi /= dum;
	q3mi /= dum;

	q0[mi] = q0mi; q1[mi] = q1mi; q2[mi] = q2mi; q3[mi] = q3mi;

	// Find rotation matrix
	R11[mi] = 2.0*(q0mi*q0mi + q1mi*q1mi) - 1.0;
	R12[mi] = 2.0*(q1mi*q2mi + q0mi*q3mi);
	R13[mi] = 2.0*(q1mi*q3mi - q0mi*q2mi);
	R21[mi] = 2.0*(q1mi*q2mi - q0mi*q3mi);
	R22[mi] = 2.0*(q0mi*q0mi + q2mi*q2mi) - 1.0;
	R23[mi] = 2.0*(q3mi*q2mi + q0mi*q1mi);
	pxmi = 2.0*(q1mi*q3mi + q0mi*q2mi);
	pymi = 2.0*(q3mi*q2mi - q0mi*q1mi);
	pzmi = 2.0*(q0mi*q0mi + q3mi*q3mi) - 1.0;
	dum = sqrtf(pxmi*pxmi + pymi*pymi + pzmi*pzmi);
	pxmi /= dum;
	pymi /= dum;
	pzmi /= dum;
	px[mi] = pxmi;
	py[mi] = pymi;
	pz[mi] = pzmi;

	// Assign body forces and torques
	fbx[mi] = 0.0;
	fby[mi] = 0.0;
	fbz[mi] = 0.0;
	tbx[mi] = elf*pzmi*pymi;
	tby[mi] = -elf*pzmi*pxmi;
	tbz[mi] = 0.0;

	// Make the equilibrium rotation matrix
	if (i > 0) {
		thetaeqmi = thetaeq[mi]; phieqmi = phieq[mi];
		R11eq[mi] = cosf(thetaeqmi)*cosf(phieqmi);
		R12eq[mi] = cosf(thetaeqmi)*sinf(phieqmi);
		R13eq[mi] = -sinf(thetaeqmi);
		R21eq[mi] = -sinf(phieqmi);
		R22eq[mi] = cosf(phieqmi);
		R23eq[mi] = 0.0;
		R31eq[mi] = sinf(thetaeqmi)*cosf(phieqmi);
		R32eq[mi] = sinf(thetaeqmi)*sinf(phieqmi);
		R33eq[mi] = cosf(thetaeqmi);
	}
}
