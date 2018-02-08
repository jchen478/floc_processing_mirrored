//
//  setVarArray.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 10/31/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__setVarArray__
#define __flexfric__setVarArray__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void setVarArray(float **var, int **intVar, int *ifiber,
	int *ncnt, int *ncpf, int *num_groups, int *overs, int *total_contacts,
	int *nfib, int *nseg, int *fac,
	int *bin, int *list,
	int *bnei, int *nxbin, int *nybin, int *nzbin, int *clist, int *status, int *lead_clist, int *nc,
	int *bdimx, int *bdimy, int *bdimz, int *maxCon, int *clist_pos, int *maxGr, int *maxBin,
	int *bnum, int *potCon, int *potConSize, int *groupId,
	int *cluster, int * clusterAccess, int *clusterCount,
	float *rcmx, float *rcmy, float *rcmz,
	float *rx, float *ry, float *rz, float *q0, float *q1, float *q2, float *q3, float *q0dot,
	float *q1dot, float *q2dot, float *q3dot, float *qe0, float *qe1, float *qe2, float *qe3,
	float *px, float *py, float *pz, float *ucmx, float *ucmy, float *ucmz, float *ucox, float *ucoy,
	float *ucoz, float *ux, float *uy, float *uz, float *uxfl, float *uyfl, float *uzfl, float *wx,
	float *wy, float *wz, float *R11, float *R12, float *R13, float *R21, float *R22, float *R23,
	float *R11eq, float *R12eq, float *R13eq, float *R21eq, float *R22eq, float *R23eq, float *R31eq,
	float *R32eq, float *R33eq, float *fx, float *fy, float *fz, float *tx, float *ty,
	float *tz, float *fcx, float *fcy, float *fcz, float *tcx, float *tcy, float *tcz, float *fbx,
	float *fby, float *fbz, float *tbx, float *tby, float *tbz, 

	float *C0, float *C1, float *C2, float *C3, float *C4, float *g, float *duxdx, float *duxdy,
	float *duxdz, float *duydx, float *duydy, float *duydz, float *duzdx, float *duzdy, float *duzdz,
	float *dx, float *dy, float *dz, float *Omega_x, float *Omega_y, float *Omega_z,
	float *rp, float *kb, float *kt, float *dt, float *over_cut,
	float *sidex, float *sidey, float *sidez, float *E11, float *E12, float *E13, float *E22,
	float *E23, float *E33, float *contact_cutoff, float *rep_cutoff, float *neighb_cutoff,
	float *Ya, float *Yc, float *Yh, float *thetaeq, float *phieq,
	float *delta_rx, float *mu_stat, float *mu_kin,
	float *fstar, float *fact, float *Astar, float *decatt, float *elf,
	float *GijxV, float *GijyV, float *GijzV, float *GjixV, float *GjiyV, float *GjizV,
	float *nxV, float *nyV, float *nzV, float *gV);
#endif /* defined(__flexfric__setVarArray__) */
