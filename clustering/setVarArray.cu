#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include "setVarArray.h"

using namespace std;

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
	float *nxV, float *nyV, float *nzV, float *gV){

	intVar[0] = ifiber;
	intVar[1] = ncnt;
	intVar[2] = ncpf;
	intVar[3] = num_groups;
	intVar[4] = overs;
	intVar[5] = total_contacts;
	intVar[6] = nfib;
	intVar[7] = nseg;
	intVar[8] = fac;
	intVar[12] = bin;
	intVar[13] = list;
	intVar[14] = bnei;
	intVar[15] = nxbin;
	intVar[16] = nybin;
	intVar[17] = nzbin;
	intVar[18] = clist;
	intVar[19] = status;
	intVar[20] = lead_clist;
	intVar[21] = nc;
	intVar[22] = bdimx;
	intVar[23] = bdimy;
	intVar[24] = bdimz;
	intVar[25] = maxCon;
	intVar[26] = clist_pos;
	intVar[27] = maxGr;
	intVar[28] = maxBin;
	intVar[31] = bnum;
	intVar[32] = potCon;
	intVar[33] = potConSize;
	intVar[34] = groupId;
	intVar[35] = cluster;
	intVar[36] = clusterAccess;
	intVar[37] = clusterCount;
	var[0] = rcmx;
	var[1] = rcmy;
	var[2] = rcmz;
	var[3] = rx;
	var[4] = ry;
	var[5] = rz;
	var[6] = q0;
	var[7] = q1;
	var[8] = q2;
	var[9] = q3;
	var[10] = q0dot;
	var[11] = q1dot;
	var[12] = q2dot;
	var[13] = q3dot;
	var[14] = qe0;
	var[15] = qe1;
	var[16] = qe2;
	var[17] = qe3;
	var[18] = px;
	var[19] = py;
	var[20] = pz;
	var[21] = ucmx;
	var[22] = ucmy;
	var[23] = ucmz;
	var[24] = ucox;
	var[25] = ucoy;
	var[26] = ucoz;
	var[27] = ux;
	var[28] = uy;
	var[29] = uz;
	var[30] = uxfl;
	var[31] = uyfl;
	var[32] = uzfl;
	var[33] = wx;
	var[34] = wy;
	var[35] = wz;
	var[36] = R11;
	var[37] = R12;
	var[38] = R13;
	var[39] = R21;
	var[40] = R22;
	var[41] = R23;
	var[42] = R11eq;
	var[43] = R12eq;
	var[44] = R13eq;
	var[45] = R21eq;
	var[46] = R22eq;
	var[47] = R23eq;
	var[48] = R31eq;
	var[49] = R32eq;
	var[50] = R33eq;
	var[60] = fx;
	var[61] = fy;
	var[62] = fz;
	var[63] = tx;
	var[64] = ty;
	var[65] = tz;
	var[66] = fcx;
	var[67] = fcy;
	var[68] = fcz;
	var[69] = tcx;
	var[70] = tcy;
	var[71] = tcz;
	var[72] = fbx;
	var[73] = fby;
	var[74] = fbz;
	var[75] = tbx;
	var[76] = tby;
	var[77] = tbz;
	var[93] = C0;
	var[94] = C1;
	var[95] = C2;
	var[96] = C3;
	var[97] = C4;
	var[98] = g;
	var[99] = duxdx;
	var[100] = duxdy;
	var[101] = duxdz;
	var[102] = duydx;
	var[103] = duydy;
	var[104] = duydz;
	var[105] = duzdx;
	var[106] = duzdy;
	var[107] = duzdz;
	var[108] = dx;
	var[109] = dy;
	var[110] = dz;
	var[111] = Omega_x;
	var[112] = Omega_y;
	var[113] = Omega_z;
	var[114] = rp;
	var[115] = kb;
	var[116] = kt;
	var[117] = dt;
	var[118] = over_cut;
	var[119] = sidex;
	var[120] = sidey;
	var[121] = sidez;
	var[122] = E11;
	var[123] = E12;
	var[124] = E13;
	var[125] = E22;
	var[126] = E23;
	var[127] = E33;
	var[128] = contact_cutoff;
	var[129] = rep_cutoff;
	var[130] = neighb_cutoff;
	var[131] = Ya;
	var[132] = Yc;
	var[133] = Yh;
	var[136] = thetaeq;
	var[137] = phieq;
	var[138] = delta_rx;
	var[139] = mu_stat;
	var[140] = mu_kin;
	var[144] = fstar;
	var[145] = fact;
	var[146] = Astar;
	var[147] = decatt;
	var[148] = elf;
	var[149] = GijxV;
	var[150] = GijyV;
	var[151] = GijzV;
	var[152] = GjixV;
	var[153] = GjiyV;
	var[154] = GjizV;
	var[155] = nxV;
	var[156] = nyV;
	var[157] = nzV;
	var[158] = gV;
}