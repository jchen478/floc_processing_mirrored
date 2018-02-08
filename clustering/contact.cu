#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "contact.h"

using namespace std;

__device__ void parallel_sort_para(int mi, int nj, float sx, float sy, float sz, float pxmi, float pymi, float pzmi,
	float pxnj, float pynj, float pznj, float pdotp, float rp,
	float *xmin, float *ymin);

__global__ void contact(float **var, int **intVar){

	const int npcn = 2000;

	int mi = threadIdx.x + blockIdx.x*blockDim.x;
	int nfib = *intVar[6]; 
	int nseg = *intVar[7]; 
	int *potConSize = intVar[33];
	int *cluster = intVar[35]; 
	int nPair = potConSize[mi];
	if (nPair == 0) return;

	int *ncpf = intVar[2];
	int *clist = intVar[18];
	int maxCon = *intVar[25];
	int *potCon = intVar[32];

	float *rx = var[3];
	float *ry = var[4];
	float *rz = var[5];
	float *px = var[18];
	float *py = var[19];
	float *pz = var[20];
	float *fcx = var[66];
	float *fcy = var[67];
	float *fcz = var[68];
	float *tcx = var[69];
	float *tcy = var[70];
	float *tcz = var[71];
	float rp = *var[114];
	float over_cut = *var[118];
	float sidex = *var[119];
	float sidey = *var[120];
	float sidez = *var[121];
	float contact_cutoff = *var[128];
	float rep_cutoff = *var[129];
	float delta_rx = *var[138];
	float fstar = *var[144];
	float fact = *var[145];
	float Astar = *var[146];
	float decatt = *var[147];
	float *GijxV = var[149];
	float *GijyV = var[150];
	float *GijzV = var[151];
	float *GjixV = var[152];
	float *GjiyV = var[153];
	float *GjizV = var[154];
	float *nxV = var[155];
	float *nyV = var[156];
	float *nzV = var[157];
	float *gV = var[158];

	float rxmi, rymi, rzmi, rxnj, rynj, rznj;
	float pxmi, pymi, pzmi, pxnj, pynj, pznj;
	float sxx, syy, szz, corx, cory, corz;
	float rxmi_shift, rymi_shift, rzmi_shift;
	float pdotp, xmin, ymin, dx, dy, dz, sep;
	float xi[9], yj[9], gij, nijx, nijy, nijz, forc;
	float Gijx, Gijy, Gijz, Gjix, Gjiy, Gjiz, sep_tmp;

	int nP, nj, ipos, ith, oldmi, oldnj, m, n;

	rxmi = rx[mi]; rymi = ry[mi]; rzmi = rz[mi];
	pxmi = px[mi]; pymi = py[mi]; pzmi = pz[mi];

	for (nP = 0; nP < nPair; nP++){

		nj = potCon[mi * npcn + nP];
		//printf("in loop mi nj %4d %4d\n", mi, nj); 
		rxnj = rx[nj]; rynj = ry[nj]; rznj = rz[nj];
		pxnj = px[nj]; pynj = py[nj]; pznj = pz[nj];

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


		//if (mi == 799 && nj == 1312){
		//	printf("rx mi nj %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f \n", rxmi, rymi, rzmi, rxnj, rynj, rznj); 
		//}
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
		//printf("mi nj pdoptp sep xmin ymin %4d %4d %20.10f %20.10f\n", mi, nj, pdotp, pdotp*pdotp);
		// Check if segments are parallel
		if (fabsf(pdotp*pdotp - 1.0) <= 1.0e-6) {
			//printf("parallel sort %4d %4d\n", mi, nj); 
			parallel_sort_para(mi, nj, sxx, syy, szz, pxmi, pymi, pzmi,
				pxnj, pynj, pznj, pdotp, rp, &xmin, &ymin);
			sep = (sxx + ymin*pxnj - xmin*pxmi)*(sxx + ymin*pxnj - xmin*pxmi) +
				(syy + ymin*pynj - xmin*pymi)*(syy + ymin*pynj - xmin*pymi) +
				(szz + ymin*pznj - xmin*pzmi)*(szz + ymin*pznj - xmin*pzmi);
			//printf("parallel: mi nj %4d %4d sep %15.8f xmin ymin %15.8f %15.8f\n", mi, nj, sep, xmin, ymin); 
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
		//printf("gij %15.8f\n", gij); 
		//if (mi == 799 && nj == 1312){
		//	printf("gij %15.10f\n", gij); 
		//}
		nijx = (sxx + ymin*pxnj - xmin*pxmi) / gij;
		nijy = (syy + ymin*pynj - xmin*pymi) / gij;
		nijz = (szz + ymin*pznj - xmin*pzmi) / gij;
		Gijx = xmin*pxmi + gij*nijx / 2.0;
		Gijy = xmin*pymi + gij*nijy / 2.0;
		Gijz = xmin*pzmi + gij*nijz / 2.0;
		Gjix = ymin*pxnj - gij*nijx / 2.0;
		Gjiy = ymin*pynj - gij*nijy / 2.0;
		Gjiz = ymin*pznj - gij*nijz / 2.0;
		if (gij < 2.0){
			atomicAdd(intVar[4], 1); // overs
			//printf("overs: %4d %4d %15.10f %15.10f %15.10f \n", mi, nj, xmin, ymin, gij); 
		}
		if (gij < over_cut){
			gij = over_cut;
			//printf("overs: %4d %4d %15.10f %15.10f %15.10f \n", mi, nj, xmin, ymin, gij);
		}
		forc = fstar*expf(-fact*(gij - 2.0)) - Astar*expf(-decatt*(gij - 2.0)*(gij - 2.0));
		if (sep < rep_cutoff){
			//printf("rep %4d %4d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n", mi, nj, nijx, nijy, nijz, forc, sep, gij); 
			atomicAdd(fcx + mi, -nijx*forc);
			atomicAdd(fcy + mi, -nijy*forc);
			atomicAdd(fcz + mi, -nijz*forc);
			atomicAdd(tcx + mi, -forc*xmin*(pymi*nijz - pzmi*nijy));
			atomicAdd(tcy + mi, -forc*xmin*(pzmi*nijx - pxmi*nijz));
			atomicAdd(tcz + mi, -forc*xmin*(pxmi*nijy - pymi*nijx));
			atomicAdd(fcx + nj, nijx*forc);
			atomicAdd(fcy + nj, nijy*forc);
			atomicAdd(fcz + nj, nijz*forc);
			atomicAdd(tcx + nj, forc*ymin*(pynj*nijz - pznj*nijy));
			atomicAdd(tcy + nj, forc*ymin*(pznj*nijx - pxnj*nijz));
			atomicAdd(tcz + nj, forc*ymin*(pxnj*nijy - pynj*nijx));
		}
		if (sep < contact_cutoff){
			oldmi = atomicAdd(ncpf + mi, 1);
			oldnj = atomicAdd(ncpf + nj, 1);
			clist[mi*maxCon + oldmi] = nj;
			clist[nj*maxCon + oldnj] = mi;
			nxV[mi*maxCon + oldmi] = nijx;
			nyV[mi*maxCon + oldmi] = nijy;
			nzV[mi*maxCon + oldmi] = nijz;
			gV[mi*maxCon + oldmi] = gij;
			GijxV[mi*maxCon + oldmi] = Gijx;
			GijyV[mi*maxCon + oldmi] = Gijy;
			GijzV[mi*maxCon + oldmi] = Gijz;
			GjixV[mi*maxCon + oldmi] = Gjix;
			GjiyV[mi*maxCon + oldmi] = Gjiy;
			GjizV[mi*maxCon + oldmi] = Gjiz;
			m = mi / nseg; 
			n = nj / nseg; 
			cluster[m*nfib + n] = 1; 
			cluster[n*nfib + m] = 1; 
			//printf("%6d %6d %15.10f %15.10f %15.10f\n", mi, nj, xmin, ymin, gij); 
			if (oldmi >= maxCon || oldnj >= maxCon){
				printf("link: mi nj %7d %7d exceeding maxCon, allocate more space\n", mi, nj);
			}			
		}
	}

}
__device__ void parallel_sort_para(int mi, int nj, float sx, float sy, float sz, float pxmi, float pymi, float pzmi,
	float pxnj, float pynj, float pznj, float pdotp, float rp,
	float *xmin, float *ymin){

	//printf("accessing parallel_sort_para\n");
	float posneg, pn2, dist, sijp, sijm, sjip, sjim;

	//printf("%4d %4d %15.10f p %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n", mi, nj, pdotp, pxmi, pymi, pzmi, pxnj, pynj, pznj); 

	// The different end point to fiber contact points
	sijp = pxmi*sx + pymi*sy + pzmi*sz + rp*pdotp;
	sijm = pxmi*sx + pymi*sy + pzmi*sz - rp*pdotp;
	sjip = -(pxnj*sx + pynj*sy + pznj*sz) + rp*pdotp;
	sjim = -(pxnj*sx + pynj*sy + pznj*sz) - rp*pdotp;

	//printf("parallel\n"); 
	//printf("%4d %4d sx sy sz %15.10f %15.10f %15.10f sijp sijm sjip sjim %15.10f %15.10f %15.10f %15.10f\n", mi, nj, sx, sy, sz, sijp, sijm, sjip, sjim); 

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
	//printf("xmin ymin in %12.8f %12.8f\n", *xmin, *ymin);
	//printf("xmin ymin out %4d %4d %16.10f %16.10f\n", mi, nj, *xmin, *ymin); 
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
		//printf("xmin ymin in %12.8f %12.8f\n", *xmin, *ymin);
	}
}
