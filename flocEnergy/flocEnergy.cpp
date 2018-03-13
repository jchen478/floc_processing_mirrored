#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

float elasticEnergy(float *px, float *py, float *pz,
	float *rx, float *ry, float *rz,
	float *q0, float *q1, float *q2, float *q3,
	float *R11, float *R12, float *R13,
	float *R21, float *R22, float *R23,
	float *R11eq, float *R12eq, float *R13eq,
	float *R21eq, float *R22eq, float *R23eq,
	float *R31eq, float *R32eq, float *R33eq,
	float kb, int nfibC_ID, int nseg);

int main(){

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////

	// Cluster related files
	FILE *ClusterInfo, *Cluster_results, *Cluster;
	FILE *flocEnergy;
	FILE *Parameters, *Equilibrium_Angles;
	FILE *rxfile, *ryfile, *rzfile;
	FILE *pxfile, *pyfile, *pzfile;
	FILE *q0file, *q1file, *q2file, *q3file;

	ClusterInfo = fopen("ClusterInfo.in", "r");
	Cluster_results = fopen("Cluster_results.txt", "r");
	Cluster = fopen("Cluster.txt", "r");
	flocEnergy = fopen("flocEnergy.txt", "w");
	Parameters = fopen("Parameters.in", "r");
	Equilibrium_Angles = fopen("../../Equilibrium_Angles.in", "r");
	rxfile = fopen("../../rx.txt", "rb");
	ryfile = fopen("../../ry.txt", "rb");
	rzfile = fopen("../../rz.txt", "rb");
	pxfile = fopen("../../px.txt", "rb");
	pyfile = fopen("../../py.txt", "rb");
	pzfile = fopen("../../pz.txt", "rb");
	q0file = fopen("../../q0.txt", "rb");
	q1file = fopen("../../q1.txt", "rb");
	q2file = fopen("../../q2.txt", "rb");
	q3file = fopen("../../q3.txt", "rb");

	int nC, nfibC_ID, nfib, nseg, config_write, *nfibC;
	float dt, strain, kb, Eelastic, dum;
	// nC - number of clusters
	// nfibC_ID - number of fiber in specified floc
	int m, i, n, j, mi, nj, ni, cID, idum, idum1, idum2;
	// dummy variables

	fscanf(ClusterInfo, "%d", &idum);
	fscanf(ClusterInfo, "%*[^\n]%d", &nC);
	fclose(ClusterInfo);

	nfibC = (int*)malloc(nC*sizeof(int));
	// number of fiber in cluster

	for (i = 0; i < nC; i++){
		fscanf(Cluster_results, "%d", &idum);
		fscanf(Cluster_results, " %d", &nfibC[i]);
	}

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &kb);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dt);
	fscanf(Parameters, "%*[^\n]%f", &strain);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, " %f", &dum);
	fscanf(Parameters, " %f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &config_write);
	fclose(Parameters);

	int step, nConfig;
	float *thetaeq, *phieq, theta, phi;
	float *R11eq, *R12eq, *R13eq, *R21eq, *R22eq, *R23eq;
	float *R11, *R12, *R13, *R21, *R22, *R23;
	float *R31eq, *R32eq, *R33eq;
	float *px, *py, *pz, *pxC, *pyC, *pzC;
	float *rx, *ry, *rz, *rxC, *ryC, *rzC;
	float *q0C, *q1C, *q2C, *q3C, *q0, *q1, *q2, *q3;
	int *fibID = (int*)malloc(nfib*sizeof(int));

	rxC = (float*)malloc(nfib*nseg*sizeof(float));
	ryC = (float*)malloc(nfib*nseg*sizeof(float));
	rzC = (float*)malloc(nfib*nseg*sizeof(float));
	pxC = (float*)malloc(nfib*nseg*sizeof(float));
	pyC = (float*)malloc(nfib*nseg*sizeof(float));
	pzC = (float*)malloc(nfib*nseg*sizeof(float));
	q0C = (float*)malloc(nfib*nseg*sizeof(float));
	q1C = (float*)malloc(nfib*nseg*sizeof(float));
	q2C = (float*)malloc(nfib*nseg*sizeof(float));
	q3C = (float*)malloc(nfib*nseg*sizeof(float));

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
	thetaeq = (float*)malloc(nfib*nseg*sizeof(float));
	phieq = (float*)malloc(nfib*nseg*sizeof(float));


	R11 = (float*)malloc(nfib*nseg*sizeof(float));
	R12 = (float*)malloc(nfib*nseg*sizeof(float));
	R13 = (float*)malloc(nfib*nseg*sizeof(float));
	R21 = (float*)malloc(nfib*nseg*sizeof(float));
	R22 = (float*)malloc(nfib*nseg*sizeof(float));
	R23 = (float*)malloc(nfib*nseg*sizeof(float));
	R11eq = (float*)malloc(nfib*nseg*sizeof(float));
	R12eq = (float*)malloc(nfib*nseg*sizeof(float));
	R13eq = (float*)malloc(nfib*nseg*sizeof(float));
	R21eq = (float*)malloc(nfib*nseg*sizeof(float));
	R22eq = (float*)malloc(nfib*nseg*sizeof(float));
	R23eq = (float*)malloc(nfib*nseg*sizeof(float));
	R31eq = (float*)malloc(nfib*nseg*sizeof(float));
	R32eq = (float*)malloc(nfib*nseg*sizeof(float));
	R33eq = (float*)malloc(nfib*nseg*sizeof(float));

	nConfig = int(strain / dt / float(config_write)) + 1;

	// Read Equilibrium Angles //
	for (m = 0; m < nfib; m++){
		for (i = 0; i < nseg; i++){
			mi = m*nseg + i;
			if (i > 0){
				fscanf(Equilibrium_Angles, "%d %d %f %f", &idum1, &idum2, thetaeq + mi, phieq + mi);
			}
		}
	}
	fclose(Equilibrium_Angles);

	// read till last frame
	for (step = 0; step < nConfig; step++){

		// r
		fread(&dum, sizeof(float), 1, rxfile);
		fread(rx, sizeof(float), nfib*nseg, rxfile);
		fread(&dum, sizeof(float), 1, ryfile);
		fread(ry, sizeof(float), nfib*nseg, ryfile);
		fread(&dum, sizeof(float), 1, rzfile);
		fread(rz, sizeof(float), nfib*nseg, rzfile);

		// p
		fread(&dum, sizeof(float), 1, pxfile);
		fread(px, sizeof(float), nfib*nseg, pxfile);
		fread(&dum, sizeof(float), 1, pyfile);
		fread(py, sizeof(float), nfib*nseg, pyfile);
		fread(&dum, sizeof(float), 1, pzfile);
		fread(pz, sizeof(float), nfib*nseg, pzfile);

		// q
		fread(&dum, sizeof(float), 1, q0file);
		fread(q0, sizeof(float), nfib*nseg, q0file);
		fread(&dum, sizeof(float), 1, q1file);
		fread(q1, sizeof(float), nfib*nseg, q1file);
		fread(&dum, sizeof(float), 1, q2file);
		fread(q2, sizeof(float), nfib*nseg, q2file);
		fread(&dum, sizeof(float), 1, q3file);
		fread(q3, sizeof(float), nfib*nseg, q3file);
	}

	// calculate elastic energy for each floc
	for (cID = 0; cID < nC; cID++){

		// number of fibers in floc
		nfibC_ID = nfibC[cID];

		fscanf(Cluster, " %d", &idum);
		fscanf(Cluster, " %d", &idum);
		for (m = 0; m < nfibC_ID; m++){
			// indices of fibers in floc
			fscanf(Cluster, " %d", &fibID[m]);
			n = fibID[m];
			// data specific for fiber
			for (i = 0; i < nseg; i++){
				mi = m*nseg + i;
				ni = n*nseg + i;
				rxC[mi] = rx[ni];
				ryC[mi] = ry[ni];
				rzC[mi] = rz[ni];
				pxC[mi] = px[ni];
				pyC[mi] = py[ni];
				pzC[mi] = pz[ni];
				q0C[mi] = q0[ni];
				q1C[mi] = q1[ni];
				q2C[mi] = q2[ni];
				q3C[mi] = q3[ni];
				if (i >= 1){
					theta = thetaeq[ni];
					phi = phieq[ni];
					R11eq[mi] = cosf(theta)*cosf(phi);
					R12eq[mi] = cosf(theta)*sinf(phi);
					R13eq[mi] = -sinf(theta);
					R21eq[mi] = -sinf(phi);
					R22eq[mi] = cosf(phi);
					R23eq[mi] = 0.0;
					R31eq[mi] = sinf(theta)*cosf(phi);
					R32eq[mi] = sinf(theta)*sinf(phi);
					R33eq[mi] = cosf(theta);
				}
			}
		}

		// calculate elastic energy
		Eelastic = elasticEnergy(pxC, pyC, pzC,
			rxC, ryC, rzC,
			q0C, q1C, q2C, q3C,
			R11, R12, R13, R21, R22, R23,
			R11eq, R12eq, R13eq, R21eq, R22eq,
			R23eq, R31eq, R32eq, R33eq,
			kb, nfibC_ID, nseg);
		fprintf(flocEnergy, "%4d %4d %10.6f\n", cID, nfibC_ID, Eelastic);
		//printf("%4d %4d %10.6f\n", cID, nfibC_ID, Eelastic);
		fflush(flocEnergy);
	}

	free(rx); free(ry); free(rz);
	free(px); free(py); free(pz);
	free(q0); free(q1); free(q2); free(q3);
	free(rxC); free(ryC); free(rzC);
	free(pxC); free(pyC); free(pzC);
	free(q0C); free(q1C); free(q2C); free(q3C);
	free(thetaeq); free(phieq);

	free(R11);
	free(R12);
	free(R13);
	free(R21);
	free(R22);
	free(R23);
	free(R11eq);
	free(R12eq);
	free(R13eq);
	free(R21eq);
	free(R22eq);
	free(R23eq);
	free(R31eq);
	free(R32eq);
	free(R33eq);
	free(nfibC);
	free(fibID);

	fclose(rxfile); fclose(ryfile); fclose(rzfile);
	fclose(pxfile); fclose(pyfile); fclose(pzfile);
	fclose(q0file); fclose(q1file); fclose(q2file); fclose(q3file);
	fclose(Cluster_results);
	fclose(Cluster);
	fclose(flocEnergy);

	return 0;
}

float elasticEnergy(float *px, float *py, float *pz,
	float *rx, float *ry, float *rz,
	float *q0, float *q1, float *q2, float *q3,
	float *R11, float *R12, float *R13,
	float *R21, float *R22, float *R23,
	float *R11eq, float *R12eq, float *R13eq,
	float *R21eq, float *R22eq, float *R23eq,
	float *R31eq, float *R32eq, float *R33eq,
	float kb, int nfibC_ID, int nseg){

	float Energy = 0.0;
	float cx, cy, cz, ang_theta, ang_phi;
	float zeqx, zeqy, zeqz, zdotz, dum;
	float yitx, yity, yitz, yieqx, yieqy, yieqz;
	int m, i, mi;

	for (m = 0; m < nfibC_ID; m++){
		for (i = 0; i < nseg; i++){
			mi = m*nseg + i;
			R11[mi] = 2.0*(q0[mi] * q0[mi] + q1[mi] * q1[mi]) - 1.0;
			R12[mi] = 2.0*(q1[mi] * q2[mi] + q0[mi] * q3[mi]);
			R13[mi] = 2.0*(q1[mi] * q3[mi] - q0[mi] * q2[mi]);
			R21[mi] = 2.0*(q1[mi] * q2[mi] - q0[mi] * q3[mi]);
			R22[mi] = 2.0*(q0[mi] * q0[mi] + q2[mi] * q2[mi]) - 1.0;
			R23[mi] = 2.0*(q3[mi] * q2[mi] + q0[mi] * q1[mi]);
		}
	}

	for (m = 0; m < nfibC_ID; m++){
		for (i = 1; i < nseg; i++){
			ang_theta = 0.0;
			ang_phi = 0.0;
			mi = m*nseg + i;
			// calculations for theta
			zeqx = R11[mi - 1] * R31eq[mi] + R21[mi - 1] * R32eq[mi] + px[mi - 1] * R33eq[mi];
			zeqy = R12[mi - 1] * R31eq[mi] + R22[mi - 1] * R32eq[mi] + py[mi - 1] * R33eq[mi];
			zeqz = R13[mi - 1] * R31eq[mi] + R23[mi - 1] * R32eq[mi] + pz[mi - 1] * R33eq[mi];
			zdotz = px[mi] * zeqx + py[mi] * zeqy + pz[mi] * zeqz;

			if (zdotz <= 1.0 - 1.0E-6){
				if (zdotz < -1.0){
					zdotz = -1.0;
				}
				ang_theta = acosf(zdotz);
			}

			// calculations for phi
			cx = rx[mi] - rx[mi - 1];
			cy = ry[mi] - ry[mi - 1];
			cz = rz[mi] - rz[mi - 1];
			dum = sqrtf(cx*cx + cy*cy + cz*cz);
			cx = cx / dum;
			cy = cy / dum;
			cz = cz / dum;
			zeqx = R11[mi - 1] * R21eq[mi] + R21[mi - 1] * R22eq[mi] + px[mi - 1] * R23eq[mi];
			zeqy = R12[mi - 1] * R21eq[mi] + R22[mi - 1] * R22eq[mi] + py[mi - 1] * R23eq[mi];
			zeqz = R13[mi - 1] * R21eq[mi] + R23[mi - 1] * R22eq[mi] + pz[mi - 1] * R23eq[mi];
			yitx = R21[mi] - cx*(cx*R21[mi] + cy*R22[mi] + cz*R23[mi]);
			yity = R22[mi] - cy*(cx*R21[mi] + cy*R22[mi] + cz*R23[mi]);
			yitz = R23[mi] - cz*(cx*R21[mi] + cy*R22[mi] + cz*R23[mi]);
			dum = sqrtf(yitx*yitx + yity*yity + yitz*yitz);

			yitx = yitx / dum;
			yity = yity / dum;
			yitz = yitz / dum;
			yieqx = zeqx - cx*(cx*zeqx + cy*zeqy + cz*zeqz);
			yieqy = zeqy - cy*(cx*zeqx + cy*zeqy + cz*zeqz);
			yieqz = zeqz - cz*(cx*zeqx + cy*zeqy + cz*zeqz);
			dum = sqrtf(yieqx*yieqx + yieqy*yieqy + yieqz*yieqz);

			yieqx = yieqx / dum;
			yieqy = yieqy / dum;
			yieqz = yieqz / dum;
			zdotz = yitx*yieqx + yity*yieqy + yitz*yieqz;

			if (zdotz <= (1.0 - 1.0E-6)){
				if (zdotz < -1.0){
					zdotz = -1.0;
				}
				ang_phi = acosf(zdotz);
			}
			Energy += ang_theta*ang_theta + 0.67*ang_phi*ang_phi;
		}
	}


	return kb*Energy / float(nfibC_ID);
}

