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
	FILE *flocElastic; 

	ClusterInfo = fopen("ClusterInfo.in","r");
	Cluster_results = fopen("../Cluster_results.txt","r");
	Cluster = fopen("../Cluster.txt","r");
	flocElastic = fopen("flocElastic.txt","w");

	// Obtain cluster info to read the files
	int cID, nC;
	fscanf(ClusterInfo, "%d", &cID);
	fscanf(ClusterInfo, "%*[^\n]%d", &nC);
	fclose(ClusterInfo); 
	
	// dummy variables
	int m, i, j, mi, idum, idum1, idum2;
	
	// allocate arrays to store cluster info
	int *nfibC = (int*)malloc(nC*sizeof(int)); 
	// number of fiber in cluster

	for (i = 0; i < nC; i++){
		fscanf(Cluster_results, "%d", &idum);
		fscanf(Cluster_results, " %d", &nfibC[i]);
	}	
	

	int nfibC_ID; // number of fiber in specified ID 

	nfibC_ID = nfibC[cID]; 

	// allocate array for specified floc
	int *fibID = (int*)malloc(nfibC_ID*sizeof(int));

	for (i = 0; i < nC; i++){

		fscanf(Cluster, " %d", &idum);
		fscanf(Cluster, " %d", &idum);
		if( i == cID){
			printf("nfibC_ID %4d\n",nfibC_ID);
			for(j = 0; j < nfibC_ID; j++){
				fscanf(Cluster, " %d", &fibID[j]);
				printf("%d\n", fibID[j]);
			}
			break; 
		}
		for(j = 0; j < nfibC[i]; j++){
			fscanf(Cluster, " %d", &idum);
		}
	}

	free(nfibC);
	fclose(Cluster_results); 
	fclose(Cluster); 

	// read files and print out specified floc
	FILE *Parameters, *Equilibrium_Angles; 
	FILE *rxfile, *ryfile, *rzfile;
	FILE *rxfileFloc, *ryfileFloc, *rzfileFloc;
	FILE *pxfile, *pyfile, *pzfile;
	FILE *pxfileFloc, *pyfileFloc, *pzfileFloc;
	FILE *uxfile, *uyfile, *uzfile;
	FILE *uxfileFloc, *uyfileFloc, *uzfileFloc;
	FILE *wxfile, *wyfile, *wzfile;
	FILE *wxfileFloc, *wyfileFloc, *wzfileFloc;
	FILE *q0file, *q1file, *q2file, *q3file;
	FILE *q0fileFloc, *q1fileFloc, *q2fileFloc, *q3fileFloc;
	FILE *center_mass, *center_massFloc; 
	
	int nfib, nseg, config_write;
	float dt, strain, kb, Eelastic, dum;
	
	Parameters = fopen("../Parameters.in", "r");
	Equilibrium_Angles = fopen("../../../Equilibrium_Angles.in", "r");
	// rcm
	center_mass = fopen("../../../center_mass.txt", "rb");
	center_massFloc = fopen("center_mass.txt", "wb");

	// r
	rxfile = fopen("../../../rx.txt", "rb");
	ryfile = fopen("../../../ry.txt", "rb");
	rzfile = fopen("../../../rz.txt", "rb");
	rxfileFloc = fopen("rx.txt", "wb");
	ryfileFloc = fopen("ry.txt", "wb");
	rzfileFloc = fopen("rz.txt", "wb");
	// p
	pxfile = fopen("../../../px.txt", "rb");
	pyfile = fopen("../../../py.txt", "rb");
	pzfile = fopen("../../../pz.txt", "rb");
	pxfileFloc = fopen("px.txt", "wb");
	pyfileFloc = fopen("py.txt", "wb");
	pzfileFloc = fopen("pz.txt", "wb");
	// u
	uxfile = fopen("../../../ux.txt", "rb");
	uyfile = fopen("../../../uy.txt", "rb");
	uzfile = fopen("../../../uz.txt", "rb");
	uxfileFloc = fopen("ux.txt", "wb");
	uyfileFloc = fopen("uy.txt", "wb");
	uzfileFloc = fopen("uz.txt", "wb");
	// w
	wxfile = fopen("../../../wx.txt", "rb");
	wyfile = fopen("../../../wy.txt", "rb");
	wzfile = fopen("../../../wz.txt", "rb");
	wxfileFloc = fopen("wx.txt", "wb");
	wyfileFloc = fopen("wy.txt", "wb");
	wzfileFloc = fopen("wz.txt", "wb");
	// q
	q0file = fopen("../../../q0.txt", "rb");
	q1file = fopen("../../../q1.txt", "rb");
	q2file = fopen("../../../q2.txt", "rb");
	q3file = fopen("../../../q3.txt", "rb");
	q0fileFloc = fopen("q0.txt", "wb");
	q1fileFloc = fopen("q1.txt", "wb");
	q2fileFloc = fopen("q2.txt", "wb");
	q3fileFloc = fopen("q3.txt", "wb");

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
	float *r, *rc; 
	float *rC, *rcC; 
	float *thetaeq, *phieq;
	float *thetaeqC, *phieqC;
	float *R11eq, *R12eq, *R13eq, *R21eq, *R22eq, *R23eq;
	float *R11, *R12, *R13, *R21, *R22, *R23;
	float *R31eq, *R32eq, *R33eq;
	float *pxC, *pyC, *pzC, *q0C, *q1C, *q2C, *q3C;
	float *rxC, *ryC, *rzC;

	rC = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	rcC = (float*)malloc(nfibC_ID*sizeof(float));
	rxC = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	ryC = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	rzC = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	pxC = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	pyC = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	pzC = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	q0C = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	q1C = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	q2C = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	q3C = (float*)malloc(nfibC_ID*nseg*sizeof(float));
	r = (float*)malloc(nfib*nseg*sizeof(float)); 
	rc = (float*)malloc(nfib*sizeof(float));
	thetaeq = (float*)malloc(nfib*nseg*sizeof(float)); 
	phieq = (float*)malloc(nfib*nseg*sizeof(float)); 
	thetaeqC = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	phieqC = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	
	R11 = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R12 = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R13 = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R21 = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R22 = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R23 = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R11eq = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R12eq = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R13eq = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R21eq = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R22eq = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R23eq = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R31eq = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R32eq = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	R33eq = (float*)malloc(nfibC_ID*nseg*sizeof(float)); 
	
	nConfig = int(strain / dt / float(config_write))+1;
	
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

	for (i = 0; i < nfibC_ID; i++){
		for (j = 1; j < nseg; j++){
			thetaeqC[i*nseg+j] = thetaeq[fibID[i]*nseg+j];				
			phieqC[i*nseg+j] = phieq[fibID[i]*nseg+j];				
		}
	}
	
	free(thetaeq); free(phieq); 

	for (i = 0; i < nfibC_ID; i++){
		for (j = 1; j < nseg; j++){
			mi = i*nseg+j;
			R11eq[mi] = cosf(thetaeqC[i*nseg+j])*cosf(phieqC[i*nseg+j]);
			R12eq[mi] = cosf(thetaeqC[i*nseg+j])*sinf(phieqC[i*nseg+j]);
			R13eq[mi] = -sinf(thetaeqC[i*nseg+j]);
			R21eq[mi] = -sinf(phieqC[i*nseg+j]);
			R22eq[mi] = cosf(phieqC[i*nseg+j]);
			R23eq[mi] = 0.0;
			R31eq[mi] = sinf(thetaeqC[i*nseg+j])*cosf(phieqC[i*nseg+j]);
			R32eq[mi] = sinf(thetaeqC[i*nseg+j])*sinf(phieqC[i*nseg+j]);
			R33eq[mi] = cosf(thetaeqC[i*nseg+j]);
		}
	}
	free(thetaeqC); free(phieqC); 

	
	// reprint configuration files
	for (step = 0; step < nConfig; step++){
		
		// rc
		fread(&dum, sizeof(float), 1, center_mass);
		fwrite(&dum, sizeof(float), 1, center_massFloc);
		// rcmx 
		fread(rc, sizeof(float), nfib, center_mass);
		for (i = 0; i < nfibC_ID; i++){
			rcC[i] = rc[fibID[i]];
		}
		fwrite(rcC, sizeof(float), nfibC_ID, center_massFloc);
		// rcmy
		fread(rc, sizeof(float), nfib, center_mass);
		for (i = 0; i < nfibC_ID; i++){
			rcC[i] = rc[fibID[i]];
		}
		fwrite(rcC, sizeof(float), nfibC_ID, center_massFloc);
		// rcmz
		fread(rc, sizeof(float), nfib, center_mass);
		for (i = 0; i < nfibC_ID; i++){
			rcC[i] = rc[fibID[i]];
		}
		fwrite(rcC, sizeof(float), nfibC_ID, center_massFloc);
		// rx
		fread(&dum, sizeof(float), 1, rxfile);
		fread(r, sizeof(float), nfib*nseg, rxfile);
		fwrite(&dum, sizeof(float), 1, rxfileFloc); 
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				rxC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(rxC, sizeof(float), nfibC_ID*nseg, rxfileFloc);
		// ry
		fread(&dum, sizeof(float), 1, ryfile);
		fread(r, sizeof(float), nfib*nseg, ryfile);
		fwrite(&dum, sizeof(float), 1, ryfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				ryC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(ryC, sizeof(float), nfibC_ID*nseg, ryfileFloc);
		// rz
		fread(&dum, sizeof(float), 1, rzfile);
		fread(r, sizeof(float), nfib*nseg, rzfile);
		fwrite(&dum, sizeof(float), 1, rzfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				rzC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(rzC, sizeof(float), nfibC_ID*nseg, rzfileFloc);
		// px 
		fread(&dum, sizeof(float), 1, pxfile);
		fread(r, sizeof(float), nfib*nseg, pxfile);
		fwrite(&dum, sizeof(float), 1, pxfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				pxC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(pxC, sizeof(float), nfibC_ID*nseg, pxfileFloc);
		// py
		fread(&dum, sizeof(float), 1, pyfile);
		fread(r, sizeof(float), nfib*nseg, pyfile);
		fwrite(&dum, sizeof(float), 1, pyfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				pyC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(pyC, sizeof(float), nfibC_ID*nseg, pyfileFloc);
		// pz
		fread(&dum, sizeof(float), 1, pzfile);
		fread(r, sizeof(float), nfib*nseg, pzfile);
		fwrite(&dum, sizeof(float), 1, pzfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				pzC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(pzC, sizeof(float), nfibC_ID*nseg, pzfileFloc);
		// ux 
		fread(&dum, sizeof(float), 1, uxfile);
		fread(r, sizeof(float), nfib*nseg, uxfile);
		fwrite(&dum, sizeof(float), 1, uxfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				rC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(rC, sizeof(float), nfibC_ID*nseg, uxfileFloc);
		// uy 
		fread(&dum, sizeof(float), 1, uyfile);
		fread(r, sizeof(float), nfib*nseg, uyfile);
		fwrite(&dum, sizeof(float), 1, uyfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				rC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(rC, sizeof(float), nfibC_ID*nseg, uyfileFloc);
		// uy 
		fread(&dum, sizeof(float), 1, uzfile);
		fread(r, sizeof(float), nfib*nseg, uzfile);
		fwrite(&dum, sizeof(float), 1, uzfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				rC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(rC, sizeof(float), nfibC_ID*nseg, uzfileFloc);
		// wx 
		fread(&dum, sizeof(float), 1, wxfile);
		fread(r, sizeof(float), nfib*nseg, wxfile);
		fwrite(&dum, sizeof(float), 1, wxfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				rC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(rC, sizeof(float), nfibC_ID*nseg, wxfileFloc);
		// wy
		fread(&dum, sizeof(float), 1, wyfile);
		fread(r, sizeof(float), nfib*nseg, wyfile);
		fwrite(&dum, sizeof(float), 1, wyfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				rC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(rC, sizeof(float), nfibC_ID*nseg, wyfileFloc);
		// wz
		fread(&dum, sizeof(float), 1, wzfile);
		fread(r, sizeof(float), nfib*nseg, wzfile);
		fwrite(&dum, sizeof(float), 1, wzfileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				rC[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(rC, sizeof(float), nfibC_ID*nseg, wzfileFloc);
		// q0
		fread(&dum, sizeof(float), 1, q0file);
		fread(r, sizeof(float), nfib*nseg, q0file);
		fwrite(&dum, sizeof(float), 1, q0fileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				q0C[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(q0C, sizeof(float), nfibC_ID*nseg, q0fileFloc);
		// q1
		fread(&dum, sizeof(float), 1, q1file);
		fread(r, sizeof(float), nfib*nseg, q1file);
		fwrite(&dum, sizeof(float), 1, q1fileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				q1C[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(q1C, sizeof(float), nfibC_ID*nseg, q1fileFloc);
		// q2
		fread(&dum, sizeof(float), 1, q2file);
		fread(r, sizeof(float), nfib*nseg, q2file);
		fwrite(&dum, sizeof(float), 1, q2fileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				q2C[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(q2C, sizeof(float), nfibC_ID*nseg, q2fileFloc);
		// q3
		fread(&dum, sizeof(float), 1, q3file);
		fread(r, sizeof(float), nfib*nseg, q3file);
		fwrite(&dum, sizeof(float), 1, q3fileFloc);
		for (i = 0; i < nfibC_ID; i++){
			for (j = 0; j < nseg; j++){
				q3C[i*nseg+j] = r[fibID[i]*nseg+j];				
			}
		}
		fwrite(q3C, sizeof(float), nfibC_ID*nseg, q3fileFloc);

		// calculate elastic energy
		Eelastic = elasticEnergy(pxC, pyC, pzC, 
					 rxC, ryC, rzC,
					 q0C, q1C, q2C, q3C,
					 R11, R12, R13, R21, R22, R23,
					 R11eq, R12eq, R13eq, R21eq, R22eq,
					 R23eq, R31eq, R32eq, R33eq,
					 kb, nfibC_ID, nseg);
		fprintf(flocElastic, "%10.5f %10.6f\n", float(step*config_write)*dt, Eelastic);	
	}

	free(r); free(rc); 
	free(rC); free(rcC); 
	free(rxC); free(ryC); free(rzC); 
	free(pxC); free(pyC); free(pzC); 
	free(q0C); free(q1C); free(q2C); free(q3C); 

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

	fclose(rxfile); fclose(ryfile); fclose(rzfile);
	fclose(rxfileFloc); fclose(ryfileFloc); fclose(rzfileFloc);

	fclose(pxfile); fclose(pyfile); fclose(pzfile);
	fclose(pxfileFloc); fclose(pyfileFloc); fclose(pzfileFloc);

	fclose(uxfile); fclose(uyfile); fclose(uzfile);
	fclose(uxfileFloc); fclose(uyfileFloc); fclose(uzfileFloc);

	fclose(wxfile); fclose(wyfile); fclose(wzfile);
	fclose(wxfileFloc); fclose(wyfileFloc); fclose(wzfileFloc);

	fclose(q0file); fclose(q1file); fclose(q2file); fclose(q3file);
	fclose(q0fileFloc); fclose(q1fileFloc); fclose(q2fileFloc); fclose(q3fileFloc);

	fclose(center_mass); fclose(center_massFloc); 

	fclose(flocElastic); 
	
	free(fibID);

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

	/*for (m = 0; m < nfibC_ID; m++){
		for (i = 0; i < nseg; i++){
			mi = m*nseg+i;
			printf("%d %d %10.6f %10.6f %10.6f\n", m, i, rx[mi], ry[mi], rz[mi]);
		}
	}*/

	for (m = 0; m < nfibC_ID; m++){
		for (i = 0; i < nseg; i++){
			mi = m*nseg+i;
			R11[mi] = 2.0*(q0[mi]*q0[mi] + q1[mi]*q1[mi]) - 1.0;
			R12[mi] = 2.0*(q1[mi]*q2[mi] + q0[mi]*q3[mi]);
			R13[mi] = 2.0*(q1[mi]*q3[mi] - q0[mi]*q2[mi]);
			R21[mi] = 2.0*(q1[mi]*q2[mi] - q0[mi]*q3[mi]);
			R22[mi] = 2.0*(q0[mi]*q0[mi] + q2[mi]*q2[mi]) - 1.0;
			R23[mi] = 2.0*(q3[mi]*q2[mi] + q0[mi]*q1[mi]);
		}
	}

	for (m = 0; m < nfibC_ID; m++){
		for (i = 1; i < nseg; i++){
			ang_theta = 0.0; 
			ang_phi = 0.0; 
			mi = m*nseg+i;
			// calculations for theta
			zeqx = R11[mi-1]*R31eq[mi]+R21[mi-1]*R32eq[mi]+px[mi-1]*R33eq[mi];		
			zeqy = R12[mi-1]*R31eq[mi]+R22[mi-1]*R32eq[mi]+py[mi-1]*R33eq[mi];		
			zeqz = R13[mi-1]*R31eq[mi]+R23[mi-1]*R32eq[mi]+pz[mi-1]*R33eq[mi];		
			zdotz = px[mi]*zeqx+py[mi]*zeqy+pz[mi]*zeqz;
		
			if (zdotz <= 1.0 - 1.0E-6){
				if (zdotz < -1.0){
					zdotz = -1.0;
				}
				ang_theta = acosf(zdotz);
			}
			
			// calculations for phi
			cx = rx[mi] - rx[mi-1];
			cy = ry[mi] - ry[mi-1];
			cz = rz[mi] - rz[mi-1];
			dum = sqrtf(cx*cx + cy*cy + cz*cz);
			cx = cx / dum;
			cy = cy / dum;
			cz = cz / dum;
			zeqx = R11[mi-1]*R21eq[mi]+R21[mi-1]*R22eq[mi]+px[mi-1]*R23eq[mi];		
			zeqy = R12[mi-1]*R21eq[mi]+R22[mi-1]*R22eq[mi]+py[mi-1]*R23eq[mi];		
			zeqz = R13[mi-1]*R21eq[mi]+R23[mi-1]*R22eq[mi]+pz[mi-1]*R23eq[mi];		
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

	
	return kb*Energy/float(nfibC_ID);
}

