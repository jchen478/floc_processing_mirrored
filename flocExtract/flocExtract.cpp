#include <stdio.h>
#include <stdlib.h>
using namespace std; 

int main(){

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////
	FILE *Parameters, *Centers_of_Mass, *Euler_Parameters;
	FILE *Cluster_fibID, *ClusterInfo; 
	FILE *center_mass;
	FILE *q0file, *q1file, *q2file, *q3file;

	Parameters = fopen("../Parameters.in", "r"); 
	Cluster_fibID = fopen("../Cluster_fibID.txt","r");
	ClusterInfo = fopen("ClusterInfo.in", "r"); 
	Centers_of_Mass = fopen("Centers_of_Mass.in", "w"); 
	Euler_Parameters = fopen("Euler_Parameters.in", "w"); 
	center_mass = fopen("../../../center_mass.txt", "rb");
	q0file = fopen("../../../q0.txt", "rb");
	q1file = fopen("../../../q1.txt", "rb");
	q2file = fopen("../../../q2.txt", "rb");
	q3file = fopen("../../../q3.txt", "rb");

	int cID; 
	fscanf(ClusterInfo, "%d ",&cID);

	int nfib, nseg, config_write; 
	float dt, strain, dum;
	
	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
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

	int m, i, mi, n, idum, nStep; 	
	nStep = int((strain/dt)/float(config_write)) + 1;

	int *fibID;
	float *rcmx, *rcmy, *rcmz;
	float *q0, *q1, *q2, *q3;

	fibID = (int*)malloc(nfib*sizeof(int));
	rcmx = (float*)malloc(nfib*sizeof(float));
	rcmy = (float*)malloc(nfib*sizeof(float));
	rcmz = (float*)malloc(nfib*sizeof(float));
	q0 = (float*)malloc(nfib*nseg*sizeof(float));
	q1 = (float*)malloc(nfib*nseg*sizeof(float));
	q2 = (float*)malloc(nfib*nseg*sizeof(float));
	q3 = (float*)malloc(nfib*nseg*sizeof(float));

	for (m = 0; m < nfib; m++){
		fscanf(Cluster_fibID, " %d %d", &idum, fibID + m);
	}

	for (n = 0; n < nStep; n++){
		fread(&dum, sizeof(float), 1, center_mass);
		fread(&dum, sizeof(float), 1, q0file);
		fread(&dum, sizeof(float), 1, q1file);
		fread(&dum, sizeof(float), 1, q2file);
		fread(&dum, sizeof(float), 1, q3file);
		fread(rcmx, sizeof(float), nfib, center_mass);
		fread(rcmy, sizeof(float), nfib, center_mass);
		fread(rcmz, sizeof(float), nfib, center_mass);
		fread(q0, sizeof(float), nfib*nseg, q0file);
		fread(q1, sizeof(float), nfib*nseg, q1file);
		fread(q2, sizeof(float), nfib*nseg, q2file);
		fread(q3, sizeof(float), nfib*nseg, q3file);
	}	

	for (m = 0; m < nfib; m++){
		if (fibID[m] != cID)
			continue;
		fprintf(Centers_of_Mass, "%4d%17.8E%17.8E%17.8E\n",
			m + 1, rcmx[m], rcmy[m], rcmz[m]);
		for (i = 0; i < nseg; i++){
			mi = m*nseg + i;
			fprintf(Euler_Parameters, "%6d%4d%17.8E%17.8E%17.8E%17.8E\n",
				m + 1, i + 1, q0[mi], q1[mi], q2[mi], q3[mi]);
		}
	} 

	free(fibID); free(rcmx); free(rcmy); free(rcmz);
	free(q0); free(q1); free(q2); free(q3);

	fclose(Centers_of_Mass); fclose(Euler_Parameters);
	fclose(ClusterInfo); fclose(Cluster_fibID);  

	return 0;
}
