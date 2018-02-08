#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){

	// Relevant files
	FILE *ClusterInfo, *Cluster_results;
	FILE *Sdist;
	ClusterInfo = fopen("ClusterInfo.in", "r");
	Cluster_results = fopen("Cluster_results.txt", "r");
	Sdist = fopen("Sdist.txt", "w");

	// Variables in input files
	int nfib, idum;
	int nFloc;
	float strain, dum;

	// Read in cluster info
	fscanf(ClusterInfo, "%d", &idum);
	fscanf(ClusterInfo, "%*[^\n]%d", &nFloc);

	// Read in number of fibers in each floc
	int *nfibFloc, *fiberId, i;
	nfibFloc = (int*)malloc(nFloc*sizeof(int)); 
	fiberId = (int*)malloc(nfib*sizeof(int));
	for (i = 0; i < nFloc; i++){
		fscanf(Cluster_results, "%d", &idum);
		fscanf(Cluster_results, " %d", nfibFloc + i);
	}

	// close files
	fclose(ClusterInfo); fclose(Cluster_results);
	
	int *Sbin, max; 
	Sbin = (int*)calloc(10000,sizeof(int));

	max = 0;
	for (i = 0; i < nFloc; i++){
		Sbin[nfibFloc[i]]++; 
		if (nfibFloc[i] > max)
			max = nfibFloc[i];
		if (max >= 10000){
			fprintf(Sdist,"Number of contacts exceed 10000, abort.\n");
			fclose(Sdist);
			free(nfibFloc); free(fiberId); free(Sbin);
			return 0; 
		}
	}
	for (i = 2; i <= max; i++){
		fprintf(Sdist,"%4d %4d\n", i, Sbin[i]);
	}
		
	// free memory
	free(nfibFloc); free(fiberId); free(Sbin);
	fclose(Sdist); 

	return 0;
}
