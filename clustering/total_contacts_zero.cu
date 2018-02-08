#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include "total_contacts_zero.h"

using namespace std;

__global__ void total_contacts_zero(float **var, int **intVar){

	int *overs = intVar[4]; 
	int *total_contacts = intVar[5];
	*total_contacts = 0; 
	*overs = 0; 
}
