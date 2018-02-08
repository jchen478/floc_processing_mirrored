//
//  clustering.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 6/19/17.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__clusteringZero__
#define __flexfric__clusteringZero__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void clusteringZero(float **var, int **intVar);

#endif /* defined(__flexfric__clustering__) */
