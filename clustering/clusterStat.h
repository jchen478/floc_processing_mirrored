//
//  clusterStat.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 6/16/17.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__clusterStat__
#define __flexfric__clusterStat__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void clusterStat(float **var, int **intVar);

#endif /* defined(__flexfric__clusterStat__) */