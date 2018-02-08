//
//  clustering.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 6/5/17.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__clustering__
#define __flexfric__clustering__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void clustering(float **var, int **intVar);

#endif /* defined(__flexfric__clustering__) */