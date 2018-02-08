//
//  initialize.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 7/12/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__initialize__
#define __flexfric__initialize__

#include <stdio.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

__global__ void initialize(float **var, int **intVar);

#endif /* defined(__flexfric__initialize__) */
