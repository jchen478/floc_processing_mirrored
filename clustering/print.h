//
//  print.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 6/16/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__print__
#define __flexfric__print__

#include <stdio.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

__global__ void print(float **var, int **intVar);

#endif /* defined(__flexfric__print__) */