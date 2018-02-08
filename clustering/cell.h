//
//  cell.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 11/23/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__cell__
#define __flexfric__cell__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void cell(float **var, int **intVar);

#endif /* defined(__flexfric__cell__) */