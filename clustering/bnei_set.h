//
//  bnei_set.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 11/28/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__bnei_set__
#define __flexfric__bnei_set__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void bnei_set(float **var, int **intVar);

#endif /* defined(__flexfric__bnei_set__) */