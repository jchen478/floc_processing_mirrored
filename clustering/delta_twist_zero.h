//
//  delta_twist_zero.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 7/16/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__delta_twist_zero__
#define __flexfric__delta_twist_zero__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void delta_twist_zero(float **var, int **intVar);

#endif /* defined(__flexfric__delta_twist_zero__) */
