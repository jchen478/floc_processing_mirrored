//
//  initialize_regrow2.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 7/29/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__initialize_regrow2__
#define __flexfric__initialize_regrow2__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void initialize_regrow2(float **var, int **intVar);

#endif /* defined(__flexfric__initialize_regrow2__) */