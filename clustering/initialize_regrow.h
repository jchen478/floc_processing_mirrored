//
//  initialize_regrow.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 7/27/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__initialize_regrow__
#define __flexfric__initialize_regrow__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void initialize_regrow(float **var, int **intVar);

#endif /* defined(__flexfric__initialize_regrow__) */