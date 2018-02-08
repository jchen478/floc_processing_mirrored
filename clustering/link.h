//
//  link.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 11/27/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__link__
#define __flexfric__link__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void link(float **var, int **intVar);

#endif /* defined(__flexfric__link__) */