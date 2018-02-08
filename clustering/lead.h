//
//  lead.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 11/30/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__lead__
#define __flexfric__lead__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void lead(float **var, int **intVar);

#endif /* defined(__flexfric__lead__) */