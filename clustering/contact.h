//
//  contact.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 3/9/17.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__contact__
#define __flexfric__contact__

#include <stdio.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

__global__ void contact(float **var, int **intVar);

#endif /* defined(__flexfric__contact__) */