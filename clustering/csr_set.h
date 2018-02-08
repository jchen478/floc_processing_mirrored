//
//  csr_set.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 6/21/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__csr_set__
#define __flexfric__csr_set__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void csr_set(int **intVar);

#endif /* defined(__flexfric__csr_set__) */
