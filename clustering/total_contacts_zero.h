//
//  total_contacts_zero.h
//  flexfric_cuda_cu
//
//  Created by Jing-Yao Chen on 7/19/16.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#ifndef __flexfric__total_contacts_zero__
#define __flexfric__total_contacts_zero__

#include <stdio.h>
#include <cuda_runtime_api.h>

__global__ void total_contacts_zero(float **var, int **intVar);

#endif /* defined(__flexfric__total_contacts_zero__) */
