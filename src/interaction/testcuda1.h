#ifndef TESTCUDA1_H_
#define TESTCUDA1_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <iostream>
#include <cuda.h>
#define BLOCK_SIZE 256
#define SOFTENING 1e-9f



void testcuda_call();

#endif
