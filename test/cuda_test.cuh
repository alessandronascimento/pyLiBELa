#ifndef CUDA_TEST_CUH
#define CUDA_TEST_CUH

__global__
void saxpy(int n, float a, float *x, float *y);
void execution();
#endif
