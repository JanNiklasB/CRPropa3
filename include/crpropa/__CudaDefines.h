#ifndef CRPROPA__CUDADEFINES_H
#define CRPROPA__CUDADEFINES_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CONSTANT __constant__
#include <cuda_runtime.h>
#include <cuda.h>

#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CONSTANT
#endif

#endif