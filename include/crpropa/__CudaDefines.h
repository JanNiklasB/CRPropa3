#ifndef CRPROPA__CUDADEFINES_H
#define CRPROPA__CUDADEFINES_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CONSTANT __constant__
#include <cuda_runtime.h>
#include <cuda.h>
// define variable namespace to cuda::std, 
// only needed in device functions, std still usable
namespace crpropa{
	namespace crstd = cuda::std;
}
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CONSTANT
// define variable namespace to std
namespace crpropa{
	namespace crstd = std;
}
#endif

#endif