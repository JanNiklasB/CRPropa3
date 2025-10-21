#ifndef CRPROPA__CUDADEFINES_H
#define CRPROPA__CUDADEFINES_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CONSTANT __constant__
#include <cuda_runtime.h>
#include <cuda.h>
#include <cuda/std/cstdlib>
#include <iostream>

// see answer of talonmies on https://stackoverflow.com/a/14038590/30769038
inline void gpuAssert(cudaError_t code, char * file, int line, bool Abort=true){
	if (code != 0) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code),file,line);
		if (Abort) exit(code);
	}
}
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CONSTANT
#endif

#endif