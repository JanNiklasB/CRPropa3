#ifndef CRPROPA_COMMON_H
#define CRPROPA_COMMON_H

#include "crpropa/__CudaDefines.h"

#include <string>
#include <vector>

/**
 @file
 @brief Common helper functions
 */

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/// Returns the full path to a CRPropa data file
std::string getDataPath(std::string filename);

/// Returns the install prefix
std::string getInstallPrefix();

/// Returns a certain digit from a given integer
CUDA_CALLABLE_MEMBER inline int digit(const int& value, const int& d) {
	return (value % (d * 10)) / d;
}

template<typename T>
CUDA_CALLABLE_MEMBER T* push_back(T* src, std::size_t& size, T value){
	T* tmp = new T[size+1];
	for (std::size_t i=0; i<size; i++)
		tmp[i] = src[i];
	tmp[size] = value;
	delete[] src;
	src = tmp;
	size++;
	return src;
}

template<typename T>
CUDA_CALLABLE_MEMBER T* insert(T* src, std::size_t& size, T value, std::size_t idx){
	T* tmp = new T[size+1];

	bool inserted=false;
	for (std::size_t i=0; i<size+1; i++){
		if(i==idx){
			tmp[i] = value;
			inserted = true;
			continue;
		}
		if(inserted) tmp[i] = src[i-1];
		else tmp[i] = src[i];
	}
	delete[] src;
	src = tmp;
	size++;
	return src;
}

template<typename T>
CUDA_CALLABLE_MEMBER T* erase(T* src, std::size_t& size, std::size_t idx){
	T* tmp = new T[size-1];

	bool removed=false;
	for (std::size_t i=0; i<size; i++){
		if(i==idx){
			removed = true;
			continue;
		}
		if(removed) tmp[i-1] = src[i];
		else tmp[i] = src[i];
	}
	delete[] src;
	src = tmp;
	size--;
	return src;
}

/// Same Behaviour as std::lower_bound, but accepts the array directly and returns the index instead
template<typename T>
CUDA_CALLABLE_MEMBER size_t lower_bound(T x, const T *X, std::size_t size) {
	std::size_t count=size-1, step;
	T it;
	size_t i1 = 0;
	while (count>0) {
		step = count/2;
		it = X[step+i1];
		if (it<x) {
			i1 = step+1;
			count -= step+1;
		}
		else count = step;
	}
	return i1;
}

/// Same Behaviour as std::upper_bound, but accepts the array directly and returns the index instead
template<typename T>
CUDA_CALLABLE_MEMBER size_t upper_bound(T x, const T *X, std::size_t size) {
	std::size_t count=size-1, step;
	T it;
	size_t i1 = 0;
	while (count>0) {
		step = count/2;
		it = X[step+i1];
		if (!(x<it)) {
			i1 = step+1;
			count -= step+1;
		}
		else count = step;
	}
	return i1;
}

template<typename T>
CUDA_CALLABLE_MEMBER size_t _quickfind(T x, const T *arr, std::size_t left, std::size_t right) {
	if(x==arr[right-1]) return right-1;
	if(x==arr[left]) return left;
	int middle = left + (right-left)/2;
	if(x<arr[middle]) return _quickfind(x, arr, left, middle);
	if(x>arr[middle]) return _quickfind(x, arr, middle+1, right);
	return middle;
}

/// Find index of a sorted array (sort with quicksort)
template<typename T>
CUDA_CALLABLE_MEMBER size_t findSorted(T x, const T *arr, std::size_t size) {
	return _quickfind(x, arr, 0, size);
}

template<typename T>
CUDA_CALLABLE_MEMBER void _swap(T *arr, std::size_t index1, std::size_t index2){
	T tmp = arr[index2];
	arr[index2] = arr[index1];
	arr[index1] = tmp;
}

template<typename T>
CUDA_CALLABLE_MEMBER int _divide(T *arr, std::size_t left, std::size_t right){
	std::size_t i = left;
	std::size_t j = right-1;
	T pivot = arr[right];

	while(i<j){
		while(i<j && arr[i]<=pivot) i++;
		while(j>i && arr[i]>pivot) j--;
		if(arr[i]>arr[j]) _swap(arr, i, j);
	}

	if(arr[i]>pivot) _swap(arr, i, right);
	else i = right;
	return i;
}
template<typename T>
CUDA_CALLABLE_MEMBER void _quicksort(T *arr, std::size_t left, std::size_t right){
	if(left < right){
		std::size_t divider = _divide(arr, left, right);
		_quicksort(arr, left, divider-1);
		_quicksort(arr, divider+1, right);
	}
}

template<typename T>
CUDA_CALLABLE_MEMBER void sort(T *arr, std::size_t size){
	quicksort(arr, 0, size);
}

// Return value xclip which is the closest to x, so that lower <= xclip <= upper
template <typename T>
CUDA_CALLABLE_MEMBER T clip(const T& x, const T& lower, const T& upper) {
	return std::max(lower, std::min(x, upper));
}

/// Perform linear interpolation on a set of n tabulated data points X[0 .. n-1] -> Y[0 .. n-1]
/// Returns Y[0] if x < X[0] and Y[n-1] if x > X[n-1]
double interpolate(double x, const std::vector<double>& X,
	const std::vector<double>& Y);

/// Perform linear interpolation on a set of n tabulated data points X[0 .. n-1] -> Y[0 .. n-1]
/// Returns Y[0] if x < X[0] and Y[n-1] if x > X[n-1]
CUDA_CALLABLE_MEMBER double interpolate(double x, const double* X,
	const double* Y, int size);


/// Perform bilinear interpolation on a set of (n,m) tabulated data points X[0 .. n-1], Y[0 .. m-1] -> Z[0.. n-1*m-1]
/// Returns 0 if x < X[0] or x > X[n-1] or y < Y[0] or y > Y[m-1]
double interpolate2d(double x, double y, const std::vector<double>& X,
		const std::vector<double>& Y, const std::vector<double>& Z);

/// Perform bilinear interpolation on a set of (n,m) tabulated data points X[0 .. n-1], Y[0 .. m-1] -> Z[0.. n-1*m-1]
/// Returns 0 if x < X[0] or x > X[n-1] or y < Y[0] or y > Y[m-1]
CUDA_CALLABLE_MEMBER double interpolate2d(double x, double y, const double* X,
		const double* Y, const double* Z, int size);

/// Perform linear interpolation on equidistant tabulated data
/// Returns Y[0] if x < lo and Y[n-1] if x > hi
double interpolateEquidistant(double x, double lo, double hi,
		const std::vector<double>& Y);

CUDA_CALLABLE_MEMBER double interpolateEquidistant(double x, double lo, double hi,
		const double* Y, int size);


/// Find index of value in a sorted vector X that is closest to x
size_t closestIndex(double x, const std::vector<double> &X);
/// Find index of value in a sorted array X that is closest to x
CUDA_CALLABLE_MEMBER size_t closestIndex(double x, const double* X, int size);

/** @}*/


/// pow implementation as template for integer exponents pow_integer<2>(x)
/// evaluates to x*x
template <unsigned int exponent>
CUDA_CALLABLE_MEMBER inline double pow_integer(double base)
{
  return pow_integer<(exponent >> 1)>(base*base) * (((exponent & 1) > 0) ? base : 1);
}

template <>
CUDA_CALLABLE_MEMBER inline double pow_integer<0>(double base)
{
  return 1;
}

/// - input:  function over which to integrate, integration limits A and B
/// - output: 8-points Gauß-Legendre integral
template<typename Integrand>
CUDA_CALLABLE_MEMBER double gaussInt(Integrand&& integrand, double A, double B) {
	const double X[8] = {.0950125098, .2816035507, .4580167776, .6178762444, .7554044083, .8656312023, .9445750230, .9894009349};
	const double W[8] = {.1894506104, .1826034150, .1691565193, .1495959888, .1246289712, .0951585116, .0622535239, .0271524594};
	const double XM = 0.5 * (B + A);
	const double XR = 0.5 * (B - A);
	double SS = 0.;
	for (int i = 0; i < 8; ++i) {
		double DX = XR * X[i];
		SS += W[i] * (integrand(XM + DX) + integrand(XM - DX));
	}
	return XR * SS;
}

} // namespace crpropa

#endif // CRPROPA_COMMON_H
