#include "crpropa/Common.h"

#include "kiss/path.h"
#include "kiss/logger.h"

#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>

#define index(i,j) ((j)+(i)*Y.size())

CUDA_CALLABLE_MEMBER size_t lower_bound(double x, const double *X, int size) {
	int count=size-1, step;
	double it;
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

CUDA_CALLABLE_MEMBER size_t upper_bound(double x, const double *X, int size){
	int count=size-1, step;
	double it;
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

namespace crpropa {

std::string getDataPath(std::string filename) {
	static std::string dataPath;
	if (dataPath.size())
		return concat_path(dataPath, filename);

	const char *env_path = getenv("CRPROPA_DATA_PATH");
	if (env_path) {
		if (is_directory(env_path)) {
			dataPath = env_path;
			KISS_LOG_INFO << "getDataPath: use environment variable, "
					<< dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}

#ifdef CRPROPA_INSTALL_PREFIX
	{
		std::string _path = CRPROPA_INSTALL_PREFIX "/share/crpropa";
		if (is_directory(_path)) {
			dataPath = _path;
			KISS_LOG_INFO
			<< "getDataPath: use install prefix, " << dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}
#endif

	{
		std::string _path = executable_path() + "../data";
		if (is_directory(_path)) {
			dataPath = _path;
			KISS_LOG_INFO << "getDataPath: use executable path, " << dataPath
					<< std::endl;
			return concat_path(dataPath, filename);
		}
	}

#ifdef CRPROPA_BINARY_DIR
	{
		std::string _path = CRPROPA_BINARY_DIR "/data";
		if (is_directory(_path)) {
			dataPath = _path;
			KISS_LOG_INFO
			<< "getDataPath: use binary path, " << dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}
#endif

	dataPath = "data";
	KISS_LOG_INFO << "getDataPath: use default, " << dataPath << std::endl;
	return concat_path(dataPath, filename);
}


std::string getInstallPrefix()
{
  std::string _path = "";
  #ifdef CRPROPA_INSTALL_PREFIX
    _path += CRPROPA_INSTALL_PREFIX;
  #endif
  return _path;
};

double interpolate(double x, const std::vector<double> &X,
		const std::vector<double> &Y) {
	return interpolate(x, X.data(), Y.data(), X.size());
}

double interpolate(double x, const double* X,
		const double* Y, int Size) {
	size_t i = upper_bound(x, X, Size);
	if (i == 0)
		return Y[0];
	if (i == Size-1)
		return Y[Size-1];
	return Y[i] + (x - X[i]) * (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
}

double interpolate2d(double x, double y, const std::vector<double> &X,
		const std::vector<double> &Y, const std::vector<double> &Z) {

	std::vector<double>::const_iterator itx = std::upper_bound(X.begin(), X.end(), x);
	std::vector<double>::const_iterator ity = std::upper_bound(Y.begin(), Y.end(), y);

	if (x > X.back() || x < X.front())
		return 0;
	if (y > Y.back() || y < Y.front())
		return 0;

	if (itx == X.begin() && ity == Y.begin())
		return Z.front();
	if (itx == X.end() && ity == Y.end())
		return Z.back();

	size_t i = itx - X.begin() - 1;
	size_t j = ity - Y.begin() - 1;

	double Q11 = Z[index(i,j)];
	double Q12 = Z[index(i,j+1)];
	double Q21 = Z[index(i+1,j)];
	double Q22 = Z[index(i+1,j+1)];

	double R1 = ((X[i+1]-x)/(X[i+1]-X[i]))*Q11+((x-X[i])/(X[i+1]-X[i]))*Q21;
	double R2 = ((X[i+1]-x)/(X[i+1]-X[i]))*Q12+((x-X[i])/(X[i+1]-X[i]))*Q22;

	return ((Y[j+1]-y)/(Y[j+1]-Y[j]))*R1+((y-Y[j])/(Y[j+1]-Y[j]))*R2;
}

double interpolateEquidistant(double x, double lo, double hi,
		const std::vector<double> &Y) {
	return interpolateEquidistant(x, lo, hi, Y.data(), Y.size());
}

double interpolateEquidistant(double x, double lo, double hi,
		const double* Y, int YSize) {
	if (x <= lo)
		return Y[0];
	if (x >= hi)
		return Y[YSize-1];

	double dx = (hi - lo) / (YSize - 1);
	double p = (x - lo) / dx;
	size_t i = floor(p);
	return Y[i] + (p - i) * (Y[i + 1] - Y[i]);
}

size_t closestIndex(double x, const std::vector<double> &X) {
	return closestIndex(x, X.data(), X.size());
}

size_t closestIndex(double x, const double *X, int size) {
	size_t i1 = lower_bound(x, X, size);
	if (i1 == 0)
		return i1;
	size_t i0 = i1 - 1;
	if (std::fabs(X[i0] - x) < std::fabs(X[i1] - x))
		return i0;
	else
		return i1;
}

} // namespace crpropa

