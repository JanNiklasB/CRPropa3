#include "crpropa/Common.h"

#include "kiss/path.h"
#include "kiss/logger.h"

#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>

#define index(i,j) ((j)+(i)*size)

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
		const double* Y, int size) {
	size_t i = upper_bound<double>(x, X, size);
	if (i == 0)
		return Y[0];
	if (i == size-1)
		return Y[size-1];
	return Y[i] + (x - X[i]) * (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
}

double interpolate2d(double x, double y, const std::vector<double> &X,
		const std::vector<double> &Y, const std::vector<double> &Z) {
	
	return interpolate2d(x, y, X.data(), Y.data(), Z.data(), Z.size());
}

double interpolate2d(double x, double y, const double* X,
		const double* Y, const double* Z, int size) {

	size_t i = upper_bound<double>(x, X, size);
	size_t j = upper_bound<double>(y, Y, size);

	if (x > X[size-1] || x < X[0])
		return 0;
	if (y > Y[size-1] || y < Y[0])
		return 0;

	if (X[i] == X[0] && Y[j] == Y[0])
		return Z[0];
	if (X[i] == X[size-1] && Y[j] == Y[size-1])
		return Z[size-1];


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
		const double* Y, int size) {
	if (x <= lo)
		return Y[0];
	if (x >= hi)
		return Y[size-1];

	double dx = (hi - lo) / (size - 1);
	double p = (x - lo) / dx;
	size_t i = floor(p);
	return Y[i] + (p - i) * (Y[i + 1] - Y[i]);
}

size_t closestIndex(double x, const std::vector<double> &X) {
	return closestIndex(x, X.data(), X.size());
}

size_t closestIndex(double x, const double *X, int size) {
	size_t i1 = lower_bound<double>(x, X, size);
	if (i1 == 0)
		return i1;
	size_t i0 = i1 - 1;
	if (crstd::fabs(X[i0] - x) < crstd::fabs(X[i1] - x))
		return i0;
	else
		return i1;
}

} // namespace crpropa

