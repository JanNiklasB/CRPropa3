#include "crpropa/Cosmology.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include <cmath>
#include <stdexcept>

namespace crpropa {

void Cosmology::update() {
	double dH = c_light / H0; // Hubble distance

	double* E = new double[n];
	E[0] = 1;

	// Relation between comoving distance r and redshift z (cf. J.A. Peacock, Cosmological physics, p. 89 eq. 3.76)
	// dr = c / H(z) dz, integration using midpoint rule
	double dlz = log10(zmax) - log10(zmin);
	for (int i = 1; i < n; i++) {
		Z[i] = zmin * pow(10, i * dlz / (n - 1)); // logarithmic even spacing
		double dz = (Z[i] - Z[i - 1]); // redshift step
		E[i] = sqrt(omegaL + omegaM * pow_integer<3>(1 + Z[i]));
		Dc[i] = Dc[i - 1] + dH * dz * (1 / E[i] + 1 / E[i - 1]) / 2;
		Dl[i] = (1 + Z[i]) * Dc[i];
		Dt[i] = Dt[i - 1]
				+ dH * dz
						* (1 / ((1 + Z[i]) * E[i])
								+ 1 / ((1 + Z[i - 1]) * E[i - 1])) / 2;
	}
	delete[] E;
}

Cosmology::Cosmology() {
	// Cosmological parameters (K.A. Olive et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014))
	H0 = 67.3 * 1000 * meter / second / Mpc; // default values
	omegaM = 0.315;
	omegaL = 1 - omegaM;

	 Z = new double[n];
	Dc = new double[n];
	Dl = new double[n];
	Dt = new double[n];

	Z[0] = 0;
	Dc[0] = 0;
	Dl[0] = 0;
	Dt[0] = 0;

	ZSize = n;
	DcSize = n;
	DlSize = n;
	DtSize = n;

	update();
}

Cosmology::~Cosmology(){
	delete[] Z ;
	delete[] Dc;
	delete[] Dl;
	delete[] Dt;
}

void Cosmology::setParameters(double hubbleParameter, double omegaMatter) {
	H0 = hubbleParameter * 1e5 / Mpc;
	omegaM = omegaMatter;
	omegaL = 1 - omegaMatter;
	update();
}

double Cosmology::hubbleRate(double z) const {
	return H0 * sqrt(omegaL + omegaM * pow(1 + z, 3));
}

double Cosmology::getOmegaL() const {
	return omegaL;
}

double Cosmology::getOmegaM() const {
	return omegaM;
}

double Cosmology::getH0() const {
	return H0;
}

double Cosmology::comovingDistance2Redshift(double d) const {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DcPtr[this->DcSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->Dc, this->Z, this->ZSize);
}

double Cosmology::redshift2ComovingDistance(double z) const {
	#ifndef __CUDACC__
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > this->zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	#endif
	return interpolate(z, this->Z, this->Dc, this->DcSize);
}

double Cosmology::luminosityDistance2Redshift(double d) const {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DlPtr[this->DlSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->Dl, this->Z, this->ZSize);
}

double Cosmology::redshift2LuminosityDistance(double z) const {
	#ifndef __CUDACC__
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > this->zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	#endif
	return interpolate(z, this->Z, this->Dl, this->DlSize);
}

double Cosmology::lightTravelDistance2Redshift(double d) const {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DtPtr[this->DtSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->Dt, this->Z, this->ZSize);
}

double Cosmology::redshift2LightTravelDistance(double z) const {
	#ifndef __CUDACC__
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > this->zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	#endif
	return interpolate(z, this->Z, this->Dt, this->DtSize);
}

double Cosmology::comoving2LightTravelDistance(double d) const {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DcPtr[this->DcSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->Dc, this->Dt, this->DtSize);
}

double Cosmology::lightTravel2ComovingDistance(double d) const {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DtPtr[this->DtSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->Dt, this->Dc, this->DcSize);
}

} // namespace crpropa
