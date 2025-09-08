#include "crpropa/Cosmology.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include <cmath>
#include <stdexcept>

namespace crpropa {

void Cosmology::update() {
	double dH = c_light / H0; // Hubble distance

	std::vector<double> E(n);
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
}

Cosmology::Cosmology() {
	// Cosmological parameters (K.A. Olive et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014))
	H0 = 67.3 * 1000 * meter / second / Mpc; // default values
	omegaM = 0.315;
	omegaL = 1 - omegaM;

	Z.resize(n);
	Dc.resize(n);
	Dl.resize(n);
	Dt.resize(n);

	Z[0] = 0;
	Dc[0] = 0;
	Dl[0] = 0;
	Dt[0] = 0;

	ZPtr = Z.data();
	DcPtr = Dc.data();
	DlPtr = Dl.data();
	DtPtr = Dt.data();

	ZSize = Z.size();
	DcSize = Dc.size();
	DlSize = Dl.size();
	DtSize = Dt.size();

	update();
}

void Cosmology::setParameters(double hubbleParameter, double omegaMatter) {
	H0 = hubbleParameter * 1e5 / Mpc;
	omegaM = omegaMatter;
	omegaL = 1 - omegaMatter;
	update();
}

double Cosmology::hubbleRate(double z) {
	return H0 * sqrt(omegaL + omegaM * pow(1 + z, 3));
}

double Cosmology::getOmegaL() {
	return omegaL;
}

double Cosmology::getOmegaM() {
	return omegaM;
}

double Cosmology::getH0() {
	return H0;
}

double Cosmology::comovingDistance2Redshift(double d) {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DcPtr[this->DcSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->DcPtr, this->ZPtr, this->ZSize);
}

double Cosmology::redshift2ComovingDistance(double z) {
	#ifndef __CUDACC__
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > this->zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	#endif
	return interpolate(z, this->ZPtr, this->DcPtr, this->DcSize);
}

double Cosmology::luminosityDistance2Redshift(double d) {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DlPtr[this->DlSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->DlPtr, this->ZPtr, this->ZSize);
}

double Cosmology::redshift2LuminosityDistance(double z) {
	#ifndef __CUDACC__
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > this->zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	#endif
	return interpolate(z, this->ZPtr, this->DlPtr, this->DlSize);
}

double Cosmology::lightTravelDistance2Redshift(double d) {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DtPtr[this->DtSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->DtPtr, this->ZPtr, this->ZSize);
}

double Cosmology::redshift2LightTravelDistance(double z) {
	#ifndef __CUDACC__
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > this->zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	#endif
	return interpolate(z, this->ZPtr, this->DtPtr, this->DtSize);
}

double Cosmology::comoving2LightTravelDistance(double d) {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DcPtr[this->DcSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->DcPtr, this->DtPtr, this->DtSize);
}

double Cosmology::lightTravel2ComovingDistance(double d) {
	#ifndef __CUDACC__
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > this->DtPtr[this->DtSize-1])
		throw std::runtime_error("Cosmology: d > dmax");
	#endif
	return interpolate(d, this->DtPtr, this->DcPtr, this->DcSize);
}

static Cosmology cosmology;

ref_ptr<Cosmology> getStaticCosmology(){
	return ref_ptr<Cosmology>(&cosmology);
}

void setCosmologyParameters(double hubbleParameter, double omegaMatter){
	cosmology.setParameters(hubbleParameter, omegaMatter);
}

double hubbleRate(double redshift){
	return cosmology.hubbleRate(redshift);
}

double omegaL(){
	return cosmology.getOmegaL();
}

double omegaM(){
	return cosmology.getOmegaM();
}

double H0(){
	return cosmology.getH0();
}

double comovingDistance2Redshift(double distance){
	return cosmology.comovingDistance2Redshift(distance);
}

double redshift2ComovingDistance(double redshift){
	return cosmology.redshift2ComovingDistance(redshift);
}

double luminosityDistance2Redshift(double distance){
	return cosmology.luminosityDistance2Redshift(distance);
}

double redshift2LuminosityDistance(double redshift){
	return cosmology.redshift2LuminosityDistance(redshift);
}

double lightTravelDistance2Redshift(double distance){
	return cosmology.lightTravelDistance2Redshift(distance);
}

double redshift2LightTravelDistance(double redshift){
	return cosmology.redshift2LightTravelDistance(redshift);
}

double comoving2LightTravelDistance(double distance){
	return cosmology.comoving2LightTravelDistance(distance);
}

double lightTravel2ComovingDistance(double distance){
	return cosmology.lightTravel2ComovingDistance(distance);
}

} // namespace crpropa
