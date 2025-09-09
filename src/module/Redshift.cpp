#include "crpropa/module/Redshift.h"
#include "crpropa/Units.h"
#include "crpropa/Cosmology.h"

#include <limits>

namespace crpropa {

Redshift::Redshift(){
	#ifdef __CUDACC__
	cosmology = new Cosmology;
	#else
	cosmology = getStaticCosmology();
	#endif
}

ref_ptr<Cosmology> Redshift::getCosmology(){
	return cosmology;
}

ref_ptr<Cosmology> Redshift::setCosmology(ref_ptr<Cosmology> cosmology){
	this->cosmology = cosmology;
}

void Redshift::process(Candidate *c) const {
	double z = c->getRedshift();

	// check if z = 0
	if (z <= std::numeric_limits<double>::min())
		return;

	// use small step approximation:  dz = H(z) / c * ds
	double dz = cosmology->hubbleRate(z) / c_light * c->getCurrentStep();

	// prevent dz > z
	dz = std::min(dz, z);

	// update redshift
	c->setRedshift(z - dz);

	// adiabatic energy loss: dE / dz = E / (1 + z)
	double E = c->current.getEnergy();
	c->current.setEnergy(E * (1 - dz / (1 + z)));
}

std::string Redshift::getDescription() const {
	std::stringstream s;
	s << "Redshift: h0 = " << cosmology->hubbleRate() / 1e5 * Mpc << ", omegaL = "
			<< cosmology->getOmegaL() << ", omegaM = " << cosmology->getOmegaM();
	return s.str();
}

FutureRedshift::FutureRedshift(){
	#ifdef __CUDACC__
	cosmology = new Cosmology;
	#else
	cosmology = getStaticCosmology();
	#endif
}

ref_ptr<Cosmology> FutureRedshift::getCosmology(){
	return cosmology;
}

ref_ptr<Cosmology> FutureRedshift::setCosmology(ref_ptr<Cosmology> cosmology){
	this->cosmology = cosmology;
}

void FutureRedshift::process(Candidate *c) const {
	double z = c->getRedshift();

	// check if z = -1
	if (z <= -1)
		return;

	// use small step approximation:  dz = H(z) / c * ds
	double dz = cosmology->hubbleRate(z) / c_light * c->getCurrentStep();

	// update redshift
	c->setRedshift(z - dz);

	// adiabatic energy loss: dE / dz = E / (1 + z)
	double E = c->current.getEnergy();
	c->current.setEnergy(E * (1 - dz / (1 + z)));
}

std::string FutureRedshift::getDescription() const {
	std::stringstream s;
	s << "FutureRedshift: h0 = " << cosmology->hubbleRate() / 1e5 * Mpc << ", omegaL = "
			<< cosmology->getOmegaL() << ", omegaM = " << cosmology->getOmegaM();
	return s.str();
}

} // namespace crpropa
