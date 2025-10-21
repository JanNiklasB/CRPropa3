#include "crpropa/ParticleState.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"

#include "HepPID/ParticleIDMethods.hh"

#include <cstdlib>
#include <sstream>

namespace crpropa {

ParticleState::ParticleState()
	: id(0), energy(0.), position(0.), direction(0.), pmass(0.), charge(0.), nuclearMassTable(NULL){
	setId(id);
	setEnergy(0);
	setPosition(Vector3d(0,0,0));
	setDirection(Vector3d(-1,0,0));
}

ParticleState::ParticleState(NuclearMassTable* nuclearMassTable, int id, double E, Vector3d pos, Vector3d dir)
	: id(0), energy(0.), position(0.), direction(0.), pmass(0.), charge(0.), nuclearMassTable(nuclearMassTable){
	setId(id);
	setEnergy(E);
	setPosition(pos);
	setDirection(dir);
}

ParticleState::ParticleState(int id, double E, Vector3d pos, Vector3d dir)
	: id(0), energy(0.), position(0.), direction(0.), pmass(0.), charge(0.), nuclearMassTable(NULL){
	setId(id);
	setEnergy(E);
	setPosition(pos);
	setDirection(dir);
}

void ParticleState::setPosition(const Vector3d &pos) {
	position = pos;
}

const Vector3d &ParticleState::getPosition() const {
	return position;
}

void ParticleState::setDirection(const Vector3d &dir) {
	direction = dir / dir.getR();
}

const Vector3d &ParticleState::getDirection() const {
	return direction;
}

void ParticleState::setEnergy(double newEnergy) {
	energy = std::max(0., newEnergy); // prevent negative energies
}

double ParticleState::getEnergy() const {
	return energy;
}

double ParticleState::getRigidity() const {
	return fabs(energy / charge);
}

void ParticleState::setId(int newId) {
	id = newId;
	if(nuclearMassTable)
		pmass = nuclearMassTable->particleMass(id);
	else
		pmass = 0;
	if (isNucleus(id)) {
		charge = chargeNumber(id) * eplus;
		if (id < 0)
			charge *= -1; // anti-nucleus
	} else {
		charge = HepPID::charge(id) * eplus;
	}
}

int ParticleState::getId() const {
	return id;
}

double ParticleState::getMass() const {
	return pmass;
}

void ParticleState::setNuclearMassTable(NuclearMassTable* nuclearMassTable){
	this->nuclearMassTable = nuclearMassTable;
}

NuclearMassTable* ParticleState::getNuclearMassTable() const{
	return nuclearMassTable;
}

double ParticleState::getCharge() const {
	return charge;
}

double ParticleState::getLorentzFactor() const {
	return energy / (pmass * c_squared);
}

void ParticleState::setLorentzFactor(double lf) {
	#ifdef __CUDACC__
	lf = cuda::std::max(0., lf);  // prevent negative Lorentz factors
	#else
	lf = std::max(0., lf); // prevent negative Lorentz factors
	#endif
	energy = lf * pmass * c_squared;
}

Vector3d ParticleState::getVelocity() const {
	return direction * c_light;
}

Vector3d ParticleState::getMomentum() const {
	return direction * (energy / c_light);
}

std::string ParticleState::getDescription() const {
	std::stringstream ss;
	ss << "Particle " << id << ", ";
	ss << "E = " << energy / EeV << " EeV, ";
	ss << "x = " << position / Mpc << " Mpc, ";
	ss << "p = " << direction;
	return ss.str();
}

} // namespace crpropa
