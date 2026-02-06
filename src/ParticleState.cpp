#include "crpropa/ParticleState.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"

#include "HepPID/ParticleIDMethods.hh"

#include <cstdlib>
#include <cmath>
#include <sstream>

namespace crpropa {

ParticleState::ParticleState(int id, double E, Vector3d pos, Vector3d dir): id(0), energy(0.), position(0.), direction(0.), pmass(0.), charge(0.)
{
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
	energy = std::max(pmass*c_squared, newEnergy); // energy has to be at least restmass
}

double ParticleState::getEnergy() const {
	return energy;
}

double ParticleState::getRigidity() const {
	return fabs(energy / charge);
}

void ParticleState::setId(int newId) {
	id = newId;
	pmass = particleMass(id);
	if (isNucleus(id)) {
		charge = chargeNumber(id) * eplus;
		if (id < 0)
			charge *= -1; // anti-nucleus
	} else {
		charge = HepPID::charge(id) * eplus;
	}

	// energy has to be at least restmass
	energy = std::max(pmass*c_squared, energy);
}

int ParticleState::getId() const {
	return id;
}

double ParticleState::getMass() const {
	return pmass;
}

void ParticleState::setMass(double Mass) {
	pmass = Mass;
}

double ParticleState::getCharge() const {
	return charge;
}

void ParticleState::setCharge(int ChargeNumber) {
	charge = ChargeNumber * eplus;
}

double ParticleState::getLorentzFactor() const {
	return energy / (pmass * c_squared);
}

void ParticleState::setLorentzFactor(double lf) {
	lf = std::max(0., lf); // prevent negative Lorentz factors
	setEnergy(lf * pmass * c_squared);
}

Vector3d ParticleState::getVelocity() const {
	if (getMass()==0) return direction*c_light;
	return direction * c_light*sqrt(1-1/pow(getLorentzFactor(), 2));
}

Vector3d ParticleState::getMomentum() const {
	return direction * sqrt( pow(energy/c_light, 2) - pow(pmass*c_light, 2) );
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
