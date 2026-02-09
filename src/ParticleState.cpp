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
	direction = dir.getUnitVector();
}

const Vector3d &ParticleState::getDirection() const {
	return direction;
}

void ParticleState::setEnergy(double newEnergy) {
	energy = std::max(0., newEnergy);
}

double ParticleState::getEnergy() const {
	return energy;
}

double ParticleState::getRigidity() const {
	return fabs(energy / charge);
}

void ParticleState::setId(int newId) {
	id = newId;
	setMass(particleMass(id));
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
	if (getMass()==0)
		return INFINITY;
	return getEnergy()/getMass()/c_squared + 1;
}

void ParticleState::setLorentzFactor(double lf) {
	lf = std::max(0., lf); // prevent negative Lorentz factors
	setEnergy((lf-1) * pmass * c_squared);
}

double ParticleState::getBeta() const {
	return getVelocity().getR()/c_squared;
}

Vector3d ParticleState::getVelocity() const {
	Vector3d velocity;
	if (getMass()==0) 
		velocity = getDirection()*c_light;
	else if (getLorentzFactor()==1)  // can happen if if gamma-1 < numericalPrecission
		velocity = getDirection() * sqrt(getEnergy()*2/getMass());  // non relativistic case
	else
		velocity = getDirection() * c_light*sqrt(1-1/pow(getLorentzFactor(), 2));

	return velocity;
}

Vector3d ParticleState::getMomentum() const {
	return getLorentzFactor()*getMass()*getVelocity();
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
