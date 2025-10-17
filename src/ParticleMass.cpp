#include "crpropa/ParticleMass.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"

#include "kiss/convert.h"
#include "kiss/logger.h"

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>

namespace crpropa {

NuclearMassTable::NuclearMassTable() {
	init();
}

NuclearMassTable::~NuclearMassTable(){
	delete[] table;
}

void NuclearMassTable::init() {
	std::string filename = getDataPath("nuclear_mass.txt");
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("crpropa: could not open file " + filename);

	int Z, N;
	double mass;
	while (infile.good()) {
		if (infile.peek() != '#') {
			infile >> Z >> N >> mass;
			push_back(table, tableSize, mass);
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	infile.close();
}

double NuclearMassTable::getMass(std::size_t idx) {
	return table[idx];
}

double NuclearMassTable::particleMass(int id) {
	if (isNucleus(id))
		return nuclearMass(id);
	if (abs(id) == 11)
		return mass_electron;
	return 0.0;
}

double NuclearMassTable::nuclearMass(int id) {
	int A = massNumber(id);
	int Z = chargeNumber(id);
	return nuclearMass(A, Z);
}

double NuclearMassTable::nuclearMass(int A, int Z) {
	if ((A < 1) or (A > 56) or (Z < 0) or (Z > 26) or (Z > A)) {
		#ifndef __CUDACC__
		KISS_LOG_WARNING <<
		"nuclearMass: nuclear mass not found in the mass table for " <<
	        "A = " << A << ", Z = " << Z << ". " <<
		"Approximated value used A * amu - Z * m_e instead.";
		#endif
		return A * amu - Z * mass_electron;
	}
	int N = A - Z;
	return getMass(Z * 31 + N);
}

static NuclearMassTable nuclearMassTable;

double particleMass(int id){
	return nuclearMassTable.particleMass(id);
}

double nuclearMass(int id){
	return nuclearMassTable.nuclearMass(id);
}

double nuclearMass(int A, int Z){
	return nuclearMassTable.nuclearMass(A, Z);
}

} // namespace crpropa
