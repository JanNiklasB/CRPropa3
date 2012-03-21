#include "mpc/module/Output.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <kiss/convert.h>

namespace mpc {

std::string getOutputString(ParticleState particle) {
	std::stringstream s;
	s << particle.getId() << ", ";
	s << particle.getEnergy() / EeV << ", ";
	Vector3 pos = particle.getPosition() / Mpc;
	s << pos.x() << ", ";
	s << pos.y() << ", ";
	s << pos.z() << ", ";
	Vector3 dir = particle.getDirection();
	s << dir.x() << ", ";
	s << dir.y() << ", ";
	s << dir.z();
	return s.str();
}

TrajectoryOutput::TrajectoryOutput(std::string name) {
	outfile.open(name.c_str());
	outfile << "# Age, HepId, E, posX, posY, posZ, dirX, dirY, dirZ, event\n";
}

TrajectoryOutput::~TrajectoryOutput() {
	outfile.close();
}

void TrajectoryOutput::process(Candidate *candidate) {
#pragma omp critical
	{
		outfile << candidate->getTrajectoryLength() / Mpc << ", "
				<< getOutputString(candidate->current) << "\n";
	}
}

std::string TrajectoryOutput::getDescription() const {
	return "Trajectory output";
}

FlaggedOutput::FlaggedOutput(std::string name, Candidate::Status flag) {
	this->flag = flag;
	outfile.open(name.c_str());
	outfile
			<< "id, x, y, z, E, phi, theta, distance, i_id, i_x, i_y, i_z, i_E, i_phi, i_theta\n";
}

FlaggedOutput::~FlaggedOutput() {
	outfile.close();
}

void FlaggedOutput::process(Candidate *candidate) const {
	if (candidate->getStatus() != flag)
		return;

#pragma omp critical
	{
		outfile << candidate->current.getId() << ", ";
		outfile << candidate->current.getPosition().x() / Mpc << ", ";
		outfile << candidate->current.getPosition().y() / Mpc << ", ";
		outfile << candidate->current.getPosition().z() / Mpc << ", ";
		outfile << candidate->current.getEnergy() / EeV << ", ";
		outfile << candidate->current.getDirection().phi() << ", ";
		outfile << candidate->current.getDirection().theta() << ", ";
		outfile << candidate->getTrajectoryLength() / Mpc << ", ";
		outfile << candidate->initial.getId() << ", ";
		outfile << candidate->initial.getPosition().x() / Mpc << ", ";
		outfile << candidate->initial.getPosition().y() / Mpc << ", ";
		outfile << candidate->initial.getPosition().z() / Mpc << ", ";
		outfile << candidate->initial.getEnergy() / EeV << ", ";
		outfile << candidate->initial.getDirection().phi() << ", ";
		outfile << candidate->initial.getDirection().theta();
		outfile << std::endl;
	}
}

std::string FlaggedOutput::getDescription() const {
	switch (flag) {
	case Candidate::Active:
		return "FlaggedOutput, Active";
	case Candidate::Detected:
		return "FlaggedOutput, Detected";
	case Candidate::OutOfBounds:
		return "FlaggedOutput, OutOfBounds";
	case Candidate::Stopped:
		return "FlaggedOutput, Stopped";
	case Candidate::UserDefined:
	default:
		return "FlaggedOutput, user defined (" + kiss::str(flag) + ")";
	}

}

void ShellOutput::process(Candidate *candidate) const {
#pragma omp critical
	{
		std::cout << std::fixed << std::showpoint << std::setprecision(2)
				<< std::setw(6);
		std::cout << candidate->getTrajectoryLength() / Mpc << " Mpc,  ";
		std::cout << candidate->current.getId() << ",  ";
		std::cout << candidate->current.getEnergy() / EeV << " EeV,  ";
		std::cout << candidate->current.getPosition() / Mpc << " Mpc, Status: ";
		std::cout << candidate->getStatus();
		std::cout << std::endl;
	}
}

std::string ShellOutput::getDescription() const {
	return "ShellOutput";
}

} // namespace mpc