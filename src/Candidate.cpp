#include "crpropa/Candidate.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include <stdexcept>

namespace crpropa {

CUDA_CALLABLE_MEMBER uint64_t nextSerialNumberGlobal = 0;

Candidate::Candidate() :
redshift(0), trajectoryLength(0), weight(1.), currentStep(0), nextStep(0), active(true), parent(0), tagOrigin("PRIM"), time(0), NuclearMassPtr(NULL) {
	ParticleState state;
	source = state;
	created = state;
	previous = state;
	current = state;

	nextSerialNumber = &nextSerialNumberGlobal;

	#if defined(OPENMP_3_1)
		#pragma omp atomic capture
		{serialNumber = (*nextSerialNumber)++;}
	#elif defined(__GNUC__)
		{serialNumber = __sync_add_and_fetch(nextSerialNumber, 1);}
	#else
		#pragma omp critical(serialNumber)
		{serialNumber = (*nextSerialNumber)++;}
	#endif
}

Candidate::Candidate(int id, double E, Vector3d pos, Vector3d dir, double z, double weight, std::string tagOrigin) :
redshift(z), trajectoryLength(0), weight(weight), currentStep(0), nextStep(0), active(true), parent(0), tagOrigin(tagOrigin), time(0) {
	NuclearMassPtr = NULL;
	ParticleState state(NuclearMassPtr, id, E, pos, dir);
	source = state;
	created = state;
	previous = state;
	current = state;

	nextSerialNumber = &nextSerialNumberGlobal;

	#if defined(OPENMP_3_1)
		#pragma omp atomic capture
		{serialNumber = (*nextSerialNumber)++;}
	#elif defined(__GNUC__)
		{serialNumber = __sync_add_and_fetch(nextSerialNumber, 1);}
	#else
		#pragma omp critical(serialNumber)
		{serialNumber = (*nextSerialNumber)++;}
	#endif
}

Candidate::Candidate(NuclearMassTable* NuclearMassTablePtr, int id, double E, Vector3d pos, Vector3d dir, double z, double weight, std::string tagOrigin) :
redshift(z), trajectoryLength(0), weight(weight), currentStep(0), nextStep(0), active(true), parent(0), tagOrigin(tagOrigin), time(0), NuclearMassPtr(NuclearMassTablePtr) {
	ParticleState state(NuclearMassTablePtr, id, E, pos, dir);
	source = state;
	created = state;
	previous = state;
	current = state;

	nextSerialNumber = &nextSerialNumberGlobal;

	#if defined(OPENMP_3_1)
		#pragma omp atomic capture
		{serialNumber = (*nextSerialNumber)++;}
	#elif defined(__GNUC__)
		{serialNumber = __sync_add_and_fetch(nextSerialNumber, 1);}
	#else
		#pragma omp critical(serialNumber)
		{serialNumber = (*nextSerialNumber)++;}
	#endif
}

Candidate::Candidate(const ParticleState &state) :
	source(state), created(state), current(state), previous(state), redshift(0), trajectoryLength(0), currentStep(0), nextStep(0), active(true), parent(0), tagOrigin ("PRIM"), time(0) {
	
	NuclearMassPtr = state.getNuclearMassTable();
	nextSerialNumber = &nextSerialNumberGlobal;

	#if defined(OPENMP_3_1)
		#pragma omp atomic capture
		{serialNumber = (*nextSerialNumber)++;}
	#elif defined(__GNUC__)
		{serialNumber = __sync_add_and_fetch(nextSerialNumber, 1);}
	#else
		#pragma omp critical(serialNumber)
		{serialNumber = (*nextSerialNumber)++;}
	#endif
}

Candidate::~Candidate(){
	delete[] secondaries;
}

bool Candidate::isActive() const {
	return active;
}

void Candidate::setActive(bool b) {
	active = b;
}

int Candidate::getSecondarySize() const {
	return secondariesSize;
}

double Candidate::getRedshift() const {
	return redshift;
}

double Candidate::getTrajectoryLength() const {
	return trajectoryLength;
}

double Candidate::getVelocity() const {
	return c_light;
}

double Candidate::getWeight() const {
	return weight;
}

double Candidate::getCurrentStep() const {
	return currentStep;
}

double Candidate::getNextStep() const {
	return nextStep;
}

void Candidate::setRedshift(double z) {
	redshift = z;
}

void Candidate::setTrajectoryLength(double a) {
	trajectoryLength = a;
}

void Candidate::setWeight(double w) {
	weight = w;
}

void Candidate::updateWeight(double w) {
	weight *= w;
}

void Candidate::setCurrentStep(double lstep) {
	currentStep = lstep;
	trajectoryLength += lstep;
	time += lstep / getVelocity();
}

void Candidate::setNextStep(double step) {
	nextStep = step;
}

void Candidate::limitNextStep(double step) {
	#ifdef __CUDACC__
	nextStep = cuda::std::min(nextStep, step);
	#else
	nextStep = std::min(nextStep, step);
	#endif
}

void Candidate::setProperty(const std::string &name, const Variant &value) {
	properties[name] = value;
}

void Candidate::setTagOrigin (std::string tagOrigin) {
	this->tagOrigin = tagOrigin;
}

std::string Candidate::getTagOrigin () const {
	return tagOrigin;
}

void Candidate::setTime(double t) {
	time = t;
}

double Candidate::getTime() const {
	return time;
}

const Variant &Candidate::getProperty(const std::string &name) const {
	PropertyMap::const_iterator i = properties.find(name);
	if (i == properties.end())
		throw std::runtime_error("Unknown candidate property: " + name);
	return i->second;
}

bool Candidate::removeProperty(const std::string& name) {
	PropertyMap::iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	properties.erase(i);
	return true;
}

bool Candidate::hasProperty(const std::string &name) const {
	PropertyMap::const_iterator i = properties.find(name);
	if (i == properties.end())
		return false;
	return true;
}

void Candidate::addSecondary(Candidate *c) {
	push_back(secondaries, secondariesSize, c);
}

void Candidate::addSecondary(int id, double energy, double w, std::string tagOrigin) {
	ref_ptr<Candidate> secondary = new Candidate;
	secondary->setRedshift(redshift);
	secondary->setTrajectoryLength(trajectoryLength);
	secondary->setTime(time);
	secondary->setWeight(weight * w);
	secondary->setTagOrigin(tagOrigin);
	secondary->setNuclearMassTable(NuclearMassPtr);
	for (PropertyMap::const_iterator it = properties.begin(); it != properties.end(); ++it) {
		secondary->setProperty(it->first, it->second);		
	}
	secondary->source = source;
	secondary->previous = previous;
	secondary->created = previous;
	secondary->current = current;
	secondary->current.setId(id);
	secondary->current.setEnergy(energy);
	secondary->parent = this;
	addSecondary(secondary);
}

void Candidate::addSecondary(int id, double energy, Vector3d position, double w, std::string tagOrigin) {
	ref_ptr<Candidate> secondary = new Candidate;
	secondary->setRedshift(redshift);
	secondary->setTrajectoryLength(trajectoryLength - (current.getPosition() - position).getR());
	secondary->setTime(time - (current.getPosition() - position).getR() / getVelocity());
	secondary->setWeight(weight * w);
	secondary->setTagOrigin(tagOrigin);
	secondary->setNuclearMassTable(NuclearMassPtr);
	for (PropertyMap::const_iterator it = properties.begin(); it != properties.end(); ++it) {
		secondary->setProperty(it->first, it->second);		
	}
	secondary->source = source;
	secondary->previous = previous;
	secondary->created = previous;
	secondary->current = current;
	secondary->current.setId(id);
	secondary->current.setEnergy(energy);
	secondary->current.setPosition(position);
	secondary->created.setPosition(position);
	secondary->parent = this;
	addSecondary(secondary);
}

void Candidate::clearSecondaries() {
	delete[] secondaries;
	secondaries = NULL;
	secondariesSize = 0;
}

std::string Candidate::getDescription() const {
	std::stringstream ss;
	ss << "CosmicRay at z = " << getRedshift() << "\n";
	ss << "  source:  " << source.getDescription() << "\n";
	ss << "  current: " << current.getDescription();
	return ss.str();
}

void Candidate::copy(const Candidate* C){
	source = C->source;
	created = C->created;
	current = C->current;
	previous = C->previous;

	properties = C->properties;
	active = C->active;
	redshift = C->redshift;
	weight = C->weight;
	trajectoryLength = C->trajectoryLength;
	time = C->time;
	currentStep = C->currentStep;
	nextStep = C->nextStep;
	parent = C->parent;
	NuclearMassPtr = C->NuclearMassPtr;
}

Candidate* Candidate::clone(bool recursive) const {
	ref_ptr<Candidate> cloned = new Candidate;
	cloned->copy(this);
	if (recursive) {
		for (int i=0; i<secondariesSize; i++){
			ref_ptr<Candidate> s = secondaries[i]->clone(recursive);
			s->parent = cloned;
			cloned->addSecondary(s);
		}
	}
	return cloned;
}

uint64_t Candidate::getSerialNumber() const {
	return serialNumber;
}

void Candidate::setSerialNumber(const uint64_t snr) {
	serialNumber = snr;
}

uint64_t Candidate::getSourceSerialNumber() const {
	if (parent)
		return parent->getSourceSerialNumber();
	else
		return serialNumber;
}

uint64_t Candidate::getCreatedSerialNumber() const {
	if (parent)
		return parent->getSerialNumber();
	else
		return serialNumber;
}

void Candidate::setNextSerialNumber(uint64_t snr) {
	nextSerialNumberGlobal = snr;
	#ifdef __CUDACC__
	// need to update global after modifing, direct modification is only done on host
	// other modifications are done over pointer to global which is synced automatically
	cudaMemcpyToSymbol(nextSerialNumberGlobal, &nextSerialNumberGlobal, sizeof(uint64_t));
	#endif
}

uint64_t Candidate::getNextSerialNumber() {
	return nextSerialNumberGlobal;
}

void Candidate::setNuclearMassTable(NuclearMassTable* NuclearMassTablePtr) {
	NuclearMassPtr = NuclearMassTablePtr;
	source.setNuclearMassTable(NuclearMassPtr);
	created.setNuclearMassTable(NuclearMassPtr);
	previous.setNuclearMassTable(NuclearMassPtr);
	current.setNuclearMassTable(NuclearMassPtr);
}

NuclearMassTable* Candidate::getNuclearMassTable() const {
	return NuclearMassPtr;
}

void Candidate::restart() {
	setActive(true);
	setTrajectoryLength(0);
	setTime(0);
	previous = source;
	current = source;
}

} // namespace crpropa
