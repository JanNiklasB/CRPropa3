#include "crpropa/module/Tools.h"
#include "crpropa/Clock.h"

#include <iostream>
#include <sstream>

using namespace std;

namespace crpropa {

PerformanceModule::~PerformanceModule() {
	double total = 0;
	for (size_t i = 0; i < modulesSize; i++) {
		_module_info &m = modules[i];
		total += m.time;
	}
	printf("Performance for %d calls:\n", calls);
	for (size_t i = 0; i < modulesSize; i++) {
		_module_info &m = modules[i];
		printf(
			" - %d % -> %s: %d\n", 
			floor((1000 * m.time / total) + 0.5)/10, 
			m.module->getDescription(), 
			(m.time / calls)
		);
	}
	delete[] modules;
	Module::~Module();
}

void PerformanceModule::add(Module *module) {
	_module_info info;
	info.module = module;
	info.time = 0;
	push_back(modules, modulesSize, info);
}

void PerformanceModule::process(Candidate *candidate) const {
	double* times = new double[modulesSize];
	for (size_t i = 0; i < modulesSize; i++) {
		_module_info &m = modules[i];
		#ifndef __CUDACC__
		double start = Clock::getInstance().getMillisecond();
		m.module->process(candidate);
		double end = Clock::getInstance().getMillisecond();
		times[i] = end - start;
		#else
		auto start = cuda::std::chrono::high_resolution_clock::now();
		m.module->process(candidate);
		auto end = cuda::std::chrono::high_resolution_clock::now();
		cuda::std::chrono::duration time = end - start;
		times[i] = cuda::std::chrono::duration_cast<cuda::std::chrono::milliseconds>(time).count();
		#endif
	}

	#pragma omp critical(PerformanceModule)
	{
		for (size_t i = 0; i < modulesSize; i++) {
			_module_info &m = modules[i];
			m.time += times[i];
		}
		calls++;
	}
	delete[] times;
}

string PerformanceModule::getDescription() const {
	stringstream sstr;
	sstr << "PerformanceModule (";
	for (size_t i = 0; i < modulesSize; i++) {
		_module_info &m = modules[i];
		if (i > 0)
			sstr << ", ";
		sstr << m.module->getDescription();
	}
	sstr << ")";
	return sstr.str();
}

// ----------------------------------------------------------------------------

ParticleFilter::ParticleFilter(const std::set<int> &ids) {
	idsSize = ids.size();
	this->ids = new int[idsSize];
	int counter=0;
	for(auto id=ids.begin(); id!=ids.end(); id++){
		this->ids[counter] = *id;
		counter++;
	}
}

ParticleFilter::~ParticleFilter(){
	delete[] ids;
	AbstractCondition::~AbstractCondition();
}

void ParticleFilter::addId(int id) {
	insert(ids, idsSize, id, lower_bound(id, ids, idsSize));
}

void ParticleFilter::removeId(int id) {
	erase(ids, idsSize, findSorted(id, ids, idsSize));
}

std::set<int> ParticleFilter::getIds() {
	return std::set<int>(ids, ids+idsSize);
}

void ParticleFilter::process(Candidate* candidate) const {
	if (findSorted(candidate->current.getId(), ids, idsSize) == (idsSize-1))
		reject(candidate);
	else
		accept(candidate);
}

string ParticleFilter::getDescription() const {
	stringstream sstr;
	sstr << "ParticleFilter: ";
	for (int i=0; i<idsSize; i++){
		sstr << ids[i] << ", ";
	}
	sstr << ")";
	return sstr.str();
}

// ----------------------------------------------------------------------------
EmissionMapFiller::EmissionMapFiller(EmissionMap *emissionMap) : emissionMap(emissionMap) {
}

void EmissionMapFiller::setEmissionMap(EmissionMap *emissionMap) {
	this->emissionMap = emissionMap;
}

void EmissionMapFiller::process(Candidate* candidate) const {
	if (emissionMap) {
		#pragma omp critical(EmissionMap)
		{
			emissionMap->fillMap(candidate->source);
		}
	}
}

string EmissionMapFiller::getDescription() const {
	return "EmissionMapFiller";
}

} // namespace crpropa
