#include "crpropa/module/Tools.h"
#include "crpropa/Clock.h"

#include <iostream>
#include <sstream>

using namespace std;

namespace crpropa {

PerformanceModule::~PerformanceModule() {
	double total = 0;
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		total += m.time;
	}
	cout << "Performance for " << calls << " calls:" << endl;
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		cout << " - " << floor((1000 * m.time / total) + 0.5) / 10 << "% -> "
				<< m.module->getDescription() << ": " << (m.time / calls)
				<< endl;
	}
}

void PerformanceModule::add(Module *module) {
	_module_info info;
	info.module = module;
	info.time = 0;
	modules.push_back(info);
	modulesPtr = modules.data();
	modulesSize = modules.size();
}

void PerformanceModule::process(Candidate *candidate) const {
	double* times(modulesSize);
	for (size_t i = 0; i < modulesSize; i++) {
		_module_info &m = modulesPtr[i];
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
			_module_info &m = modulesPtr[i];
			m.time += times[i];
		}
		calls++;
	}
}

string PerformanceModule::getDescription() const {
	stringstream sstr;
	sstr << "PerformanceModule (";
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		if (i > 0)
			sstr << ", ";
		sstr << m.module->getDescription();
	}
	sstr << ")";
	return sstr.str();
}

// ----------------------------------------------------------------------------

ParticleFilter::ParticleFilter() {
}

ParticleFilter::ParticleFilter(const std::set<int> &ids) : ids(ids) {
	#ifdef __CUDACC__
	updatePtr();
	#endif
}

ParticleFilter::~ParticleFilter(){
	#ifdef __CUDACC__
	delete[] idsPtr;
	#endif
}

#ifdef __CUDACC__
void ParticleFilter::updatePtr(){
	idsSize = ids.size();
	delete[] idsPtr;
	idsPtr = new int[idsSize];
	int counter = 0;
	for (auto i=ids.begin(); i!=ids.end(); i++){
		idsPtr[counter] = *i;
		counter++;
	}
}
#endif

void ParticleFilter::addId(int id) {
	ids.insert(id);
	#ifdef __CUDACC__
	updatePtr();
	#endif
}

void ParticleFilter::removeId(int id) {
	ids.erase(id);
	#ifdef __CUDACC__
	updatePtr();
	#endif
}

std::set<int> &ParticleFilter::getIds() {
	return ids;
}

void ParticleFilter::process(Candidate* candidate) const {
	#ifndef __CUDACC__
	if (ids.find(candidate->current.getId()) == ids.end())
	#else
	if (*cuda::std::find(&idsPtr[0], &idsPtr[idsSize-1], candidate->current.getId()) == idsPtr[idsSize-1])
	#endif
		reject(candidate);
	else
		accept(candidate);
}

string ParticleFilter::getDescription() const {
	stringstream sstr;
	sstr << "ParticleFilter: ";
	for (std::set<int>::const_iterator i = ids.begin(); i != ids.end(); i++) {
		sstr << *i << ", ";
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
