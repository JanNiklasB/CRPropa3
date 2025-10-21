#include "crpropa/module/ParticleCollector.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/Units.h"

namespace crpropa {

ParticleCollector::ParticleCollector()
	: nBuffer(10e6), clone(false), recursive(false)  {
	container = new Candidate*[nBuffer]; // for 1e6 candidates ~ 500MB of RAM
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer)
	: nBuffer(nBuffer), clone(false), recursive(false)  {
	container = new Candidate*[nBuffer];
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer, const bool clone)
	: nBuffer(nBuffer), clone(clone), recursive(false) {
	container = new Candidate*[nBuffer];
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer, const bool clone, const bool recursive)
	: nBuffer(nBuffer), clone(clone), recursive(recursive) {
	container = new Candidate*[nBuffer];
}

void ParticleCollector::process(Candidate *c) const {
	#pragma omp critical(ModifyContainer)
	{
		// countermeasure against overflow, but slows down simulation significantly
		if(containerSize>=nBuffer){
			Candidate** tmp = new Candidate*[containerSize+1];
			for(int i=0; i<containerSize; i++)
				tmp[i] = container[i];
			delete[] container;
			container = tmp;
		}
		if(clone)
				container[containerSize] = c->clone(recursive);
		else
			container[containerSize] = c;
		containerSize++;
	}
}

void ParticleCollector::process(ref_ptr<Candidate> c) const {
	ParticleCollector::process(c.get());
}

void ParticleCollector::reprocess(Module *action) const {
	for (int itr=0; itr<containerSize; itr++){
		if (clone)
			action->process(container[itr]->clone(false));
		else
			action->process(container[itr]);
	}
}

void ParticleCollector::dump(const std::string &filename) const {
	TextOutput output(filename.c_str(), Output::Everything);
	reprocess(&output);
	output.close();
}

void ParticleCollector::load(const std::string &filename){
	TextOutput::load(filename.c_str(), this);
}

ParticleCollector::~ParticleCollector() {
	clearContainer();
}

std::size_t ParticleCollector::size() const {
	return containerSize;
}

ref_ptr<Candidate> ParticleCollector::operator[](const std::size_t i) const {
	return container[i];
}

void ParticleCollector::clearContainer() {
	delete[] container;
	container = NULL;
	containerSize = 0;
}

std::vector<Candidate* > ParticleCollector::getContainer() {
	return std::vector<Candidate* >(container, container+containerSize);
}

void ParticleCollector::setClone(bool b) {
	clone = b;
}

bool ParticleCollector::getClone() const {
	return clone;
}

std::string ParticleCollector::getDescription() const {
	return "ParticleCollector";
}

void ParticleCollector::getTrajectory(ModuleList* mlist, std::size_t i, Module *output) const {
	ref_ptr<Candidate> c_tmp = container[i]->clone();

	c_tmp->restart();

	mlist->add(output);
	mlist->run(c_tmp);
	mlist->remove(mlist->size()-1);
}

void ParticleCollector::getTrajectory(ref_ptr<ModuleList> mlist, std::size_t i, ref_ptr<Module> output) const {
	ParticleCollector::getTrajectory((ModuleList*) mlist, i, (Module*) output);
}

} // namespace crpropa
