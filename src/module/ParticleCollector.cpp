#include "crpropa/module/ParticleCollector.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/Units.h"

namespace crpropa {

ParticleCollector::ParticleCollector()
	: nBuffer(10e6), clone(false), recursive(false)  {
	container.reserve(nBuffer); // for 1e6 candidates ~ 500MB of RAM
	containerPtr = container.data();
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer)
	: nBuffer(nBuffer), clone(false), recursive(false)  {
	container.reserve(nBuffer);
	containerPtr = container.data();
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer, const bool clone)
	: nBuffer(nBuffer), clone(clone), recursive(false) {
	container.reserve(nBuffer);
	containerPtr = container.data();
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer, const bool clone, const bool recursive)
	: nBuffer(nBuffer), clone(clone), recursive(recursive) {
	container.reserve(nBuffer);
	containerPtr = container.data();
}

void ParticleCollector::process(Candidate *c) const {
	#pragma omp critical(ModifyContainer)
	{
		// countermeasure against overflow, but slows down simulation significantly
		if(containerSize>=nBuffer){
			ref_ptr<Candidate>* tmp = new ref_ptr<Candidate>[containerSize+1];
			for(int i=0; i<containerSize; i++)
				tmp[i] = containerPtr[i];
			delete[] containerPtr;
			containerPtr = tmp;
		}
		if(clone)
				containerPtr[containerSize] = c->clone(recursive);
		else
			containerPtr[containerSize] = c;
		containerSize++;
	}
}

void ParticleCollector::process(ref_ptr<Candidate> c) const {
	ParticleCollector::process(c.get());
}

void ParticleCollector::reprocess(Module *action) const {
	checkVector();
	for (ParticleCollector::iterator itr = container.begin(); itr != container.end(); ++itr){
		if (clone)
			action->process((*(itr->get())).clone(false));
		else
			action->process(itr->get());
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

void ParticleCollector::checkVector() const{
	if(containerSize>nBuffer){
		container = std::vector<ref_ptr<Candidate>>(&containerPtr[0], &containerPtr[containerSize-1]);
	}
}

std::size_t ParticleCollector::size() const {
	return containerSize;
}

ref_ptr<Candidate> ParticleCollector::operator[](const std::size_t i) const {
	return containerPtr[i];
}

void ParticleCollector::clearContainer() {
	container.clear();
	containerPtr = NULL;
	containerSize = 0;
}

std::vector<ref_ptr<Candidate> >& ParticleCollector::getContainer() const {
	checkVector();
	return container;
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

ParticleCollector::iterator ParticleCollector::begin() {
	checkVector();
	return container.begin();
}

ParticleCollector::const_iterator ParticleCollector::begin() const {
	checkVector();
	return container.begin();
}

ParticleCollector::iterator ParticleCollector::end() {
	checkVector();
	return container.end();
}

ParticleCollector::const_iterator ParticleCollector::end() const {
	checkVector();
	return container.end();
}

void ParticleCollector::getTrajectory(ModuleList* mlist, std::size_t i, Module *output) const {
	ref_ptr<Candidate> c_tmp = containerPtr[i]->clone();

	c_tmp->restart();

	mlist->add(output);
	mlist->run(c_tmp);
	mlist->remove(mlist->size()-1);
}

void ParticleCollector::getTrajectory(ref_ptr<ModuleList> mlist, std::size_t i, ref_ptr<Module> output) const {
	ParticleCollector::getTrajectory((ModuleList*) mlist, i, (Module*) output);
}

} // namespace crpropa
