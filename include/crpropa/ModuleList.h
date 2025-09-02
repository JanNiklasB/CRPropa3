#ifndef CRPROPA_MODULE_LIST_H
#define CRPROPA_MODULE_LIST_H

#include <algorithm>
#include <csignal>
#include <iostream>
#include <vector>
#include <exception>
#include <sstream>

#include "crpropa/Candidate.h"
#include "crpropa/Module.h"
#include "crpropa/Source.h"
#include "crpropa/module/Output.h"

namespace crpropa {

/**
 @class ModuleList
 @brief The simulation itself: A list of simulation modules
 */
class ModuleList: public Module {
public:
	typedef std::vector<ref_ptr<Module> > module_list_t;
	typedef std::vector<ref_ptr<Candidate> > candidate_vector_t;

	#ifdef __CUDACC__
	__shared__ int maxThreadsPerBlock;
	#endif

	ModuleList();
	CUDA_CALLABLE_MEMBER virtual ~ModuleList();
	void setShowProgress(bool show = true); ///< activate a progress bar

	void add(Module* module);
	void remove(std::size_t i);
	std::size_t size() const;
	ref_ptr<Module> operator[](const std::size_t i);

	CUDA_CALLABLE_MEMBER void process(Candidate* candidate) const; ///< call process in all modules
	CUDA_CALLABLE_MEMBER void process(ref_ptr<Candidate> candidate) const; ///< call process in all modules

	void run(Candidate* candidate, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a single candidate
	void run(ref_ptr<Candidate> candidate, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a single candidate
	void run(const candidate_vector_t *candidates, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a candidate vector
	void run(SourceInterface* source, size_t count, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a number of candidates from the given source
	#ifdef __CUDACC__
	void copySecondaries(ref_ptr<Candidate> candidate);
	void cudarun(SourceInterface* source, size_t count, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a number of candidates from the given source
	#endif

	std::string getDescription() const;
	void showModules() const;
	
	/** iterator goodies */
	typedef module_list_t::iterator iterator;
	typedef module_list_t::const_iterator const_iterator;
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;

	void setInterruptAction(Output* action);
	void dumpCandidate(ref_ptr<Candidate> cand) const;

private:
	module_list_t modules;
	bool showProgress;
	Output* interruptAction;
	bool haveInterruptAction = false;
	std::vector<int> notFinished; // list with not finished numbers of candidates
	#ifdef __CUDACC__
	ref_ptr<Module>* modulesPtr=NULL;
	int modulesSize=0;
	#endif
};

/**
 @class ModuleListRunner
 @brief Run the provided ModuleList when process is called.
 */
class ModuleListRunner: public Module {
private:
	ref_ptr<ModuleList> mlist;
public:

	ModuleListRunner(ModuleList *mlist);
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const; ///< call run of wrapped ModuleList
	std::string getDescription() const;
};

#ifdef __CUDACC__
__global__ void _cudarun(
	ref_ptr<Candidate>* candidates,
	int candidatesSize,
	const ModuleList* ModuleList,
	bool recursive = true, 
	bool secondariesFirst = false
);
#endif

} // namespace crpropa

#endif // CRPROPA_MODULE_LIST_H
