#ifndef CRPROPA_MODULE_LIST_H
#define CRPROPA_MODULE_LIST_H

#include "crpropa/__CudaDefines.h"

#include "thrust/universal_vector.h"

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

template<class T, class U>
__global__ void createInstance(T** ptr, U* original){
	(*ptr) = new U(*original);
}

/**
 @class ModuleList
 @brief The simulation itself: A list of simulation modules
 */
class ModuleList : public Referenced {
public:
	typedef std::vector<ref_ptr<Module> > module_list_t;
	typedef std::vector<ref_ptr<Candidate> > candidate_vector_t;

	#ifdef __CUDACC__
	__shared__ int maxThreadsPerBlock;
	#endif

	CUDA_CALLABLE_MEMBER ModuleList() : showProgress(false){}
	virtual ~ModuleList();
	void setShowProgress(bool show = true); ///< activate a progress bar

	template<class Type> void add(Type* module);
	void remove(std::size_t i);
	std::size_t size() const;
	ref_ptr<Module> operator[](const std::size_t i);

	CUDA_CALLABLE_MEMBER void process(Candidate* candidate) const; ///< call process in all modules
	CUDA_CALLABLE_MEMBER void process(ref_ptr<Candidate> candidate) const; ///< call process in all modules

	void run(Candidate* candidate, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a single candidate
	void run(ref_ptr<Candidate> candidate, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a single candidate
	void run(const candidate_vector_t *candidates, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a candidate vector
	void run(SourceInterface* source, std::size_t count, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a number of candidates from the given source
	#ifdef __CUDACC__
	void copySecondaries(ref_ptr<Candidate> candidate);
	void cudarun(SourceInterface* source, std::size_t count, bool recursive = true, bool secondariesFirst = false); ///< run simulation for a number of candidates from the given source
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
	#ifdef __CUDACC__
	std::vector<Module**> dmodules;
	Module*** dmodulesPtr=NULL;
	int dmodulesSize=0;
	bool isInKernel = false;
	#endif

	bool showProgress;
	Output* interruptAction;
	bool haveInterruptAction = false;
	std::vector<int> notFinished; // list with not finished numbers of candidates
};

template<class Type>
void ModuleList::add(Type *module) {
	modules.push_back(module);

	#ifdef __CUDACC__
	Module** _ptr;
	gpuErrchk(cudaMalloc(&_ptr, sizeof(Type**)));
	
	createInstance<<<1,1>>>(_ptr, module);
	cudaDeviceSynchronize();
	gpuErrchk(cudaGetLastError());
	
	dmodules.push_back(_ptr);
	dmodulesPtr = dmodules.data();
	dmodulesSize = dmodules.size();
	#endif
}

/**
 @class ModuleListRunner
 @brief Run the provided ModuleList when process is called.
 */
class ModuleListRunner: public Module {
private:
	ModuleList* mlist;
public:
	CUDA_CALLABLE_MEMBER ModuleListRunner(){};
	ModuleListRunner(ModuleList *mlist);
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const; ///< call run of wrapped ModuleList
	std::string getDescription() const;
};

#ifdef __CUDACC__
__global__ void _cudarun(
	Candidate** candidates,
	int candidatesSize,
	const ModuleList* MLIST,
	bool recursive = true, 
	bool secondariesFirst = false
);
#endif

} // namespace crpropa

#endif // CRPROPA_MODULE_LIST_H
