#ifndef CRPROPA_MODULETOOLS_H
#define CRPROPA_MODULETOOLS_H

#include "crpropa/__CudaDefines.h"
#include "crpropa/Module.h"
#include "crpropa/EmissionMap.h"

#include <vector>
#include <set>
#include <chrono>

namespace crpropa {
/**
 * \addtogroup Tools
 * @{
 */

/**
 @class PerformanceModule
 @brief Module to monitor the simulation performance

 Add modules under investigation to this module instead of the ModuleList.
 */
class PerformanceModule: public Module {
private:
	struct _module_info {
		double time;
		Module* module;
	};

	mutable _module_info* modules;
	mutable int modulesSize;
	mutable size_t calls;

public:
	CUDA_CALLABLE_MEMBER PerformanceModule(){}
	~PerformanceModule();
	void add(Module* module);
	CUDA_CALLABLE_MEMBER void process(Candidate* candidate) const;
	std::string getDescription() const;
};

/**
  @class ParticleFilter
  @brief Reject Particles not listed in filter.
*/
class ParticleFilter: public AbstractCondition {
	int* ids, idsSize;

public:
	CUDA_CALLABLE_MEMBER ParticleFilter(){};
	ParticleFilter(const std::set<int> &ids);
	~ParticleFilter();

	void addId(int id);
	void removeId(int remove);
	std::set<int> getIds();

	CUDA_CALLABLE_MEMBER void process(Candidate* candidate) const;
	std::string getDescription() const;
};


/**
  @class EmissionMapFiller
  @brief Fill EmissionMap with source particle state
*/
class EmissionMapFiller: public Module {
	EmissionMap* emissionMap;
public:
	CUDA_CALLABLE_MEMBER EmissionMapFiller(){};
	EmissionMapFiller(EmissionMap *emissionMap);
	void setEmissionMap(EmissionMap *emissionMap);
	CUDA_CALLABLE_MEMBER void process(Candidate* candidate) const;
	std::string getDescription() const;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_MODULETOOLS_H
