#ifndef CRPROPA_MODULE_H
#define CRPROPA_MODULE_H

#include "crpropa/__CudaDefines.h"
#include "crpropa/Candidate.h"
#include "crpropa/Referenced.h"
#include "crpropa/Common.h"
#include "cuda.h"

#include <string>

namespace crpropa {

class Candidate;

/**
 @class Module
 @brief Abstract base class for modules
 */
class Module: public Referenced {
	const char* description="";
public:
	Module();
	virtual ~Module();

	virtual std::string getDescription() const;
	void setDescription(const std::string &description);
	CUDA_CALLABLE_MEMBER virtual void process(Candidate *candidate) const {};
	CUDA_CALLABLE_MEMBER inline void process(ref_ptr<Candidate> candidate) const {
		process(candidate.get());
	}
	// this function needs to be created manually on device since it is virtual!
	CUDA_CALLABLE_MEMBER virtual void test() const{
		printf("Module::test(pointer)\n");
	}
};


/**
 @class AbstractCondition
 @brief Abstract Module providing common features for conditional modules.
 */
class AbstractCondition: public Module {
protected:
	Module *rejectAction, *acceptAction;
	bool makeRejectedInactive, makeAcceptedInactive;
	char *rejectFlagKey, *rejectFlagValue;
	int rejectFlagKeySize, rejectFlagValueSize;
	char *acceptFlagKey, *acceptFlagValue;
	int acceptFlagKeySize, acceptFlagValueSize;

	CUDA_CALLABLE_MEMBER void reject(Candidate *candidate) const;
	CUDA_CALLABLE_MEMBER inline void reject(ref_ptr<Candidate> candidate) const {
		reject(candidate.get());
	}

	CUDA_CALLABLE_MEMBER void accept(Candidate *candidate) const;
	CUDA_CALLABLE_MEMBER inline void accept(ref_ptr<Candidate> candidate) const {
		accept(candidate.get());
	}

public:
	AbstractCondition();
	virtual ~AbstractCondition();
	CUDA_CALLABLE_MEMBER void onReject(Module *rejectAction);
	CUDA_CALLABLE_MEMBER void onAccept(Module *acceptAction);
	CUDA_CALLABLE_MEMBER void setMakeRejectedInactive(bool makeInactive);
	CUDA_CALLABLE_MEMBER void setMakeAcceptedInactive(bool makeInactive);
	void setRejectFlag(std::string key, std::string value);
	void setAcceptFlag(std::string key, std::string value);

	// return the reject flag (key & value), delimiter is the "&".
	std::string getRejectFlag();

	// return the accept flag (key & value), delimiter is the "&"
	std::string getAcceptFlag();
};

/**
 @class Deactivation
 @brief Direct deactivation of the candidate. Can be used for debuging.
*/
class Deactivation: public AbstractCondition {
	public: 
		CUDA_CALLABLE_MEMBER void process(Candidate *cand) const { reject(cand); }
};


} // namespace crpropa

#endif /* CRPROPA_MODULE_H */
