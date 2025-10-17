#ifndef SIMPLEPROPAGATION_H
#define SIMPLEPROPAGATION_H

#include "crpropa/__CudaDefines.h"
#include "crpropa/Module.h"
#include "crpropa/Units.h"

namespace crpropa {
/**
 * \addtogroup Propagation 
 * @{
 */

/**
 @class SimplePropagation
 @brief Simple rectilinear propagation in absence of magnetic fields.

 This module implements rectilinear propagation.
 The step size is guaranteed to be larger than minStep and smaller than maxStep.
 It always proposes a next step size of maxStep.
 */
class SimplePropagation: public Module {
private:
	double minStep, maxStep;

public:
	CUDA_CALLABLE_MEMBER SimplePropagation(){}
	SimplePropagation(double minStep = (0.1 * kpc), double maxStep = (1 * Gpc));

	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
	CUDA_CALLABLE_MEMBER void test() const{
		printf("SimplePropagation::test()\n");
	}
	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);
	double getMinimumStep() const;
	double getMaximumStep() const;
	std::string getDescription() const;
};
/** @}*/

} // namespace crpropa

#endif // SIMPLEPROPAGATION_H

