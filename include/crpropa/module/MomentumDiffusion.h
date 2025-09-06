#ifndef CRPROPA_MOMENTUMDIFFUSION_H
#define CRPROPA_MOMENTUMDIFFUSION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <stdexcept>

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include "kiss/logger.h"

namespace crpropa {

/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class ConstantMomentumDiffusion
 * Simplest model for diffusion in momentum space
 */

class ConstantMomentumDiffusion: public Module {

private:
	double Dpp; // Diffusion coefficient
	double limit; // maximal fractional energy loss

public:
	/** Constructor
	@param Dpp 	momentum diffusion coefficient
	*/
	ConstantMomentumDiffusion(double Dpp);

	/** Constructor
	@param Dpp 		momentum diffusion coefficient
	@param limit 	maximal fractional energy loss
	*/
	ConstantMomentumDiffusion(double Dpp, double limit);

	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
	CUDA_CALLABLE_MEMBER double calculateAScalar(double p) const;
	CUDA_CALLABLE_MEMBER double calculateBScalar() const;

	void setLimit(double l);
	void setDpp(double Dpp);

	double getLimit() const;
	double getDpp() const;

	std::string getDescription() const;

};

/** @}*/

}; //end namespace crpropa

#endif // CRPROPA_MOMENTUMDIFFUSION_H
