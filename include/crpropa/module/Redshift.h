#ifndef CRPROPA_REDSHIFT_H
#define CRPROPA_REDSHIFT_H

#include "crpropa/__CudaDefines.h"
#include "crpropa/Module.h"
#include "crpropa/Cosmology.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class Redshift
 @brief Updates redshift and applies adiabatic energy loss according to the traveled distance.
 */
class Redshift: public Module {
private:
	ref_ptr<Cosmology> cosmology;
	
public:
	Redshift();
	Redshift(ref_ptr<Cosmology> cosmology) : cosmology(cosmology) {}

	ref_ptr<Cosmology> getCosmology();
	ref_ptr<Cosmology> setCosmology(ref_ptr<Cosmology> cosmology);

	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class FutureRedshift
 @brief Updates redshift and applies adiabatic energy loss according to the traveled distance. Extends to negative redshift values to allow for symmetric time windows around z=0.
 */
class FutureRedshift: public Module {
private:
	ref_ptr<Cosmology> cosmology;

public:
	FutureRedshift();
	FutureRedshift(ref_ptr<Cosmology> cosmology) : cosmology(cosmology) {}

	ref_ptr<Cosmology> getCosmology();
	ref_ptr<Cosmology> setCosmology(ref_ptr<Cosmology> cosmology);

	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_REDSHIFT_H
