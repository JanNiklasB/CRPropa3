#ifndef CRPROPA_ELASTICSCATTERING_H
#define CRPROPA_ELASTICSCATTERING_H

#include "crpropa/__CudaDefines.h"
#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Units.h"

#include <vector>
#include <cmath>

namespace crpropa {

/**
 @class ElasticScattering
 @brief Elastic scattering of background photons on cosmic-ray nuclei.
 */
class ElasticScattering: public Module {
private:
	ref_ptr<PhotonField> photonField;

	std::vector<double> tabRate; // elastic scattering rate
	std::vector<std::vector<double> > tabCDF; // CDF as function of background photon energy
	std::string interactionTag = "ES";

	const double lgmin = 6.;  // minimum log10(Lorentz-factor)
	const double lgmax = 14.; // maximum log10(Lorentz-factor)
	const size_t nlg = 201;   // number of Lorentz-factor steps
	const double epsmin = log10(2 * eV) + 3;    // log10 minimum photon background energy in nucleus rest frame for elastic scattering
	const double epsmax = log10(2 * eV) + 8.12; // log10 maximum photon background energy in nucleus rest frame for elastic scattering
	const size_t neps = 513; // number of photon background energies in nucleus rest frame

public:
	/** Constructor
	 @param photonField		target photon field
	 */
	ElasticScattering(ref_ptr<PhotonField> photonField);
	void initRate(std::string filename);
	void initCDF(std::string filename);
	void setPhotonField(ref_ptr<PhotonField> photonField);
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
	
	std::string getInteractionTag() const;
	void setInteractionTag(std::string tag);
};

} // namespace crpropa

#endif // CRPROPA_ELASTICSCATTERING_H
