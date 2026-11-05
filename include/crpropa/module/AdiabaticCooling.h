#ifndef CRPROPA_ADIABATICCOOLING_H
#define CRPROPA_ADIABATICCOOLING_H

#include "crpropa/Module.h"
#include "crpropa/advectionField/AdvectionField.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
@class AdiabaticCooling
@brief Implements adiabatic cooling/heating due to advection.
*/

class AdiabaticCooling: public Module {
private:
	ref_ptr<AdvectionField> advectionField;
	double limit;

public:
	/** Default constructor.
	 @param advectionField 	The advection field used for the adiabatic energy change
	 */
	AdiabaticCooling(ref_ptr<AdvectionField> advectionField);
	/** Constructor
	 @param advectionField 	The advection field used for the adiabatic energy change
	 @param limit 			Maximum relative energy change allowed
	 */
	AdiabaticCooling(ref_ptr<AdvectionField> advectionField, double limit);
	void process(Candidate *c) const;

	void setLimit(double l);

	double getLimit() const;

};
/** @}*/



}; // end namesspace crpropa
#endif // CRPROPA_ADIABATICCOOLING_H
