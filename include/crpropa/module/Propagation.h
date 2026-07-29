#ifndef CRPROPA_PROPAGATION_H
#define CRPROPA_PROPAGATION_H

#include "crpropa/Module.h"

namespace crpropa {
/**
 * \addtogroup Propagation
 * @{
 */

/**
 @class Propagation
 @brief Base class for all kinetic propagators
 */
class Propagation: public Module {

public:
	class Y {
	public:
		Vector3d x, u; /*< phase-point: position and direction */

		Y() {
		}

		Y(const Vector3d &x, const Vector3d &u) :
				x(x), u(u) {
		}

		Y(double f) :
				x(Vector3d(f, f, f)), u(Vector3d(f, f, f)) {
		}

		Y operator *(double f) const {
			return Y(x * f, u * f);
		}

		Y &operator +=(const Y &y) {
			x += y.x;
			u += y.u;
			return *this;
		}
	};

	/** Calculates the new position and direction of the particle based on the solution of the Lorentz force
	 * @param pos	current position of the candidate
	 * @param dir	current direction of the candidate
	 * @param dt	current timestep/lengthstep size of the candidate
	 * @param z		current redshift is needed to calculate the magnetic field
	 * @param current current particle state
	 * @return	  return the new calculated position and direction of the candidate 
	 */
	virtual Y dY(Vector3d pos, Vector3d dir, double dt, double z, const ParticleState &current) const = 0;

	/** Adapt step size if required and calculates the new position and direction of the particle with the usage of the function dY
	 * @param y		 current position and direction of candidate
	 * @param out	 position and direction of candidate after the step
	 * @param error	 error for the current step
	 * @param dt	 current timestep/lengthstep size of the candidate
	 * @param z		 current red shift
	 * @param p		 current particle state
	 */
	virtual void tryStep(const Y &y, Y &out, Y &error, double dt, double z, const ParticleState &p) const = 0;
};
/** @}*/

} // namespace crpropa

#endif // PROPAGATION_H
