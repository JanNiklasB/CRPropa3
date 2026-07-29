#ifndef CRPROPA_PROPAGATIONCK_H
#define CRPROPA_PROPAGATIONCK_H

#include "crpropa/module/Propagation.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"
#include "kiss/logger.h"

namespace crpropa {
/**
 * \addtogroup Propagation 
 * @{
 */

/**
 @class PropagationCK
 @brief Rectilinear propagation through magnetic fields using the Cash-Karp method.

 This module solves the equations of motion of a relativistic charged particle when propagating through a magnetic field.\n
 It uses the Runge-Kutta integration method with Cash-Karp coefficients.\n
 The step size control tries to keep the relative error close to, but smaller than the designated tolerance.
 Additionally a minimum and maximum size for the steps can be set.
 For neutral particles a rectilinear propagation is applied and a next step of the maximum step size proposed.
 */
class PropagationCK: public Propagation {
private:
	std::vector<double> a, b, bs; /*< Cash-Karp coefficients */
	ref_ptr<MagneticField> field;
	double tolerance; /*< target relative error of the numerical integration */
	double minStep; /*< minimum step size of the propagation */
	double maxStep; /*< maximum step size of the propagation */

public:
	/** Constructor for the adaptive Kash Carp.
	 * @param field
	 * @param tolerance	 tolerance is criterion for step adjustment. Step adjustment takes place only if minStep < maxStep
	 * @param minStep	   minStep/c_light is the minimum integration time step
	 * @param maxStep	   maxStep/c_light is the maximum integration time step. 
	 */
    PropagationCK(ref_ptr<MagneticField> field = NULL, double tolerance = 1e-4,
			double minStep = (0.1 * kpc), double maxStep = (1 * Gpc));

	void process(Candidate *candidate) const;

	/** Calculates the new position and direction of the particle based on the solution of the Lorentz force
	 * @param pos	current position of the candidate
	 * @param dir	current direction of the candidate
	 * @param dt	current timestep size of the candidate
	 * @param z		current redshift is needed to calculate the magnetic field
	 * @param current current particle state
	 * @return	  return the new calculated position and direction of the candidate 
	 */
	virtual Y dY(Vector3d pos, Vector3d dir, double dt, double z, const ParticleState &current) const override;

	/** Adapt step size if required and calculates the new position and direction of the particle with the usage of the function dY
	 * @param y		 current position and direction of candidate
	 * @param out	 position and direction of candidate after the step
	 * @param error	 error for the current step
	 * @param h		 current step size
	 * @param z		 current red shift
	 * @param p		 current particle state
	 */
	void tryStep(const Y &y, Y &out, Y &error, double h, double z, const ParticleState &p) const override;

	void setField(ref_ptr<MagneticField> field);
	void setTolerance(double tolerance);
	void setMinimumStep(double minStep);
	void setMaximumStep(double maxStep);

	 /** get functions for the parameters of the class PropagationCK, similar to the set functions */
	ref_ptr<MagneticField> getField() const;
	
	/** get magnetic field vector at current candidate position
	 * @param pos   current position of the candidate
	 * @param z	 current redshift is needed to calculate the magnetic field
	 * @return	  magnetic field vector at the position pos */
	Vector3d getFieldAtPosition(Vector3d pos, double z) const;

	double getTolerance() const;
	double getMinimumStep() const;
	double getMaximumStep() const;
	std::string getDescription() const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PROPAGATIONCK_H
