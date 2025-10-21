#ifndef CRPROPA_COSMOLOGY_H
#define CRPROPA_COSMOLOGY_H

#include "crpropa/__CudaDefines.h"
#include "Referenced.h"
#include <vector>

namespace crpropa {
/**
 * \addtogroup PhysicsDefinitions
 * @{
 */

/**
 @class Cosmology
 @brief Cosmology calculations
*/
struct Cosmology : public Referenced {
	private:
	
	double H0; // Hubble parameter at z=0
	double omegaM; // matter density parameter
	double omegaL; // vacuum energy parameter

	const size_t n = 1000;
	const double zmin = 0.0001;
	const double zmax = 100;

	double* Z =NULL;  // redshift
	double* Dc=NULL; // comoving distance [m]
	double* Dl=NULL; // luminosity distance [m]
	double* Dt=NULL; // light travel distance [m]

	size_t ZSize=0, DcSize=0, DlSize=0, DtSize=0;

	public:

	CUDA_CALLABLE_MEMBER Cosmology();
	CUDA_CALLABLE_MEMBER ~Cosmology();

	CUDA_CALLABLE_MEMBER void update();

	/**
	 Set the cosmological parameters for a flat universe. To ensure flatness omegaL is set to 1 - omegaMatter
	@param hubbleParameter	dimensionless Hubble parameter, default = 0.673
	@param omegaMatter		matter parameter, default = 0.315
	*/
	void setParameters(double hubbleParameter, double omegaMatter);

	/**
	 Hubble rate at given redshift
	H(z) = H0 * sqrt(omegaM * (1 + z)^3 + omegaL)
	*/
	CUDA_CALLABLE_MEMBER double hubbleRate(double redshift = 0) const;

	// Returns the dark energy density parameter
	CUDA_CALLABLE_MEMBER double getOmegaL() const;

	// Returns the matter density parameter
	CUDA_CALLABLE_MEMBER double getOmegaM() const;

	// Returns the hubble parameter
	CUDA_CALLABLE_MEMBER double getH0() const;

	/**
	 Redshift of a comoving object at a given comoving distance to an observer at z = 0.
	d_comoving(z) = c/H0 * int_0^z dz' / E(z')
	*/
	CUDA_CALLABLE_MEMBER double comovingDistance2Redshift(double distance) const;

	/**
	 Comoving distance between an observer at z = 0 and a comoving object at z.
	d_comoving(z) = c/H0 * int_0^z dz' / E(z')
	*/
	CUDA_CALLABLE_MEMBER double redshift2ComovingDistance(double redshift) const;

	/**
	 Redshift of a comoving object at a given luminosity distance to an observer at z = 0.
	d_luminosity(z) = (1 + z) * d_comoving(z)
	*/
	CUDA_CALLABLE_MEMBER double luminosityDistance2Redshift(double distance) const;

	/**
	 Luminosity distance between an observer at z = 0 and a comoving object at z.
	d_luminosity(z) = (1 + z) * d_comoving(z)
	*/
	CUDA_CALLABLE_MEMBER double redshift2LuminosityDistance(double redshift) const;

	/**
	 Redshift of a comoving object at a given light travel distance to an observer at z = 0.
	d_lighttravel(z) = c/H0 * int_0^z dz' / ((1 + z')  *  E(z'))
	*/
	CUDA_CALLABLE_MEMBER double lightTravelDistance2Redshift(double distance) const;

	/**
	 Light travel distance between an observer at z = 0 and a comoving object at z.
	d_lighttravel(z) = c/H0 * int_0^z dz' / ((1 + z')  *  E(z'))
	*/
	CUDA_CALLABLE_MEMBER double redshift2LightTravelDistance(double redshift) const;

	// Conversion from comoving distance to light travel distance.
	CUDA_CALLABLE_MEMBER double comoving2LightTravelDistance(double distance) const;

	// Conversion from light travel distance to comoving distance.
	CUDA_CALLABLE_MEMBER double lightTravel2ComovingDistance(double distance) const;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_COSMOLOGY_H
