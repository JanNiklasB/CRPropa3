#ifndef CRPROPA_UNITS_H
#define CRPROPA_UNITS_H

#include "__CudaDefines.h"

#include <cmath>

namespace crpropa {
/**
 * \addtogroup PhysicsDefinitions
 * @{
 */

/**
 @file
 @brief Definition of SI base units and constants

 Definition of SI base units and constants used elsewhere in the code
 Based on:
 - CODATA recommended values of the fundamental physical constants: 2006
 	doi:10.1103/RevModPhys.80.633
 - IAU 2012 Resolution B2, IAU 2015 Resolution B2
 	https://www.iau.org/administration/resolutions/
 */

// SI units
CUDA_CONSTANT static const double meter = 1;
CUDA_CONSTANT static const double second = 1;
CUDA_CONSTANT static const double kilogram = 1;
CUDA_CONSTANT static const double ampere = 1;
CUDA_CONSTANT static const double mol = 1;
CUDA_CONSTANT static const double kelvin = 1;

// derived units
CUDA_CONSTANT static const double newton = 1 * kilogram * meter / second / second;
CUDA_CONSTANT static const double pascal = 1 * newton / meter / meter;
CUDA_CONSTANT static const double joule = 1 * newton * meter;
CUDA_CONSTANT static const double tesla = 1 * newton / ampere / meter;
CUDA_CONSTANT static const double volt = 1 * kilogram * meter * meter / ampere / second / second / second;
CUDA_CONSTANT static const double coulomb = 1 * ampere * second;
CUDA_CONSTANT static const double hertz = 1 / second;
CUDA_CONSTANT static const double rad = 1;
CUDA_CONSTANT static const double deg = M_PI / 180.;

// SI Prefixes
CUDA_CONSTANT static const double yocto = 1E-24;
CUDA_CONSTANT static const double zepto = 1E-21;
CUDA_CONSTANT static const double atto = 1E-18;
CUDA_CONSTANT static const double femto = 1E-15;
CUDA_CONSTANT static const double pico = 1E-12;
CUDA_CONSTANT static const double nano = 1E-9;
CUDA_CONSTANT static const double micro = 1E-6;
CUDA_CONSTANT static const double milli = 1E-3;

CUDA_CONSTANT static const double kilo = 1E3;
CUDA_CONSTANT static const double mega = 1E6;
CUDA_CONSTANT static const double giga = 1E9;
CUDA_CONSTANT static const double tera = 1E12;
CUDA_CONSTANT static const double peta = 1E15;
CUDA_CONSTANT static const double exa = 1E18;
CUDA_CONSTANT static const double zetta = 1E21;
CUDA_CONSTANT static const double yotta = 1E24;


// physical constants
CUDA_CONSTANT static const double eplus = 1.602176487e-19 * ampere * second;
CUDA_CONSTANT static const double c_light = 2.99792458e8 * meter / second;
CUDA_CONSTANT static const double c_squared = c_light * c_light;
CUDA_CONSTANT static const double amu = 1.660538921e-27 * kilogram;
CUDA_CONSTANT static const double mass_proton = 1.67262158e-27 * kilogram;
CUDA_CONSTANT static const double mass_neutron = 1.67492735e-27 * kilogram;
CUDA_CONSTANT static const double mass_electron = 9.10938291e-31 * kilogram;
CUDA_CONSTANT static const double h_planck = 6.62606957e-34 * joule * second;
CUDA_CONSTANT static const double k_boltzmann = 1.3806488e-23 * joule / kelvin;
CUDA_CONSTANT static const double mu0 = 4 * M_PI * 1e-7 * newton / ampere / ampere;
CUDA_CONSTANT static const double epsilon0 = 1.0 / mu0 / c_squared * ampere * second / volt / meter;
CUDA_CONSTANT static const double alpha_finestructure = eplus * eplus / 2. / epsilon0 / h_planck / c_light;
CUDA_CONSTANT static const double radius_electron = eplus * eplus / 4. / M_PI / epsilon0 / mass_electron / c_squared;
CUDA_CONSTANT static const double sigma_thomson = 8. * M_PI / 3. * radius_electron * radius_electron;

// gauss
CUDA_CONSTANT static const double gauss = 1e-4 * tesla;
CUDA_CONSTANT static const double microgauss = 1e-6 * gauss;
CUDA_CONSTANT static const double nanogauss = 1e-9 * gauss;
CUDA_CONSTANT static const double muG = microgauss;
CUDA_CONSTANT static const double nG = nanogauss;

CUDA_CONSTANT static const double erg = 1E-7 * joule;

// electron volt
CUDA_CONSTANT static const double electronvolt = eplus * volt;
CUDA_CONSTANT static const double kiloelectronvolt = 1e3 * electronvolt;
CUDA_CONSTANT static const double megaelectronvolt = 1e6 * electronvolt;
CUDA_CONSTANT static const double gigaelectronvolt = 1e9 * electronvolt;
CUDA_CONSTANT static const double teraelectronvolt = 1e12 * electronvolt;
CUDA_CONSTANT static const double petaelectronvolt = 1e15 * electronvolt;
CUDA_CONSTANT static const double exaelectronvolt = 1e18 * electronvolt;
CUDA_CONSTANT static const double eV = electronvolt;
CUDA_CONSTANT static const double keV = kiloelectronvolt;
CUDA_CONSTANT static const double MeV = megaelectronvolt;
CUDA_CONSTANT static const double GeV = gigaelectronvolt;
CUDA_CONSTANT static const double TeV = teraelectronvolt;
CUDA_CONSTANT static const double PeV = petaelectronvolt;
CUDA_CONSTANT static const double EeV = exaelectronvolt;

CUDA_CONSTANT static const double barn = 1E-28 * meter * meter;

// astronomical distances
CUDA_CONSTANT static const double au = 149597870700 * meter;
CUDA_CONSTANT static const double ly = 365.25 * 24 * 3600 * second * c_light;
CUDA_CONSTANT static const double parsec = 648000 / M_PI * au;
CUDA_CONSTANT static const double kiloparsec = 1e3 * parsec;
CUDA_CONSTANT static const double megaparsec = 1e6 * parsec;
CUDA_CONSTANT static const double gigaparsec = 1e9 * parsec;
CUDA_CONSTANT static const double pc = parsec;
CUDA_CONSTANT static const double kpc = kiloparsec;
CUDA_CONSTANT static const double Mpc = megaparsec;
CUDA_CONSTANT static const double Gpc = gigaparsec;

// meter
CUDA_CONSTANT static const double kilometer = 1000 * meter;
CUDA_CONSTANT static const double centimeter = 0.01 * meter;
CUDA_CONSTANT static const double km = kilometer;
CUDA_CONSTANT static const double cm = centimeter;

// second
CUDA_CONSTANT static const double nanosecond = 1e-9 * second;
CUDA_CONSTANT static const double microsecond = 1e-6 * second;
CUDA_CONSTANT static const double millisecond = 1e-3 * second;
CUDA_CONSTANT static const double minute = 60 * second;
CUDA_CONSTANT static const double hour = 3600 * second;
CUDA_CONSTANT static const double day = 24 * hour;
CUDA_CONSTANT static const double year = 365.25 * 24 * hour;
CUDA_CONSTANT static const double kiloyear = 1e3 * year;
CUDA_CONSTANT static const double Megayear = 1e6 * year;
CUDA_CONSTANT static const double Gigayear = 1e9 * year;
CUDA_CONSTANT static const double ns = nanosecond;
CUDA_CONSTANT static const double mus = microsecond;
CUDA_CONSTANT static const double ms = millisecond;
CUDA_CONSTANT static const double sec = second;
CUDA_CONSTANT static const double yr = year;
CUDA_CONSTANT static const double kyr = kiloyear;
CUDA_CONSTANT static const double Myr = Megayear;
CUDA_CONSTANT static const double Gyr = Gigayear;

// volume
CUDA_CONSTANT static const double ccm = cm*cm*cm;

/** @}*/

} // namespace crpropa

#endif // CRPROPA_UNITS_H
