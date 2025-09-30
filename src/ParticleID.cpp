#include "crpropa/ParticleID.h"

#include "HepPID/ParticleName.hh"
#include "kiss/convert.h"

#include <string>
#include <cmath>
#include <iostream>

namespace CudaHepPID{

// Author:  Lynn Garren
// the following functions in namespace CudaHepPID were directly taken from the
// HepPID library by Lynn Garren and then slightly modified for cuda use
CUDA_CALLABLE_MEMBER enum location { nj=1, nq3, nq2, nq1, nl, nr, n, n8, n9, n10 };
CUDA_CALLABLE_MEMBER unsigned short digit( location loc, const int & pid );
CUDA_CALLABLE_MEMBER int A(const int & pid );
CUDA_CALLABLE_MEMBER int Z(const int & pid );
CUDA_CALLABLE_MEMBER bool isNucleus( const int & pid );

// Ion numbers are +/- 10LZZZAAAI. 
int Z( const int & pid )
{
	// a proton can also be a Hydrogen nucleus
	if( abs(pid) == 2212 ) { return 1; }
	if( isNucleus(pid) ) return (abs(pid)/10000)%1000;
	return 0;
}

// Ion numbers are +/- 10LZZZAAAI. 
int A( const int & pid )
{
	// a proton can also be a Hydrogen nucleus
	if( abs(pid) == 2212 ) { return 1; }
	if( isNucleus(pid) ) return (abs(pid)/10)%1000;
	return 0;
}

bool isNucleus( const int & pid )
{
	// a proton can also be a Hydrogen nucleus
	if( abs(pid) == 2212 ) { return true; }
	// new standard: +/- 10LZZZAAAI
	if( ( digit(n10,pid) == 1 ) && ( digit(n9,pid) == 0 ) ) {
	// charge should always be less than or equal to baryon number
	// the following line is A >= Z
	if( (abs(pid)/10)%1000 >= (abs(pid)/10000)%1000 ) { return true; }
	}
	return false;
}

//  split the PID into constituent integers
unsigned short digit( location loc, const int & pid )
{
	//  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj
	//  the location enum provides a convenient index into the PID
		//
		//  Modified for CRPropa: use precalculated values isntead of pow for
		//  performance
		static unsigned int p10[] = { 1, 10, 100, 1000,  10000, 100000, 1000000,
			10000000, 100000000, 1000000000};
	return (abs(pid)/ p10[loc-1])%10;
	//    int numerator = (int) std::pow(10.0,(loc-1));
	//    return (abspid(pid)/numerator)%10;
}

}

namespace crpropa {

int nucleusId(int a, int z) {
	#ifndef __CUDACC__
	if (z < 0)
		throw std::runtime_error(
				"crpropa::Nucleus: no nucleus with Z < 0, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	if (a < 1)
		throw std::runtime_error(
				"crpropa::Nucleus: no nucleus with A < 1, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	if (a < z)
		throw std::runtime_error(
				"crpropa::Nucleus: no nucleus with A < Z, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	#endif
	return 1000000000 + z * 10000 + a * 10;
}

int chargeNumber(int id) {
	return CudaHepPID::Z(id);
}

int massNumber(int id) {
	if (id == 2112)
		return 1;
	return CudaHepPID::A(id);
}

bool isNucleus(int id) {
	if (id == 2112)
		return true; // consider neutron as nucleus
	return CudaHepPID::isNucleus(id);
}

std::string convertIdToName(int id) {
	// handle a few extra cases that HepPID doesn't like
	if (id == 1000000010) // neutron
		id = 2112;
	if (id == -1000000010) // anti-neutron
		id = -2112;
	if (id == -1000010010) // anti-proton
		id = -2212;
	return HepPID::particleName(id);
}

}
