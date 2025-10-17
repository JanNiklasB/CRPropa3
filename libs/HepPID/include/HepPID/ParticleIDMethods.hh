// ----------------------------------------------------------------------
//
// ParticleIDMethods.hh
// Author:  Lynn Garren
//
//  various utilities to extract information from the particle ID
//
//  In the standard numbering scheme, the PID digits (base 10) are:
//            +/- n nr nl nq1 nq2 nq3 nj
//  It is expected that any 7 digit number used as a PID will adhere to 
//  the Monte Carlo numbering scheme documented by the PDG.
//  Note that many "new" particles not explicitly defined already 
//  can be expressed within this numbering scheme. 
//
//  These are the same methods that can be found in HepPDT::ParticleID
// ----------------------------------------------------------------------
#ifndef PARTICLE_ID_METHODS_HH
#define PARTICLE_ID_METHODS_HH

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

namespace HepPID {

///  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj
///  The location enum provides a convenient index into the PID.
enum location { nj=1, nq3, nq2, nq1, nl, nr, n, n8, n9, n10 };

/// return the digit at a named location in the PID
CUDA_CALLABLE_MEMBER unsigned short digit( location loc, const int & pid );

/// if this is a nucleus (ion), get A
/// Ion numbers are +/- 10LZZZAAAI. 
CUDA_CALLABLE_MEMBER int A(const int & pid );

/// if this is a nucleus (ion), get Z
/// Ion numbers are +/- 10LZZZAAAI. 
CUDA_CALLABLE_MEMBER int Z(const int & pid );

/// if this is a nucleus (ion), get nLambda
/// Ion numbers are +/- 10LZZZAAAI. 
CUDA_CALLABLE_MEMBER int lambda( const int & pid );

/// absolute value of particle ID
CUDA_CALLABLE_MEMBER int abspid( const int & pid );

/// extract fundamental ID (1-100) if this is a "fundamental" particle
CUDA_CALLABLE_MEMBER int fundamentalID( const int & pid );
/// if this is a fundamental particle, does it have a valid antiparticle?
bool hasFundamentalAnti( const int & pid );

/// returns everything beyond the 7th digit 
/// (e.g. outside the standard numbering scheme)
CUDA_CALLABLE_MEMBER int extraBits( const int & pid );

// ---  boolean methods:
//
/// is this a valid ID?
bool isValid( const int & pid );
/// is this a valid meson ID?
CUDA_CALLABLE_MEMBER bool isMeson( const int & pid );
/// is this a valid baryon ID?
CUDA_CALLABLE_MEMBER bool isBaryon( const int & pid );
/// is this a valid diquark ID?
CUDA_CALLABLE_MEMBER bool isDiQuark( const int & pid );
/// is this a valid hadron ID?
CUDA_CALLABLE_MEMBER bool isHadron( const int & pid );
/// is this a valid lepton ID?
CUDA_CALLABLE_MEMBER bool isLepton( const int & pid );
/// is this a valid ion ID?
CUDA_CALLABLE_MEMBER bool isNucleus( const int & pid );
/// is this a valid pentaquark ID?
CUDA_CALLABLE_MEMBER bool isPentaquark( const int & pid );
/// is this a valid SUSY ID?
CUDA_CALLABLE_MEMBER bool isSUSY( const int & pid );
/// is this a valid R-hadron ID?
CUDA_CALLABLE_MEMBER bool isRhadron( const int & pid );
/// is this a valid Dyon (magnetic monopole) ID?
CUDA_CALLABLE_MEMBER bool isDyon( const int & pid );
/// Check for QBall or any exotic particle with electric charge beyond the qqq scheme
/// Ad-hoc numbering for such particles is 100xxxx0, where xxxx is the charge in tenths. 
CUDA_CALLABLE_MEMBER bool isQBall( const int & pid );

/// does this particle contain an up quark?
CUDA_CALLABLE_MEMBER bool hasUp( const int & pid );
/// does this particle contain a down quark?
CUDA_CALLABLE_MEMBER bool hasDown( const int & pid );
/// does this particle contain a strange quark?
CUDA_CALLABLE_MEMBER bool hasStrange( const int & pid );
/// does this particle contain a charm quark?
CUDA_CALLABLE_MEMBER bool hasCharm( const int & pid );
/// does this particle contain a bottom quark?
CUDA_CALLABLE_MEMBER bool hasBottom( const int & pid );
/// does this particle contain a top quark?
CUDA_CALLABLE_MEMBER bool hasTop( const int & pid );

// ---  other information
//
/// jSpin returns 2J+1, where J is the total spin
CUDA_CALLABLE_MEMBER int  jSpin( const int & pid );
/// sSpin returns 2S+1, where S is the spin
CUDA_CALLABLE_MEMBER int  sSpin( const int & pid );
/// lSpin returns 2L+1, where L is the orbital angular momentum
CUDA_CALLABLE_MEMBER int  lSpin( const int & pid );
/// return 3 times the charge (3 x quark charge is an int)
/// If this is a Q-ball, return 30 times the charge.
CUDA_CALLABLE_MEMBER int threeCharge( const int & pid );
/// return the actual charge
CUDA_CALLABLE_MEMBER double charge( const int & pid );


} // HepPID

#endif // PARTICLE_ID_METHODS_HH
