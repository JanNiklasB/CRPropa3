#ifndef CRPROPA_CANDIDATE_H
#define CRPROPA_CANDIDATE_H

#include "crpropa/__CudaDefines.h"
#include "crpropa/ParticleState.h"
#include "crpropa/Referenced.h"
#include "crpropa/AssocVector.h"
#include "crpropa/Variant.h"

#include <vector>
#include <map>
#include <sstream>
#include <stdint.h>

#ifdef __CUDACC__
#include "thrust/device_vector.h"
#endif

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/**
 @class Candidate Candidate.h include/crpropa/Candidate.h
 @brief All information about the cosmic ray.

 The Candidate is a passive object, that holds the information about the state
 of the cosmic ray and the simulation itself.
 */
class Candidate: public Referenced {
public:
	ParticleState source; /**< Particle state at the source */
	ParticleState created; /**< Particle state of parent particle at the time of creation */
	ParticleState current; /**< Current particle state */
	ParticleState previous; /**< Particle state at the end of the previous step */
	NuclearMassTable* NuclearMassPtr=NULL;  ///< initialized in ModuleList.h to save storage and enable GPU support

	Candidate** secondaries=NULL; /**< Secondary particles from interactions */
	std::size_t secondariesSize=0;

	#ifdef __CUDACC__
	typedef AssocVector<const char*, Variant> PropertyMap;
	#else
	typedef Loki::AssocVector<const char*, Variant> PropertyMap;
	#endif
	PropertyMap properties; /**< Map of property names and their values. */

	/** Parent candidate. 0 if no parent (initial particle). Must not be a ref_ptr to prevent circular referencing. */
	Candidate *parent=NULL;	
private:
	bool active; /**< Active status */
	double weight; /**< Weight of the candidate */
	double redshift; /**< Current simulation time-point in terms of redshift z */
	double trajectoryLength; /**< Comoving distance [m] the candidate has traveled so far */
	double currentStep; /**< Size of the currently performed step in [m] comoving units */
	double nextStep; /**< Proposed size of the next propagation step in [m] comoving units */
	const char* tagOrigin; /**< Name of interaction/source process which created this candidate*/
	double time; /**< Time [s] that has passed in the laboratory frame of reference */

	uint64_t* nextSerialNumber=NULL;
	uint64_t serialNumber;

public:
	CUDA_CALLABLE_MEMBER Candidate();

	/** Creates a candidate
	 * \param id  Particle ID
	 * \param energy Particle Energy
	 * \param position Initial position
	 * \param direction Initial direction
	 * \param z Current simulation time-point in terms of redshift z
	 * \param weight Weight of candidate
	 * \param tagOrigin Name of interaction/source process which created this candidate
	 */
	Candidate(
		int id,
		double energy = 0,
		Vector3d position = Vector3d(0, 0, 0),
		Vector3d direction = Vector3d(-1, 0, 0),
		double z = 0,
		double weight = 1., 
		const char* tagOrigin = "PRIM"
	);

	/** Creates a candidate with predefined NuclearMassTable
	 * \param NuclearMassTable Ptr to an instance of NuclearMassTable
	 * \param id  Particle ID
	 * \param energy Particle Energy
	 * \param position Initial position
	 * \param direction Initial direction
	 * \param z Current simulation time-point in terms of redshift z
	 * \param weight Weight of candidate
	 * \param tagOrigin Name of interaction/source process which created this candidate
	 */
	CUDA_CALLABLE_MEMBER Candidate(
		NuclearMassTable* NuclearMassTable,
		int id = 0,
		double energy = 0,
		Vector3d position = Vector3d(0, 0, 0),
		Vector3d direction = Vector3d(-1, 0, 0),
		double z = 0,
		double weight = 1., 
		const char* tagOrigin = "PRIM"
	);

	/**
	 Creates a candidate, initializing the Candidate::source, Candidate::created,
	 Candidate::previous and Candidate::current state with the argument.
	 */
	CUDA_CALLABLE_MEMBER Candidate(const ParticleState &state);

	CUDA_CALLABLE_MEMBER ~Candidate();

	CUDA_CALLABLE_MEMBER bool isActive() const;
	CUDA_CALLABLE_MEMBER void setActive(bool b);

	CUDA_CALLABLE_MEMBER int getSecondarySize() const;

	CUDA_CALLABLE_MEMBER void setTrajectoryLength(double length);
	CUDA_CALLABLE_MEMBER double getTrajectoryLength() const;
	
	CUDA_CALLABLE_MEMBER double getVelocity() const;

	CUDA_CALLABLE_MEMBER void setRedshift(double z);
	CUDA_CALLABLE_MEMBER double getRedshift() const;

	/**
	 Sets weight of each candidate.
	 Weights are calculated for each tracked secondary.
	 */
	CUDA_CALLABLE_MEMBER void setWeight(double weight);
	CUDA_CALLABLE_MEMBER void updateWeight(double weight);
	CUDA_CALLABLE_MEMBER double getWeight() const;

	/**
	 Sets the current step and increases the trajectory length accordingly.
	 Only the propagation module should use this.
	 */
	CUDA_CALLABLE_MEMBER void setCurrentStep(double step);
	CUDA_CALLABLE_MEMBER double getCurrentStep() const;

	/**
	 Sets the proposed next step.
	 Only the propagation module should use this.
	 */
	CUDA_CALLABLE_MEMBER void setNextStep(double step);
	CUDA_CALLABLE_MEMBER double getNextStep() const;

	/**
	 Sets the tagOrigin of the candidate. Can be used to trace back the interactions
	 */
	CUDA_CALLABLE_MEMBER void setTagOrigin(const char* tagOrigin);
	void setTagOrigin (const std::string& tagOrigin);
	CUDA_CALLABLE_MEMBER const char* getTagOrigin() const;

	/**
	 Sets the time of the candidate.
	 */
	CUDA_CALLABLE_MEMBER void setTime(double t);
	CUDA_CALLABLE_MEMBER double getTime() const;

	/**
	 Make a bid for the next step size: the lowest wins.
	 */
	CUDA_CALLABLE_MEMBER void limitNextStep(double step);

	/// This function sets a property and creates one if it does not exist
	void setProperty(const std::string &name, const Variant &value);
	/// This functions sets a property only if it exists
	CUDA_CALLABLE_MEMBER void setProperty(const char* name, const Variant &value);
	CUDA_CALLABLE_MEMBER const Variant &getProperty(const char* name) const;
	const Variant &getProperty(const std::string &name) const;
	bool removeProperty(const std::string &name);
	bool hasProperty(const std::string &name) const;
	CUDA_CALLABLE_MEMBER bool hasProperty(const char* name) const;

	/**
	 Add a new candidate to the list of secondaries.
	 @param c Candidate

	 Adds a new candidate to the list of secondaries of this candidate.
	 The secondaries Candidate::source and Candidate::previous state are set to the _source_ and _previous_ state of its parent.
	 The secondaries Candidate::created and Candidate::current state are set to the _current_ state of its parent, except for the secondaries current energy and particle id.
	 Trajectory length and redshift are copied from the parent.
	 */
	CUDA_CALLABLE_MEMBER void addSecondary(Candidate *c);
	CUDA_CALLABLE_MEMBER inline void addSecondary(ref_ptr<Candidate> c) { addSecondary(c.get()); };
	/**
	 Add a new candidate to the list of secondaries.
	 @param id			particle ID of the secondary
	 @param energy		energy of the secondary
	 @param w			weight of the secondary
	 @param tagOrigin 	tag of the secondary
	 */
	CUDA_CALLABLE_MEMBER void addSecondary(int id, double energy, double w = 1., const char* tagOrigin = "SEC");
	/**
	 Add a new candidate to the list of secondaries.
	 @param id			particle ID of the secondary
	 @param energy		energy of the secondary
	 @param position	start position of the secondary
	 @param w			weight of the secondary
	 @param tagOrigin 	tag of the secondary
	 */
	CUDA_CALLABLE_MEMBER void addSecondary(int id, double energy, Vector3d position, double w = 1., const char* tagOrigin = "SEC");
	CUDA_CALLABLE_MEMBER void clearSecondaries();

	std::string getDescription() const;

	/** Unique (inside process) serial number (id) of candidate */
	CUDA_CALLABLE_MEMBER uint64_t getSerialNumber() const;
	CUDA_CALLABLE_MEMBER void setSerialNumber(const uint64_t snr);

	/** Serial number of candidate at source*/
	CUDA_CALLABLE_MEMBER uint64_t getSourceSerialNumber() const;

	/** Serial number of candidate at creation */
	CUDA_CALLABLE_MEMBER uint64_t getCreatedSerialNumber() const;

	/** Set the next serial number to use */
	CUDA_CALLABLE_MEMBER static void setNextSerialNumber(uint64_t snr);

	/** Get the next serial number that will be assigned */
	CUDA_CALLABLE_MEMBER static uint64_t getNextSerialNumber();

	/// Set NuclearMassTablePtr
	CUDA_CALLABLE_MEMBER void setNuclearMassTable(NuclearMassTable* NuclearMassTablePtr);

	/// Get NuclearMassTablePtr
	CUDA_CALLABLE_MEMBER NuclearMassTable* getNuclearMassTable() const;

	CUDA_CALLABLE_MEMBER void copy(const Candidate* Secondary);

	/**
	 Create an exact clone of candidate
	 @param recursive	recursively clone and add the secondaries
	 */
	CUDA_CALLABLE_MEMBER Candidate* clone(bool recursive = false) const;

	/**
	 Copy the source particle state to the current state
	 and activate it if inactive, e.g. restart it
	*/
	CUDA_CALLABLE_MEMBER void restart();
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_CANDIDATE_H
