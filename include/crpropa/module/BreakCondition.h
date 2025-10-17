#ifndef CRPROPA_BREAKCONDITION_H
#define CRPROPA_BREAKCONDITION_H

#include "crpropa/__CudaDefines.h"
#include "crpropa/Module.h"

namespace crpropa {
/**
 * \addtogroup Condition 
 * @{
 */

/**
 @class MaximumTrajectoryLength
 @brief Deactivates the candidate beyond a maximum trajectory length

 This module deactivates the candidate at a given maximum trajectory length.
 In that case the property ("Deactivated", module::description) is set.
 It also limits the candidates next step size to ensure the maximum trajectory length is not exceeded.
 */
class MaximumTrajectoryLength: public AbstractCondition {
	double maxLength;
	Vector3d* observerPositions;
	int observerPositionsSize;
public:
	CUDA_CALLABLE_MEMBER MaximumTrajectoryLength();
	MaximumTrajectoryLength(double length);
	~MaximumTrajectoryLength();
	void setMaximumTrajectoryLength(double length);
	double getMaximumTrajectoryLength() const;
	void addObserverPosition(const Vector3d &position);
	const std::vector<Vector3d>& getObserverPositions() const;
	std::string getDescription() const;
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
};

/**
 @class MinimumEnergy
 @brief Deactivates the candidate below a minimum energy

 This module deactivates the candidate below a given minimum energy.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumEnergy: public AbstractCondition {
	double minEnergy;
public:
	CUDA_CALLABLE_MEMBER MinimumEnergy();
	MinimumEnergy(double minEnergy);
	void setMinimumEnergy(double energy);
	double getMinimumEnergy() const;
	std::string getDescription() const;
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
};


/**
 @class MinimumRigidity
 @brief Deactivates the candidate below a minimum rigidity

 This module deactivates the candidate below a given minimum rigidity (E/Z in EeV).
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumRigidity: public AbstractCondition {
	double minRigidity;
public:
	CUDA_CALLABLE_MEMBER MinimumRigidity();
	MinimumRigidity(double minRigidity);
	void setMinimumRigidity(double minRigidity);
	double getMinimumRigidity() const;
	std::string getDescription() const;
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
};

/**
 @class MinimumRedshift
 @brief Deactivates the candidate below a minimum redshift

 This module deactivates the candidate below a given minimum redshift.
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumRedshift: public AbstractCondition {
	double zmin;
public:
	CUDA_CALLABLE_MEMBER MinimumRedshift();
	MinimumRedshift(double zmin);
	void setMinimumRedshift(double z);
	double getMinimumRedshift();
	std::string getDescription() const;
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
};

/**
 @class MinimumChargeNumber
 @brief Deactivates the candidate below a minimum number

 This module deactivates the candidate below a given minimum charge number.
 A minimum charge number of 26 deactivates all (anti-) isotopes which 
 are ranked in the periodic table before iron (Fe). 
 In that case the property ("Deactivated", module::description) is set.
 */
class MinimumChargeNumber: public AbstractCondition {
	int minChargeNumber;
public:
	CUDA_CALLABLE_MEMBER MinimumChargeNumber();
	MinimumChargeNumber(int minChargeNumber = 0);
	void setMinimumChargeNumber(int chargeNumber);
	int getMinimumChargeNumber() const;
	std::string getDescription() const;
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
};

/**
 @class MinimumEnergyPerParticleId
 @brief Deactivates the candidate below a minimum energy for specific particle Ids.

 This module deactivates the candidate below a given minimum energy for specific particle types.
 In that case the property ("Deactivated", module::description) is set.
 All particles whose minimum energy is not specified follow the more general minEnergyOthers condition.
 */
class MinimumEnergyPerParticleId: public AbstractCondition {
	double* minEnergies;
	int minEnergiesSize;
	int* particleIds;
	int particleIdsSize;
	double minEnergyOthers;
public:
	CUDA_CALLABLE_MEMBER MinimumEnergyPerParticleId();
	MinimumEnergyPerParticleId(double minEnergyOthers);
	~MinimumEnergyPerParticleId();
	void setMinimumEnergyOthers(double energy);
	double getMinimumEnergyOthers() const;
	void add(int id, double energy);
	std::string getDescription() const;
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
};


/**
 @class DetectionLength
 @brief Detects the candidate at a given trajectoryLength
 
 This break condition can be used for non-regular time observation of the particle density. See also ObserverTimeEvolution.
 */
class DetectionLength: public AbstractCondition {
	double detLength;
public:
	CUDA_CALLABLE_MEMBER DetectionLength();
	DetectionLength(double length);
	void setDetectionLength(double length);
	double getDetectionLength() const;
	std::string getDescription() const;
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_BREAKCONDITION_H
