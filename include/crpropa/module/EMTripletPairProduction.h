#ifndef CRPROPA_EMTRIPLETPAIRPRODUCTION_H
#define CRPROPA_EMTRIPLETPAIRPRODUCTION_H

#include <fstream>
#include <cmath>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Geometry.h"
#include "crpropa/InteractionRates.h"

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class EMTripletPairProduction
 @brief Electron triplet pair production of electrons with background photons.

 This module simulates electron triplet pair production of electrons with background photons for several photon fields.
 The secondary electrons from this interaction are optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 The surface is defined to include the nodes of the grid contained within.
*/
class EMTripletPairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField;
	bool haveElectrons;
	double limit;
	double thinning;
  ref_ptr<Surface> surface;
	std::string interactionTag = "EMTP";
  ref_ptr<InteractionRates> interactionRates;

  // tabulated interaction rate 1/lambda(E)
	std::vector<double> tabEnergy;  //!< electron energy in [J]
	std::vector<double> tabRate;  //!< interaction rate in [1/m]
	
	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE;  //!< electron energy in [J]
	std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
	std::vector< std::vector<double> > tabCDF;  //!< cumulative interaction rate

public:
	/** Constructor
   The object used to load, store and access to the interaction rates of the process is the interactionRates pointer.
	 @param photonField		target photon field
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param thinning		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			step size limit as fraction of mean free path
   @param surface   suface to enclose the grid nodes to be loaded
	 */
	EMTripletPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, double thinning = 0, double limit = 0.1, ref_ptr<Surface> surface = nullptr);

	// set the target photon field
	void setPhotonField(ref_ptr<PhotonField> photonField);

	// decide if secondary electrons are added to the simulation	
	void setHaveElectrons(bool haveElectrons);

	/** limit the step to a fraction of the mean free path
	 @param limit	fraction of the mean free path, should be between 0 and 1
	*/
	void setLimit(double limit);

	/** Apply thinning with a given thinning factor
	 * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
	 */
	void setThinning(double thinning);
    
  /** Apply a surface that confine the position dependent photon field
   * @param surface closed surface to confine the grid to be  uploaded
   */
  void setSurface(ref_ptr<Surface> surface);
  ref_ptr<Surface> getSurface() const;
    
	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;
	
  /** Loads the interaction rate
   * (THIS FUNCTION WILL BE DEPRICATED SOON, use initRate(filename intRates) instead)
   * @param filename The name of the file containing the interaction rates
   */
  void initRate(std::string filename);
  /** Loads the cumultative interaction rate
   * (THIS FUNCTION WILL BE DEPRICATED SOON, use initRate(filename intRates) instead)
   * @param filename The name of the file containing the interaction rates
   */
  void initCumulativeRate(std::string filename);

    /** Loads the interaction rate in InteractionRates class
   * This function loads the interaction rate, in the proper object 
   * of the InteractionRates class, for the homogenouos background photon fields.
   * @param filename The name of the file containing the interaction rates
   * @param intRatesHom TODO
   */
  void initRate(std::string filename, ref_ptr<InteractionRatesHomogeneous> intRatesHom);
  /** Loads the cumultative interaction rate in InteractionRates class
   * This function is used to load the rates, in the dedicated object of the InteractionRates
   * class, for spatial dependent photon fields in the interaction module constructor.
   * @param filename The name of the file containing the interaction rates
   * @param intRatesPosDep TODO
   */
  void initCumulativeRate(std::string filename, ref_ptr<InteractionRatesHomogeneous> intRatesHom);
  
  /** Loads the interaction rate in InteractionRates class
   * This function is used to load the rates, in the dedicated object of the InteractionRates
   * class, for spatial dependent photon fields in the interaction module constructor.
   * @param filename The name of the file containing the interaction rates
   * @param intRatesPosDep TODO
   */
  void initRatePositionDependentPhotonField(std::string filename, ref_ptr<InteractionRatesPositionDependent> intRatesPosDep);
  /** Loads the cumultative interaction rate in InteractionRates class
   * This function is used to load the rates, in the dedicated object of the InteractionRates
   * class, for spatial dependent photon fields in the interaction module constructor.
   * @param filename The name of the file containing the interaction rates
   * @param intRatesPosDep TODO
   */
  void initCumulativeRatePositionDependentPhotonField(std::string filepath, ref_ptr<InteractionRatesPositionDependent> intRatesPosDep);
    
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate) const;

};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_EMTRIPLETPAIRPRODUCTION_H
