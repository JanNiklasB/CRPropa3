#include "crpropa/Candidate.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/ElasticScattering.h"
#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/module/Redshift.h"
#include "crpropa/module/EMPairProduction.h"
#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/module/SynchrotronRadiation.h"

#include <fstream>
#include <iostream>

namespace crpropa {

void testSynchrotronPhotonEnergy(){
	double brms = 1 * muG; 
	std::cout << "test before SynchrotronRadiationConstructor\n";
	SynchrotronRadiation sync(brms, true);
	sync.setSecondaryThreshold(0.); // allow all secondaries for testing

	std::cout << "before Candidate constructor\n";
	double E = 1 * TeV;
	Candidate c(11, E);
	c.setCurrentStep(10 * pc); 
	c.setNextStep(10 * pc);
	
	double lf = c.current.getLorentzFactor();
	double Rg = E / eplus / c_light / (brms * sqrt(2. / 3) ); // factor 2/3 for avg magnetic field direction. 
	double Ecrit = 3. / 4 * h_planck / M_PI * c_light * pow(lf, 3) / Rg;

	std::cout << "test before process\n";
	sync.process(c);

	// check avg energy of the secondary photons 
	double Esec = 0; 
	std::cout << "test before getEnergy loop\n";
	for (size_t i = 0; i < c.secondaries.size(); i++) {
		std::cout << i << std::endl;
		Esec += c.secondaries[i] -> current.getEnergy();
	}
	Esec /= c.secondaries.size();
}


} // namespace crpropa

int main(int argc, char **argv) {
	crpropa::testSynchrotronPhotonEnergy();
	// ::testing::InitGoogleTest(&argc, argv);
	// return RUN_ALL_TESTS();
	return 0;
}