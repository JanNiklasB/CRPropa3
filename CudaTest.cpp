#include "CRPropa.h"
#include <vector>
#include <iostream>

using namespace crpropa;
using namespace std;

void SimpleCudaSim(){
	// int numCandidates = 100;

	// ModuleList SIM;
	// ref_ptr<MagneticFieldList> Field = new MagneticFieldList();
	
	// ref_ptr<TurbulenceSpectrum> Spectrum = new TurbulenceSpectrum(1.*muG, 0.1*pc, 100*pc);
	// ref_ptr<PlaneWaveTurbulence> Turbulence = new PlaneWaveTurbulence(*Spectrum, 100, 1);
	// Field->addField(Turbulence);

	// SIM.add(new PropagationBP(Field, 0.001, 0.001*pc, 20.*pc));
	// SIM.add(new MaximumTrajectoryLength(1*kpc));
	
	// ref_ptr<Observer> Obs = new Observer();
	// ref_ptr<TextOutput> output = new TextOutput("../test.txt", Output::Event1D);
	// Obs->add(new ObserverTimeEvolution(0.001*pc, 1*kpc, 10000, false));
	// Obs->setDeactivateOnDetection(false);
	// Obs->onDetection(output);
	// SIM.add(Obs);

	// ref_ptr<Source> source = new Source();
	// source->add( new SourceUniform1D(1 * Mpc, 1000 * Mpc) );
	// source->add( new SourceRedshift1D() );
	// SIM.setShowProgress(true);
	
	// SIM.run(source.get(), 1000);

	printf("test\n");
	ModuleList sim;
	// sim.add(new SimplePropagation(1*kpc, 10*Mpc));
	// sim.add(new Redshift());
	// sim.add(new PhotoPionProduction(new CMB()));
	// sim.add(new PhotoPionProduction(new IRB_Dominguez11()));
	// sim.add(new PhotoDisintegration(new CMB()));
	// sim.add(new PhotoDisintegration(new IRB_Dominguez11()));
	// sim.add(new NuclearDecay());
	// sim.add(new ElectronPairProduction(new CMB()));
	// sim.add(new ElectronPairProduction(new IRB_Dominguez11()));
	// sim.add(new MinimumEnergy(1*EeV));

	// printf("test\n");
	// ref_ptr<Observer> obs = new Observer();
	// obs->add(new Observer1D());
	// obs->onDetection(new TextOutput("../events.txt", Output::Event1D));
	// sim.add(obs);

	// printf("test\n");
	// ref_ptr<Source> source = new Source();
	// source->add(new SourceUniform1D(1*Mpc, 1000*Mpc));
	// source->add(new SourceRedshift1D());

	// printf("test\n");
	// ref_ptr<SourceComposition> composition =
	// 	new SourceComposition(1*EeV, 100*EeV, -1);
	// composition->add(1,  1,  1);
	// composition->add(4,  2,  1);
	// composition->add(14, 7,  1);
	// composition->add(56, 26, 1);
	// source->add(composition);

	// sim.setShowProgress(true);
	// sim.run(source.get(), 2000, true);
}

// #include "crpropa/Logging.h"

int main(){
	// int numCandidates = 100;
	// cudaDeviceProp prop;
	// cudaGetDeviceProperties(&prop, 0);

	// int threadsPerBlock = prop.maxThreadsDim[0];
	// int blocksPerGrid   = (numCandidates + threadsPerBlock - 1) / threadsPerBlock;
	// testKernal<<<blocksPerGrid, threadsPerBlock>>>();
	printf("test\n");
	SimpleCudaSim();
	// logInfo("Test");

	return 0;
}