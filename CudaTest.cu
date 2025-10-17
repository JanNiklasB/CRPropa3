#include "CRPropa.h"
#include "cuda_runtime.h"
#include <vector>

using namespace crpropa;
using namespace std;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

// just start a single simple simulation in a kernal and initialize everything here, see if it works...
__global__ void testKernal(){
}

void SimpleCudaSim(){
	int numCandidates = 1;

	ModuleList SIM;
	ref_ptr<MagneticFieldList> Field = new MagneticFieldList();
	
	ref_ptr<TurbulenceSpectrum> Spectrum = new TurbulenceSpectrum(1.*muG, 0.1*pc, 100*pc);
	ref_ptr<PlaneWaveTurbulence> Turbulence = new PlaneWaveTurbulence(*Spectrum, 100, 1);
	Field->addField(Turbulence);

	SIM.add(new SimplePropagation(0.1*pc, 0.1*pc));
	// SIM.add(new MaximumTrajectoryLength(10*pc));

	// ref_ptr<ParticleCollector> Collector = new ParticleCollector(100, false);
	// SIM.add(Collector);
	
	ref_ptr<Observer> Obs = new Observer();
	// ref_ptr<TextOutput> output = new TextOutput("../test.txt", Output::Event1D);
	ref_ptr<ShellOutput> output = new ShellOutput();
	Obs->add(new ObserverTimeEvolution(0.001*pc, 1*kpc, 10000, false));
	Obs->setDeactivateOnDetection(false);
	Obs->onDetection(output);
	// SIM.add(Obs);

	ref_ptr<Source> source = new Source();
	source->add( new SourceUniform1D(1 * Mpc, 1000 * Mpc) );
	
	SIM.cudarun(source.get(), numCandidates);

	// Collector->dump("../testCandidates.txt");
	printf("test\n");
}


int main(){
	// int numCandidates = 100;
	// cudaDeviceProp prop;
	// cudaGetDeviceProperties(&prop, 0);

	// int threadsPerBlock = prop.maxThreadsDim[0];
	// int blocksPerGrid   = (numCandidates + threadsPerBlock - 1) / threadsPerBlock;
	// testKernal<<<blocksPerGrid, threadsPerBlock>>>();
	SimpleCudaSim();

	// cudaDeviceSynchronize();

	

	return 0;
}