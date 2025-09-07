#include "crpropa/module/EMInverseComptonScattering.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMInverseComptonScattering::EMInverseComptonScattering(ref_ptr<PhotonField> photonField, bool havePhotons, double thinning, double limit) {
	setPhotonField(photonField);
	setHavePhotons(havePhotons);
	setLimit(limit);
	setThinning(thinning);
}

void EMInverseComptonScattering::setPhotonField(ref_ptr<PhotonField> photonField) {
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMInverseComptonScattering: " + fname);
	initRate(getDataPath("EMInverseComptonScattering/rate_" + fname + ".txt"));
	initCumulativeRate(getDataPath("EMInverseComptonScattering/cdf_" + fname + ".txt"));
}

void EMInverseComptonScattering::setHavePhotons(bool havePhotons) {
	this->havePhotons = havePhotons;
}

void EMInverseComptonScattering::setLimit(double limit) {
	this->limit = limit;
}

void EMInverseComptonScattering::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMInverseComptonScattering::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMInverseComptonScattering: could not open file " + filename);

	// clear previously loaded tables
	tabEnergy.clear();
	tabRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEnergy.push_back(pow(10, a) * eV);
				tabRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
	
	tabEnergyPtr = tabEnergy.data();
	tabEnergySize = tabEnergy.size();
	tabRatePtr = tabRate.data();
	tabRateSize = tabRate.size();
}

void EMInverseComptonScattering::initCumulativeRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMInverseComptonScattering: could not open file " + filename);

	// clear previously loaded tables
	tabE.clear();
	tabs.clear();
	tabCDF.clear();
	if (tabCDFPtr) delete[] tabCDFPtr, tabCDFInnerSizes;
	
	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');

	// read s values in first line
	double a;
	infile >> a; // skip first value
	while (infile.good() and (infile.peek() != '\n')) {
		infile >> a;
		tabs.push_back(pow(10, a) * eV * eV);
	}

	// read all following lines: E, cdf values
	while (infile.good()) {
		infile >> a;
		if (!infile)
			break;  // end of file
		tabE.push_back(pow(10, a) * eV);
		std::vector<double> cdf;
		for (int i = 0; i < tabs.size(); i++) {
			infile >> a;
			cdf.push_back(a / Mpc);
		}
		tabCDF.push_back(cdf);
	}
	infile.close();

	tabEPtr = tabE.data();
	tabESize = tabE.size();
	tabsPtr = tabs.data();
	tabsSize = tabs.size();

	tabCDFPtr = new double*[tabCDF.size()];
	tabCDFInnerSizes = new int[tabCDF.size()];
	for (int i=0; i<tabCDF.size(); i++){
		tabCDFPtr[i] = tabCDF[i].data();
		tabCDFInnerSizes[i] = tabCDF[i].size();
	}
	tabCDFSize = tabCDF.size();

}

// Class to calculate the energy distribution of the ICS photon and to sample from it
class ICSSecondariesEnergyDistribution {
	private:
		double** data=NULL;
		int* dataInnerSize=NULL;
		int dataSize=0;
		double* s_values=NULL;
		int s_valuesSize=0;
		size_t Ns;
		size_t Nrer;
		double s_min;
		double s_max;
		double dls;

	public:
		// differential cross-section, see Lee '96 (arXiv:9604098), eq. 23 for x = Ee'/Ee
		CUDA_CALLABLE_MEMBER double dSigmadE(double x, double beta) {
			double q = ((1 - beta) / beta) * (1 - 1./x);
			return ((1 + beta) / beta) * (x + 1./x + 2 * q + q * q);
		}

		// create the cumulative energy distribution of the up-scattered photon
		CUDA_CALLABLE_MEMBER ICSSecondariesEnergyDistribution() {
			Ns = 1000;
			Nrer = 1000;
			s_min = mec2 * mec2;
			s_max = 2e23 * eV * eV;
			dls = (log(s_max) - log(s_min)) / Ns;

			data = new double*[Ns];
			dataInnerSize = new int[Ns];
			dataSize = Ns;
			for (int i=0; i<Ns; i++){
				data[i] = new double[Nrer];
				dataInnerSize[i] = Nrer;
			}

			// tabulate s bin borders
			s_values = new double[Ns+1];
			for (size_t i = 0; i < Ns + 1; ++i)
				s_values[i] = s_min * exp(i*dls);


			// for each s tabulate cumulative differential cross section
			for (size_t i = 0; i < Ns; i++) {
				double s = s_min * exp((i+0.5) * dls);
				double beta = (s - s_min) / (s + s_min);
				double x0 = (1 - beta) / (1 + beta);
				double dlx = -log(x0) / Nrer;

				// cumulative midpoint integration
				data[i][0] = dSigmadE(x0, beta) * expm1(dlx);
				for (size_t j = 1; j < Nrer; j++) {
					double x = x0 * exp((j+0.5) * dlx);
					double dx = exp((j+1) * dlx) - exp(j * dlx);
					data[i][j] = dSigmadE(x, beta) * dx;
					data[i][j] += data[i][j-1];
				}
			}
		}

		CUDA_CALLABLE_MEMBER ~ICSSecondariesEnergyDistribution(){
			delete[] data;
			delete[] dataInnerSize;
			delete[] s_values;
		}

		// draw random energy for the up-scattered photon Ep(Ee, s)
		CUDA_CALLABLE_MEMBER double sample(double Ee, double s) {
			size_t idx = lower_bound<double>(s, s_values, s_valuesSize);
			double *s0 = data[idx];
			Random &random = Random::instance();
			size_t j = random.randBin(s0, dataInnerSize[idx]) + 1; // draw random bin (upper bin boundary returned)
			double beta = (s - s_min) / (s + s_min);
			double x0 = (1 - beta) / (1 + beta);
			double dlx = -log(x0) / Nrer;
			double binWidth = x0 * (exp(j * dlx) - exp((j-1) * dlx));
			double Ep = (x0 * exp((j-1) * dlx) + binWidth) * Ee;
			return crstd::min(Ee, Ep); // prevent Ep > Ee from numerical inaccuracies
		}
};

void EMInverseComptonScattering::performInteraction(Candidate *candidate) const {
	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	if (E < tabEPtr[0] or E > tabEPtr[tabESize-1])
		return;

	// sample the value of s
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabEPtr, tabESize);
	size_t j = random.randBin(tabCDFPtr[i], tabCDFInnerSizes[i]);
	double s_kin = pow(10, log10(tabsPtr[j]) + (random.rand() - 0.5) * 0.1);
	double s = s_kin + mec2 * mec2;

	// sample electron energy after scattering
	static ICSSecondariesEnergyDistribution distribution;
	double Enew = distribution.sample(E, s);

	// add up-scattered photon
	if (havePhotons) {
		double Esecondary = E - Enew;
		double f = Enew / E;
		if (random.rand() < pow(1 - f, thinning)) {
			double w = 1. / pow(1 - f, thinning);
			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			candidate->addSecondary(22, Esecondary / (1 + z), pos, w, interactionTag);
		}
	}

	// update the primary particle energy; do this after adding the secondary to correctly set the secondary's parent
	candidate->current.setEnergy(Enew / (1 + z));
}

void EMInverseComptonScattering::process(Candidate *candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;

	// scale the particle energy instead of background photons
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	if (E < tabEnergyPtr[0] or (E > tabEnergyPtr[tabEnergySize-1]))
		return;

	// interaction rate
	double rate = interpolate(E, tabEnergyPtr, tabRatePtr, tabRateSize);
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

	// run this loop at least once to limit the step size
	double step = candidate->getCurrentStep();
	Random &random = Random::instance();
	do {
		double randDistance = -log(random.rand()) / rate;

		// check for interaction; if it doesn't ocurr, limit next step
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}
		performInteraction(candidate);

		// repeat with remaining step
		step -= randDistance;
	} while (step > 0);
}

void EMInverseComptonScattering::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMInverseComptonScattering::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
