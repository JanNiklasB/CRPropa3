#include "crpropa/module/EMPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>


namespace crpropa {

static const double mec2 = mass_electron * c_squared;

EMPairProduction::EMPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons, double thinning, double limit) {
	setPhotonField(photonField);
	setThinning(thinning);
	setLimit(limit);
	setHaveElectrons(haveElectrons);
}

void EMPairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
	this->photonField = photonField;
	std::string fname = photonField->getFieldName();
	setDescription("EMPairProduction: " + fname);
	initRate(getDataPath("EMPairProduction/rate_" + fname + ".txt"));
	initCumulativeRate(getDataPath("EMPairProduction/cdf_" + fname + ".txt"));
}

void EMPairProduction::setHaveElectrons(bool haveElectrons) {
	this->haveElectrons = haveElectrons;
}

void EMPairProduction::setLimit(double limit) {
	this->limit = limit;
}

void EMPairProduction::setThinning(double thinning) {
	this->thinning = thinning;
}

void EMPairProduction::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMPairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
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

void EMPairProduction::initCumulativeRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMPairProduction: could not open file " + filename);

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

// Hold an data array to interpolate the energy distribution on
class PPSecondariesEnergyDistribution {
	private:
		double* tab_s=NULL;
		int tab_sSize=0;
		double** data=NULL;
		int dataSize=0;
		int* dataInnerSize=NULL;
		size_t N;

	public:
		// differential cross section for pair production for x = Epositron/Egamma, compare Lee 96 arXiv:9604098
		CUDA_CALLABLE_MEMBER double dSigmadE_PPx(double x, double beta) {
			double A = (x / (1. - x) + (1. - x) / x );
			double B =  (1. / x + 1. / (1. - x) );
			double y = (1 - beta * beta);
			return A + y * B - y * y / 4 * B * B;
		}

		CUDA_CALLABLE_MEMBER PPSecondariesEnergyDistribution() {
			N = 1000;
			size_t Ns = 1000;
			double s_min = 4 * mec2 * mec2;
			double s_max = 1e23 * eV * eV;
			double dls = log(s_max / s_min) / Ns;

			data = new double*[Ns];
			dataSize = Ns;
			dataInnerSize = new int[Ns];
			for (int i=0; i<Ns; i++){
				data[i] = new double[N];
				dataInnerSize[i] = N;
			}
			tab_s = new double[Ns+1];

			for (size_t i = 0; i < Ns + 1; ++i)
				tab_s[i] = s_min * exp(i*dls); // tabulate s bin borders

			for (size_t i = 0; i < Ns; i++) {
				double s = s_min * exp(i*dls + 0.5*dls);
				double beta = sqrt(1 - s_min/s);
				double x0 = (1 - beta) / 2;
				double dx = log((1 + beta) / (1 - beta)) / N;

				// cumulative midpoint integration
				data[i][0] = dSigmadE_PPx(x0, beta) * expm1(dx);
				for (size_t j = 1; j < N; j++) {
					double x = x0 * exp(j*dx + 0.5*dx);
					double binWidth = exp((j+1)*dx)-exp(j*dx);
					data[i][j] = dSigmadE_PPx(x, beta) * binWidth + data[i][j-1];
				}
			}
		}

		CUDA_CALLABLE_MEMBER ~PPSecondariesEnergyDistribution() {
			delete[] tab_s, data, dataInnerSize;
		}

		// sample positron energy from cdf(E, s_kin)
		CUDA_CALLABLE_MEMBER double sample(double E0, double s) {
			// get distribution for given s
			size_t idx = lower_bound<double>(s, tab_s, tab_sSize);
			if (idx > dataSize)
				return NAN;
				
			double* s0 = data[idx];

			// draw random bin
			Random &random = Random::instance();
			size_t j = random.randBin(s0, dataInnerSize[idx]) + 1;

			double s_min = 4. * mec2 * mec2;
			double beta = sqrtl(1. - s_min / s);
			double x0 = (1. - beta) / 2.;
			double dx = log((1 + beta) / (1 - beta)) / N;
			double binWidth = x0 * (exp(j*dx) - exp((j-1)*dx));
			if (random.rand() < 0.5)
				return E0 * (x0 * exp((j-1) * dx) + binWidth);
			else
				return E0 * (1 - (x0 * exp((j-1) * dx) + binWidth));
		}
};

void EMPairProduction::performInteraction(Candidate *candidate) const {
	
	// photon is lost after interacting
	candidate->setActive(false);

	// check if secondary electron pair needs to be produced
	if (not haveElectrons)
		return;
		
	// scale particle energy instead of background photon energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// check if in tabulated energy range
	if (E < tabEPtr[0] or (E > tabEPtr[tabESize-1]))
		return;

	// sample the value of s
	Random &random = Random::instance();
	size_t i = closestIndex(E, tabEPtr, tabESize);  // find closest tabulation point
	size_t j = random.randBin(tabCDFPtr[i], tabCDFInnerSizes[i]);
	double lo = std::max(4 * mec2 * mec2, tabsPtr[j-1]);  // first s-tabulation point below min(s_kin) = (2 me c^2)^2; ensure physical value
	double hi = tabsPtr[j];
	double s = lo + random.rand() * (hi - lo);

	// sample electron / positron energy
	static PPSecondariesEnergyDistribution interpolation;
	double Ee = interpolation.sample(E, s);
	double Ep = E - Ee;
	double f = Ep / E;

	// for some backgrounds Ee=nan due to precision limitations.
	if (not std::isfinite(Ee) || not std::isfinite(Ep))
		return;

	// sample random position along current step
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	// apply sampling
	if (random.rand() < pow(f, thinning)) {
		double w = 1. / pow(f, thinning);
		candidate->addSecondary(11, Ep / (1 + z), pos, w, interactionTag);
	}
	if (random.rand() < pow(1 - f, thinning)){
		double w = 1. / pow(1 - f, thinning);
		candidate->addSecondary(-11, Ee / (1 + z), pos, w, interactionTag);	
	}
}

void EMPairProduction::process(Candidate *candidate) const {
	// check if photon
	if (candidate->current.getId() != 22)
		return;

	// scale particle energy instead of background photon energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// check if in tabulated energy range
	if ((E < tabEnergyPtr[0]) or (E > tabEnergyPtr[tabEnergySize-1]))
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
		} else {
			performInteraction(candidate);
			return;
		}
		step -= randDistance; 
	} while (step > 0.);

}

void EMPairProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string EMPairProduction::getInteractionTag() const {
	return interactionTag;
}

} // namespace crpropa
