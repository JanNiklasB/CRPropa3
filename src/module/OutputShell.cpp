#include "crpropa/module/OutputShell.h"
#include "crpropa/Units.h"

#include <iomanip>

namespace crpropa {

void ShellOutput::process(Candidate* c) const {
	#pragma omp critical(ShellOutput)
	{
		// showpoint is set by "#" width of 6 is set by 6, precision is set by .3
		printf("%#6.3f Mpc,  ", c->getTrajectoryLength() / Mpc);
		printf("%#6.3f, ", c->getRedshift());
		printf("%#6.3f, ", c->current.getId());
		printf("%#6.3f EeV, ", c->current.getEnergy() / EeV);
		printf("%#6.3f Mpc, ", c->current.getPosition() / Mpc);
		printf("%#6.3f\n", c->current.getDirection());
	}
}

std::string ShellOutput::getDescription() const {
	return "Shell output";
}

void ShellOutput1D::process(Candidate* c) const {
	#pragma omp critical(ShellOutput)
	{
		// showpoint is set by "#" width of 6 is set by 6, precision is set by .3
		printf("%#6.3f Mpc, ", c->current.getPosition() / Mpc);
		printf("%#6.3f, ", c->getRedshift());
		printf("%#6.3f, ", c->current.getId());
		printf("%#6.3f EeV\n", c->current.getEnergy() / EeV);
	}
}

std::string ShellOutput1D::getDescription() const {
	return "Shell output for 1D";
}

void ShellPropertyOutput::process(Candidate* c) const {
	#ifdef __CUDACC__
	#pragma omp critical(ShellOutput)
	{
		for (size_t i=0 ; i <= c->properties.size(); i++) {
			printf("  %f, %f\n", c->properties[i].first, c->properties[i].second);
		}
	}

	#else
	Candidate::PropertyMap::const_iterator i = c->properties.begin();
	#pragma omp critical(ShellOutput)
	{
		for ( ; i != c->properties.end(); i++) {
			printf("  %f, %f\n", i->first, i->second);
		}
	}
	#endif
}

std::string ShellPropertyOutput::getDescription() const {
	return "Shell property output";
}

} // namespace crpropa
