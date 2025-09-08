#ifndef CRPROPA_OUTPUTSHELL_H
#define CRPROPA_OUTPUTSHELL_H

#include "crpropa/__CudaDefines.h"
#include "crpropa/Module.h"
#include "crpropa/AssocVector.h"
#include "crpropa/Variant.h"

namespace crpropa {
/**
 * \addtogroup Output
 * @{
 */

/**
 @class ShellOutput
 @brief Show the trajectory in the shell.
 */
class ShellOutput: public Module {
public:
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class ShellOutput1D
 @brief Show the trajectory in the shell.
 */
class ShellOutput1D: public Module {
public:
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class ShellPropertyOutput
 @brief Show the candidate properties in the shell.
 */
class ShellPropertyOutput: public Module {
public:
	typedef Loki::AssocVector<std::string, Variant> PropertyMap;
	CUDA_CALLABLE_MEMBER void process(Candidate *candidate) const;
	std::string getDescription() const;
};
/** @}*/

} // namespace cprpropa

#endif // CRPROPA_OUTPUTSHELL_H
