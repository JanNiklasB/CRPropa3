/// Example plugin for CRPropa.
///
/// Please consider sharing the awesome plugin with you fellow researchers by
/// creating a seperate repository for your project. We maintain a list of
/// plugins to CRPropa on our webpage and are happy to add a link to your
/// project, just send us: (name of the plugin, short description, url)

#include <crpropa/Module.h>
#include <crpropa/Source.h>


/// A custom C++ module
class MyModule : public crpropa::Module
{
public:
	/// The parent's constructor need to be called on initialization!
	MyModule();
	// the process function needs the CUDA_CALLABLE_MEMBER tag, this enables usage on cuda cards,
	// any function used in process also needs the CUDA_CALLABLE_MEMBER (__device__ __host__) tag
	CUDA_CALLABLE_MEMBER void process(crpropa::Candidate *candidate) const;
};


/// A custom source feature
class AddMyProperty: public crpropa::SourceFeature
{
public:
	/// The parent's constructor need to be called on initialization!
	AddMyProperty();
	void prepareCandidate(crpropa::Candidate &candidate) const;
};
