#include "myPlugin.h"

#include <string>

const std::string myPropertyName = "counter";

// The parent's constructor need to be called on initialization!
MyModule::MyModule() : crpropa::Module()
{
	setDescription("MyPlugin::MyModule");
}

void MyModule::process(crpropa::Candidate *candidate) const
{
	// To enable parallelization, the modules have to be stateless - the
	// process method should thus not modify internal variables!
	// Furthermore, to enable parallelization per cuda the process method is only allowed to
	// use cuda friendly functions and methods, so no vector or string, etc.
	// Cuda provides many std function by default (see documentation if not sure)
	// CRPropa provides crstd which either resolves to cuda::std or std
	std::cout << "MyPlugin::MyModule::process() called\n";
	if(candidate->hasProperty(myPropertyName)) {
		uint32_t v = candidate->getProperty(myPropertyName);
		std::cout << " My property: " << myPropertyName << " = " <<  v << std::endl;
		candidate->setProperty(myPropertyName, crpropa::Variant::fromUInt32(v + 1));
	}
}

// ------------------------------------------------------------------
// The parent's constructor need to be called on initialization!
AddMyProperty::AddMyProperty() : crpropa::SourceFeature()
{
	description = "AddMyProperty: Adds my property to the candidate";
}
void AddMyProperty::prepareCandidate(crpropa::Candidate &candidate) const
{
	candidate.setProperty(myPropertyName, crpropa::Variant::fromUInt32(0));
}
