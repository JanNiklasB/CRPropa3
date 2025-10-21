#include "crpropa/Module.h"

#include <typeinfo>

namespace crpropa {

Module::Module() {
	const std::type_info &info = typeid(*this);
	setDescription(info.name());
}

Module::~Module() {}

std::string Module::getDescription() const {
	return std::string(description);
}

void Module::setDescription(const std::string &d) {
	description = d.c_str();
}

AbstractCondition::AbstractCondition() :
		makeRejectedInactive(true), makeAcceptedInactive(false) {
	setRejectFlag("Rejected", typeid(*this).name());
}

AbstractCondition::~AbstractCondition(){
	delete[] rejectFlagKey;
	delete[] rejectFlagValue;
	delete[] acceptFlagKey;
	delete[] acceptFlagValue;
	Module::~Module();
}

void AbstractCondition::reject(Candidate *candidate) const {
	if (!candidate)
		return;

	if (rejectAction)
		rejectAction->process(candidate);

	if (rejectFlagKeySize)
		candidate->setProperty(rejectFlagKey, rejectFlagValue);

	if (makeRejectedInactive)
		candidate->setActive(false);
}

void AbstractCondition::accept(Candidate *candidate) const {
	if (!candidate)
		return;

	if (acceptAction)
		acceptAction->process(candidate);

	if (acceptFlagKeySize)
		candidate->setProperty(acceptFlagKey, acceptFlagValue);

	if (makeAcceptedInactive)
		candidate->setActive(false);
}

void AbstractCondition::setMakeRejectedInactive(bool deactivate) {
	makeRejectedInactive = deactivate;
}

void AbstractCondition::setMakeAcceptedInactive(bool deactivate) {
	makeAcceptedInactive = deactivate;
}

void AbstractCondition::onReject(Module *action) {
	rejectAction = action;
}

void AbstractCondition::onAccept(Module *action) {
	acceptAction = action;
}

void AbstractCondition::setRejectFlag(std::string key, std::string value) {
	delete[] rejectFlagKey;
	delete[] rejectFlagValue;

	rejectFlagKeySize = key.size();
	rejectFlagValueSize = value.size();

	rejectFlagKey = new char[rejectFlagKeySize];
	rejectFlagValue = new char[rejectFlagValueSize];
	key.copy(rejectFlagKey, rejectFlagKeySize);
	value.copy(rejectFlagValue, rejectFlagValueSize);
}

void AbstractCondition::setAcceptFlag(std::string key, std::string value) {
	delete[] acceptFlagKey;
	delete[] acceptFlagValue;

	acceptFlagKeySize = key.size();
	acceptFlagValueSize = value.size();

	acceptFlagKey = new char[acceptFlagKeySize];
	acceptFlagValue = new char[acceptFlagValueSize];
	key.copy(acceptFlagKey, acceptFlagKeySize);
	value.copy(acceptFlagValue, acceptFlagValueSize);
}

std::string AbstractCondition::getRejectFlag() {
	std::string out = std::string(rejectFlagKey, rejectFlagKeySize)
		+ "&" + std::string(rejectFlagValue, rejectFlagKeySize);
	return out;
}

std::string AbstractCondition::getAcceptFlag() {
	std::string out = std::string(acceptFlagKey, acceptFlagKeySize)
		+ "&" + std::string(acceptFlagValue, acceptFlagValueSize);
	return out;
}


} // namespace crpropa
