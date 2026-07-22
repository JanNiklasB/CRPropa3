#include "crpropa/electricField/ElectricField.h"

namespace crpropa{

void ElectricFieldList::addField(ref_ptr<ElectricField> field) {
	fields.push_back(field);
}

Vector3d ElectricFieldList::getField(const Vector3d &position) const {
	Vector3d b;
	for (int i = 0; i < fields.size(); i++)
		b += fields[i]->getField(position);
	return b;
}

}