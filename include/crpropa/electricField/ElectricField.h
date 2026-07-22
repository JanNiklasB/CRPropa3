#ifndef CRPROPA_ELECTRICFIELD_H
#define CRPROPA_ELECTRICFIELD_H

#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"

namespace crpropa{

/**
 @class ElectricField
 @brief Abstract base class for electric fields.
 */
class ElectricField: public Referenced {
	public:
		virtual ~ElectricField() {
		}
		/// A default virtual function to retrieve the field at a given position
		virtual Vector3d getField(const Vector3d &position) const {
			return Vector3d(0,0,0);
		};
		/// A default virtual function to retrieve the field at a given postion and redshift
		virtual Vector3d getField(const Vector3d &position, double z) const {
			return getField(position);
		};
};

/**
 @class ElectricFieldList
 @brief Electric field decorator implementing a superposition of fields.
 */
class ElectricFieldList: public ElectricField {
	private:
		std::vector<ref_ptr<ElectricField> > fields;
	public:
		/** A method to add another electric field
		 * @param field Another electric field
		 */
		void addField(ref_ptr<ElectricField> field);
		/** A method to retrieve the electric field at a given position
		 * @param position Vector3d with the cartesian position data
		 */
		Vector3d getField(const Vector3d &position) const;
};

/**
 @class UniformElectricField
 @brief Electric field with one E-field vector.
 */
class UniformElectricField: public ElectricField {
	private:
		Vector3d value;
	public:
		/**
		 * Constructor
		 * @param value electric field strength
		*/
		UniformElectricField(const Vector3d &value) :
				value(value) {
		}
		/** A method to retrieve the electric field at a given position
		 * @param position Vector3d with the cartesian position data
		 */
		Vector3d getField(const Vector3d &position) const {
			return value;
		}
};

}


#endif