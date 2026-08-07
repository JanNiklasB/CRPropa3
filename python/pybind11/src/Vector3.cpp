#include <ios> // std::_Ios_Seekdir
#include <iterator> // __gnu_cxx::__normal_iterator
#include <locale> // std::locale
#include <memory> // std::allocator
#include <ostream> // std::basic_ostream
#include <sstream> // __str__
#include <streambuf> // std::basic_streambuf
#include <string> // std::basic_string
#include <string> // std::char_traits

#include <functional>
#include <pybind11/pybind11.h>
#include <string>
#include <crpropa/Vector3.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>, false)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*, false)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

void bind_unknown_unknown(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // crpropa::Vector3 file: line:36
		pybind11::class_<crpropa::Vector3<double>, std::shared_ptr<crpropa::Vector3<double>>> cl(M("crpropa"), "Vector3_double_t", "");
		cl.def( pybind11::init( [](){ return new crpropa::Vector3<double>(); } ) );
		cl.def( pybind11::init( [](crpropa::Vector3<double> const &o){ return new crpropa::Vector3<double>(o); } ) );
		cl.def( pybind11::init<const double *>(), pybind11::arg("v") );

		cl.def( pybind11::init<const float *>(), pybind11::arg("v") );

		cl.def( pybind11::init<const double &, const double &, const double &>(), pybind11::arg("X"), pybind11::arg("Y"), pybind11::arg("Z") );

		cl.def( pybind11::init<double>(), pybind11::arg("t") );

		cl.def("setX", (void (crpropa::Vector3<double>::*)(const double)) &crpropa::Vector3<double>::setX, "C++: crpropa::Vector3<double>::setX(const double) --> void", pybind11::arg("X"));
		cl.def("setY", (void (crpropa::Vector3<double>::*)(const double)) &crpropa::Vector3<double>::setY, "C++: crpropa::Vector3<double>::setY(const double) --> void", pybind11::arg("Y"));
		cl.def("setZ", (void (crpropa::Vector3<double>::*)(const double)) &crpropa::Vector3<double>::setZ, "C++: crpropa::Vector3<double>::setZ(const double) --> void", pybind11::arg("Z"));
		cl.def("setXYZ", (void (crpropa::Vector3<double>::*)(const double, const double, const double)) &crpropa::Vector3<double>::setXYZ, "C++: crpropa::Vector3<double>::setXYZ(const double, const double, const double) --> void", pybind11::arg("X"), pybind11::arg("Y"), pybind11::arg("Z"));
		cl.def("setR", (void (crpropa::Vector3<double>::*)(const double)) &crpropa::Vector3<double>::setR, "C++: crpropa::Vector3<double>::setR(const double) --> void", pybind11::arg("r"));
		cl.def("setRThetaPhi", (void (crpropa::Vector3<double>::*)(const double, const double, const double)) &crpropa::Vector3<double>::setRThetaPhi, "C++: crpropa::Vector3<double>::setRThetaPhi(const double, const double, const double) --> void", pybind11::arg("r"), pybind11::arg("theta"), pybind11::arg("phi"));
		cl.def("getX", (double (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getX, "C++: crpropa::Vector3<double>::getX() const --> double");
		cl.def("getY", (double (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getY, "C++: crpropa::Vector3<double>::getY() const --> double");
		cl.def("getZ", (double (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getZ, "C++: crpropa::Vector3<double>::getZ() const --> double");
		cl.def("getR", (double (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getR, "C++: crpropa::Vector3<double>::getR() const --> double");
		cl.def("getR2", (double (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getR2, "C++: crpropa::Vector3<double>::getR2() const --> double");
		cl.def("getPhi", (double (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getPhi, "C++: crpropa::Vector3<double>::getPhi() const --> double");
		cl.def("getTheta", (double (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getTheta, "C++: crpropa::Vector3<double>::getTheta() const --> double");
		cl.def("getUnitVector", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getUnitVector, "C++: crpropa::Vector3<double>::getUnitVector() const --> class crpropa::Vector3<double>");
		cl.def("getUnitVectorTheta", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getUnitVectorTheta, "C++: crpropa::Vector3<double>::getUnitVectorTheta() const --> class crpropa::Vector3<double>");
		cl.def("getUnitVectorPhi", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::getUnitVectorPhi, "C++: crpropa::Vector3<double>::getUnitVectorPhi() const --> class crpropa::Vector3<double>");
		cl.def("getAngleTo", (double (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::getAngleTo, "C++: crpropa::Vector3<double>::getAngleTo(const class crpropa::Vector3<double> &) const --> double", pybind11::arg("v"));
		cl.def("isParallelTo", (bool (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &, double) const) &crpropa::Vector3<double>::isParallelTo, "C++: crpropa::Vector3<double>::isParallelTo(const class crpropa::Vector3<double> &, double) const --> bool", pybind11::arg("v"), pybind11::arg("maxAngle"));
		cl.def("getDistanceTo", (double (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::getDistanceTo, "C++: crpropa::Vector3<double>::getDistanceTo(const class crpropa::Vector3<double> &) const --> double", pybind11::arg("point"));
		cl.def("getParallelTo", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::getParallelTo, "C++: crpropa::Vector3<double>::getParallelTo(const class crpropa::Vector3<double> &) const --> class crpropa::Vector3<double>", pybind11::arg("v"));
		cl.def("getPerpendicularTo", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::getPerpendicularTo, "C++: crpropa::Vector3<double>::getPerpendicularTo(const class crpropa::Vector3<double> &) const --> class crpropa::Vector3<double>", pybind11::arg("v"));
		cl.def("getRotated", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &, double) const) &crpropa::Vector3<double>::getRotated, "C++: crpropa::Vector3<double>::getRotated(const class crpropa::Vector3<double> &, double) const --> class crpropa::Vector3<double>", pybind11::arg("axis"), pybind11::arg("angle"));
		cl.def("clip", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(double, double) const) &crpropa::Vector3<double>::clip, "C++: crpropa::Vector3<double>::clip(double, double) const --> class crpropa::Vector3<double>", pybind11::arg("lower"), pybind11::arg("upper"));
		cl.def("abs", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::abs, "C++: crpropa::Vector3<double>::abs() const --> class crpropa::Vector3<double>");
		cl.def("floor", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::floor, "C++: crpropa::Vector3<double>::floor() const --> class crpropa::Vector3<double>");
		cl.def("ceil", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::ceil, "C++: crpropa::Vector3<double>::ceil() const --> class crpropa::Vector3<double>");
		cl.def("min", (double (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::min, "C++: crpropa::Vector3<double>::min() const --> double");
		cl.def("max", (double (crpropa::Vector3<double>::*)() const) &crpropa::Vector3<double>::max, "C++: crpropa::Vector3<double>::max() const --> double");
		cl.def("dot", (double (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::dot, "C++: crpropa::Vector3<double>::dot(const class crpropa::Vector3<double> &) const --> double", pybind11::arg("v"));
		cl.def("cross", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::cross, "C++: crpropa::Vector3<double>::cross(const class crpropa::Vector3<double> &) const --> class crpropa::Vector3<double>", pybind11::arg("v"));
		cl.def("__eq__", (bool (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::operator==, "C++: crpropa::Vector3<double>::operator==(const class crpropa::Vector3<double> &) const --> bool", pybind11::arg("v"));
		cl.def("__add__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::operator+, "C++: crpropa::Vector3<double>::operator+(const class crpropa::Vector3<double> &) const --> class crpropa::Vector3<double>", pybind11::arg("v"));
		cl.def("__add__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const double &) const) &crpropa::Vector3<double>::operator+, "C++: crpropa::Vector3<double>::operator+(const double &) const --> class crpropa::Vector3<double>", pybind11::arg("f"));
		cl.def("__sub__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::operator-, "C++: crpropa::Vector3<double>::operator-(const class crpropa::Vector3<double> &) const --> class crpropa::Vector3<double>", pybind11::arg("v"));
		cl.def("__sub__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const double &) const) &crpropa::Vector3<double>::operator-, "C++: crpropa::Vector3<double>::operator-(const double &) const --> class crpropa::Vector3<double>", pybind11::arg("f"));
		cl.def("__mul__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::operator*, "C++: crpropa::Vector3<double>::operator*(const class crpropa::Vector3<double> &) const --> class crpropa::Vector3<double>", pybind11::arg("v"));
		cl.def("__mul__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(double) const) &crpropa::Vector3<double>::operator*, "C++: crpropa::Vector3<double>::operator*(double) const --> class crpropa::Vector3<double>", pybind11::arg("v"));
		cl.def("__truediv__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::operator/, "C++: crpropa::Vector3<double>::operator/(const class crpropa::Vector3<double> &) const --> class crpropa::Vector3<double>", pybind11::arg("v"));
		cl.def("__truediv__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const double &) const) &crpropa::Vector3<double>::operator/, "C++: crpropa::Vector3<double>::operator/(const double &) const --> class crpropa::Vector3<double>", pybind11::arg("f"));
		cl.def("__mod__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &) const) &crpropa::Vector3<double>::operator%, "C++: crpropa::Vector3<double>::operator%(const class crpropa::Vector3<double> &) const --> class crpropa::Vector3<double>", pybind11::arg("v"));
		cl.def("__mod__", (class crpropa::Vector3<double> (crpropa::Vector3<double>::*)(const double &) const) &crpropa::Vector3<double>::operator%, "C++: crpropa::Vector3<double>::operator%(const double &) const --> class crpropa::Vector3<double>", pybind11::arg("f"));
		cl.def("__isub__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &)) &crpropa::Vector3<double>::operator-=, "C++: crpropa::Vector3<double>::operator-=(const class crpropa::Vector3<double> &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__isub__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const double &)) &crpropa::Vector3<double>::operator-=, "C++: crpropa::Vector3<double>::operator-=(const double &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("__iadd__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &)) &crpropa::Vector3<double>::operator+=, "C++: crpropa::Vector3<double>::operator+=(const class crpropa::Vector3<double> &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__iadd__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const double &)) &crpropa::Vector3<double>::operator+=, "C++: crpropa::Vector3<double>::operator+=(const double &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("__imul__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &)) &crpropa::Vector3<double>::operator*=, "C++: crpropa::Vector3<double>::operator*=(const class crpropa::Vector3<double> &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__imul__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const double &)) &crpropa::Vector3<double>::operator*=, "C++: crpropa::Vector3<double>::operator*=(const double &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("__itruediv__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &)) &crpropa::Vector3<double>::operator/=, "C++: crpropa::Vector3<double>::operator/=(const class crpropa::Vector3<double> &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__itruediv__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const double &)) &crpropa::Vector3<double>::operator/=, "C++: crpropa::Vector3<double>::operator/=(const double &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("__imod__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &)) &crpropa::Vector3<double>::operator%=, "C++: crpropa::Vector3<double>::operator%=(const class crpropa::Vector3<double> &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__imod__", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const double &)) &crpropa::Vector3<double>::operator%=, "C++: crpropa::Vector3<double>::operator%=(const double &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("assign", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const class crpropa::Vector3<double> &)) &crpropa::Vector3<double>::operator=, "C++: crpropa::Vector3<double>::operator=(const class crpropa::Vector3<double> &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("assign", (class crpropa::Vector3<double> & (crpropa::Vector3<double>::*)(const double &)) &crpropa::Vector3<double>::operator=, "C++: crpropa::Vector3<double>::operator=(const double &) --> class crpropa::Vector3<double> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));

		cl.def("__str__", [](crpropa::Vector3<double> const &o) -> std::string { std::ostringstream s; using namespace crpropa; s << o; return s.str(); } );

		{ // crpropa::Vector3<double>::(anonymous union at /home/jnb/Repos/CRPDEV/include/crpropa/Vector3.h:38:2) file: line:38

			{ // crpropa::Vector3<double>::(anonymous union)::(anonymous struct at /home/jnb/Repos/CRPDEV/include/crpropa/Vector3.h:39:3) file: line:39
				cl.def_readwrite("x", &crpropa::Vector3<double>::x);
				cl.def_readwrite("y", &crpropa::Vector3<double>::y);
				cl.def_readwrite("z", &crpropa::Vector3<double>::z);
			}

		}

	}
	{ // crpropa::Vector3 file: line:36
		pybind11::class_<crpropa::Vector3<float>, std::shared_ptr<crpropa::Vector3<float>>> cl(M("crpropa"), "Vector3_float_t", "");
		cl.def( pybind11::init( [](){ return new crpropa::Vector3<float>(); } ) );
		cl.def( pybind11::init( [](crpropa::Vector3<float> const &o){ return new crpropa::Vector3<float>(o); } ) );
		cl.def( pybind11::init<const double *>(), pybind11::arg("v") );

		cl.def( pybind11::init<const float *>(), pybind11::arg("v") );

		cl.def( pybind11::init<const float &, const float &, const float &>(), pybind11::arg("X"), pybind11::arg("Y"), pybind11::arg("Z") );

		cl.def( pybind11::init<float>(), pybind11::arg("t") );

		cl.def("setX", (void (crpropa::Vector3<float>::*)(const float)) &crpropa::Vector3<float>::setX, "C++: crpropa::Vector3<float>::setX(const float) --> void", pybind11::arg("X"));
		cl.def("setY", (void (crpropa::Vector3<float>::*)(const float)) &crpropa::Vector3<float>::setY, "C++: crpropa::Vector3<float>::setY(const float) --> void", pybind11::arg("Y"));
		cl.def("setZ", (void (crpropa::Vector3<float>::*)(const float)) &crpropa::Vector3<float>::setZ, "C++: crpropa::Vector3<float>::setZ(const float) --> void", pybind11::arg("Z"));
		cl.def("setXYZ", (void (crpropa::Vector3<float>::*)(const float, const float, const float)) &crpropa::Vector3<float>::setXYZ, "C++: crpropa::Vector3<float>::setXYZ(const float, const float, const float) --> void", pybind11::arg("X"), pybind11::arg("Y"), pybind11::arg("Z"));
		cl.def("setR", (void (crpropa::Vector3<float>::*)(const float)) &crpropa::Vector3<float>::setR, "C++: crpropa::Vector3<float>::setR(const float) --> void", pybind11::arg("r"));
		cl.def("setRThetaPhi", (void (crpropa::Vector3<float>::*)(const float, const float, const float)) &crpropa::Vector3<float>::setRThetaPhi, "C++: crpropa::Vector3<float>::setRThetaPhi(const float, const float, const float) --> void", pybind11::arg("r"), pybind11::arg("theta"), pybind11::arg("phi"));
		cl.def("getX", (float (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getX, "C++: crpropa::Vector3<float>::getX() const --> float");
		cl.def("getY", (float (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getY, "C++: crpropa::Vector3<float>::getY() const --> float");
		cl.def("getZ", (float (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getZ, "C++: crpropa::Vector3<float>::getZ() const --> float");
		cl.def("getR", (float (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getR, "C++: crpropa::Vector3<float>::getR() const --> float");
		cl.def("getR2", (float (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getR2, "C++: crpropa::Vector3<float>::getR2() const --> float");
		cl.def("getPhi", (float (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getPhi, "C++: crpropa::Vector3<float>::getPhi() const --> float");
		cl.def("getTheta", (float (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getTheta, "C++: crpropa::Vector3<float>::getTheta() const --> float");
		cl.def("getUnitVector", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getUnitVector, "C++: crpropa::Vector3<float>::getUnitVector() const --> class crpropa::Vector3<float>");
		cl.def("getUnitVectorTheta", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getUnitVectorTheta, "C++: crpropa::Vector3<float>::getUnitVectorTheta() const --> class crpropa::Vector3<float>");
		cl.def("getUnitVectorPhi", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::getUnitVectorPhi, "C++: crpropa::Vector3<float>::getUnitVectorPhi() const --> class crpropa::Vector3<float>");
		cl.def("getAngleTo", (float (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::getAngleTo, "C++: crpropa::Vector3<float>::getAngleTo(const class crpropa::Vector3<float> &) const --> float", pybind11::arg("v"));
		cl.def("isParallelTo", (bool (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &, float) const) &crpropa::Vector3<float>::isParallelTo, "C++: crpropa::Vector3<float>::isParallelTo(const class crpropa::Vector3<float> &, float) const --> bool", pybind11::arg("v"), pybind11::arg("maxAngle"));
		cl.def("getDistanceTo", (float (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::getDistanceTo, "C++: crpropa::Vector3<float>::getDistanceTo(const class crpropa::Vector3<float> &) const --> float", pybind11::arg("point"));
		cl.def("getParallelTo", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::getParallelTo, "C++: crpropa::Vector3<float>::getParallelTo(const class crpropa::Vector3<float> &) const --> class crpropa::Vector3<float>", pybind11::arg("v"));
		cl.def("getPerpendicularTo", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::getPerpendicularTo, "C++: crpropa::Vector3<float>::getPerpendicularTo(const class crpropa::Vector3<float> &) const --> class crpropa::Vector3<float>", pybind11::arg("v"));
		cl.def("getRotated", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &, float) const) &crpropa::Vector3<float>::getRotated, "C++: crpropa::Vector3<float>::getRotated(const class crpropa::Vector3<float> &, float) const --> class crpropa::Vector3<float>", pybind11::arg("axis"), pybind11::arg("angle"));
		cl.def("clip", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(float, float) const) &crpropa::Vector3<float>::clip, "C++: crpropa::Vector3<float>::clip(float, float) const --> class crpropa::Vector3<float>", pybind11::arg("lower"), pybind11::arg("upper"));
		cl.def("abs", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::abs, "C++: crpropa::Vector3<float>::abs() const --> class crpropa::Vector3<float>");
		cl.def("floor", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::floor, "C++: crpropa::Vector3<float>::floor() const --> class crpropa::Vector3<float>");
		cl.def("ceil", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::ceil, "C++: crpropa::Vector3<float>::ceil() const --> class crpropa::Vector3<float>");
		cl.def("min", (float (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::min, "C++: crpropa::Vector3<float>::min() const --> float");
		cl.def("max", (float (crpropa::Vector3<float>::*)() const) &crpropa::Vector3<float>::max, "C++: crpropa::Vector3<float>::max() const --> float");
		cl.def("dot", (float (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::dot, "C++: crpropa::Vector3<float>::dot(const class crpropa::Vector3<float> &) const --> float", pybind11::arg("v"));
		cl.def("cross", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::cross, "C++: crpropa::Vector3<float>::cross(const class crpropa::Vector3<float> &) const --> class crpropa::Vector3<float>", pybind11::arg("v"));
		cl.def("__eq__", (bool (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::operator==, "C++: crpropa::Vector3<float>::operator==(const class crpropa::Vector3<float> &) const --> bool", pybind11::arg("v"));
		cl.def("__add__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::operator+, "C++: crpropa::Vector3<float>::operator+(const class crpropa::Vector3<float> &) const --> class crpropa::Vector3<float>", pybind11::arg("v"));
		cl.def("__add__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const float &) const) &crpropa::Vector3<float>::operator+, "C++: crpropa::Vector3<float>::operator+(const float &) const --> class crpropa::Vector3<float>", pybind11::arg("f"));
		cl.def("__sub__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::operator-, "C++: crpropa::Vector3<float>::operator-(const class crpropa::Vector3<float> &) const --> class crpropa::Vector3<float>", pybind11::arg("v"));
		cl.def("__sub__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const float &) const) &crpropa::Vector3<float>::operator-, "C++: crpropa::Vector3<float>::operator-(const float &) const --> class crpropa::Vector3<float>", pybind11::arg("f"));
		cl.def("__mul__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::operator*, "C++: crpropa::Vector3<float>::operator*(const class crpropa::Vector3<float> &) const --> class crpropa::Vector3<float>", pybind11::arg("v"));
		cl.def("__mul__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(float) const) &crpropa::Vector3<float>::operator*, "C++: crpropa::Vector3<float>::operator*(float) const --> class crpropa::Vector3<float>", pybind11::arg("v"));
		cl.def("__truediv__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::operator/, "C++: crpropa::Vector3<float>::operator/(const class crpropa::Vector3<float> &) const --> class crpropa::Vector3<float>", pybind11::arg("v"));
		cl.def("__truediv__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const float &) const) &crpropa::Vector3<float>::operator/, "C++: crpropa::Vector3<float>::operator/(const float &) const --> class crpropa::Vector3<float>", pybind11::arg("f"));
		cl.def("__mod__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &) const) &crpropa::Vector3<float>::operator%, "C++: crpropa::Vector3<float>::operator%(const class crpropa::Vector3<float> &) const --> class crpropa::Vector3<float>", pybind11::arg("v"));
		cl.def("__mod__", (class crpropa::Vector3<float> (crpropa::Vector3<float>::*)(const float &) const) &crpropa::Vector3<float>::operator%, "C++: crpropa::Vector3<float>::operator%(const float &) const --> class crpropa::Vector3<float>", pybind11::arg("f"));
		cl.def("__isub__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &)) &crpropa::Vector3<float>::operator-=, "C++: crpropa::Vector3<float>::operator-=(const class crpropa::Vector3<float> &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__isub__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const float &)) &crpropa::Vector3<float>::operator-=, "C++: crpropa::Vector3<float>::operator-=(const float &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("__iadd__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &)) &crpropa::Vector3<float>::operator+=, "C++: crpropa::Vector3<float>::operator+=(const class crpropa::Vector3<float> &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__iadd__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const float &)) &crpropa::Vector3<float>::operator+=, "C++: crpropa::Vector3<float>::operator+=(const float &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("__imul__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &)) &crpropa::Vector3<float>::operator*=, "C++: crpropa::Vector3<float>::operator*=(const class crpropa::Vector3<float> &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__imul__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const float &)) &crpropa::Vector3<float>::operator*=, "C++: crpropa::Vector3<float>::operator*=(const float &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("__itruediv__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &)) &crpropa::Vector3<float>::operator/=, "C++: crpropa::Vector3<float>::operator/=(const class crpropa::Vector3<float> &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__itruediv__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const float &)) &crpropa::Vector3<float>::operator/=, "C++: crpropa::Vector3<float>::operator/=(const float &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("__imod__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &)) &crpropa::Vector3<float>::operator%=, "C++: crpropa::Vector3<float>::operator%=(const class crpropa::Vector3<float> &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("__imod__", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const float &)) &crpropa::Vector3<float>::operator%=, "C++: crpropa::Vector3<float>::operator%=(const float &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));
		cl.def("assign", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const class crpropa::Vector3<float> &)) &crpropa::Vector3<float>::operator=, "C++: crpropa::Vector3<float>::operator=(const class crpropa::Vector3<float> &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("v"));
		cl.def("assign", (class crpropa::Vector3<float> & (crpropa::Vector3<float>::*)(const float &)) &crpropa::Vector3<float>::operator=, "C++: crpropa::Vector3<float>::operator=(const float &) --> class crpropa::Vector3<float> &", pybind11::return_value_policy::automatic, pybind11::arg("f"));

		cl.def("__str__", [](crpropa::Vector3<float> const &o) -> std::string { std::ostringstream s; using namespace crpropa; s << o; return s.str(); } );

		{ // crpropa::Vector3<float>::(anonymous union at /home/jnb/Repos/CRPDEV/include/crpropa/Vector3.h:38:2) file: line:38

			{ // crpropa::Vector3<float>::(anonymous union)::(anonymous struct at /home/jnb/Repos/CRPDEV/include/crpropa/Vector3.h:39:3) file: line:39
				cl.def_readwrite("x", &crpropa::Vector3<float>::x);
				cl.def_readwrite("y", &crpropa::Vector3<float>::y);
				cl.def_readwrite("z", &crpropa::Vector3<float>::z);
			}

		}

	}
}
