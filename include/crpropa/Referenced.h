#ifndef CRPROPA_REFERENCED_H
#define CRPROPA_REFERENCED_H

#include <cstddef>
#include <type_traits>

#ifdef DEBUG
#include <iostream>
#include <typeinfo>
#endif

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/**
 @class Referenced
 @brief Base class for reference counting

 A form of memory management is needed to prevent memory leaks when using MPC in Python via SWIG.
 This base class enables reference counting.
 Every reference increases the reference counter, every dereference decreases it.
 When the counter is decreased to 0, the object is deleted.
 Candidate, Module, MagneticField and Source inherit from this class
 */
class Referenced {
public:

	inline Referenced() :
			_referenceCount(0) {
	}

	inline Referenced(const Referenced&) :
			_referenceCount(0) {
	}

	inline Referenced& operator =(const Referenced&) {
		return *this;
	}

	inline size_t addReference() const {
		int newRef;
#if defined(OPENMP_3_1)
		#pragma omp atomic capture
		{newRef = _referenceCount++;}
#elif defined(__GNUC__)
		newRef = __sync_add_and_fetch(&_referenceCount, 1);
#else
		#pragma omp critical(newRef)
		{newRef = _referenceCount++;}
#endif
		return newRef;
	}

	inline size_t removeReference() const {
#ifdef DEBUG
		if (_referenceCount == 0)
			std::cerr
					<< "WARNING: Remove reference from Object with NO references: "
					<< typeid(*this).name() << std::endl;
#endif
		int newRef;
#if defined(OPENMP_3_1)
		#pragma omp atomic capture
		{newRef = _referenceCount--;}
#elif defined(__GNUC__)
		newRef = __sync_sub_and_fetch(&_referenceCount, 1);
#else
		#pragma omp critical(newRef)
		{newRef = _referenceCount--;}
#endif

		if (newRef == 0) {
			delete this;
		}
		return newRef;
	}

	int removeReferenceNoDelete() const {
		return --_referenceCount;
	}

	inline size_t getReferenceCount() const {
		return _referenceCount;
	}

protected:

	virtual inline ~Referenced() {
#ifdef DEBUG
		if (_referenceCount)
			std::cerr << "WARNING: Deleting Object with references: "
					<< typeid(*this).name() << std::endl;
#endif
	}

	mutable size_t _referenceCount;
};

inline void intrusive_ptr_add_ref(Referenced* p) {
	p->addReference();
}
inline void intrusive_ptr_release(Referenced* p) {
	p->removeReference();
}

// helper Reference class functions, which function to use is decided at compile time!
template<class T>
typename std::enable_if_t<std::is_base_of<Referenced, T>::value>
inline addReferenceIf(T* ptr) noexcept {
	if (ptr) ptr->addReference();
}
template<class T>
typename std::enable_if_t<!std::is_base_of<Referenced, T>::value>
inline addReferenceIf(T* ptr) noexcept {}

template<class T>
typename std::enable_if_t<std::is_base_of<Referenced, T>::value>
inline removeReferenceIf(T* ptr) noexcept {
	if (ptr) ptr->removeReference();
}
template<class T>
typename std::enable_if_t<!std::is_base_of<Referenced, T>::value>
inline removeReferenceIf(T* ptr) noexcept {}

template<class T>
typename std::enable_if_t<std::is_base_of<Referenced, T>::value>
inline removeReferenceNoDeleteIf(T* ptr) noexcept {
	if (ptr) ptr->removeReferenceNoDelete();
}
template<class T>
typename std::enable_if_t<!std::is_base_of<Referenced, T>::value>
inline removeReferenceNoDeleteIf(T* ptr) noexcept {}

/**
 @class ref_ptr
 @brief Referenced pointer

 This class manages the memory of objects that inherit the class Referenced.
 Every other object is not managed but instead treated as raw pointer and should
 be managed by the original variable owner.
 */
template<class T>
class ref_ptr {
public:
	typedef T element_type;

	/** Default Constructor */
	ref_ptr() : _ptr(0) {}
	#ifndef SWIG 
	/** Constructor from reference
	 * Use this constructor if you want to use a stack pointer,
	 * so every pointer not created with a new operator.
	 */
	ref_ptr(T& ref) : _ptr(&ref) {}
	#endif
	/** Constructor from pointer
	 * Should not be used with stack pointer,
	 * so every pointer not created with a new operator.
	 */
	ref_ptr(T* ptr) : _ptr(ptr) {
		addReferenceIf(_ptr);
	}
	/** Constructor from ref_ptr with same template type */
	ref_ptr(const ref_ptr& rp) :
			_ptr(rp._ptr) {
		addReferenceIf(_ptr);
	}
	/** Constructor from ref_ptr with other template type, must be convertable */
	template<class Other> ref_ptr(const ref_ptr<Other>& rp) :
			_ptr(rp._ptr) {
		addReferenceIf(_ptr);
	}

	/** destructor */
	~ref_ptr() {
		removeReferenceIf(_ptr);
		_ptr = 0;
	}

	/** Assign operator for ref_ptr with same template type */
	ref_ptr& operator =(const ref_ptr& rp) {
		assign(rp);
		return *this;
	}
	/** Assign operator for ref_ptr with other template type, must be convertable */
	template<class Other> ref_ptr& operator =(const ref_ptr<Other>& rp) {
		assign(rp);
		return *this;
	}
	/** Assign operator for pointer, do not assign stack pointer */
	inline ref_ptr& operator =(T* ptr) {
		if (_ptr == ptr)
			return *this;
		T* tmp_ptr = _ptr;
		_ptr = ptr;
		addReferenceIf(_ptr);
		removeReferenceIf(tmp_ptr);
		return *this;
	}
	#ifndef SWIG
	/** Assign operator for stack references */
	inline ref_ptr& operator =(T& ref){
		if (_ptr == &ref)
			return *this;
		_ptr = &ref;
		return *this;
	}
	#endif

	operator T*() const {
		return _ptr;
	}

	T& operator*() const {
		return *_ptr;
	}
	T* operator->() const {
		return _ptr;
	}
	T* get() const {
		return _ptr;
	}

	bool operator!() const {
		return _ptr == 0;
	} // not required
	bool valid() const {
		return _ptr != 0;
	}

	T* release() {
		T* tmp = _ptr;
		removeReferenceNoDeleteIf(_ptr);
		_ptr = 0;
		return tmp;
	}

	void swap(ref_ptr& rp) {
		T* tmp = _ptr;
		_ptr = rp._ptr;
		rp._ptr = tmp;
	}

private:

	template<class Other> void assign(const ref_ptr<Other>& rp) {
		if (_ptr == rp._ptr)
			return;
		T* tmp_ptr = _ptr;
		_ptr = rp._ptr;
		addReferenceIf(_ptr);
		removeReferenceIf(tmp_ptr);
	}

	template<class Other> friend class ref_ptr;

	T* _ptr;
};

template<class T> inline
void swap(ref_ptr<T>& rp1, ref_ptr<T>& rp2) {
	rp1.swap(rp2);
}

template<class T> inline T* get_pointer(const ref_ptr<T>& rp) {
	return rp.get();
}

template<class T, class Y> inline ref_ptr<T> static_pointer_cast(
		const ref_ptr<Y>& rp) {
	return static_cast<T*>(rp.get());
}

template<class T, class Y> inline ref_ptr<T> dynamic_pointer_cast(
		const ref_ptr<Y>& rp) {
	return dynamic_cast<T*>(rp.get());
}

template<class T, class Y> inline ref_ptr<T> const_pointer_cast(
		const ref_ptr<Y>& rp) {
	return const_cast<T*>(rp.get());
}

/** @}*/
} // namespace crpropa

#endif // CRPROPA_REFERENCED_H
