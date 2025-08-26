#ifndef CRPROPA_REFERENCED_H
#define CRPROPA_REFERENCED_H

#include <cstddef>

#ifdef DEBUG
#include <iostream>
#include <typeinfo>
#endif

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#include <cuda_runtime.h>
#else
#define CUDA_CALLABLE_MEMBER
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

	/** Default Constructor
	 * just sets _referenceCount to 0
	*/
	inline Referenced() :
			_referenceCount(0) {
	}

	inline Referenced(const Referenced&) :
			_referenceCount(0) {
	}

	/// Equal operator, does only return this class
	inline Referenced& operator =(const Referenced&) {
		return *this;
	}

	/// Increases reference count by one (atomic) and returns new value
	CUDA_CALLABLE_MEMBER inline size_t addReference() const {
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

	/** Decreases reference count by one (atomic) and returns new value
	 * When the reference count reaches 0 then delete this class
	 * Exception when reference equals zero before substracting
	 */
	CUDA_CALLABLE_MEMBER inline size_t removeReference() const {
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

	/// Decreases reference count without eventually deleting class
	CUDA_CALLABLE_MEMBER int removeReferenceNoDelete() const {
		return --_referenceCount;
	}

	/// Getter function for reference count
	CUDA_CALLABLE_MEMBER inline size_t getReferenceCount() const {
		return _referenceCount;
	}

protected:

	/// Destructor
	CUDA_CALLABLE_MEMBER virtual inline ~Referenced() {
		#ifdef DEBUG
			if (_referenceCount)
				std::cerr << "WARNING: Deleting Object with references: "
						<< typeid(*this).name() << std::endl;
		#endif
	}

	mutable size_t _referenceCount;
};

/** Calls addReference function for given Referenced pointer */
CUDA_CALLABLE_MEMBER inline void intrusive_ptr_add_ref(Referenced* p) {
	p->addReference();
}
/** Calls removeReference function for given Referenced pointer */
CUDA_CALLABLE_MEMBER inline void intrusive_ptr_release(Referenced* p) {
	p->removeReference();
}

/**
 @class ref_ptr
 @brief Referenced pointer
 * The standard smart pointer class for crpropa objects.
 */
template<class T>
class ref_ptr {
public:
	typedef T element_type; ///< custom name for template type T

	/** Default constructor
	 * sets the pointer to NULL
	 */
	CUDA_CALLABLE_MEMBER ref_ptr() : _ptr(0) {}
	/** Constructor
	 * Takes a raw pointer and assigns it to _ptr, then calls pointer->addReference()
	 * \param ptr A pointer to any Referenced object or object which inherits Referenced
	 */
	CUDA_CALLABLE_MEMBER ref_ptr(T* ptr) :	_ptr(ptr) {
		if (_ptr)
			_ptr->addReference();
	}
	/** Constructor
	 * Takes a ref_ptr object and assigns its _ptr to this _ptr,
	 * then calls pointer->addReference()
 	 * \param rp ref_ptr object
	 */
	CUDA_CALLABLE_MEMBER ref_ptr(const ref_ptr& rp) : _ptr(rp._ptr) {
		if (_ptr)
			_ptr->addReference();
	}
	/** Constructor
	 * Takes any template instance of ref_ptr and assigns its _ptr to his _ptr,
	 * then calls pointer->addReference()
	 * \param rp Any ref_ptr<typename> object
	 */
	template<class Other> CUDA_CALLABLE_MEMBER ref_ptr(const ref_ptr<Other>& rp) : _ptr(rp._ptr) {
		if (_ptr)
			_ptr->addReference();
	}

	/** Destructor
	 * calls _ptr->removeReference() which deletes referenced object if it sets ref counter to 0
	 * and then sets _ptr to NULL
	 */
	CUDA_CALLABLE_MEMBER ~ref_ptr() {
		if (_ptr)
			_ptr->removeReference();
		_ptr = 0;
	}

	/** Equal operator
	 * Calls assign function with given ref_ptr
	 */
	CUDA_CALLABLE_MEMBER ref_ptr& operator =(const ref_ptr& rp) {
		assign(rp);
		return *this;
	}
	/** Equal operator
	 * Calls assign function with given ref_ptr<typename>
	 */
	template<class Other> CUDA_CALLABLE_MEMBER ref_ptr& operator =(const ref_ptr<Other>& rp) {
		assign(rp);
		return *this;
	}
	/** Equal operator
	 * Calls assign function with given raw pointer
	 * First checks if the stored pointer is allready the given pointer
	 * After that does the same as assign
	 */
	CUDA_CALLABLE_MEMBER inline ref_ptr& operator =(T* ptr) {
		if (_ptr == ptr)
			return *this;
		T* tmp_ptr = _ptr;
		_ptr = ptr;
		if (_ptr)
			_ptr->addReference();
		if (tmp_ptr)
			tmp_ptr->removeReference();
		return *this;
	}

	CUDA_CALLABLE_MEMBER operator T*() const {
		return _ptr;
	}

	/// Dereference operator
	CUDA_CALLABLE_MEMBER T& operator*() const {
		return *_ptr;
	}
	/// Arror operator
	CUDA_CALLABLE_MEMBER T* operator->() const {
		return _ptr;
	}
	/// Getter function for the stored pointer
	CUDA_CALLABLE_MEMBER T* get() const {
		return _ptr;
	}

	/// Checks if _ptr == NULL
	CUDA_CALLABLE_MEMBER bool operator!() const {
		return _ptr == 0;
	} // not required
	/// Checks if _ptr != NULL
	CUDA_CALLABLE_MEMBER bool valid() const {
		return _ptr != 0;
	}

	/// Saves _ptr to tmp_ptr, calls _ptr->removeReferenceNoDelete() and returns tmp_ptr
	CUDA_CALLABLE_MEMBER T* release() {
		T* tmp = _ptr;
		if (_ptr)
			_ptr->removeReferenceNoDelete();
		_ptr = 0;
		return tmp;
	}

	/// swaps the stored _ptr with pointer stored in rp and vice versa
	CUDA_CALLABLE_MEMBER void swap(ref_ptr& rp) {
		T* tmp = _ptr;
		_ptr = rp._ptr;
		rp._ptr = tmp;
	}

	// // not so trivial like this, classes do not exist in cuda, only functions, base types, structs etc.
	// // probably keep standard Referenced.h and implement cuda HostToDevice and cuda DeviceToHost for each class
	// #ifdef __CUDACC__
	// CUDA_CALLABLE_MEMBER void hostToDevice(){
	// 	T* tmp = NULL;  // create temporary pointer
	// 	cudaMalloc((void**)&tmp, sizeof(T));  // reserve device memory
	// 	cudaMemcpy(tmp, _ptr, sizeof(T), cudaMemcpyHostToDevice); // copy data from _ptr to tmp
	// 	_ptr->removeReference();  // remove reference since now on device
	// 	_ptr = tmp;  // set _ptr to tmp
	// }

	// CUDA_CALLABLE_MEMBER void deviceToHost(){
	// 	T* tmp = NULL;  // assign host memory, do not need malloc since only one object expected
	// 	cudaMemcpy(tmp, _ptr, sizeof(T), cudaMemcpyDeviceToHost);
	// 	cudaFree(_ptr);  // free memory befor freeing _ptr
	// 	_ptr->removeReferenceNoDelete();  // remove reference since now on device but do not delete
	// 	_ptr = tmp;  // set _ptr to tmp
	// }
	// #endif

private:

	/** Function to assign other ref_ptr
	 * Takes any other ref_ptr<typename> and assigns its pointer to this _ptr
	 * if it is not the same already
	 * \param rp Any ref_ptr<typename> object
	 */
	template<class Other> CUDA_CALLABLE_MEMBER void assign(const ref_ptr<Other>& rp) {
		if (_ptr == rp._ptr)
			return;
		T* tmp_ptr = _ptr;
		_ptr = rp._ptr;
		if (_ptr)
			_ptr->addReference();
		if (tmp_ptr)
			tmp_ptr->removeReference();
	}

	template<class Other> friend class ref_ptr;  ///< other instances of ref_ptr are allowed to access privates

	T* _ptr;  ///< stored pointer
};

template<class T> CUDA_CALLABLE_MEMBER inline
void swap(ref_ptr<T>& rp1, ref_ptr<T>& rp2) {
	rp1.swap(rp2);
}

template<class T> CUDA_CALLABLE_MEMBER inline T* get_pointer(const ref_ptr<T>& rp) {
	return rp.get();
}

template<class T, class Y> CUDA_CALLABLE_MEMBER inline ref_ptr<T> static_pointer_cast(
		const ref_ptr<Y>& rp) {
	return static_cast<T*>(rp.get());
}

template<class T, class Y> CUDA_CALLABLE_MEMBER inline ref_ptr<T> dynamic_pointer_cast(
		const ref_ptr<Y>& rp) {
	return dynamic_cast<T*>(rp.get());
}

template<class T, class Y> CUDA_CALLABLE_MEMBER inline ref_ptr<T> const_pointer_cast(
		const ref_ptr<Y>& rp) {
	return const_cast<T*>(rp.get());
}

/** @}*/
} // namespace crpropa

#endif // CRPROPA_REFERENCED_H
