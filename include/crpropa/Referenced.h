#ifndef CRPROPA_REFERENCED_H
#define CRPROPA_REFERENCED_H

#include <cstddef>
#include <memory>

namespace crpropa {
/**
 * \addtogroup Core
 * @{
 */

/**
 @class ref_ptr
@brief Referenced pointer
*/
template<class T>
class ref_ptr {
	public:
	typedef T element_type;

	ref_ptr() : _ptr(0), _shared_ptr(0) {}
	ref_ptr(T* ptr) : _ptr(ptr), _shared_ptr(0), _is_shared(false) {}
	ref_ptr(const ref_ptr& rp) {
		_shared_ptr = rp._shared_ptr;
		_is_shared = rp._is_shared;
		_ptr = rp._ptr;
	}
	template<class Other> ref_ptr(const ref_ptr<Other>& rp) {
		_shared_ptr = rp._shared_ptr;
		_is_shared = rp._is_shared;
		_ptr = rp._ptr;
	}
	template<class Other> ref_ptr(const std::shared_ptr<Other>& shared_ptr) {
		_shared_ptr = shared_ptr;
		_is_shared = true;
		_ptr = 0;
	}

	~ref_ptr() {
		_ptr = 0;
		_shared_ptr = 0;
	}

	ref_ptr& operator =(const ref_ptr& rp) {
		assign(rp);
		return *this;
	}

	template<class Other> ref_ptr& operator =(const ref_ptr<Other>& rp) {
		assign(rp);
		return *this;
	}


	inline ref_ptr& operator =(long int ptr) {
		_shared_ptr = ptr;
		_ptr = ptr;
		return *this;
	}

	// operator T*() const {
	// 	return get();
	// }

	T& operator*() const {
		if (_is_shared)
			return *_shared_ptr;
		else
			return *_ptr;
	}
	T* operator->() const {
		if (_is_shared)
			return _shared_ptr.get();
		else
			return _ptr;
	}

	bool operator==(const ref_ptr& rp) const {
		return get()==rp.get();
	}

	T* get() const {
		if (_is_shared)
			return _shared_ptr.get();
		else
			return _ptr;
	}

	/** Returns stored shared_ptr
	 * In case of a stored normal pointer this only return NULL;
	 * User need to check if returned pointer is valid
	 */
	std::shared_ptr<T> get_shared() const {
		return _shared_ptr;
	}

	bool valid() const {
		if (_is_shared)
			return _shared_ptr != 0;
		else
			return _ptr != 0;
	}

	/** swaps two reference pointer
	 * Returns zero if unsuccessfull
	 */
	int swap(ref_ptr& rp) {
		if(rp._is_shared!=_is_shared)
			return 0;
		if(_is_shared){
			_shared_ptr.swap(rp._shared_ptr);
		} else {
			T* tmp = _ptr;
			_ptr = rp._ptr;
			rp._ptr = tmp;
		}
		return 1;
	}

	bool is_shared() const{
		return _is_shared;
	}

	private:

	template<class Other> void assign(const ref_ptr<Other>& rp) {
		_shared_ptr = rp._shared_ptr;
		_is_shared = rp._is_shared;
		_ptr = rp._ptr;
	}

	template<class Other> friend class ref_ptr;

	T* _ptr = 0;
	std::shared_ptr<T> _shared_ptr = 0;
	bool _is_shared=false;
};

template<class T> inline
void swap(ref_ptr<T>& rp1, ref_ptr<T>& rp2) {
	rp1.swap(rp2);
}

template<class T, class Y> 
inline ref_ptr<T> static_pointer_cast(const ref_ptr<Y>& rp) {
	if (rp.is_shared())
		return std::static_pointer_cast<T>(rp.get_shared());
	else
		return static_cast<T*>(rp.get());
}

template<class T, class Y> 
inline ref_ptr<T> dynamic_pointer_cast(const ref_ptr<Y>& rp) {
	if (rp.is_shared())
		return std::dynamic_pointer_cast<T>(rp.get_shared());
	else
	return dynamic_cast<T*>(rp.get());
}

template<class T, class Y> 
inline ref_ptr<T> const_pointer_cast(const ref_ptr<Y>& rp) {
	if (rp.is_shared())
		return std::const_pointer_cast<T>(rp);
	else
	return const_cast<T*>(rp.get());
}

/** @}*/
} // namespace crpropa

#endif // CRPROPA_REFERENCED_H
