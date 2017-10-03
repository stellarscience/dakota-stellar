#include "OptppSmartPtr.h"

#ifndef SMARTPTR_C
#define SMARTPTR_C

namespace OPTPP {

template <class T> 
inline SmartPtr<T>::~SmartPtr() {
	if (refCount_ != 0 && --(*refCount_) == 0) {
		delete ptr_;
		ptr_ = 0;
		delete refCount_;
		refCount_ = 0;
	}
}

template <class T> 
inline const SmartPtr<T>& SmartPtr<T>::operator =(const SmartPtr<T>& rhs) {
	if (ptr_ != rhs.ptr_) {
		if (refCount_ != 0 && --(*refCount_) == 0) {
			delete ptr_;
			delete refCount_;
		}
		ptr_ = rhs.ptr_;
		refCount_ = rhs.refCount_;
		if (refCount_ != 0)
			++(*refCount_);
	}
	return *this;
}

} // namespace OPTPP

#endif

