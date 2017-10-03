#ifndef SMARTPTR_H
#define SMARTPTR_H


//#include "BoolHack.h"
#include <string>
#include <stdexcept>

#ifndef inline
#define cond_inline 
#else
#define cond_inline inline
#endif



/**
 * \ingroup General
 * Templated reference-counting pointer class. 
 * @author Kevin Long
 */

template <class T>
class SmartPtr {
public:
	cond_inline SmartPtr(T* ptr = 0);
	cond_inline SmartPtr(const SmartPtr<T>& other);
	~SmartPtr();
	const SmartPtr<T>& operator =(const SmartPtr<T>& rhs);
	cond_inline T* operator ->();
	cond_inline bool isNull() const;
	cond_inline const T* operator ->() const;
	/** Read/write dereferencing */
	cond_inline T& operator *();
	/** Read-only dereferencing */
	cond_inline const T& operator *() const;
	cond_inline operator const T* () const;

protected:
	T* ptr_;
	int* refCount_;
};

template <class T> cond_inline
SmartPtr<T>::SmartPtr(T* ptr)
	: ptr_(ptr),
	refCount_(0)
{
	if (ptr_) {
		refCount_ = new int;
		if (refCount_ == 0)
			throw std::bad_alloc();
		*refCount_ = 1;
	}

}

template <class T> cond_inline
SmartPtr<T>::SmartPtr(const SmartPtr<T>& other)
	: ptr_(other.ptr_),
	refCount_(other.refCount_)
{
	if (refCount_ != 0)
		++(*refCount_);

}

template <class T>
SmartPtr<T>::~SmartPtr()
{
        if (refCount_ != 0 && --(*refCount_) == 0)
                {
                        delete ptr_;
                        ptr_ = 0;
                        delete refCount_;
                        refCount_ = 0;
                }
}

template <class T>
const SmartPtr<T>& SmartPtr<T>::operator =(const SmartPtr<T>& rhs)
{
        if (ptr_ != rhs.ptr_)
        {
                if (refCount_ != 0 && --(*refCount_) == 0)
                {
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

template <class T> cond_inline
bool SmartPtr<T>::isNull() const {
  return (ptr_ == 0);
}

template <class T> cond_inline
T* SmartPtr<T>::operator ->() {
	if (ptr_ == 0)
	//	throw std::runtime_error("SmartPtr<T>::operator ->() on null pointer");	
		throw std::bad_alloc();
	return ptr_;
}

template <class T> cond_inline
const T* SmartPtr<T>::operator ->() const {

	if (ptr_ == 0 ) {
	//	throw std::runtime_error("SmartPtr<T>::operator ->() on null pointer");
		throw std::bad_alloc();
//		throw std::runtime_error();
	}

	return ptr_;
}

template <class T> cond_inline
T& SmartPtr<T>::operator *() {
	if (ptr_ == 0)
		throw std::runtime_error("SmartPtr<T>::operator *() on null pointer");
	return *ptr_;
}

template <class T> cond_inline
const T& SmartPtr<T>::operator *() const {
	if (ptr_ == 0)
		throw std::runtime_error("SmartPtr<T>::operator *() on null pointer");
	return *ptr_;
}

template <class T> cond_inline
SmartPtr<T>::operator const T* () const {
	return ptr_;
}


#undef cond_inline
#endif
