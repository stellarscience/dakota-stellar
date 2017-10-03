#ifndef MEMORYEXCEPTION_H
#define MEMORYEXCEPTION_H

#include "ExceptionBase.h"
#include <iostream>


/**
 * \ingroup Exception
 * Thrown when a new or malloc fails.
 */
class MemoryException : public ExceptionBase
{
 public:
	
	/**
	 * Panic: print a message and terminate.
	 */
	static void raise(const std::string& msg);
};

#endif
