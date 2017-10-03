#ifndef EXCEPTIONBASE_H
#define EXCEPTIONBASE_H
#include <string>
#include <iostream>
/**#ifdef ANSI_HEADERS
#include <cstdlib>
#else
#include <stdlib.h>
#endif
*/



/**
 * \ingroup Exception
 * Base class for exceptions.
 */
class ExceptionBase
{
 public:
	ExceptionBase() {;}
	ExceptionBase(const std::string& mesg);
	ExceptionBase(const char* mesg);

	/**
	 * Print a message to stderr (through the IO server)
	 */
	virtual void print() const ;


	virtual const std::string message() const {return mesg_;}
	virtual const std::string toString() const {return mesg_;}


	/**
	 * Throw an ExceptionBase with the given message
	 */
	static void raise(const std::string& msg)  ;
	static void raise(const char* msg)  ;

	/**
	 * Throw a TraceBack, allowing the exception to be traced back through
	 * the call stack.
	 */

	//virtual void trace(const std::string& msg) const ;


	/**
	 * print a message, shut down the PMachine, and exit
	 */
	virtual void bailout() const ;

	/**
	 * exit immediately. Call in extreme circumstances such as
	 * out-of-memory errors.
	 */
	static void panic(const std::string msg);

 protected:
	std::string mesg_;
};



#endif
