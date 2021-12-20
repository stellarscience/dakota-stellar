/*  _________________________________________________________________________
 *
 *  UTILIB: A utility library for developing portable C++ codes.
 *  Copyright (c) 2008 Sandia Corporation.
 *  This software is distributed under the BSD License.
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *  For more information, see the README file in the top UTILIB directory.
 *  _________________________________________________________________________
 */

/**
 * \file _generic.h
 *
 * Defines and typedefs used everywhere.
 *
 * \author William E. Hart
 */

#ifndef utilib_generic_h
#define utilib_generic_h

#include <utilib/std_headers.h>
#include <sys/types.h>

#if defined(__cplusplus)
namespace utilib {
#endif

/**
 * \enum EnumDataOwned
 *
 * Ownership categories for objects with reference counts.
 */
#ifndef _ENUMDATAOWNED_
#define _ENUMDATAOWNED_
enum EnumDataOwned
{
  DataNotOwned=0,	/**< Data owned by some other object */
  DataOwned=1,          /**< Memory allocated by object itself */
  AcquireOwnership=1,	/**< Synonym for DataOwned */
  AssumeOwnership=2	/**< We own it now, but it comes from elsewhere */
                        /**< Once the object has been made this is      */
                        /**< identical to DataOwned                     */
};
#endif

/**
 * \enum OrderSense
 *
 * Bias used for comparisons in a dynamic ADT.
 */
#ifndef ORDERSENSE
#define ORDERSENSE
enum OrderSense
{
  increasing=1, /**< Order from least to greatest */
  decreasing=-1, /**< Order from greatest to least */
  minimal=1,
  maximal=-1
};
#endif

/**
 * \def BUF_SIZE
 *
 * A default size for buffers
 */
enum {BUF_SIZE=256 };

/**
 * \def ERR
 *
 * The default value of error values.
 */
enum { ERR = -999 };

#if defined(__cplusplus)
} // namespace end
#endif

/**
 * \def OK
 *
 * Value used to indicate that an operation worked.
 */
#if defined(__cplusplus)
namespace utilib {
constexpr int OK = 0;
}
#else
#ifndef OK
#define OK		0
#endif
#endif


/**
 * \def TRUE
 *
 * The boolean value for true.
 */
#if defined(__cplusplus)
namespace utilib {
constexpr int TRUE = 1;
}
#else
#ifndef TRUE
#define TRUE		1
#endif
#endif

/**
 * \def FALSE
 *
 * The boolean value for false.
 */
#if defined(__cplusplus)
namespace utilib {
constexpr int FALSE = 0;
}
#else
#ifndef FALSE
#define FALSE		0
#endif
#endif

/**
 * \def ON
 *
 * Used to incidate the on state.
 */
#if defined(__cplusplus)
namespace utilib {
constexpr int ON = 0;
}
#else
#ifndef ON
#define ON		1
#endif
#endif

/**
 * \def OFF
 *
 * Used to incidate the off state.
 */
#if defined(__cplusplus)
namespace utilib {
constexpr int OFF = 0;
}
#else
#ifndef OFF
#define OFF		0
#endif
#endif

/**
 * \def YES
 *
 * Used to incidate a yes response.
 */
#if defined(__cplusplus)
namespace utilib {
constexpr int YES = 1;
}
#else
#ifndef YES
#define YES		1
#endif
#endif

/**
 * \def NO
 *
 * Used to incidate a no response.
 */
#if defined(__cplusplus)
namespace utilib {
constexpr int NO = 0;
}
#else
#ifndef NO
#define NO		0
#endif
#endif

/**
 * \def NULL
 *
 * Defines the value of empty pointers.
 */
//#if defined(__cplusplus)
//namespace utilib {
//constexpr int NULL = 0;
//}
//#else
//#ifdef NULL
//#undef NULL		/* Always override the definition of NULL */
//#endif
//#define NULL		0
//#endif

/**
 * \def EOF
 *
 * The end-of-file value.
 */
//#ifndef EOF
//#define EOF		(-1)
//#endif

/**
 * \def PAUSE
 *
 * A macro that waits for the user to hit a key.
 */
#define PAUSE()	fflush(stdout); while(fgetc(stdin) == EOF);

/**
 * \def _(args)
 *
 * Used to provide a consistent definition for non-ansi C and ansi C.
 */
#if defined(UTILIB_STDC_HEADERS) || defined(__cplusplus)
#define _(args) args
#else
#define _(args) ()
#endif


#ifdef DEBUG			/* Debug defines to see if conflicts exist */
#define TRUE	1
#define FALSE	0
#define OK	0
#define ON	1
#define OFF	0
#define YES	1
#define NO	0
#define NULL	0
#define EOF	(-1)
#define ERR	-999
#endif


#if !defined(UTILIB_HAVE_STD) && !defined(__cplusplus)
/**
 * \typedef size_t
 *
 * The typedef for \a size_t arguments.
 */
#ifndef _SIZE_T
#define _SIZE_T
typedef unsigned size_t;
#endif
#endif

/**
 * \typedef VOID
 *
 * The void type is a \a char in standard C.
 */
#ifndef VOID
typedef char VOID;
#endif

/**
 * \typedef size_type
 *
 * Used to provide a consistent definition of the size_t type.
 */
typedef size_t size_type;

#endif
