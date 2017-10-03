#ifndef _MAINEFFECTS3ANALYZERH_
#define _MAINEFFECTS3ANALYZERH_

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
#include "ddace_config.h"
#endif /* HAVE_CONFIG_H */

#include <cstdlib>
#include <string>
#include <sstream>

/**
 * Changes:	Added consts to make functions safer.
 *		moved ctors from .cpp to .h
**/

class ColumnHeader {
public:

	// Added const. Also set variables in .h file instead of .cpp
        static std::string FACTOR;
        static std::string RESPONSE;


	ColumnHeader() {};
        ColumnHeader(std::string n, std::string u, std::string f) :
		title(n), 
		abbreviatedTitle(computeAbbreviatedTitle(n)), 
		factorOrResponse(f),
		name(n), 
		units(u)
	 {};

        std::string getTitle() const { return title; };
        void setTitle(const std::string t) { title = t; };

        std::string getAbbreviatedTitle() const { return abbreviatedTitle; };

        void setAbbreviatedTitle(const std::string t) { abbreviatedTitle = t;};

        std::string getFactorOrResponse() const { return factorOrResponse; };
        void setFactorOrResponse(const std::string f) { factorOrResponse = f; }

        std::string computeAbbreviatedTitle(std::string title);

private:
        std::string title;
        std::string abbreviatedTitle;
        std::string factorOrResponse;

        std::string name;
        std::string units;
};


#endif
