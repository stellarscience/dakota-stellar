#ifndef _DATAVLAUEH_
#define _DATAVALUEH_

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
#include "ddace_config.h"
#endif /* HAVE_CONFIG_H */

#include <string>
#include <sstream>


/**
 *  Changes: 	moved a lot of code from .cpp file to .h file
 *  		changed "get" functions to const
 *		changed to more proper class ctors, using : instead of =
*/

class DataValue {
public:
        static std::string EMPTY;
        static std::string STRING;
        static std::string INTEGER;
        static std::string DOUBLE;
        

	DataValue() {};
        DataValue(std::string s) : dataType(STRING), stringValue(s) {};
        DataValue(int i) : dataType(INTEGER), intValue(i) {};
        DataValue(double d) : dataType(DOUBLE), doubleValue(d) {};        

        std::string getDataType() const { return dataType; };
        void setDataType(std::string dt) { dataType = dt; };

        std::string getStringValue() const { return stringValue; };

        int getIntValue() const { return intValue; };

        double getDoubleValue() const { return doubleValue; };

        bool equals(DataValue dataValue);

        std::string toString();

private:
        std::string dataType;
        std::string stringValue;
        int intValue;
        double doubleValue;
	

};



#endif
