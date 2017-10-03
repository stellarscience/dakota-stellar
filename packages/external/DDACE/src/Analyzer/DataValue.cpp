#include "DataValue.h"

std::string DataValue::EMPTY = "empty";
std::string DataValue::STRING = "string";
std::string DataValue::INTEGER = "integer";
std::string DataValue::DOUBLE = "double";

/**
int DataValue::getIntValue() {
    return intValue;
}

double DataValue::getDoubleValue() {
    return doubleValue;
}
*/

bool DataValue::equals(DataValue other){
    if (getDataType() != other.getDataType())
        return false;
    if (getDataType() == STRING)
        return (getStringValue() == other.getStringValue());
    if (getDataType() == INTEGER)
        return (getIntValue() == other.getIntValue());
    if (getDataType() ==DOUBLE)
        return (getDoubleValue() == other.getDoubleValue());
    return false;
}

std::string DataValue::toString() {
    std::ostringstream ss;
    ss << "DataValue:";
    ss << "dataType=" << getDataType() << " ";
    if (getDataType()==STRING)
        ss << "value=" << getStringValue();
    if (getDataType()==INTEGER){
        ss << "value=" << getIntValue();
    }
    if (getDataType()==DOUBLE)
        ss << "value=" << getDoubleValue();
    return(ss.str());
}
