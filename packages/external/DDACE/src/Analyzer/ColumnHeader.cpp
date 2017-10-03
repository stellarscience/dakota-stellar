#include "ColumnHeader.h"

std::string ColumnHeader::FACTOR = "factor";
std::string ColumnHeader::RESPONSE = "response";
/**
std::string ColumnHeader::getTitle(){
    return(title);
}

void ColumnHeader::setTitle(std::string t) {
    title = t;
}

std::string ColumnHeader::getAbbreviatedTitle() {
    return(abbreviatedTitle);
}


void ColumnHeader::setAbbreviatedTitle(std::string abbrTitle) {
    abbreviatedTitle = abbrTitle;
}

std::string ColumnHeader::getFactorOrResponse(){
    return(factorOrResponse);
}


void ColumnHeader::setFactorOrResponse(std::string f){
    factorOrResponse = f;
}
*/
std::string ColumnHeader::computeAbbreviatedTitle(std::string title){
    std::ostringstream ss;
    char *ptr = (char *)title.c_str();
    while (*ptr != '\0'){
         if (isupper(*ptr)) {
	     ss << std::string(ptr,1);
         }
         ptr++;
    }
    return(ss.str());
}


