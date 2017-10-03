#ifndef DDACEXMLREADER_H
#define DDACEXMLREADER_H


#include "DDaceXMLHandler.h"

/**
 *
 * DDaceXMLReader is the class for interfacing to the XML readers.  
 * The constructor takes one parameter: the name of the XML file to read
 */
class DDaceXMLReader { public: DDaceXMLReader(const std::string& filename);

	DDace createObject() {return handler_.createObject();}
	std::vector<std::string> getVarNames() const {return handler_.getVarNames();}
	std::vector<std::string> getOutputNames() const 
	  {return handler_.getOutputNames();}
	std::string getArchiveFilename()const
	  {return handler_.getArchiveFilename();}
	
 protected:
	void init(const std::string& filename);

	DDaceXMLHandler handler_;
};

#endif

