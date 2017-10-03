#include "DDaceXMLReader.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
using namespace std;
#include "DDace.h"
#include "SmartPtr.inl"
#include "XMLUtils.h"
#include "FileInputSource.h"
#include "XMLParser.h"

DDaceXMLReader::DDaceXMLReader(const String& filename)
	: handler_()
{
	try
		{
			XMLUtils::init();
			init(filename);
		}
	catch(ExceptionBase& e)
		{
			e.trace("in DDaceXMLReader ctor");
		}
}

void DDaceXMLReader::init(const String& filename)
{
	try
		{
			if (DDace::getComm().getRank() != 0) return;
			
			XMLParser parser;
			FileInputSource src(filename);
			
			parser.setDocumentHandler(&handler_);
			
			parser.parse(src);
		}
	catch(ExceptionBase& e)
		{
			e.trace("in DDaceXMLReader::init()");
		}
}

	
	

	
	
