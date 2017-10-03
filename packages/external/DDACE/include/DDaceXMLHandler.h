#ifndef DDACEXMLHANDLER_H
#define DDACEXMLHANDLER_H

#include "DDaceSampler.h"
#include "DDace.h"
#include "XMLHandlerBase.h"
#include "DDaceBool.h"
#include "DDaceSamplePoint.h"


typedef hash_set<std::string, hash<std::string>, eqstr>		Hashtable;

/**
 * DDaceXMLHandler defines the interface that is used to handle the
 * XML interface.  This class is derived from a CPPUtilities class:
 * XMLHandlerBase
 * */

class DDaceXMLHandler : public XMLHandlerBase
{
 public:
  DDaceXMLHandler();
  virtual ~DDaceXMLHandler(){;}
  
  DDace createObject();
  
  std::vector<std::string> getVarNames() const ;
  std::vector<std::string> getOutputNames() const ;
  std::string getArchiveFilename() const;
  virtual void startElement(const std::string& name, 
			    const Hashtable<String, String>& attributes);
  
  virtual void endElement(const std::string& name);
  
  virtual void characters(const std::string& chars, 
			  const unsigned int length);
 protected:
  
  DDaceSampler createSampler() const ;

  void init(const std::string& filename);
  static Hashtable samplerTypes_;
  std::string samplerName_;
  Hashtable samplerParams_;
  std::vector<std::string> varNames_;
  std::vector<Distribution> varDist_;
  std::vector<std::string> outputNames_;
  std::string archiveFilename_;
  int currentTag_;
  std::string charStr_;
  std::vector<DDaceSamplePoint> pts_;
  std::vector<std::vector<double> > results_;
  std::vector<DDaceRunStatus> status_;
  

	// status variables
  bool isSamplePt_;
  bool isSampleResult_;
  /** 
   * if a "SampleResult" is in the XML input, the XML file is an archive
   * file.  If it's an archive file, we need to read in the status values
   * for each run.  Otherwise, we need to leave the status values to 
   * "Run not started" for the creation of the DDace object.  Thus we
   * need to set a flag for whether or not the file is an archive.
   */
  bool isArchiveFile_;    
  bool getPtsComplete_;
};

struct eqstr
{
  bool operator()(const std::string s1, const string s2) const
  {
    return s1.compare(s2) == 0;
  }
};


#endif
