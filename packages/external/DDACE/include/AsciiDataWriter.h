#ifndef ASCIIDATAWRITER_H
#define ASCIIDATAWRITER_H

#include "DataWriter.h"

/** 
 * AsciiDataWriter writes selected archived data or data in a DDace object to
 * an ascii formatted text file.
 */  

class AsciiDataWriter : public DataWriterBase
{
 public:

  /** construct data from the DDace object containing the archive info */
  AsciiDataWriter(const DDace& ddaceObj)  
    : DataWriterBase(ddaceObj)
    {;}
  
  /** construct data from the DDace  archive file */
  AsciiDataWriter(const String& archiveName)
    : DataWriterBase(archiveName)
    {;}
  
  virtual ~AsciiDataWriter(){;}
  
  DataWriterBase* clone() const;

  /** write the plot data to a file */
  void writeToFile(const String& filename) const;
  
  /** write the plot data to a stream */
  void write(ostream& os) const;
  
 private:

};
#endif
