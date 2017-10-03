#include "DDaceInputReader.h"
#include "StrUtils.h"

DDaceInputReader::DDaceInputReader(const Array<String>& cmds)
	: DDaceReaderBase(cmds)
{;}

DDaceInputReader::DDaceInputReader(const String& filename)
	: DDaceReaderBase(filename)
{;}

void DDaceInputReader::getOutputNames(Array<String>& outputNames) const 
{
  outputNames.resize(0);
  for (int i=0; i<tokens_.length(); i++)
    {
      if (tokens_[i][0].allCaps() == "RETURN")
				{
					String big = reassembleFromTokens(tokens_[i], 1);
					splitList(big, outputNames);
					break;
				}
			if (tokens_[i][0].allCaps() == "OUTPUT")
				{
					if (tokens_[i].length() == 2)
						{
							outputNames.append(tokens_[i][1]);
						}
					else
						{
							fatalError("DDaceInputReader::getOutputNames bad OUTPUT line");
						}
				}
		}
}

void DDaceInputReader::getArchiveFilename(String& archiveFilename) const 
{
  if (!lookup("ARCHIVE", archiveFilename))
    {
      archiveFilename = "ddaceArchive";
    }
}



#ifdef BLAP
void DDaceInputReader::splitList(const String& big, Array<String>& list) const
{
	if (big.subString(0,1)!="[") 
		{
			list.resize(1);
			list[0] = big;
			return;
		}
	
	int parenDepth = 0;
	int localCount = 0;
	String tmp(big);
	list.resize(0);

	// start at 1 to ignore '[';
	
	for (int i=1; i<big.length(); i++)
		{
			if (big[i]=='(') parenDepth++;
			if (big[i]==')') parenDepth--;
			if (big[i]==']') 
				{
					tmp[localCount]='\0'; 
					list.append(tmp);
					break;
				}
			if (big[i]==',' && parenDepth==0)
				{
					tmp[localCount]='\0';
					list.append(tmp);
					tmp = big;
					localCount = 0;
					continue;
				}
			tmp[localCount] = big[i];
			localCount++;
		}
}
#endif						
	



