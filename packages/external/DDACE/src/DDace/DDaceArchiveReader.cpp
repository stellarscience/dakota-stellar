#include "DDaceArchiveReader.h"
#include <cctype>
#include "StrUtils.h"

DDaceArchiveReader::DDaceArchiveReader(const Array<String>& cmds)
	: DDaceReaderBase(cmds), nSamples_(0), nOutputs_(0)
{
	DDaceSampler sampler;
	getSampler(sampler);
	nSamples_ = sampler.nSamples();

	Array<String> outputNames;
	getOutputNames(outputNames);
	nOutputs_ = outputNames.length();
}

DDaceArchiveReader::DDaceArchiveReader(const String& filename)
	: DDaceReaderBase(filename), nSamples_(0), nOutputs_(0)
{
	DDaceSampler sampler;
	getSampler(sampler);
	nSamples_ = sampler.nSamples();

	Array<String> outputNames;
	getOutputNames(outputNames);
	nOutputs_ = outputNames.length();
}


void DDaceArchiveReader::getOutputNames(Array<String>& outputNames) const 
{
	outputNames.resize(0);
	for (int i=0; i<tokens_.length(); i++)
    {
      if (tokens_[i].length() < 2) continue;
      if (tokens_[i][0].allCaps() == "OUTPUT")
				{
					if (tokens_[i].length()==2)
						{
							outputNames.append(tokens_[i][1]);
						}
					else
						{
							outputNames.append(reassembleFromTokens(tokens_[i], 1));
						}
				}
    }
}


void DDaceArchiveReader::getArchivedData(Array<DDaceSamplePoint>& pts,
																				 Array<Array<double> >& results,
																				 Array<DDaceRunStatus>& status) const
{
	int i = 0;
	while (i < tokens_.length())
    {
			const Array<String>& toks = tokens_[i++];
      if (toks.length() < 2) continue;
			if (toks[0].allCaps() == "BEGIN" 
					&& toks[1].allCaps() == "DATA")
				{
					readResults(pts, results, status, i);
					break;
				}
		}
}

void DDaceArchiveReader::readResults(Array<DDaceSamplePoint>& pts,
																		 Array<Array<double> >& results,
																		 Array<DDaceRunStatus>& status,
																		 int startLine) const 
{
	int line = startLine;
	Array<double> x(nInputs_);

	pts.resize(nSamples_);
	results.resize(nSamples_);
	status.resize(nSamples_);
	
	for (int i=0; i<nSamples_; i++)
		{
			int index = atoi(tokens_[line++][0].cString());
			if (index != i) fatalError("DDaceArchiveReader::readResults index mismatch");
			int j;
			for (j=0; j<nInputs_; j++)
				{
					x[j] = atof(tokens_[line++][0].cString());
				}
			pts[i] = DDaceSamplePoint(i, x);
			
			results[i].resize(nOutputs_);
			for (j=0; j<nOutputs_; j++)
				{
					const Array<String>& toks = tokens_[line++];
					if (isalpha(toks[0][0]))
						{
							results[i][j] = 0.0;
							status[i] = (DDaceRunStatus) atoi(toks[1].cString());
						}
					else
						{
							results[i][j] = atof(toks[0].cString());
							status[i] = DDaceRunOK;
						}
				}
		}
}
					





