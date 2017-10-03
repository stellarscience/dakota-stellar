#ifndef DDACEOALHSAMPLER_H
#define DDACEOALHSAMPLER_H

/*
 * Orthogonal Array based Latin Hypercube sampler.
 */

#include "DDaceSampler.h"
#include "UniformDistribution.h"
#include <vector>
#include <string>
#include <algorithm>
#include <ostream>


class DDaceOALHSampler : public DDaceSamplerBase
{

  class triple
  {
      unsigned int first, second, third;

    public:
      triple();
      triple( unsigned int f, unsigned int s, unsigned int t );

      bool operator()( const triple& t1, const triple& t2 ) const;
  };

  public:
    /**
     * Constructor, generates an orthogonal array
     *   and a permutation matrix used for 
     *   generating sample points.
     * @param nSamples    number of samples to generate
     * @param nInputs     number of inputs per sample
     * @param nSymbols    number of input symbols
     *                    example: a binary input x so
     *                      nSymbols of x is 2
     * @param Strength    strength of the orthogonal array
     * @param randomize   whether or not to randomize the
     *                    OA before computing the U-design
     * @param lower       lower bound of sample space
     * @param upper       upper bound of sample space
     */
    DDaceOALHSampler( int nSamples, int nInputs,
                      int Strength, bool randomize,
                      double lower, double upper );

    DDaceOALHSampler( int nSamples, int nInputs,
                      int Strength, bool randomize,
                      std::vector<Distribution>& dist);

    virtual ~DDaceOALHSampler() {;}

    virtual std::vector<DDaceSamplePoint>& getSamples( std::vector<DDaceSamplePoint>& samplePoints ) const;
    virtual std::vector<std::vector<int> > getP() const {return A_;}

    virtual DDaceSamplerBase* clone() const ;
    virtual void print(std::ostream& os) const;
    virtual const std::string& typeName() const { return typeName_; }
       
    virtual int    nSymbols() const { return nSymbols_; }
    virtual double lowerBound() const { return lower_; }
    virtual double upperBound() const { return upper_; }

    virtual int getParameter( const std::string& parameterName ) const;

    std::vector<std::vector<int> > getDesign() const { return U_; }
    std::vector<std::vector<int> > getOA() const { return A_; }

  protected:
    void initPattern();
    void createPMatrix();
    void createUDesign();
    void randomizeOA();

    std::vector<std::vector<int> >        A_;      // orthogonal array
    std::vector<std::vector<int> >        P_;      // permutation matrix
    std::vector<std::vector<int> >        U_;      // U design array

    int nSymbols_;
    int Strength_;
    int lambda_;

    bool randomize_;

    double lower_;
    double upper_;

    static std::string typeName_;

};

#endif
    
