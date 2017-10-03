#include "DDaceOALHSampler.h"
#include "DDaceSampler.h"
#include "UniformDistribution.h"
#include <math.h>
#include <string.h>
#include <ostream>
#include <iostream>
#include <map>

using namespace std;

string DDaceOALHSampler::typeName_ = "DDaceOALHSampler";

extern "C" {
  int bose_link( int n, int ninputs, int str, int ***AA );
  int OA_strength( int q,int nrow,int ncol,int** A,int *str,int verbose );
}


DDaceOALHSampler::triple::triple()
{
  first = second = third = 0;
}

DDaceOALHSampler::triple::triple( unsigned int f, unsigned int s, unsigned int t )
{
  first  = f;
  second = s;
  third  = t;
}

bool DDaceOALHSampler::triple::operator()( const triple& t1, const triple& t2 ) const
{
  bool ret = false;

  if( t1.first < t2.first ) { ret = true; }
  else if( t1.first == t2.first ) {
    if( t1.second < t2.second ) { ret = true; }
    else if( t1.second == t2.second ) {
      if( t1.third < t2.third ) { ret = true; }
    }
  }
  
  return ret;
}

DDaceOALHSampler::DDaceOALHSampler( int nSamples, int nInputs,
                                    int Strength, bool randomize,
                                    double lower, double upper )
:     DDaceSamplerBase(nSamples, nInputs, false),
      A_(0), P_(0), U_(0),
      Strength_(Strength), randomize_(randomize), lower_(lower), upper_(upper) 
{
	// The following code is adapted from Charles Tong's original
	// implementation

	// compute number of symbols
	double tmp = (double) nSamples_;
	tmp = pow(tmp, 0.5000001);
	nSymbols_ = (int) tmp;
	int nSamples1 = nSymbols_ * nSymbols_;
	if( nSamples1 < nSamples_ ) {
	   int nSamples2 = (nSymbols_+1)*(nSymbols_+1);
	   if( (nSamples_ - nSamples1) < (nSamples2 - nSamples_)) {
	      nSamples_ = nSamples1;
	   }
	   else {
	      nSamples_ = nSamples2;
	      nSymbols_++;
	   }
	}

   
   lambda_ = (int)(nSamples_*(1/(pow((double)nSymbols_,(double)Strength_))));

   initPattern();
}

DDaceOALHSampler::DDaceOALHSampler( int nSamples, int nInputs,
                                    int Strength, bool randomize,
                                    std::vector<Distribution>& dist )
:     DDaceSamplerBase(nSamples, nInputs, false, dist),
      A_(0), P_(0), U_(0),
      Strength_(Strength), randomize_(randomize), lower_(0), upper_(0) 
{
	// The following code is adapted from Charles Tong's original
	// implementation

	// compute number of symbols
	double tmp = (double) nSamples_;
	tmp = pow(tmp, 0.5000001);
	nSymbols_ = (int) tmp;
	int nSamples1 = nSymbols_ * nSymbols_;
	if( nSamples1 < nSamples_ ) {
	   int nSamples2 = (nSymbols_+1)*(nSymbols_+1);
	   if( (nSamples_ - nSamples1) < (nSamples2 - nSamples_)) {
	      nSamples_ = nSamples1;
	   }
	   else {
	      nSamples_ = nSamples2;
	      nSymbols_++;
	   }
	}

   
   lambda_ = (int)(nSamples_*(1/(pow((double)nSymbols_,(double)Strength_))));

   initPattern();
}

void DDaceOALHSampler::initPattern()
{

    //--------------------------------
    //  Generate orthogonal array
    //--------------------------------
    int i,j ;
    // we need to use C arrays here to connect to the C code. 
    int** pTmp;

    int status = bose_link( nSamples_, nInputs_, Strength_, &pTmp );
    if( pTmp == 0 ) throw bad_alloc();

    if( status >= 0 ) {
      if( status != nSamples_ ) {
	cerr << "DDaceOASampler: number samples adjusted to " << status << endl;
        nSamples_ = status;
      }
    }
    else {
      throw runtime_error( "DDaceOALHSampler::initPattern: bose cannot generate points" );
    }
	
	
    // randomize the pattern
    vector<int> ivec( nSymbols_ );
    for( i = 0; i < nInputs_; i++) {
      ivec = randomIVector( nSymbols_ );
      for( j = 0; j < nSamples_; j++) {
         pTmp[j][i] = ivec[pTmp[j][i]];
      }
    }
	
    int strength;
    OA_strength( nSymbols_, nSamples_, nInputs_, pTmp, &strength, 0 );
	
    if( strength < Strength_ ) {
      throw runtime_error( "DDaceOALHSampler::initPattern: failed strength test" );
    }
	
    // copy into safe array and throw out C array
    A_.resize( nSamples_ );
    for( i = 0; i < nSamples_; i++ ) {
       A_[i].resize(nInputs_);
       for( j = 0; j < nInputs_; j++ ) {
          A_[i][j] = pTmp[i][j];
       }
       // the compiler dislikes when you 'delete' something that
       // was allocated with 'malloc'
       free( pTmp[i] );    
    }
    free( pTmp );

    // randomize OA?
    if( randomize_ ) {
       randomizeOA();
    }

    // create the permutation matrix
    // then the U-design
    createPMatrix();
    createUDesign();
}

void DDaceOALHSampler::createPMatrix()
{
   // Basically this creates an array of integers
   // which enumerate from 1 to lambda*(s^r)
   // for example: OA(4,2,2,2) i.e. A is 4x2,
   //   2 symbols(s), and strength(r) 2.
   //
   //   A = [1,1]    the frequency lambda is
   //       [1,2]    calculated as 4*(2^(-2))=1
   //       [2,1]
   //       [2,2]
   //
   //   then the P matrix would contain the
   //   integers 1..(1*(2^2))=4
   //   P = [1,3]
   //       [2,4]
  

   int i;

    int lam_t = nSamples_/nSymbols_;

   // resize P to lam_t by nSymbols 
   P_.resize( lam_t );
   for( i = 0; i < lam_t; i++ ) {
      P_[i].resize( nSymbols_ );
   }

   // populate the permutation matrix
   // using the permutation k*L*s^(r-1)+1,
   // k*L*s^(r-1)+2, ..., k*L*s^(r-1)+L*s^(r-1)
   // (Note: L=lambda)
   for( int k = 0; k < nSymbols_; k++ ) {
      for( int i = 0; i < lam_t; i++ ) {
         P_[i][k] = k*lam_t + (i+1);
      }
   }
}

void DDaceOALHSampler::createUDesign()
{
    // Basically the U design comes from replacing
    // each element of A with a randomly chosen
    // element from the column of P indexed by the
    // element of A which is being replaced.
    // *confusing*
    //
    // example: OA(4,2,2,2) same as above
    //   A = [1,1]    P = [1,3]
    //       [1,2]        [2,4]
    //       [2,1]
    //       [2,2]
    //
    //   (Note: all indicies given below are not in C
    //          notation i.e., they start at 1 not 0
    //          and are presented in A[row][column] order)
    //   
    //   First, consider element A[1][1] = 1, this
    //   element is to be replaced by a randomly
    //   selected element from the first (i.e. 1) column
    //   of P, so A[1][1] will be replaced by
    //   either 1 or 2.
    //
    //   Next, consider element A[4][2] = 2, this
    //   element is to be replaced by a randomly
    //   selected element from the second (i.e. 2) column
    //   of P, so A[4][2] will be replaced by
    //   either 3 or 4.
    //   and so on ...

    // A more mathematically rigorous explanation
    // is given below.

    unsigned int i;
    unsigned int row;
    unsigned int col;
    unsigned int tmp;
    unsigned int rnd;

    unsigned int lam_t = nSamples_/nSymbols_;
   // unsigned int count  = lambda_*(int)pow((double)nSymbols_, (double)Strength_);
 
    std::map< triple, bool, triple >              used;
    std::map< triple, bool, triple >::iterator    end = used.end();

    // resize U to nSamples by nInputs
    U_.resize( nSamples_ );
    for( i = 0; i < nSamples_; i++ ) {
       U_[i].resize( nInputs_ );
    }

    // For each column replace each element k by a
    // permutation (k-1)Ls^(r-1)+1, (k-1)Ls^(r-1)+2,...,
    // (k-1)Ls^(r-1)+Ls^(r-1) = Ls^(r-1) where each
    // permutation has an equal probability of being
    // selected (Note: L = lambda)
    for( col = 0; col < nInputs_; col++ ) {
       for( row = 0; row < nSamples_; row++ ) {
          // store the value at A[row][col] to
          // use as an index
          tmp = A_[row][col];

          // find a random index into the P matrix
          // which hasn't yet been used
          do {
             rnd = (unsigned int)floor( DistributionBase::uniformUnitDeviate() * lam_t );
             
          } while( used.find( triple( tmp, col, rnd ) ) != end || rnd == lam_t );

          // add the index to the list so it
          // won't be used again
          used[ triple( tmp, col, rnd ) ] = true;
	  end = used.end();

          // copy the value from the P matrix
          // to the U design
          U_[row][col] = P_[rnd][tmp];
       }
    }
}

void DDaceOALHSampler::randomizeOA()
{
   unsigned int i,j;
   unsigned int ind1;
   unsigned int ind2;
   unsigned int itmp;
   
   vector<int> rows( nSamples_ );
   vector<int> cols( nInputs_ );
   vector<int> syms( nSymbols_ );
   vector<int> map( nSymbols_ );
   vector<int> tmp;

   //--------------------------
   //     randomize rows
   //--------------------------
   for( i = 0; i < nSamples_; i++ ) {
      rows[i] = i;
   }
  
   for( i = 0; i < nSamples_; i += 2 ) {
      // choose two indices randomly,
      // then flag the indices so they won't be
      // used again.
     cout << " i " << i << '\n';
      do {
        ind1 = (unsigned int)floor( DistributionBase::uniformUnitDeviate() * nSamples_ );
	cout << "Index 1 from DDACE OALHS " << ind1  << '\n';
      } while( ind1 == nSamples_ || rows[ind1] < 0 );
      (rows[ind1] *= -1) -= 1;
    
      // IF statement below inserted by Laura Swiler, Dec. 2009
      // It appears that this code is erroring on the case of nInputs or 
      // nSamples being odd.  This fix checks for those cases, 
      // and the resulting samples are both OAs and LHS, but 
      // this has not been extensively tested.
      if (i != nSamples_-1) {
	do {
	  ind2 = (unsigned int)floor( DistributionBase::uniformUnitDeviate() * nSamples_ );
	  cout << "Index 2 from DDACE OALHS " << ind2  << '\n';
	} while( ind2 == nSamples_ || rows[ind2] < 0 );
	(rows[ind2] *= -1) -= 1;
      }
      // swap the two rows
      tmp      = A_[ind1];
      A_[ind1] = A_[ind2];
      A_[ind2] = tmp;
   }

   // if nSamples is odd then there
   // is one row which didn't get swapped
   if( nSamples_ % 2 != 0 ) {
      // find the row indice which hasn't been selected
      for( i = 0; i < nSamples_; i++ ) {
         if( rows[i] >= 0 ) {
            ind1 = rows[i];
            break;
         }
      }

      // randomly choose another column
      do {
        ind2 = (int)floor( DistributionBase::uniformUnitDeviate() * nSamples_ );
      } while( ind2 == nSamples_ || ind2 == ind1 );

      // swap the two rows
      tmp      = A_[ind1];
      A_[ind1] = A_[ind2];
      A_[ind2] = tmp;
   }

   //--------------------------
   //     randomize columns
   //--------------------------
   for( i = 0; i < nInputs_; i++ ) {
      cols[i] = i;
   }

   for( i = 0; i < nInputs_; i += 2 ) {
      // choose two indices randomly,
      // then flag the indices so they won't be
      // used again.
      do {
        ind1 = (int)floor( DistributionBase::uniformUnitDeviate() * nInputs_ );
      } while( ind1 == nInputs_ || cols[ind1] < 0 );
      (cols[ind1] *= -1) -= 1;
    
      // IF statement below inserted by Laura Swiler, Dec. 2009
      // It appears that this code is erroring on the case of nInputs or 
      // nSamples being odd.  This fix checks for those cases, 
      // and the resulting samples are both OAs and LHS, but 
      // this has not been extensively tested.
      if (i != nInputs_-1) {
	do {
	  ind2 = (int)floor( DistributionBase::uniformUnitDeviate() * nInputs_ );
	} while( ind2 == nInputs_ || cols[ind2] < 0 );
	(cols[ind2] *= -1) -= 1;
      }

      // swap the two columns
      for( j = 0; j < nSamples_; j++ ) {
         itmp         =  A_[j][ind1];
         A_[j][ind1]  =  A_[j][ind2];
         A_[j][ind2]  =  itmp;
      }
   }

   // if nSamples is odd then there
   // is one column which didn't get swapped
   if( nInputs_ % 2 != 0 ) {
      // find the row indice which hasn't been selected
      for( i = 0; i < nInputs_; i++ ) {
         if( cols[i] >= 0 ) {
            ind1 = cols[i];
            break;
         }
      }

      // randomly choose another column
      do {
        ind2 = (int)floor( DistributionBase::uniformUnitDeviate() * nInputs_ );
      } while( ind2 == nInputs_ || ind2 == ind1 );

      // swap the two coumns
      for( j = 0; j < nSamples_; j++ ) {
         itmp         =  A_[j][ind1];
         A_[j][ind1]  =  A_[j][ind2];
         A_[j][ind2]  =  itmp;
      }
   }

   //---------------------------
   //    randomize symbols
   //---------------------------
   for( i = 0; i < nSymbols_; i++ ) {
      syms[i] = i;
   }
   
   // create random mapping for symbols
   for( i = 0; i < nSymbols_; i++ ) {
      // randomly select a symbol then
      // set the flag so this symbols isn't
      // used again
      do {
       itmp = (unsigned int)floor( DistributionBase::uniformUnitDeviate() * nSymbols_ );
      } while( itmp == nSymbols_ || syms[itmp] < 0 );

      map[i] = syms[itmp];
      (syms[itmp] *= -1 ) -= 1;
   }
 
   // change each symbol according to the random map
   for( i = 0; i < nSamples_; i++ ) {
      for( j = 0; j < nInputs_; j++ ) {
         A_[i][j] = map[ A_[i][j] ];
      }
   }
}
      

vector<DDaceSamplePoint>& DDaceOALHSampler::getSamples( vector<DDaceSamplePoint>& samplePoints ) const
{
    unsigned int i;
    unsigned int row;
    unsigned int col;
    vector<double> s( nSamples_ );

    std::vector<double> lower_bounds(nInputs_);
    std::vector<double> upper_bounds(nInputs_);
    for(i=0;i<nInputs_;i++){
     lower_bounds[i] = dist_[i].lowerBound();
     upper_bounds[i] = dist_[i].upperBound();
   }

    // s[i][j] ~ U(0,1)
    for( i = 0; i < nSamples_; i++ ) {
       s[i] = i/((double)nSamples_) + 
              DistributionBase::uniformUnitDeviate()/((double)nSamples_);
    }

    // X[i] = (x[ u[i][1] ], x[ u[i][2] ], ..., x[ u[i][m] ])
    // and calculate the sample points
    vector<vector<double> > X_( nSamples_ );
    vector<vector<double> > S_( nSamples_ );
    for( row = 0; row < nSamples_; row++ ) {
       X_[row].resize( nInputs_ );
       S_[row].resize( nInputs_ );
       for( col = 0; col < nInputs_; col++ ) { 
          X_[row][col] = s[ U_[row][col]-1 ];
          //S_[row][col] = lower_ + X_[row][col] * (upper_-lower_);
          S_[row][col] = lower_bounds[col] + 
                      X_[row][col] * (upper_bounds[col]-lower_bounds[col]);
       }
    }

    // copy S to samplePoints
    samplePoints.resize( nSamples_ );
    for( i = 0; i < nSamples_; i++ ) {
       samplePoints[i] = DDaceSamplePoint( 0, S_[i] );
    }

    return samplePoints;
}

DDaceSamplerBase* DDaceOALHSampler::clone() const
{
    DDaceSamplerBase* rtn = new DDaceOALHSampler( *this );
    if( rtn == NULL ) {
       throw runtime_error( "DDaceOALHSampler::clone()" );
    }
    return rtn;
}

void DDaceOALHSampler::print( ostream& os ) const
{
    os << "<OrthogonalArrayLatinHypercube "
       << "samples=\""   << nSamples_ << "\" "
       << "inputs=\""    << nInputs_  << "\" "
       << "symbols=\""   << nSymbols_ << "\" "
       << "strength=\""  << Strength_ << "\" "
       << "frequency=\"" << lambda_   << "\" "
       << "randomize=\"" << (randomize_?"true":"false") << "\" "
       << "seed=\"" << DistributionBase::seed() << "\"/>";
}

int DDaceOALHSampler::getParameter(const string& parameterName) const
{
    string tmp(parameterName);
    std::transform(tmp.begin(),tmp.end(),tmp.begin(),(int(*) (int)) toupper);

    if( tmp == "SAMPLES" ) {
      return nSamples_;
    }
    else if( tmp == "INPUTS" ) {
      return nInputs_;
    }
    else if( tmp == "SYMBOLS" ) {
      return nSymbols_;
    }
    else if( tmp == "STRENGTH" ) {
      return Strength_;
    }
    else if( tmp == "FREQUENCY" ) {
      return lambda_;
    }
    else if( tmp == "RANDOMIZED" ) {
      if( randomize_ ) {
        return 1;
      }
      return 0;
    }
    
    throw runtime_error( "DDaceOALHSampler::getParameter(): Unknown parameter name." );
}

    


