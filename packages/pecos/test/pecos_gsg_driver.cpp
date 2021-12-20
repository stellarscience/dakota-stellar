/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2008, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

/** \file pecos_gsg_driver.cpp
    \brief A driver program for PECOS */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "IncrementalSparseGridDriver.hpp"
#include "SharedProjectOrthogPolyApproxData.hpp"
#include "ProjectOrthogPolyApproximation.hpp"
#include "TensorProductDriver.hpp"
#include "CubatureDriver.hpp"
#include "pecos_data_types.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

#include "TestFunctions.hpp"

using namespace std;

#define POLYTYPE           LEGENDRE_ORTHOG
#define QUADRULE           GAUSS_PATTERSON
#define MAX_CHARS_PER_LINE 1000
#define BTYPE              GLOBAL_PROJECTION_ORTHOGONAL_POLYNOMIAL
#define NUMVARS            2
#define STARTLEV           1
#define NITER              0
#define MAXORD             5
#define NQOI               1
#define VERBOSE            1
#define FCNTYPE            "gerstner-iso1"
#define VARTHRLD           1.e-2

#ifdef GSGREST
inline bool fileExist (const char *fname) {
  std::ifstream in(fname);
  return in.good();
}

void saveData(const char *fbase, const RealMatrix &evalGrid, const IntVector &evalFlag, const RealMatrix &evalData) ;
void loadData(const char *fbase, RealMatrix &evalGrid, IntVector &evalFlag, RealMatrix &evalData) ;
RealMatrix getFeval(RealMatrix &evalGrid, IntVector  &evalFlag,
                    RealMatrix &evalData, const RealMatrix &var_sets,
	            bool &foundFlag);
#else
RealMatrix feval(const RealMatrix &dataMat, const int nQoI, String ftype);
#endif

int usage(){
  printf("usage: pecos_gsg_driver [options]\n");
  printf(" -h         : print out this help message \n");
  printf(" -d <nvar>  : dimensionality of parameter space (default=%d) \n",NUMVARS);
  printf(" -e <veps>  : tolerance for incremental variance threshold  (default=%lg) \n",VARTHRLD);
  printf(" -i <niter> : maximum no. of iterations (default=%d) \n",NITER);
  printf(" -l <stlev> : starting quadrature level (default=%d) \n",STARTLEV);
  printf(" -m <mord>  : not used - maximum order for shared-poly-data (default=%d) \n",MAXORD);
  printf(" -n <nqoi>  : no. of outputs (default=%d) \n",NQOI);
  printf(" -p <ptype> : type of basis (default=LEGENDRE) \n");
  printf(" -t <ftype> : test function (default=%s) \n",FCNTYPE);
  printf(" -v <verb>  : verbosity level (default=%d) \n",VERBOSE);
  exit(0);
  return (0);
}

/// A driver program for PECOS.
int main(int argc, char* argv[])
{

  using namespace Pecos;

  // Define defaults
  int            polyType = POLYTYPE ;
  int            quadRule = QUADRULE ;
  int            nQoI     = NQOI     ;  /* no. of Quantities of Interest (model outputs) */
  size_t         nvar     = NUMVARS  ;  /* dimensionality of input space */
  unsigned short mOrd     = MAXORD   ;  /* maximum order */
  unsigned short nIter    = NITER    ;  /* maximum no. of iterations */
  unsigned short strtlev  = STARTLEV ;  /* starting quadrature level */
  unsigned short verb     = VERBOSE  ;  /* verbosity  */
  double         varEps   = VARTHRLD ;
  String         ftype    = String(FCNTYPE);
  bool           extFunc  = false ;
  short btype = (short) BTYPE;

  String pstring, qstring;
  // Command-line arguments: read user input
  int c; 
  while ((c=getopt(argc,(char **)argv,"hfd:e:i:l:m:n:p:t:v:"))!=-1){
     switch (c) {
     case 'h':
       usage();
       break;
     case 'f':
       extFunc = true;
       break;
     case 'd':
       nvar    = strtol(optarg, (char **)NULL,0);  
       break;
     case 'e':
       varEps  = strtod(optarg, (char **)NULL);
       break;
     case 'i':
       nIter   = strtol(optarg, (char **)NULL,0);
       break;
     case 'l':
       strtlev = strtol(optarg, (char **)NULL,0);
       break;
     case 'm':
       mOrd    = strtol(optarg, (char **)NULL,0);  
       break;
     case 'n':
       nQoI    = strtol(optarg, (char **)NULL,0);  
       break;
     case 'p':
       pstring = String(optarg);
       break;
     case 't':
       ftype   = String(optarg);
       break;
     case 'v':
       verb    = strtol(optarg, (char **)NULL,0);
       break;
    default :
      break;
     }
  }

  // Retrieve command-line setup
  if (pstring == String("LEGENDRE")) {
    polyType = LEGENDRE_ORTHOG ;
    quadRule = GAUSS_PATTERSON ;
    //quadRule = GAUSS_LEGENDRE ;
  }
  else if (pstring == String("HERMITE")) {
    polyType = HERMITE_ORTHOG ;
    quadRule = GENZ_KEISTER   ;
  }

  if (verb > 2)
    PCout << "Instantiating IncrementalSparseGridDriver:\n";

  RealVector dimension_pref;        // empty -> isotropic
  short growth_rate = UNRESTRICTED_GROWTH;
  short refine_cntl = DIMENSION_ADAPTIVE_CONTROL_GENERALIZED;

  // Start
  // Can either use IntegrationDriver(driver_type) and then assign data or
  // use IntegrationDriver() and then assign letter.
  IntegrationDriver int_driver; // empty envelope
  // assign letter using assign_rep()
  std::shared_ptr<IncrementalSparseGridDriver> csg_driver =
    std::make_shared<IncrementalSparseGridDriver>
    (strtlev, dimension_pref, growth_rate, refine_cntl);
  int_driver.assign_rep(csg_driver);

  if (verb > 2)
    PCout << "Instantiating basis...\n";

  std::vector<BasisPolynomial> poly_basis(nvar); // array of envelopes
  for (int i=0; i<nvar; ++i) {
    poly_basis[i] = BasisPolynomial(polyType); // envelope operator=
    poly_basis[i].collocation_rule(quadRule);
  }
  csg_driver->initialize_grid(poly_basis);

  if (verb > 2)
    PCout << "  - done\n";

  // Instantiate Pecos Objects
  if (verb > 2)
    PCout << "Instantiating pecos objects...\n";
  ExpansionConfigOptions expcfgopt(INCREMENTAL_SPARSE_GRID, // expsolnapp
                                   DEFAULT_BASIS,        // expbassus
				   NO_COMBINE,           // exp combine type
				   NO_DISCREPANCY,       // discrepancy type
                                   SILENT_OUTPUT,        // output level
                                   true,                 // vbd flag
                                   2,                    // vbd order
                                   refine_cntl,          // refinement control
				   COVARIANCE_METRIC,    // refine metric
				   ACTIVE_EXPANSION_STATS,// refine stats type
                                   100,                  // max refine iter
                                   100,                  // max solver iter
                                   1.e-5,                // conv tol
                                   2);                   // soft conv limit
  BasisConfigOptions bcopt;
  UShortArray aord(mOrd,nvar);
  SharedBasisApproxData shared_data;                          // Envelope
  std::shared_ptr<SharedProjectOrthogPolyApproxData> shared_poly_data =
    std::make_shared<SharedProjectOrthogPolyApproxData>
    (BTYPE,aord,nvar,expcfgopt,bcopt);  // Letter
  shared_data.assign_rep(shared_poly_data); // don't increment ref count
  shared_poly_data->integration_driver_rep(csg_driver);
  shared_poly_data->polynomial_basis(poly_basis);

  // Instantiate Project poly approx
  std::vector<BasisApproximation> poly_approx(nQoI); // array of envelopes
  for ( int iQoI=0; iQoI<nQoI; iQoI++)
    poly_approx[iQoI].assign_rep
      (std::make_shared<ProjectOrthogPolyApproximation>(shared_data));

  if (verb > 2)
    PCout << "  - done\n";
 
#ifdef GSGREST
  /* Define saved data*/
  RealMatrix evalGrid;
  IntVector  evalFlag;
  RealMatrix evalData(1,nQoI);
  loadData((char *)"func",evalGrid,evalFlag,evalData);
  PCout<<evalGrid;
  PCout<<evalFlag;
  PCout<<evalData;
#endif

  // initial grid and compute reference approximation
  RealMatrix var_sets;
  csg_driver->compute_grid(var_sets);
  int numPts = var_sets.numCols();
  assert(nvar==var_sets.numRows());
  if (verb > 1) { 
    PCout<<var_sets<<endl; 
    if (verb > 2) {
      PCout << "Evaluate function on reference grid, ";
      PCout << "instantiate SurrogateData and compute coefficients ...\n"; 
    }
  }

#ifdef GSGREST
  bool foundFlag ;
  RealMatrix fev = getFeval(evalGrid,evalFlag,evalData,var_sets,foundFlag);
  if ( !foundFlag ) {
    saveData((char *)"func",evalGrid,evalFlag,evalData);
    return (0);
  } 
#else
  RealMatrix fev = feval(var_sets,nQoI,ftype);
#endif

  // Create SurrogateData instances and assign to 
  // ProjectOrthogPolyApproximation instances
  bool handle = true;
  for ( int iQoI=0; iQoI<nQoI; iQoI++) {
    SurrogateData sdi(handle);
    for( int jCol = 0; jCol < numPts; jCol++) {
      SurrogateDataVars sdv(nvar,0,0);
      SurrogateDataResp sdr(1,nvar); // no gradient or hessian
      sdv.continuous_variables(Teuchos::getCol<int,double>(Teuchos::Copy,var_sets,jCol));
      sdr.response_function(fev(jCol,iQoI));
      sdi.push_back(sdv,sdr);
    }
    poly_approx[iQoI].surrogate_data(sdi);
  } 

  shared_poly_data->allocate_data();    
  for ( int iQoI=0; iQoI<nQoI; iQoI++) {
    poly_approx[iQoI].compute_coefficients();
    poly_approx[iQoI].print_coefficients(PCout,false);
  }
  
  if (verb > 2)
    PCout << "  - done\n";

  // start refinement
  csg_driver->initialize_sets();
  UShortArraySet a;
  bool conv_within_tol = false;
  for ( size_t iter = 0; iter < nIter; iter++) {

    /* Compute base variance */
    RealVector respVariance(nQoI,0.0);  
    for ( int iQoI=0; iQoI<nQoI; iQoI++) {
      auto poly_approx_rep = std::static_pointer_cast<PolynomialApproximation>
	(poly_approx[iQoI].approx_rep());
      respVariance[iQoI] = poly_approx_rep->variance() ;
    }
    Real deltaVar = 0.0;
    std::vector<short unsigned int> asave;

    a = csg_driver->active_multi_index();
    if (verb > 1)
      PCout<<"Refine, iteration: "<<iter+1
	   <<"\n  ... starting variance:\n"<<respVariance<<'\n';
 
#ifdef GSGREST
    /* Iterate through all proposed sets and save/exit if some vals
    are not available */
    bool foundAllSets = true ;
    for (UShortArraySet::iterator it=a.begin(); it!=a.end(); ++it) {
      csg_driver->push_trial_set(*it);
      if (!shared_poly_data->push_available()) {
    	csg_driver->compute_trial_grid(var_sets);
    	fev = getFeval(evalGrid,evalFlag,evalData,var_sets,foundFlag);
        if ( !foundFlag ) foundAllSets = false ;
      } else {
	csg_driver->push_set();
      }
      csg_driver->pop_set();
    }
    if ( !foundAllSets ) {
      saveData((char *)"func",evalGrid,evalFlag,evalData);
      return (0);
    } 
#endif    

    int choose = 0;
    for (UShortArraySet::iterator it=a.begin(); it!=a.end(); ++it) {

      csg_driver->increment_smolyak_multi_index(*it);

      // Update surrogate data
      numPts = 0;
      if (shared_poly_data->push_available()) {

        if (verb > 1)
          PCout<<"Restoring existing index set:\n"<<*it<<endl;

        // Set available -> restore in csg and the rest
	csg_driver->push_set();

        size_t idxRestore = shared_poly_data->restore_index();
	shared_poly_data->pre_push_data();
	for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	  poly_approx[iQoI].push_coefficients();
          // Also restore the corresponding surrogate data
	  poly_approx[iQoI].surrogate_data().push(idxRestore,true);
	}
	shared_poly_data->post_push_data();

      }
      else {

	// New set -> compute
	// Create SurrogateData instances and assign to ProjectOrthogPolyApproximation instances
	csg_driver->compute_trial_grid(var_sets);
        numPts = var_sets.numCols();
        if (verb > 1) {
          PCout<<"Computing new index set:\n"<<*it<<endl;
          //PCout<<RealMatrix(var_sets,Teuchos::TRANS)<<endl;
	}
#ifdef GSGREST
	fev = getFeval(evalGrid,evalFlag,evalData,var_sets,foundFlag);
	if ( !foundFlag ) {
	  saveData((char *)"func",evalGrid,evalFlag,evalData);
	  return (0);
	} 
#else
        fev = feval(var_sets,nQoI,ftype);
#endif

	for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	  SurrogateData& sdi = poly_approx[iQoI].surrogate_data();
	  for( int jCol = 0; jCol < numPts; jCol++) {
  	    SurrogateDataVars sdv(nvar,0,0);
	    SurrogateDataResp sdr(1,nvar); // no gradient or hessian
	    sdv.continuous_variables(Teuchos::getCol<int,double>(Teuchos::Copy,var_sets,jCol));
	    sdr.response_function(fev(jCol,iQoI));
	    sdi.push_back(sdv,sdr);
	  } // done loop over number of points
	  sdi.pop_count(numPts);
	} // done loop over QoIs
        
	shared_poly_data->increment_data();
	for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	  poly_approx[iQoI].increment_coefficients();
	}
      }

      /* Compute (normalized) change in variance */
      RealVector respVarianceNew(nQoI,0.0);  
      for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	auto poly_approx_rep = std::static_pointer_cast<PolynomialApproximation>
	  (poly_approx[iQoI].approx_rep());
	respVarianceNew[iQoI] = poly_approx_rep->variance() ;
      }
      if (verb > 1) {
        PCout<<"  ... new variance:\n";
	for ( int iQoI=0; iQoI<nQoI; iQoI++) 
	  PCout<<respVarianceNew[iQoI]<<" ";
	PCout<<"\n";
      }
      respVarianceNew -= respVariance;
      if (verb > 1) {
        PCout<<"  ... delta variance:\n";
	for ( int iQoI=0; iQoI<nQoI; iQoI++) 
	  PCout<<respVarianceNew[iQoI]<<" ";
	PCout<<"\n";
      }
      Real normChange = respVarianceNew.normFrobenius()/csg_driver->unique_trial_points();
      if (normChange > deltaVar) {
        deltaVar = normChange;
        asave = *it;
      }

      shared_poly_data->decrement_data();
      for ( int iQoI=0; iQoI<nQoI; iQoI++) {
	poly_approx[iQoI].pop_coefficients(true);
	// Also restore the corresponding surrogate data
	poly_approx[iQoI].surrogate_data().pop(true);
      }
      csg_driver->pop_set(); // reverse order from increment

    } /* End iteration over proposed sets */

    if (verb > 1)
      PCout<<"Choosing :\n"<<asave
	   <<"\n  ... with relative variance: "<<deltaVar<<std::endl ;
    
    if ( asave.size() > 0 ) {
      csg_driver->update_sets(asave);

      //need to restore the data
      size_t idxRestore = shared_poly_data->restore_index();
      shared_poly_data->pre_push_data();
      for ( int iQoI=0; iQoI<nQoI; iQoI++) {
        poly_approx[iQoI].push_coefficients();
        poly_approx[iQoI].surrogate_data().push(idxRestore,true);
      }
      shared_poly_data->post_push_data();
      csg_driver->update_reference();
    }

    if ( deltaVar < varEps )
      { conv_within_tol = true; break; }

  } /* end iteration loop */

  // output sets; implementation above applies best refinement (no revert)
  csg_driver->finalize_sets(true, conv_within_tol, false);

  // sequence from ApproximationInterface::finalize_approximation():

  // shared pre-finalize:
  shared_poly_data->pre_finalize_data();

  // per-approximation finalize:
  for ( int iQoI=0; iQoI<nQoI; iQoI++) {
    // from Approximation::finalize() called from PecosApproximation::finalize()
    SurrogateData& sdi = poly_approx[iQoI].surrogate_data();
    size_t i, num_restore = sdi.popped_sets(); // # of popped trial sets
    for (i=0; i<num_restore; ++i)
      sdi.push(shared_poly_data->finalize_index(i),false);
    sdi.clear_popped();
    // from PecosApproximation::finalize()
    poly_approx[iQoI].finalize_coefficients();
  }

  // shared post-finalize:
  shared_poly_data->post_finalize_data();

  for ( int iQoI=0; iQoI<nQoI; iQoI++)
    poly_approx[iQoI].print_coefficients(PCout,false);

  return (0);

}

#ifdef GSGREST
/* Retrieve feval from saved data */
RealMatrix getFeval(RealMatrix &evalGrid, IntVector  &evalFlag,
                    RealMatrix &evalData, const RealMatrix &var_sets,
	            bool &foundFlag) {

  int numDim = var_sets.numRows(); // Dimensionality
  int numPts = var_sets.numCols(); // Number of pts
  int nouts  = evalData.numCols(); // Dimensionality of dependent variables

  RealMatrix fev(numPts,nouts);
  if (evalGrid.numCols() == 0) {
    /* Empty data, probably first call */
    evalGrid = var_sets ;
    evalFlag.size(numPts);
    evalData.shape(numPts,nouts);
    foundFlag = false;
    return fev ;
  }

  /* check consistency */
  assert( numDim == evalGrid.numRows());
  assert( evalGrid.numCols() == evalFlag.length()  );
  assert( evalGrid.numCols() == evalData.numRows() );

  foundFlag = true;
  /* Loop over all grid points requested */
  for (size_t i=0; i<numPts; i++) {
    bool foundPt = false;

    /* Loop over all saved grid points */
    for (size_t j=0; j<evalGrid.numCols(); j++) {

      /* compute L1 norm */
      double dsum = 0.0;
      for (size_t k = 0; k<numDim ; k++ ) 
        dsum += fabs(var_sets(k,i)-evalGrid(k,j));

      if ( dsum < 1.e-10 ) {

        if (evalFlag(j) == 1 ) {
          /* found saved point */
          foundPt = true ;
          for (size_t k = 0; k<nouts ; k++ ) 
            fev(i,k) = evalData(j,k);
	}

      }

      if (foundPt) break ;

    } /* end loop over stored points */

    if ( !foundPt ) {

      /* missing point, add it to the list */
      foundFlag = false;
      /* Increment Grid */
      RealMatrix evalGridSave = evalGrid;
      evalGrid.shapeUninitialized(evalGridSave.numRows(),evalGridSave.numCols()+1);
      for ( int j1 = 0; j1<evalGridSave.numCols(); j1++ )
        for ( int i1 = 0; i1<evalGridSave.numRows(); i1++ )
          evalGrid(i1,j1) = evalGridSave(i1,j1);
      for ( int i1 = 0; i1<evalGridSave.numRows(); i1++ )
        evalGrid(i1,evalGridSave.numCols()) = var_sets(i1,i);
      /* Increment eval flag */
      IntVector evalFlagSave = evalFlag;
      evalFlag.size(evalFlagSave.length()+1);
      for ( int i1 = 0; i1<evalFlagSave.length(); i1++ )
        evalFlag(i1) = evalFlagSave(i1);
      evalFlag(evalFlagSave.length()) = 0;
      /* Increment eval func dimensions */
      RealMatrix evalDataSave = evalData;
      evalData.shape(evalDataSave.numRows()+1,evalDataSave.numCols());
      for ( int j1 = 0; j1<evalDataSave.numCols(); j1++ )
        for ( int i1 = 0; i1<evalDataSave.numRows(); i1++ )
          evalData(i1,j1) = evalDataSave(i1,j1);
    }

  } /* end loop over stored points */

  
  return fev;

}

void getMatrixSize(const char *fname, int &nRows, int &nCols) {


  std::ifstream in(fname);
  
  if(!in){
    printf("getMatrixSize() : the requested file %s does not exist !\n",fname) ;
    exit(1) ;
  }
  
  String theLine="";

  // figure out number of lines and columns
  int ix = 0 ;
  while(in.good()) {

    getline(in,theLine);
    
    if ( theLine == "" ) break;
    if ( theLine.compare(0,1,"#") == 0 ) continue ;

    istringstream s(theLine);
    int    iy = 0 ;
    double tmp    ;
    while( s >> tmp ) iy++ ;

    if ( ( ix > 0 ) && ( iy != nCols ) )
    {
      printf("getMatrixSize() : Error at line %d !!!\n",ix+1) ;
      printf("                  no. of columns should be %d instead of %d\n",nCols,iy) ;
      exit(1) ;
    }
    
    nCols = iy ;

    ix++ ;

  }

  nRows = ix ;
  in.close();

  return ;

}

void getMatrixTrans(const char *fname, RealMatrix &dataIn) {

  int nRows = dataIn.numRows();
  int nCols = dataIn.numCols();  

  if (nRows==0 || nCols==0){
    printf("getMatrixTrans() : the requested data matrix is empty\n") ;
    exit(1) ;
  }

  std::ifstream in(fname);
  
  if(!in){
    printf("getMatrixTrans() : the requested file %s does not exist\n",fname) ;
    exit(1) ;
  }
  
  String theLine="";
  int ix = 0;

  while(in.good()){

    getline(in,theLine);

    if (theLine=="") break;
    if ( theLine.compare(0, 1, "#") == 0 ) continue ;

    istringstream s(theLine);
    int  iy = 0;
    Real tmp;
    while(s >> tmp){
      dataIn(iy,ix)=tmp;
      iy++;
    }
    if ( iy != nRows ) {
      printf("getMatrixTrans(): Error at line %d while reading %s; number of columns should be %d\n", 
             ix+1, fname,nRows); 
      exit(1);
    }
    ix++;
  }
  if ( ix != nCols ) {
    printf("getMatrixTrans(): Error while reading %s; number of rows should be %d\n",fname,nCols); 
    exit(1);
  }
  in.close();

  return;

}

void getMatrix(const char *fname, RealMatrix &dataIn) {

  int nRows = dataIn.numRows();
  int nCols = dataIn.numCols();  

  if (nRows==0 || nCols==0){
    printf("getMatrixTrans() : the requested data matrix is empty\n") ;
    exit(1) ;
  }

  std::ifstream in(fname);
  
  if(!in){
    printf("getMatrixTrans() : the requested file %s does not exist\n",fname) ;
    exit(1) ;
  }
  
  String theLine="";
  int ix = 0;

  while(in.good()){

    getline(in,theLine);

    if (theLine=="") break;
    if ( theLine.compare(0, 1, "#") == 0 ) continue ;

    istringstream s(theLine);
    int  iy = 0;
    Real tmp;
    while(s >> tmp){
      dataIn(ix,iy)=tmp;
      iy++;
    }
    if ( iy != nCols ) {
      printf("getMatrixTrans(): Error at line %d while reading %s; number of columns should be %d\n", 
             ix+1, fname,nRows); 
      exit(1);
    }
    ix++;
  }
  if ( ix != nRows ) {
    printf("getMatrixTrans(): Error while reading %s; number of rows should be %d\n",fname,nCols); 
    exit(1);
  }
  in.close();

  return;

}

void getIntVector(const char *fname, IntVector &dataIn) {

  int nvals=dataIn.length();

  if (nvals==0){
    printf("getIntVector() : the requested data array is empty\n") ;
    exit(1) ;
  }

  int nCols=1;

  std::ifstream in(fname);
  
  if(!in){
    printf("getIntVector() : the requested file %s does not exist\n",fname) ;
    exit(1) ;
  }
  
  string theLine="";
  int ix=0;
  
  while( in.good() ){

    getline(in,theLine);
    
    if (theLine=="") break;

    istringstream s(theLine);
    int iy=0;
    int tmp;
    s >> tmp;
    dataIn(ix) = tmp;
    iy++;
    if (s>>tmp) {
      printf("getIntVector(): Error at line %d while reading %s; number of columns should be %d\n", 
              ix+1, fname,nCols); 
      exit(1);
    }
    ix++;
  }
  if ( ix != nvals ) {
    printf("getIntVector(): Error while reading %s; number of rows should be %d\n",
           fname,nvals); 
    exit(1);
  }
  in.close();

  return ;

}

/* Save feval to files */
void saveData(const char *fbase, const RealMatrix &evalGrid,
               const IntVector &evalFlag, const RealMatrix &evalData) {

  char fname[100];

  std::ofstream fout;

  /* save grid */
  sprintf(fname,"%s_grid.dat",fbase);  
  fout.open(fname);
  for (size_t i=0; i < evalGrid.numCols(); i++ ) {
    for (size_t j=0; j < evalGrid.numRows(); j++ ) {
      fout.precision(16);
      fout << std::scientific << std::setw(24) << evalGrid(j,i) << " ";
    }
    fout << std::endl ;
  }
  fout.close() ;

  /* save evalFlag */
  sprintf(fname,"%s_flag.dat",fbase);  
  fout.open(fname);
  for (size_t i=0; i < evalFlag.length(); i++ ) 
    fout << evalFlag(i) << std::endl;
  fout.close() ;

  /* save output matrix */
  sprintf(fname,"%s_data.dat",fbase);  
  fout.open(fname);
  fout.precision(16);
  for (size_t j=0; j < evalData.numRows(); j++ ) {
    for (size_t i=0; i < evalData.numCols(); i++ ) {
      fout << std::scientific << std::setw(24) << evalData(j,i) << " ";
    }
    fout << std::endl ;
  }
  fout.close() ;

  return ;

}

/* Load feval from files */
void loadData(const char *fbase, RealMatrix &evalGrid,
               IntVector &evalFlag, RealMatrix &evalData) {

  char fname[100];

  /* Check if all files exist */
  sprintf(fname,"%s_grid.dat",fbase);  
  if ( !fileExist(fname) ) return ;
  sprintf(fname,"%s_flag.dat",fbase);  
  if ( !fileExist(fname) ) return ;
  sprintf(fname,"%s_data.dat",fbase);  
  if ( !fileExist(fname) ) return ;


  /* load grid */
  sprintf(fname,"%s_grid.dat",fbase);  
  int nDims, nPts ;
  getMatrixSize(fname, nPts, nDims) ;
  evalGrid.shapeUninitialized(nDims,nPts);
  getMatrixTrans(fname, evalGrid);

  /* load eval flag */
  evalFlag.sizeUninitialized(nPts);
  sprintf(fname,"%s_flag.dat",fbase);  
  getIntVector(fname, evalFlag);

  /* load function outputs */
  int nOuts ;
  sprintf(fname,"%s_data.dat",fbase);  
  getMatrixSize(fname, nPts, nOuts) ;
  evalData.shapeUninitialized(nPts,nOuts);
  assert(nPts == evalGrid.numCols());
  getMatrix(fname, evalData);

  bool allEvals = true ;
  for (size_t i=0;i<nPts;i++) 
    if (evalFlag(i) != 1) {
      PCout<<" Point "<<i<<" is not evaluated!"<<std::endl;
      allEvals = false;
    }

  if (!allEvals) {
      PCout<<"  -> Exit !"<<std::endl;    
      std::terminate() ;
  }

  return ;

}

#else

RealMatrix feval(const RealMatrix &dataMat, const int nQoI, String ftype) 
{

  assert(nQoI==1);

  int i, j ;

  int numDim = dataMat.numRows(); // Dimensionality
  int numPts = dataMat.numCols(); // Number of pts

  /* Count the number of function evaluations; */
  RealMatrix fev(numPts,nQoI);
  for (i=0; i<numPts; ++i) {
    RealVector xIn(numDim);
    for (j=0; j<numDim; ++j) xIn[j] = dataMat(j,i);
    //fev(ieval,0)=genz(String("cp1"), xIn);
    if ( ftype == String("pol2") )
      fev(i,0) = custPol(ftype, xIn);
    else if ( ftype == String("gerstner-iso1") )
      fev(i,0) = gerstner(String("iso1"), xIn);
    else if ( ftype == String("gerstner-iso2") )
      fev(i,0) = gerstner(String("iso2"), xIn);
    else if ( ftype == String("gerstner-iso3") )
      fev(i,0) = gerstner(String("iso3"), xIn);
    else if ( ftype == String("gerstner-aniso1") )
      fev(i,0) = gerstner(String("aniso1"), xIn);
    else if ( ftype == String("gerstner-aniso2") )
      fev(i,0) = gerstner(String("aniso2"), xIn);
    else if ( ftype == String("gerstner-aniso3") )
      fev(i,0) = gerstner(String("aniso3"), xIn);
    else
      throw(std::runtime_error("Pecos::feval() unknown ftype\n"));
  }
  
  return fev;

}

#endif
