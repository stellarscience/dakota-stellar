#include "Random.h"

// lapack's uniform random number generator.
#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
  #define DLARAN_F77 F77_FUNC(dlaran,DLARAN)
#else
  #define DLARAN_F77 DDACE_FORTRAN_GLOBAL(dlaran,DLARAN)
#endif
#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif
void DLARAN_F77(int&);

int Random::globSeed[4] = {0, 0, 0, 0};
bool Random::globSeedReady = false;

int Random::timeSeed()
{
  struct timeval tp;
  struct timezone tzp;
  
#if defined(_WIN32) || defined(_WIN64)
  union {
    __int64 ns100; /*time since 1 Jan 1601 in 100ns units */
    FILETIME ft;
  } now;
  
  GetSystemTimeAsFileTime (&now.ft);
  tp.tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
  tp.tv_sec = (long) ((now.ns100 - 116444736000000000LL) / 10000000LL);
#else
  gettimeofday(&tp, &tzp);
#endif
  
  // use number of microseconds as seed
  int t = tp.tv_usec;

  cerr << "time seed is " << t << endl;

  return t;
}

void Random::initRandom(int seed)
{
	int r = 4096;
	int r2 = r*r;
	globSeed[0] = 0;
	globSeed[1] = seed/r2;
	globSeed[2] = (seed - globSeed[1]*r2)/r;
	globSeed[3] = (seed - globSeed[1]*r2 - globSeed[2]*r);
	if ((globSeed[3] % 2) == 0)
		{
			globSeed[3]++;
		}
	globSeedReady = true;
}

double Random::uniformDeviate(double xmin, double xmax)
{
	double u;
	
	if (!globSeedReady) initRandom(timeSeed());

	u = DLARAN_F77(globSeed);
	return (xmin + (xmax-xmin)*u);
}	

Array<int> Random::randomIVector(int leng)
{
	Array<int> rtn(leng);
  int    k, iran1, iran2, itmp;

  for (k=0; k<leng; k++) rtn[k] = k;
  for (k=0; k<3*leng; k++) {
    iran1 = (int) (leng * uniformDeviate(0,1));
    iran2 = (int) (leng * uniformDeviate(0,1));
    iran1 = (iran1 == leng) ? 0 : iran1;
    iran2 = (iran2 == leng) ? 0 : iran2;
    itmp  = rtn[iran2];
    rtn[iran2] = rtn[iran1];
    rtn[iran1] = itmp;
  }
	return rtn;
}
