/* Declarations for dynamic allocation */

int *ivector(), **imatrix();

/** 
 * Definitions for Galois Fields used by OA Sampler
 */

struct GF {
  int n,p,q;	/* hold the sizes of the dynamic
		 * memory pointed to below */
  int *xton;
  int **plus;
  int **times;
  int *inv;
  int *neg;
  int *root;
  int **poly;
};

