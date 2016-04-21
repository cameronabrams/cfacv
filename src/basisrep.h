#ifndef _BASISREP_H_
#define _BASISREP_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct BASIS_FUNCTION_CENTER {
  int D;       // dimensionality of CV space
  double * r;  // location (dim)
  double a;    // amplitude (kcal/mol);
  double s;    // sigma (dim)
  int i;       // unique index
} bCntr;

typedef struct LOCAL_RECONSTRUCTION {
  int D;       // dimensionality of CV space
  double * r;  // location (dim)
  double * f;  // force (kcal/mol/dim)
  double e;    // energy (kcal/mol)
  int i;       // unique index
} lRecon;

char * next_word ( char * p );
int sscanf_stringList ( char * str, char ** slist, int n );
int sscanf_doubleList ( char * str, double * dlist, int n ); 
lRecon * lRecon_new ( int D );
void kernel ( double z, double sig, double * f, double * g );
int vecdiff ( double * c, double * a, double * b, int D, int * per, double * dom );
double vecdot ( double * a, double * b, int D );
double norm ( double * x, int D );
int localReconstruct ( lRecon * z, bCntr * bsc, int nC, int D, int * per, double * dom );
bCntr * readBasisSet ( int * nC, int D, char * fn );
int writeBasisSet ( bCntr * bsc, int nC, int D, char * fn );
int writeReconstruction ( lRecon * recon, int nR, int D, char * fn );

#endif
