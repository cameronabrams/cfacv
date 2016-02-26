#ifndef _CENTERS_H_
#define _CENTERS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct ATOMCENTERSTRUCT {
  int n;
  int * ind;
  double * m;
  double M;
} atomCenterStruct;


atomCenterStruct * New_atomCenterStruct ( int n );

typedef struct CENTERSTRUCT * pcenterStruct;
typedef struct CENTERSTRUCT { 
  int id;
  int maxN;
  int iN;
  int * mList;
  double rg;
  double cm[3];
  pcenterStruct left; /* only used in residue blocking */
  pcenterStruct right;
} centerStruct;

centerStruct * New_centerStruct ( int id, int maxN );
void centerStuct_addMember ( centerStruct * c, int i );

void center_rg ( centerStruct * c, double * x, double * y, double * z );

int rgyr_sort ( centerStruct * c, int * bin, double * x, double * y, double * z, 
		int nAtom, int minAtom, double * rg, unsigned int Seed  );
centerStruct * Null_centerStruct ( void );
int bin_sort ( int * bin, double * x, double * y, double * z, int nAtoms, int nCenters, int nCycles, 
	       unsigned int Seed );



#endif
