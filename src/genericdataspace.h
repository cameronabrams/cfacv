#ifndef _DATASPACE_H_
#define _DATASPACE_H_
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// genericdataspace.h
typedef struct GENDATAMATRIXSTRUCT {
  int n;
  int m;
  double ** A;
} GenDataSpace;

typedef struct IGENDATAMATRIXSTRUCT {
  int n;
  int m;
  int ** I;
} IGenDataSpace;

GenDataSpace * NewGenDataSpace ( int N, int M );
void FreeGenDataSpace ( GenDataSpace * gds );
IGenDataSpace * NewIGenDataSpace ( int N, int M );
double * GenDataSpace_getAddr ( GenDataSpace * gds, int i );
int * IGenDataSpace_getAddr ( IGenDataSpace * igds, int i );
void GenDataSpace_Set ( GenDataSpace * gds, int i, int j, double val );
double GenDataSpace_Get ( GenDataSpace * gds, int i, int j );
void GenDataSpace_Scale ( GenDataSpace * gds, double x );
void GenDataSpace_Diff_RunningAverage ( GenDataSpace * gds );

void GenDataSpace_DistMap ( GenDataSpace * map, GenDataSpace * r1, GenDataSpace * r2 );
void GenDataSpace_DiffMap ( GenDataSpace * map1, GenDataSpace * map0);
void GenDataSpace_AddToMap ( GenDataSpace * map1, GenDataSpace * map0);

void GenDataSpace_Accumulate ( GenDataSpace * gds, int i, int j, double val );
void GenDataSpace_WriteToFile ( GenDataSpace * gds, char * filename );
void GenDataSpace_WriteToFile_Gnuplot ( GenDataSpace * gds, char * filename, int outerpad, char * directive );
void GenDataSpace_WriteToFile_Gnuplot_IntIndices ( GenDataSpace * gds, int * I, int * J, char * filename, int outerpad, char * directive );
void GenDataSpace_CorrCoeff (GenDataSpace * timedat, GenDataSpace * corrdat);
#endif
