%module genericdataspace
%include carrays.i
%array_functions(double, array);
%array_functions(int, arrayint);
%inline %{
double get_double(double *a, int index) {
	return a[index];
}
%}
%{
#include <stdlib.h>
#include "genericdataspace.h"
%}
extern GenDataSpace * NewGenDataSpace ( int N, int M );
extern void FreeGenDataSpace ( GenDataSpace * gds );
extern IGenDataSpace * NewIGenDataSpace ( int N, int M );
extern double * GenDataSpace_getAddr ( GenDataSpace * gds, int i );
extern int * IGenDataSpace_getAddr ( IGenDataSpace * igds, int i );
extern void GenDataSpace_Set ( GenDataSpace * gds, int i, int j, double val );
extern double GenDataSpace_Get ( GenDataSpace * gds, int i, int j );
extern void GenDataSpace_Scale ( GenDataSpace * gds, double x );
extern void GenDataSpace_Diff_RunningAverage ( GenDataSpace * gds );
extern void GenDataSpace_Accumulate ( GenDataSpace * gds, int i, int j, double val );
extern void GenDataSpace_WriteToFile ( GenDataSpace * gds, char * filename );
extern void GenDataSpace_WriteToFile_Gnuplot ( GenDataSpace * gds, char * filename, int outerpad, char * directive );
extern void GenDataSpace_WriteToFile_Gnuplot_IntIndices ( GenDataSpace * gds, int * I, int * J, char * filename, int outerpad, char * directive );
extern void GenDataSpace_DistMap ( GenDataSpace * map, GenDataSpace * r1, GenDataSpace * r2 );
extern void GenDataSpace_DiffMap ( GenDataSpace * map1, GenDataSpace * map0);
extern void GenDataSpace_AddToMap ( GenDataSpace * map1, GenDataSpace * map0);
extern void GenDataSpace_CorrCoeff (GenDataSpace * timedat, GenDataSpace * corrdat);
