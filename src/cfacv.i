%module cfa_cvlibc
%include carrays.i
%array_functions(double, array);
%array_functions(int, arrayint);
%inline %{
double get_double(double *a, int index) {
  return a[index];
}
%}
%{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cfacv.h"
%}
extern void cfacvBanner ( void );
extern FILE * my_fopen ( char * name, char * code );
extern DataSpace * NewDataSpace ( int N, int M, int K, long int seed ) ;
extern SMDataSpace* New_stringMethod_Dataspace ( int ni, int nz, int outputlevel, double nu, int evolve_ends, int dual );
extern int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m );
extern int DataSpace_AddCenter_MassOnly ( DataSpace * ds, double m );
extern int DataSpace_AddCV ( DataSpace * ds, char * typ, int nind, int * ind ) ;
extern int DataSpace_AddRestr ( DataSpace * ds, double k, double targ, int nCV, double * cvc, char * rftypstr, double zmin, double zmax  );
//extern int DataSpace_updateRestraintValue ( DataSpace * ds, int i, double rv ) ;
extern int DataSpace_AddSphericalBoundaryToRestr ( DataSpace * ds, int ir, double * ctr, double rad );
extern int DataSpace_AddTamdOpt ( DataSpace * ds, int ir, double g, double kt, double dt );
extern int DataSpace_AddTmdOpt  ( DataSpace * ds, int ir, double target, int t0, int t1 );
extern int DataSpace_getN ( DataSpace * ds );
extern void DataSpace_ReportAll ( DataSpace * ds );
extern void DataSpace_Report_Z_and_F ( DataSpace * ds, int step, FILE * fp );
extern void DataSpace_WriteCVtoArray ( DataSpace * ds, int * active, double * res );
extern void DataSpace_ReportCVs ( DataSpace * ds, int step );
extern double * SMDataSpace_image_z ( SMDataSpace * sm, int i );
extern double * SMDataSpace_image_oldz ( SMDataSpace * sm, int i );
extern double * SMDataSpace_image_g ( SMDataSpace * sm, int i );
extern double * SMDataSpace_image_M ( SMDataSpace * sm, int i );
extern int SMDataSpace_reparameterize ( SMDataSpace * sm );
extern int SMDataSpace_MoveString ( SMDataSpace * sm, double stepsize );
extern int SMDataSpace_climb ( SMDataSpace * sm );
extern int SMDataSpace_set_reparam_tol ( SMDataSpace * sm, double reparam_tol, int maxiter );
extern double * DataSpace_centerPos ( DataSpace * ds, int i );
extern double * DataSpace_z ( DataSpace * ds );
extern double * DataSpace_g ( DataSpace * ds );
extern double * DataSpace_MM ( DataSpace * ds );
extern int DataSpace_nz ( DataSpace * ds );
extern int DataSpace_ComputeCVs ( DataSpace * ds );
extern int DataSpace_AssignZsFromCVs ( DataSpace * ds );
extern int DataSpace_SetRestraints ( DataSpace * ds, double * rval );
extern int DataSpace_metricTensor_Setup ( DataSpace * ds );
extern int DataSpace_metricTensor_Reset ( DataSpace * ds );
extern int DataSpace_Tally ( DataSpace * ds );
extern int DataSpace_forceAccumulators_Reset ( DataSpace * ds );

extern int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep );
extern int DataSpace_VoronoiCellCheck ( DataSpace * ds, int home_cell );
extern double DataSpace_RestraintEnergy ( DataSpace * ds );
extern void DataSpace_ReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );
extern int DataSpace_SetRestraints ( DataSpace * ds, double * rval );
extern int DataSpace_checkdata ( DataSpace * ds );
extern int DataSpace_dump ( DataSpace * ds ); 
extern FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds );
extern void DataSpace_BinaryReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );
extern int bin_sort ( int * bin, double * x, double * y, double * z, int nAtoms, int nCenters, int nCycles, 
		  unsigned int Seed );
