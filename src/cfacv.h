#ifndef _CFACV_H_
#define _CFACV_H_

#define CFACVC_VERSION 0.10

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "measurements.h"
#include "centers.h"
#include "wrapcoords.h"

typedef struct FORCEACCUMSTRUCT {
  double f; // accumulator
  int n; // number of entries accumulated
} forceAccumStruct;

enum {BOND, ANGLE, DIHED, CARTESIAN_X, CARTESIAN_Y, CARTESIAN_Z, NULL_CV};

typedef struct CVSTRUCT {
  int typ;
  int nC;
  int * ind; // indices of atomblocks that contribute to this CV
  double val;
  double ** gr;
} cvStruct;

typedef struct MTMAPPAIR { 
  int a;  // global index of shared atom
  int i; // index in first CV's ind[] array
  int j; // index in second CV's ind[] array
} mtmappair;

typedef struct METRICTENSORSTRUCT {
  double * MM;  // tensor is stored as a 1-D array; M[i][j] same as M[i*m+j]
  int m;           // dimension of matrix (m x m) == dimension of CV space
  int n;           // number of accumulated updates to M
  int nn;          // number of non-zero elements
  int * cva, * cvb;    // [i] collective variables particpating in the i'th element of M
  /* to access a non-zero element of M, we use the single index i [0..nn-1] as
     M[cva[i]*m+cvb[i]] */
  int * nca;       // number of common atoms participating in the i'th element of M
  mtmappair ** ca;       // array of common atom indices
  /*
    i = [0,nn-1]          // for each non-zero element
      c = [0,nca[i]]     //    -> for each common atom 
        at = ca[i][c]   //        -> id common atom and update M tally
        M[cva[i]*m+cvb[i]] += 1.0/mass(at) * grad(cva[i],at)*grad(cvb[i],at);
   */
} metricTensorStruct;

typedef struct TAMDOPTSTRUCT {
  double kT;
  double gamma;
  double dt;
  double ginv;
  double noise;
  int periodic;
  double half_domain;
} tamdOptStruct;


typedef double (*tmdUpdateFunc) (double, double, int );

typedef struct TMDOPTSTRUCT {
  int t0;
  int t1;
  double invinterval;
  double initval;
  double increment;
  double target;
  tmdUpdateFunc update;
} tmdOptStruct;


enum {HARMONIC, PERIODIC, VORONOICENTER, SPHERECENTER, NULL_RF};
typedef double (*restrForceFunc) ( double, double, double, double);
typedef double (*restrEnergyFunc) ( double, double, double, double);

typedef struct RESTRSTRUCT {
  double k;
  double z;
  double val;
  double min;
  double max;
  double half_domain;
  double f;
  double u;
  int nCV;
  double * cvc;
  /* spherical boundary for milestoning */
  double * sbc; // center
  double sbr;   // radius

  tamdOptStruct * tamdOpt;
  double tamd_noise;
  double tamd_restraint;
  tmdOptStruct * tmdOpt;
  int rfityp;
  restrForceFunc forceFunc;
  restrEnergyFunc energyFunc;

  forceAccumStruct * facc;

} restrStruct;


void cfacvBanner ( void );

cvStruct * New_cvStruct ( int typ, int nC, int * ind );

tamdOptStruct * New_tamdOptStruct ( double g, double kt, double dt, int riftyp, double half_domain );

tmdOptStruct * New_tmdOptStruct ( double target, int t0, int t1, int periodic );

restrStruct * New_restrStruct ( double k, double z, int nCV, double * cvc, char * rftypstr, double zmin, double zmax  );

forceAccumStruct * New_forceAccumStruct ( void );

metricTensorStruct * New_metricTensorStruct ( int M ) ;

typedef struct DATASPACESTRUCT {
  int N;                    // number of centers
  int M;                    // number of CVs
  int K;                    // number of restraints
  int iN;
  int iM;
  int iK;
  unsigned short * Xi;      // random number generator seed
  double ** R;              // [N][3] Cartesian positions of centers
  atomCenterStruct ** ac;   // [N] array of atom-group centers
  cvStruct ** cv;           // [M] array of CV's (location of system in CV space)
  metricTensorStruct * mt;  // metric tensor structure; implements MxM matrix
  restrStruct ** restr;     // [K] array of restraints applied to CV's or linear combinations of CV's
  double * z;               // shortcut copy of restraint target values copied to/from restr[]->z
  double * oldz;            // shortcut copy of previous step's target values; copied from *->z before 
                            // *->z overwrittien from restr[]->z
  double * f;               // shortcut copy of restraint force vaules copied from restr[]->f
  double * MM;              // shortcut copy of metric tensor copied from mt->MM
} DataSpace;

typedef struct SMDATASPACESTRUCT {
  int outputlevel;
  int nz; // dimensionality of CV space
  int ni; // number of images on string
  int evolve_ends; // flag indicating whether ends evolve
  double ** MM; // array of metric tensors over all images
  double ** z; // array of z-positions over all images
  double ** oldz; // saving previous step's location
  double ** g; // dF/dz for over all images
  double ** zn; // after reparameterization
  int * ztyp; // type of CV
  double * L; // running arclength
  double * s; // equally-spaced running arclength
  double reparam_tol;
  int reparam_maxiter;
  double nu; // factor for climbing string method
} SMDataSpace;

FILE * my_fopen ( char * name, char * code ) ;

DataSpace * NewDataSpace ( int N, int M, int K, long int seed );

SMDataSpace* New_stringMethod_Dataspace ( int ni, int nz, int outputlevel, double nu, int evolve_ends );

int DataSpace_getN ( DataSpace * ds );
double * SMDataSpace_image_z ( SMDataSpace * sm, int i );
double * SMDataSpace_image_oldz ( SMDataSpace * sm, int i );
double * SMDataSpace_image_g ( SMDataSpace * sm, int i );
double * SMDataSpace_image_M ( SMDataSpace * sm, int i );
int SMDataSpace_climb ( SMDataSpace * sm );
int SMDataSpace_MoveString ( SMDataSpace * sm, double stepsize );
int SMDataSpace_reparameterize ( SMDataSpace * sm );
int SMDataSpace_set_reparam_tol ( SMDataSpace * sm, double reparam_tol, int maxiter );
double * DataSpace_centerPos ( DataSpace * ds, int i );
double * DataSpace_z ( DataSpace * ds );
double * DataSpace_g ( DataSpace * ds );
double * DataSpace_MM ( DataSpace * ds );
int DataSpace_nz ( DataSpace * ds );
int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m );
int DataSpace_AddCenter_MassOnly ( DataSpace * ds, double m );
int DataSpace_AddCV ( DataSpace * ds, char * typ, int nind, int * ind ) ;
int DataSpace_AddRestr  ( DataSpace * ds, double k, double targ, int nCV, double * cvc, char * rftypstr, double zmin, double zmax );
//int DataSpace_updateRestraintValue ( DataSpace * ds, int i, double rv );
int DataSpace_AddSphericalBoundaryToRestr ( DataSpace * ds, int ir, double * ctr, double rad );

int DataSpace_AddForceAccum ( DataSpace * ds, int ir );
int DataSpace_AddTamdOpt ( DataSpace * ds, int ir, double g, double kt, double dt );
int DataSpace_AddTmdOpt  ( DataSpace * ds, int ir, double target, int t0, int t1 );
int DataSpace_ComputeCVs ( DataSpace * ds );
int DataSpace_AssignZsFromCVs ( DataSpace * ds );
int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep );
int DataSpace_VoronoiCellCheck ( DataSpace * ds, int home_cell );
double DataSpace_RestraintEnergy ( DataSpace * ds );
void DataSpace_ReportAll ( DataSpace * ds );
void DataSpace_WriteCVtoArray ( DataSpace * ds, int * active, double * res );
void DataSpace_ReportCVs ( DataSpace * ds, int step );
void DataSpace_ReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );
void DataSpace_Report_Z_and_F ( DataSpace * ds, int step, FILE * fp );
void DataSpace_metricTensor_fprintf ( DataSpace * ds, FILE * fp );
int DataSpace_metricTensor_Reset ( DataSpace * ds );
int DataSpace_forceAccumulators_Reset ( DataSpace * ds );

int DataSpace_SetRestraints ( DataSpace * ds, double * rval );
int DataSpace_metricTensor_Setup ( DataSpace * ds );
int DataSpace_Tally ( DataSpace * ds );

int DataSpace_checkdata ( DataSpace * ds );
int DataSpace_dump ( DataSpace * ds ); 
FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds );
void DataSpace_BinaryReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp );
int bin_sort ( int * bin, double * x, double * y, double * z, int nAtoms, int nCenters, int nCycles, unsigned int Seed );
#endif
