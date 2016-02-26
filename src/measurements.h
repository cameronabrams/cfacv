#ifndef _MEASUREMENTS_H_
#define _MEASUREMENTS_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void pt ( double * angles, double * hx, double * hy, double * hz, double * ox, double * oy, double * oz, int n);

#ifndef M_SQRT3
#define M_SQRT3 1.73205080756888
#endif
double my_whitenoise ( unsigned short xsubi[3]  );

double my_getbond  ( double p0[3], double p1[3], double g0[3], double g1[3] );
double my_getangle ( double p0[3], double p1[3], double p2[3], double g0[3], double g1[3], double g2[3] );
double my_getdihed ( double p1[3], double p2[3], double p3[3], double p4[3],
		     double g1[3], double g2[3], double g3[3], double g4[3] );
#endif
