%module libmeasurements
%include carrays.i
%array_functions(double, array);
%array_functions(int, arrayint);
%inline %{
double get_double(double *a, int index) {
	return a[index];
}
%}
%{
#include "measurements.h"
%}
extern void pt ( double * angles, double * hx, double * hy, double * hz, double * ox, double * oy, double * oz, int n);
extern double my_getbond  ( double p0[3], double p1[3], double g0[3], double g1[3] );
extern double my_getangle ( double p0[3], double p1[3], double p2[3], double g0[3], double g1[3], double g2[3] );
extern double my_getdihed ( double p1[3], double p2[3], double p3[3], double p4[3],
		     double g1[3], double g2[3], double g3[3], double g4[3] );
