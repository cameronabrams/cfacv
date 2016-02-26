#ifndef _WRAPCOORDS_H_
#define _WRAPCOORDS_H_

void wrapCoords_raw ( double * x, double * y, double * z, int n, 
		  double Lx, double Ly, double Lz,
		  double ox, double oy, double oz );
void wrapCoords_water ( double * x, double * y, double * z, int n, 
		  double Lx, double Ly, double Lz,
		  double ox, double oy, double oz );

void wrapCoords ( double * x, double * y, double * z, double * m, int n,
		  int * mi, int * mn, double nMol,
		  double Lx, double Ly, double Lz,
		  double ox, double oy, double oz );

#endif
