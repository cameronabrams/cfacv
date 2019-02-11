#include "measurements.h"

//#ifndef _PARANOIA_
//#define _PARANOIA_ 1
//#endif

void pt ( double * angles, double * hx, double * hy, double * hz, double * ox, double * oy, double * oz, int n) {
  int i;
  double p0[3], p1[3], p2[3], p3[3];
  double g0[3], g1[3], g2[3], g3[3];
  for (i=0;i<(n-1);i++) {
    p0[0]=hx[i]; p0[1]=hy[i]; p0[2]=hz[i];
    p1[0]=ox[i]; p1[1]=oy[i]; p1[2]=oz[i];
    p2[0]=ox[i+1]; p2[1]=oy[i+1]; p2[2]=oz[i+1];
    p3[0]=hx[i+1]; p3[1]=hy[i+1]; p3[2]=hz[i+1];
    angles[i]=180.0/M_PI*my_getdihed(p0,p1,p2,p3,g0,g1,g2,g3);
  }
}

/* my_getbond:  accepts 2 3-component vector cartesian positions and
 * computes (i) the distance between them, and (ii) the gradients of
 * this distance with respect to each coordinate.
 */
double my_getbond ( double p0[3], double p1[3], double g0[3], double g1[3] ) {

  double b[3],dr, d;

  b[0]=p0[0]-p1[0];
  b[1]=p0[1]-p1[1];
  b[2]=p0[2]-p1[2];
  
  d=sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
  dr=1.0/d;
  g0[0]=b[0]*dr;
  g0[1]=b[1]*dr;
  g0[2]=b[2]*dr;
  g1[0]=-g0[0];
  g1[1]=-g0[1];
  g1[2]=-g0[2];

  return d;

}
/* my_getrgyr:  accepts a 2-D array of {x y z} coordinates of centers, a 2-D array of {x y z}
 * gradients, an integer with the number of centers, and a 1-D array of masses of centers
 * and computes (i) the radius of gyration between centers, and (ii) the gradients of
 * this distance with respect to each center.
 * (c) 2018 Jasmine Gardner
 */
double my_getrgyr (double ** R, double ** g0, int n, double * m, int * ind) {

  double mtot=0.0, COMx=0.0, COMy=0.0, COMz=0.0, b=0.0, rg=0.0;
  int i, ind1;

  for (i=0; i<n; i++) {
    ind1=ind[i];
    COMx+=R[ind1][0]*m[ind1];
    COMy+=R[ind1][1]*m[ind1];
    COMz+=R[ind1][2]*m[ind1];
    mtot+=m[ind1];
  }
  COMx/=(mtot);
  COMy/=(mtot);
  COMz/=(mtot);
  for (i=0; i<n; i++) {
      ind1=ind[i];
      b+=m[ind1]*(pow((R[ind1][0]-COMx),2)+pow((R[ind1][1]-COMy),2)+pow((R[ind1][2]-COMz),2));
  }
  rg=sqrt(b/mtot);

  for (i=0; i<n; i++) {
    ind1=ind[i];
    g0[i][0]=(m[ind1]/(rg*mtot))*(R[ind1][0]-COMx);
    g0[i][1]=(m[ind1]/(rg*mtot))*(R[ind1][1]-COMy);
    g0[i][2]=(m[ind1]/(rg*mtot))*(R[ind1][2]-COMz);
  }
  return rg;

}
/* my_getangle: accepts three 3-component vector cartesian positions
 * and computes (i) the angle about the second point in RADIANS and
 * (ii) the gradients of that angle with respect to each cartesian
 * coordinate.  To conform to NAMD "anglegrad", the gradients 
 * are returned in units of radians/unit-length.
 */
double my_getangle ( double p0[3], double p1[3], double p2[3], double g0[3], double g1[3], double g2[3] ) {
  double b01[3], b12[3], d01, d12, dp, ct, t, mrst, rdd, ctrd1, ctrd2;
  double sd01, sd12;

  b01[0]=p0[0]-p1[0];
  b01[1]=p0[1]-p1[1];
  b01[2]=p0[2]-p1[2];
  sd01=b01[0]*b01[0]+b01[1]*b01[1]+b01[2]*b01[2];
  d01=sqrt(sd01);
  b12[0]=p1[0]-p2[0];
  b12[1]=p1[1]-p2[1];
  b12[2]=p1[2]-p2[2];
  sd12=b12[0]*b12[0]+b12[1]*b12[1]+b12[2]*b12[2];
  d12=sqrt(sd12);
  
  dp=b01[0]*b12[0]+b01[1]*b12[1]+b01[2]*b12[2];

  rdd=1.0/(d01*d12);
  ct=-dp*rdd;
  ctrd1=ct/sd01;
  ctrd2=ct/sd12;

  t=acos(ct);

  mrst=1.0/sin(t);

  /* gradients */

  g0[0]=mrst*(b12[0]*rdd+b01[0]*ctrd1);
  g0[1]=mrst*(b12[1]*rdd+b01[1]*ctrd1);
  g0[2]=mrst*(b12[2]*rdd+b01[2]*ctrd1);

  g2[0]=mrst*(-b01[0]*rdd-ctrd2*b12[0]);
  g2[1]=mrst*(-b01[1]*rdd-ctrd2*b12[1]);
  g2[2]=mrst*(-b01[2]*rdd-ctrd2*b12[2]);

  g1[0]=-g0[0]-g2[0];
  g1[1]=-g0[1]-g2[1];
  g1[2]=-g0[2]-g2[2];

  return t;
}

/* my_getdihed: accepts four 3-component vector cartesian positions
 * and computes (i) the dihedral angle about the vector connecting the
 * second and third points defined by (1)-(2)-(3)-(4) in RADIANS, and
 * (ii) the gradients in the torsion with respect to each of the four
 * positions.  Gradients are returned in the g arrays in units
 * of radians/unit-length.
 *
 * (c) 2008 cameron abrams
 */

int mycross ( double c[3], double a[3], double b[3] ) {
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
  return 0;
}

double mydot ( double a[3], double b[3] ) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

int mysum ( double c[3], double a[3], double b[3] ) {
  c[0]=a[0]+b[0];
  c[1]=a[1]+b[1];
  c[2]=a[2]+b[2];
  return 0;
}

int mydiff ( double c[3], double a[3], double b[3] ) {
  c[0]=a[0]-b[0];
  c[1]=a[1]-b[1];
  c[2]=a[2]-b[2];
  return 0;
}

double mynorm ( double a[3] ) {
  return sqrt(mydot(a,a));
}

int myscale ( double x, double a[3] ) {
  a[0]*=x; a[1]*=x; a[2]*=x;
  return 0;
}

int mycopy ( double a[3], double b[3] ) {
  a[0]=b[0]; a[1]=b[1]; a[2]=b[2];
  return 0;
}

int myddot ( double A[][3], double a[3] ) {
  A[0][0]=a[1]*a[1]+a[2]*a[2];
  A[1][1]=a[0]*a[0]+a[2]*a[2];
  A[2][2]=a[0]*a[0]+a[1]*a[1];
  A[0][1]=A[1][0]=-a[0]*a[1];
  A[0][2]=A[2][0]=-a[0]*a[2];
  A[1][2]=A[2][1]=-a[1]*a[2];
  return 0;
}

int mytdot ( double A[][3], double a[3], double b[3] ) {
  A[0][0]=-a[1]*b[1]-a[2]*b[2];
  A[1][1]=-a[2]*b[2]-a[0]*b[0];
  A[2][2]=-a[0]*b[0]-a[1]*b[1];
  A[0][1]=2*a[1]*b[0]-a[0]*b[1];
  A[1][0]=2*a[0]*b[1]-a[1]*b[0];
  A[0][2]=2*a[2]*b[0]-a[0]*b[2];
  A[2][0]=2*a[0]*b[2]-a[2]*b[0];
  A[1][2]=2*a[2]*b[1]-a[1]*b[2];
  A[2][1]=2*a[1]*b[2]-a[2]*b[1];
  return 0;
}

int mymatvec ( double z[3], double A[][3], double a[3] ) {
  int i, j;
  for (i=0;i<3;i++) {
    z[i]=0.0;
    for (j=0;j<3;j++) {
      z[i]+=A[i][j]*a[j];
    }
  }
  return 0;
}

#ifndef SIN_THRESH
#define SIN_THRESH 0.1
#endif

double my_getdihed ( double p1[3], double p2[3], double p3[3], double p4[3],
		     double g1[3], double g2[3], double g3[3], double g4[3] ) {
  double R12[3], R23[3], R34[3], r12, r23, r34;
  double C12[3], C23[3], C123[3], c12, c23, c123;
  double cos_phi, sin_phi, phi;

  mydiff(R12,p2,p1);     r12=mynorm(R12);
  mydiff(R23,p3,p2);     r23=mynorm(R23);
  mydiff(R34,p4,p3);     r34=mynorm(R34);
  mycross(C12,R12,R23);  c12=mynorm(C12);
  mycross(C23,R23,R34);  c23=mynorm(C23);
  mycross(C123,R23,C12); c123=mynorm(C123);

  cos_phi=mydot(C12,C23)/(c12*c23);
  sin_phi=mydot(C23,C123)/(c23*c123);

  if (fabs(cos_phi+1.0)<1.e-6) phi=M_PI;
  else if (fabs(cos_phi-1.0)<1.e-6) phi=0.0;
  else phi=acos(cos_phi);
  if (sin_phi<0) phi*=-1;
  
#ifdef _PARANOIA_
  if (_PARANOIA_) {
    if (phi!=phi) {
      fprintf(stderr,"CFACV/C/PARANOIA) Tripped in my_getdihed: phi %.5lf\n",phi);
      fprintf(stderr,"CFACV/C/PARANOIA) cos_phi %.5lf\n",cos_phi);
      fprintf(stderr,"Program exits.\n");
      exit(-1);
    }
  }
#endif

  //  printf("COS) cos_phi % 10.4lf phi % 10.4lf = % 10.4lf deg\n",cos_phi,phi,180/M_PI*phi);

  myscale(1.0/c23,C23);
  
  if (fabs(sin_phi)>SIN_THRESH) {
    double rst = -1.0/sin_phi;
    double dc1[3], dc2[3], dc3[3], tmp[3];
    double g12[3],g23[3];
    
    myscale(1.0/c12,C12);

    mycopy(g12,C12);
    myscale(cos_phi,g12);
    mydiff(g12,g12,C23);
    myscale(1.0/c12,g12);
    mycopy(g23,C23);
    myscale(cos_phi,g23);
    mydiff(g23,g23,C12);
    myscale(1.0/c23,g23);

    mycross(dc1,g12,R23);
    mycross(dc3,R23,g23);
    mycross(tmp,R12,g12);
    mycross(dc2,g23,R34);
    mysum(dc2,tmp,dc2);

    myscale(rst,dc1);
    myscale(rst,dc2);
    myscale(rst,dc3);

    g1[0]=-dc1[0];           g1[1]=-dc1[1];           g1[2]=-dc1[2];
    g2[0]= dc1[0]-dc2[0];    g2[1]= dc1[1]-dc2[1];    g2[2]= dc1[2]-dc2[2];
    g3[0]= dc2[0]-dc3[0];    g3[1]= dc2[1]-dc3[1];    g3[2]= dc2[2]-dc3[2];
    g4[0]= dc3[0];           g4[1]= dc3[1];           g4[2]= dc3[2];

/*     printf("COS) % 10.4lf deg g0 % 10.4lf% 10.4lf% 10.4lf g1 % 10.4lf% 10.4lf% 10.4lf g2 % 10.4lf% 10.4lf% 10.4lf g3 % 10.4lf% 10.4lf% 10.4lf\n", */
/* 	   phi*180/M_PI,g1[0],g1[1],g1[2],g2[0],g2[1],g2[2],g3[0],g3[1],g3[2],g4[0],g4[1],g4[2]); */
  }
  else {
    double rct = 1.0/cos_phi;
    double ds1[3], ds2[3], ds3[3], tmp[3], tmp2[3];
    double DD[3][3];
    double c123_r=1.0/c123;

    if (c123_r!=c123_r) {
      fprintf(stderr,"ERROR: nan in sin-dihed @ c123 %.5lf\n",c123);
      exit(-1);
    }
    if (rct!=rct) {
      fprintf(stderr,"ERROR: nan in sin-dihed @ rct %.5lf\n",rct);
      exit(-1);
    }

    myscale(1.0/c123,C123);

    mycopy(tmp,C123);
    myscale(sin_phi,tmp);
    mydiff(ds1,C23,tmp);
    myscale(1.0/c123,ds1);
    mycopy(tmp,ds1);
    myddot(DD,R23);
    mymatvec(ds1,DD,tmp);

    mycopy(tmp,C123);
    myscale(sin_phi,tmp);
    mydiff(ds2,C23,tmp);
    myscale(1.0/c123,ds2);
    mytdot(DD,R12,R23);
    mycopy(tmp,ds2);
    mymatvec(ds2,DD,tmp);

    mycopy(tmp,C23);
    myscale(sin_phi,tmp);
    mydiff(tmp,C123,tmp);
    myscale(1.0/c23,tmp);
    mycross(tmp2,R34,tmp);
    mysum(ds2,ds2,tmp2);

    mycross(ds3,tmp,R23);

    myscale(rct,ds1);
    myscale(rct,ds2);
    myscale(rct,ds3);

    g1[0]=-ds1[0];           g1[1]=-ds1[1];           g1[2]=-ds1[2];
    g2[0]= ds1[0]-ds2[0];    g2[1]= ds1[1]-ds2[1];    g2[2]= ds1[2]-ds2[2];
    g3[0]= ds2[0]-ds3[0];    g3[1]= ds2[1]-ds3[1];    g3[2]= ds2[2]-ds3[2];
    g4[0]= ds3[0];           g4[1]= ds3[1];           g4[2]= ds3[2];

/*     printf("SIN) % 10.4lf deg g0 % 10.4lf% 10.4lf% 10.4lf g1 % 10.4lf% 10.4lf% 10.4lf g2 % 10.4lf% 10.4lf% 10.4lf g3 % 10.4lf% 10.4lf% 10.4lf\n", */
/* 	   phi*180/M_PI,g1[0],g1[1],g1[2],g2[0],g2[1],g2[2],g3[0],g3[1],g3[2],g4[0],g4[1],g4[2]); */
  }

  return phi;

}

double my_whitenoise ( unsigned short xsubi[3]  ) {
  //  return (1-2*erand48(xsubi))*M_SQRT3;
  return ((erand48(xsubi)+erand48(xsubi)+erand48(xsubi)+erand48(xsubi)+erand48(xsubi)+erand48(xsubi))-3.0)*M_SQRT2;
}

