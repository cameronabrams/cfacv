#include "basisrep.h"

/* Gaussian packet definition
 *
 * \phi_\sigma(z) = \phi(z/\sigma) = exp(-(1/2)(z/\sigma)^2)
 *
 * 
 */
void kernel ( double z, double sig, double * f, double * g ) {
  double is2=1.0/(sig*sig);
  *f=exp(-0.5*z*z*is2);
  *g=-z*is2*(*f);
}

/* returns the vector difference a-b in c */
int vecdiff ( double * c, double * a, double * b, int D, int * per, double * doms ) {
  int i;
  for (i=0;i<D;i++) {
    c[i]=a[i]-b[i];
    if (per[i]) { // dimension is periodic
       if (c[i] < -0.5*doms[i]) c[i]+=doms[i];
       if (c[i] > 0.5*doms[i])  c[i]-=doms[i];
    }
  }
  return 0;
}

/* returns the dot product of vectors a and b */
double vecdot ( double * a, double * b, int D ) {
  double dot = 0.0;
  int i;
  for (i=0;i<D;i++) dot+=a[i]*b[i];
  return dot;
}

/* returns the Euclidean norm of vector x */
double norm ( double * x, int D ) {
  return sqrt(vecdot(x,x,D));
}

int localReconstruct ( lRecon * z, bCntr * bsc, int nC ) {
  double E = 0.0;
  int i=0,j;
  bCntr * ip;
  double Zi[SDIM];
  double zi;
  double phi, gphi;

  z->e=0.0;
  for (j=0;j<SDIM;j++) {
    z->f[j]=0.0;
  }
  for (i=0;i<nC;i++) {
    ip=&bsc[i];
    vecdiff(Zi,z->r,ip->r,z->nCV,z->per,z->rdom);
    zi=norm(Zi,SDIM);
    kernel(zi,ip->s,&phi,&gphi);
    z->e+=ip->a*phi;
    for (j=0;j<SDIM;j++) {
      z->f[j]-=(zi>0)?ip->a*Zi[j]/zi*gphi:0.0;
    }
  }
  return 0;
}

#ifndef NDATA
#define NDATA 1000
#endif
// to do -- implement variable CV dimensionality
bCntr * readBasisSet ( int * nC, char * fn ) {
  int i;
  FILE * fp = fopen(fn,"r");
  if (fp) {
    bCntr * bp;
    bCntr * B = malloc(NDATA*sizeof(bCntr));
    char ln[255];
    i=0;
    while (fgets(ln,255,fp)) {
      bp=&B[i];
      sscanf(ln,"%i %lf %lf %lf %lf %lf",
	     &bp->i,&bp->r[0],&bp->r[1],&bp->r[2],&bp->a,&bp->s);
      i++;
    }
    *nC=i;
    B=realloc(B,i*sizeof(bCntr));
    return B;
  }
  *nC=0;
  return NULL;
}

int writeBasisSet ( bCntr * bsc, int nC, char * fn ) {
  FILE * fp = fopen(fn,"w");
  int i,j;
  bCntr * ip;
  for (i=0;i<nC;i++) {
    ip=&bsc[i];
    fprintf(fp,"%i ",ip->i);
    for (j=0;j<ip->nCV;j++) {
      fprintf(fp,"% 10.4lf",ip->r[j]);
    }
    fprintf(fp,"% 10.4lf%10.4lf\n",ip->a,ip->s);
  }
  fclose(fp);
}

int writeReconstruction ( lRecon * recon, int nR, char * fn ) {
  FILE * fp = fopen(fn,"w");
  int i,j;
  lRecon * ip;
  for (i=0;i<nR;i++) {
    ip=&recon[i];
    fprintf(fp,"%i ",ip->i);
    for (j=0;j<ip->nCV;j++) {
      fprintf(fp,"%.5lf ",ip->r[j]);
    }
    for (j=0;j<ip->nCV;j++) {
      fprintf(fp,"%.5lf ",ip->f[j]);
    }
    fprintf(fp,"%.5lf\n",ip->e);
  }
  fclose(fp);
}
