#include "basisrep.h"

/* below are some input routines I should really have in a general library */

char * next_word ( char * p ) {
  if (p&&*p) {
    while (*p&&!isspace(*(p++)));
    while (*p&&isspace(*(p++)));
    if (*p) return --p;
  }
  return NULL;
}

int sscanf_stringList ( char * str, char ** slist, int n ) {
  if (str) {
    char * p=str;
    int i=0;
    sscanf(p,"%s",slist[i++]);
    while (p=next_word(p)) sscanf(p,"%s",slist[i++]);
    return !(i==n);
  }
}

int sscanf_doubleList ( char * str, double * dlist, int n ) {
  if (str) {
    char * p=str;
    int i=0;
    sscanf(p,"%lf",&dlist[i++]);
    while (p=next_word(p)) sscanf(p,"%lf",&dlist[i++]);
    return !(i==n);
  }
}

lRecon * lRecon_new ( int D ) {
   lRecon * p = (lRecon*)malloc(sizeof(lRecon)); 
   p->r=(double*)malloc(D*sizeof(double));
   p->f=(double*)malloc(D*sizeof(double));
   p->D=D;
   p->e=0.0;
   p->i=0;
   return p;
}

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
int vecdiff ( double * c, double * a, double * b, int D, int * per, double * dom ) {
  int i;
  for (i=0;i<D;i++) {
    c[i]=a[i]-b[i];
    if (per[i]) { // dimension is periodic
       if (c[i] < -0.5*dom[i]) c[i]+=dom[i];
       if (c[i] > 0.5*dom[i])  c[i]-=dom[i];
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

int localReconstruct ( lRecon * z, bCntr * bsc, int nC, int D, int * per, double * dom ) {
  double E = 0.0;
  int i=0,j;
  bCntr * ip;
  double Zi[D];
  double zi;
  double phi, gphi;

  z->e=0.0;
  for (j=0;j<D;j++) {
    z->f[j]=0.0;
  }
  for (i=0;i<nC;i++) {
    ip=&bsc[i];
    vecdiff(Zi,z->r,ip->r,D,per,dom);
    zi=norm(Zi,D);
    kernel(zi,ip->s,&phi,&gphi);
    z->e+=ip->a*phi;
    for (j=0;j<D;j++) {
      z->f[j]-=(zi>0)?ip->a*Zi[j]/zi*gphi:0.0;
    }
  }
  return 0;
}

#ifndef NDATA
#define NDATA 1000
#endif
// to do -- implement variable CV dimensionality
bCntr * readBasisSet ( int * nC, int D, char * fn ) {
  int i;
  FILE * fp = fopen(fn,"r");
  if (fp) {
    bCntr * bp;
    bCntr * B = malloc(NDATA*sizeof(bCntr));
    char ln[255];
    i=0;
    while (fgets(ln,255,fp)) {
      bp=&B[i];
      bp->i=i;
      bp->r=(double*)malloc(D*sizeof(double));
      if (sscanf_doubleList(ln,bp->r,D)) {
        fprintf(stderr,"ERROR: could not read basis set from %s.\n",fn);
	exit(-1);
      }
      i++;
    }
    *nC=i;
    B=realloc(B,i*sizeof(bCntr));
    return B;
  }
  *nC=0;
  return NULL;
}

int writeBasisSet ( bCntr * bsc, int nC, int D, char * fn ) {
  FILE * fp = fopen(fn,"w");
  int i,j;
  bCntr * ip;
  for (i=0;i<nC;i++) {
    ip=&bsc[i];
    ip->i=i;
    for (j=0;j<D;j++) {
      fprintf(fp,"% 10.4lf",ip->r[j]);
    }
    fprintf(fp,"% 10.4lf%10.4lf\n",ip->a,ip->s);
  }
  fclose(fp);
}

int writeReconstruction ( lRecon * recon, int nR, int D, char * fn ) {
  FILE * fp = fopen(fn,"w");
  int i,j;
  lRecon * ip;
  for (i=0;i<nR;i++) {
    ip=&recon[i];
    for (j=0;j<D;j++) {
      fprintf(fp,"%.5lf ",ip->r[j]);
    }
    for (j=0;j<D;j++) {
      fprintf(fp,"%.5lf ",ip->f[j]);
    }
    fprintf(fp,"%.5lf\n",ip->e);
  }
  fclose(fp);
}
