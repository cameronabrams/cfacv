#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef CVDIM
#define CVDIM 3
#endif

typedef struct POINT {
  int i;
  int d;
  double * r;
} pt;

pt * new_pt ( int d ) { 
  pt * p=(pt*)malloc(sizeof(pt));
  p->r=(double*)malloc(d*sizeof(pt));
  p->i=0;
  p->d=d;
  return p;
}

pt * pt_malloc ( pt * p, int i, int d ) {
  if (p) {
    p->r=(double*)malloc(d*sizeof(double));
    p->i=i;
    p->d=d;
  }
  return p;
}

#ifndef NDATA
#define NDATA 100000
#endif

void vecdiff ( pt * c, pt * a, pt * b, int dihed) {
  int i;
  if (a->d!=b->d) {
    fprintf(stderr,"WARNING: attempting binary operation on vectors with different dimensions is not allowed.\n");
    return;
  }
  for (i=0;i<a->d;i++) {
    c->r[i]=a->r[i]-b->r[i];
    if (dihed) {
      if (c->r[i] < -M_PI) c->r[i]+=2*M_PI;
      if (c->r[i] > M_PI)  c->r[i]-=2*M_PI;
    }
  }
}

void vecadd ( pt * c, pt * a, pt * b ) {
  int i;
  for (i=0;i<a->d;i++) {
    c->r[i]=a->r[i]+b->r[i];
  }
}

void vecscale ( pt * c, pt * a, double x ) {
  int i;
  for (i=0;i<a->d;i++) {
    c->r[i]=a->r[i]*x;
  }
}

char * next_word ( char * p ) {
  if (p&&*p) {
    while (*p&&!isspace(*p)) {
      //printf("{%s}\n",p);fflush(stdout);
      p++;
    }
    while (*p&&isspace(*p)) {
      //printf(" {%s}\n",p);fflush(stdout);
      p++;
    }
    if (*p) {
      //     printf(" r{%s}\n",p);
      return p;
    }
  }
  return NULL;
}

int main ( int argc, char * argv[] ) {

  char * lfn="run.log";
  char * ofn="fra.dat";
  char ln[2000];
  FILE * fp;
  char testStr[9], *p;

  int i,line,j;
  int dim=CVDIM; // default

  double k = 100.0; // kcal/mol/A^2

  pt * z=(pt*)malloc(NDATA*sizeof(pt)), * zp;
  pt * th=(pt*)malloc(NDATA*sizeof(pt)), * thp;
  pt * fd=(pt*)malloc(NDATA*sizeof(pt)), * fdp;
  pt * rsum, * ravg;
  pt * rrsum, * rravg;
  pt Disp_th_z, * d=&Disp_th_z;

  int nz, nth, nfd;
  int dihed = 0;

  int quiet=0;

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-f")) lfn=argv[++i];
    if (!strcmp(argv[i],"-k")) k=atof(argv[++i]);
    if (!strcmp(argv[i],"-dim")) dim=atoi(argv[++i]);
    if (!strcmp(argv[i],"-ofn")) ofn=argv[++i];
    if (!strcmp(argv[i],"-dihed")) dihed=1;
    if (!strcmp(argv[i],"-quiet")) quiet=1;
  }

  for (i=0;i<NDATA;i++) {
    pt_malloc(&z[i],i,dim);
    pt_malloc(&th[i],i,dim);
    pt_malloc(&fd[i],i,dim);
  }
  pt_malloc(d,0,dim);

  fp=fopen(lfn,"r");
  if (fp) {
    line=1;
    nz=nth=nfd=0;
    while (fgets(ln,1200,fp)) {
      sscanf(ln,"%8s",testStr);
      //     fprintf(stderr,"TEST: [%s]\n",testStr);
      if (!strcmp(testStr,"CFACV/C)")) {
	p=ln+8;
	sscanf(p,"%3s",testStr);
	if (!strcmp(testStr,"Z")) {
//		fprintf(stderr,"TEST: [%s]\n",testStr);
	  p+=3;
	  zp=&z[nz];
	  p=next_word(p);
	  sscanf(p,"%i",&zp->i);
	  j=0;
	  while (p=next_word(p)) {
	    if (j==dim) {
	      exit(-1);
	    }
	    sscanf(p,"%lf",&zp->r[j++]);
	    
	    //if (!zp->i) fprintf(stdout,"   [%i][% 10.4lf]\n",j,zp->r[j-1]);
	  }
	  nz++;
	}
	else if (!strcmp(testStr,"Th")) {
	  p+=3;
	  thp=&th[nth];
	  p=next_word(p);
	  sscanf(p,"%i",&thp->i);
	  j=0;
	  while (p=next_word(p))
	    sscanf(p,"%lf",&thp->r[j++]);
	  nth++;
        }
        else if (!strcmp(testStr,"FD")) {
	  p+=3;
	  fdp=&fd[nfd];
	  p=next_word(p);
	  sscanf(p,"%i",&fdp->i);
	  j=0;
	  while (p=next_word(p))
	    sscanf(p,"%lf",&fdp->r[j++]);
	  nfd++;
        }
	else if (!strcmp(testStr,"Ver")) {
	}
        else if (!strcmp(testStr,"Set")) {
        }
        else if (!strcmp(testStr,"mt")) {
        }
        else if (!strcmp(testStr,"Met")) {
        }
	else {
	  fprintf(stderr,"ERROR: bad designation at line %i in %s: %s\n",line,lfn,testStr);
	  exit(-1);
	}
      }
      else if (!strcmp(testStr,"CFACV)")) {
        p=ln+6;
        //fprintf(stderr,"test line %s\n",p);
      }
      line++;
    }
  }
  else {
    fprintf(stderr,"ERROR: Could not open %s.\n",lfn);
    exit(-1);
  }

  if (!quiet) fprintf(stderr,"INFO) Read %i z-values, %i theta-values, and %i force-values from %s.\n",nz,nth,nfd,lfn);

  if (nz==nth && nz==nfd) {
    rsum=(pt*)calloc(nz,sizeof(pt));
    for (i=0;i<nz;i++) pt_malloc(&rsum[i],i,dim);
    ravg=(pt*)calloc(nz,sizeof(pt));
    for (i=0;i<nz;i++) pt_malloc(&ravg[i],i,dim);
    rrsum=(pt*)calloc(nfd,sizeof(pt));
    for (i=0;i<nz;i++) pt_malloc(&rrsum[i],i,dim);
    rravg=(pt*)calloc(nfd,sizeof(pt));
    for (i=0;i<nz;i++) pt_malloc(&rravg[i],i,dim);
    //fprintf(stderr,"INFO) running average allocation complete\n");fflush(stderr);
    vecdiff(d,&th[0],&z[0],dihed);
    memcpy(rsum[0].r,d->r,dim*sizeof(double));
    memcpy(ravg[0].r,d->r,dim*sizeof(double));
    memcpy(rrsum[0].r,fd[0].r,dim*sizeof(double));
    memcpy(rravg[0].r,fd[0].r,dim*sizeof(double));

    for (i=1;i<nz;i++) {
      thp=&th[i]; zp=&z[i]; fdp=&fd[i];
      vecdiff(d,thp,zp,dihed);
      vecadd(&rsum[i],&rsum[i-1],d);
      vecscale(&ravg[i],&rsum[i],k/(1.0+i));
      vecadd(&rrsum[i],&rrsum[i-1],fdp);
      vecscale(&rravg[i],&rrsum[i],1.0/(1.0+i));
    }

    fp=fopen(ofn,"w");
    fprintf(fp,"# k = %.5lf\n",k);
    for (i=0;i<nz;i++) {
      fprintf(fp,"%i ",i);
      for (j=0;j<dim;j++) fprintf(fp,"%.8lf ",ravg[i].r[j]);
      fprintf(fp," - ");
      for (j=0;j<dim;j++) fprintf(fp,"%.8lf ",rravg[i].r[j]);
      fprintf(fp,"\n");
    }
    fclose(fp);
  //  fprintf(stderr,"Generated %s.\n",ofn);
    i--;
    //thp=&th[i];
    for (j=0;j<dim;j++) fprintf(stdout,"%.8lf ",zp->r[j]);
    for (j=0;j<dim;j++) fprintf(stdout,"%.8lf ",ravg[i].r[j]);
    for (j=0;j<dim;j++) fprintf(stdout,"%.8lf ",rravg[i].r[j]);
    fprintf(stdout,"\n");
  }

//  fprintf(stderr,"Program ends.\n");

}
