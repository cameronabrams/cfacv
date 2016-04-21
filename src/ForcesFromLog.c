#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef SDIM
#define SDIM 3
#endif

typedef struct POINT {
  int i;
  double r[SDIM];
} pt;

#ifndef NDATA
#define NDATA 100000
#endif

void vecdiff ( pt * c, pt * a, pt * b, int dihed) {
  int i;
  for (i=0;i<SDIM;i++) {
    c->r[i]=a->r[i]-b->r[i];
    if (dihed) {
      if (c->r[i] < -M_PI) c->r[i]+=2*M_PI;
      if (c->r[i] > M_PI)  c->r[i]-=2*M_PI;
    }
  }
}

void vecadd ( pt * c, pt * a, pt * b ) {
  int i;
  for (i=0;i<SDIM;i++) {
    c->r[i]=a->r[i]+b->r[i];
  }
}

void vecscale ( pt * c, pt * a, double x ) {
  int i;
  for (i=0;i<SDIM;i++) {
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
  int dim=SDIM;

  double k = 100.0; // kcal/mol/A^2

  pt * z=(pt*)malloc(NDATA*sizeof(pt)), * zp;
  pt * th=(pt*)malloc(NDATA*sizeof(pt)), * thp;
  pt * rsum, * ravg;
  pt Disp_th_z, * d=&Disp_th_z;

  int nz, nth;
  int dihed = 0;

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-f")) lfn=argv[++i];
    if (!strcmp(argv[i],"-k")) k=atof(argv[++i]);
    if (!strcmp(argv[i],"-dim")) dim=atoi(argv[++i]);
    if (!strcmp(argv[i],"-ofn")) ofn=argv[++i];
    if (!strcmp(argv[i],"-dihed")) dihed=1;
  }

  fp=fopen(lfn,"r");
  if (fp) {
    line=1;
    nz=nth=0;
    while (fgets(ln,1200,fp)) {
      sscanf(ln,"%8s",testStr);
      //     fprintf(stderr,"TEST: [%s]\n",testStr);
      if (!strcmp(testStr,"CFACV/C)")) {
	p=ln+8;
	sscanf(p,"%3s",testStr);
	//	fprintf(stderr,"TEST: [%s]\n",testStr);
	if (!strcmp(testStr,"Z")) {
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
        }
	else if (!strcmp(testStr,"Ver")) {
	}
	else {
	  fprintf(stderr,"ERROR: bad designation at line %i in %s: %s\n",line,lfn,testStr);
	  exit(-1);
	}
      }
      line++;
    }
  }
  else {
    fprintf(stderr,"ERROR: Could not open %s.\n",lfn);
    exit(-1);
  }

  //fprintf(stderr,"INFO) Read %i z-values and %i theta-values from %s.\n",nz,nth,lfn);

  if (nz==nth) {
    rsum=(pt*)calloc(nz,sizeof(pt));
    ravg=(pt*)calloc(nz,sizeof(pt));
    vecdiff(d,&th[0],&z[0],dihed);
    memcpy(rsum[0].r,d->r,dim*sizeof(double));
    memcpy(ravg[0].r,d->r,dim*sizeof(double));
    for (i=1;i<nz;i++) {
      thp=&th[i]; zp=&z[i];
      vecdiff(d,thp,zp,dihed);
      vecadd(&rsum[i],&rsum[i-1],d);
      vecscale(&ravg[i],&rsum[i],k/(1.0+i));
    }

    fp=fopen(ofn,"w");
    fprintf(fp,"# k = %.5lf\n",k);
    for (i=0;i<nz;i++) {
      fprintf(fp,"%i ",i);
      for (j=0;j<dim;j++) fprintf(fp,"%.8lf ",ravg[i].r[j]);
      fprintf(fp,"\n");
    }
    fclose(fp);
  //  fprintf(stderr,"Generated %s.\n",ofn);
    i--;
    //thp=&th[i];
    for (j=0;j<dim;j++) fprintf(stdout,"%.8lf ",zp->r[j]);
    for (j=0;j<dim;j++) fprintf(stdout,"%.8lf ",ravg[i].r[j]);
    fprintf(stdout,"\n");
  }

//  fprintf(stderr,"Program ends.\n");

}
