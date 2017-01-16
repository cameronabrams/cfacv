/* Free-energy profile calculation from image forces along string
 (c) 2016 cameron f abrams
 drexel university
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct DATUM * pdat;
typedef struct DATUM {
  int dim;
  double * r;  // image location in Z-space (input)
  double * f;  // mean force vector (input)
  double * tv; // tangent vector along string (calculated)
  double tf; // tangent mean force * interimage distance(calculated)
  double dr;   // distance to next point (calculated)
  pdat next;
  pdat prev;
} dat;

int buf2wordlist ( char * buf, char ** words, int MAXWORDS ) {
  char * p=buf;
  int n=0;
  if (p) {
    while (p&&*p) {
      while (*p&&isspace(*p)) p++;
      words[n++]=p;
      if (n==MAXWORDS) return n;
      while (*p&&!isspace(*p)) p++;
      if (*p) {*p='\0'; p++;}
    }
  }
  return n;
}

double my_dot ( double * a, double * b, int n ) {
  int i;
  double dp=0.0;
  for (i=0;i<n;i++) dp+=a[i]*b[i];
  return dp;
}


dat * new_datum ( char * buf, int dim ) {
  
   int t, i;
   char * words[40];
//   fprintf(stdout,"INFO: new datum of dim %d\n",dim);fflush(stdout);

   dat * d=malloc(sizeof(dat));
   d->dim=dim;
   d->r=malloc(dim*sizeof(double));
   d->f=malloc(dim*sizeof(double));
   d->tv=malloc(dim*sizeof(double));

//   fprintf(stdout,"INFO: buf %s 2*dim %d\n",buf,2*dim); fflush(stdout);

   t=buf2wordlist(buf,words,2*dim);

//   fprintf(stdout,"INFO: new datum: buf2wordlist detects %d words\n",t);
//   fflush(stdout);

   for (i=0;i<dim;i++) {
     d->r[i]=atof(words[i]);
     d->f[i]=atof(words[i+dim]);
   }
   d->next=NULL;
   d->prev=NULL;
   return d;
}

void do_tan ( dat * d ) {
   int i;
   double l;
   if (d) {
     if (d->next && d->prev) for (i=0;i<d->dim;i++) d->tv[i]=0.5*(d->prev->r[i]-d->next->r[i]);
     else if (d->next && !d->prev) for (i=0;i<d->dim;i++) d->tv[i]=d->r[i]-d->next->r[i];
     else if (!d->next && d->prev) for (i=0;i<d->dim;i++) d->tv[i]=d->prev->r[i]-d->r[i];
     else { fprintf(stderr,"ERROR: list is not formed correctly.  This is a bug.\n"); exit(-1);}
     l=0.0;
     for (i=0;i<d->dim;i++) l+=d->tv[i]*d->tv[i];
     d->dr=sqrt(l);
     for (i=0;i<d->dim;i++) d->tv[i]/=d->dr;
     d->tf=my_dot(d->f,d->tv,d->dim);
   }
}

void compute_tangents ( dat * D ) {
   dat * d;
   for (d=D;d;d=d->next) {
      do_tan(d);
   }
}

void write_dat ( dat * D, FILE * ofp ) {
  int i=0,j=0;
  double s=0.0;
  double fep=0.0;
  dat * d;

  for (d=D;d;d=d->next) {
    fprintf(ofp,"%i %.5lf %.5le | ",j,s,fep);
    fep+=d->tf*d->dr;
    s+=d->dr;
    j++;
    for (i=0;i<d->dim;i++) fprintf(ofp,"%.6lf ",d->r[i]);
    for (i=0;i<d->dim;i++) fprintf(ofp,"%.6lf ",d->f[i]);
    for (i=0;i<d->dim;i++) fprintf(ofp,"%.6lf ",d->tv[i]);
    fprintf(ofp,"\n");
  }
}

int main ( int argc, char * argv[] ) {

  char * ifn=NULL;
  char * ofn=NULL;

  FILE * ifp=NULL;
  FILE * ofp=NULL;

  int dim=-1;

  int i;

  int ni=0;
  char ln[255];

  dat * D, * d;


  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-i")) ifn=argv[++i];
    else if (!strcmp(argv[i],"-o")) ofn=argv[++i];
    else if (!strcmp(argv[i],"-dim")) dim=atoi(argv[++i]);
  }

  if (ifn) ifp=fopen(ifn,"r");
  if (ofn) ofp=fopen(ofn,"w");
  if (dim==-1) {
    fprintf(stderr,"ERROR: must specify CV-space dimensionality with -d\n");
    exit(-1);
  }

//  fprintf(stdout,"INFO: Reading from %s\n",ifn);fflush(stdout);

  while (fgets(ln,255,ifp)) {
    if (ln[0]!='#') {
      if (!ni) { D=new_datum(ln,dim); d=D; ni++; }
      else { d->next=new_datum(ln,dim); d->next->prev=d; d=d->next; ni++; }
    }
  }
  fclose(ifp);
 
  fprintf(stdout,"INFO: %d-dimensional data read from %d images\n",dim,ni);


  compute_tangents(D);
  write_dat(D,ofp);  

  fclose(ofp);
}
