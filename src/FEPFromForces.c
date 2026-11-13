/* Free-energy profile along a string from the mean force at each image.

   Input, one line per image in string order: z_1 ... z_dim f_1 ... f_dim,
   where f is the mean restraint force kappa (theta - z) at the image, i.e.
   -grad F. Lines starting with # are skipped.

   The profile is integrated by the trapezoid rule on the chords between
   neighboring images:

     F(j+1) = F(j) + 1/2 (f_j + f_{j+1}) . (z_j - z_{j+1}),   F(0) = 0

   and s, the arc length, is the running sum of chord lengths. With -dihed,
   differences in z are wrapped into [-pi, pi].

   Output, one line per image: j s F | z_1..z_dim f_1..f_dim t_1..t_dim,
   where t is the unit tangent (central difference inside the string,
   one-sided at the ends), pointing from image j toward image j-1.

 (c) 2016-2026 cameron f abrams
 drexel university
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

typedef struct DATUM * pdat;
typedef struct DATUM {
  int dim;
  double * r;  // image location in Z-space (input)
  double * f;  // mean force vector (input)
  double * tv; // tangent vector along string (calculated)
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

static int dihed=0;

/* a - b, wrapped into [-pi, pi] when the CVs are angles */
double cvdiff ( double a, double b ) {
  double d=a-b;
  if (dihed) {
    while (d>M_PI) d-=2*M_PI;
    while (d<-M_PI) d+=2*M_PI;
  }
  return d;
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
     if (d->next && d->prev) for (i=0;i<d->dim;i++) d->tv[i]=0.5*cvdiff(d->prev->r[i],d->next->r[i]);
     else if (d->next && !d->prev) for (i=0;i<d->dim;i++) d->tv[i]=cvdiff(d->r[i],d->next->r[i]);
     else if (!d->next && d->prev) for (i=0;i<d->dim;i++) d->tv[i]=cvdiff(d->prev->r[i],d->r[i]);
     else { fprintf(stderr,"ERROR: list is not formed correctly.  This is a bug.\n"); exit(-1);}
     l=0.0;
     for (i=0;i<d->dim;i++) l+=d->tv[i]*d->tv[i];
     l=sqrt(l);
     for (i=0;i<d->dim;i++) d->tv[i]/=l;
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
    for (i=0;i<d->dim;i++) fprintf(ofp,"%.6lf ",d->r[i]);
    for (i=0;i<d->dim;i++) fprintf(ofp,"%.6lf ",d->f[i]);
    for (i=0;i<d->dim;i++) fprintf(ofp,"%.6lf ",d->tv[i]);
    fprintf(ofp,"\n");
    j++;
    /* trapezoid across the chord to the next image */
    if (d->next) {
      double l=0.0, dz;
      for (i=0;i<d->dim;i++) {
        dz=cvdiff(d->r[i],d->next->r[i]);
        l+=dz*dz;
        fep+=0.5*(d->f[i]+d->next->f[i])*dz;
      }
      s+=sqrt(l);
    }
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

  dat * D=NULL, * d=NULL;


  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-i")) ifn=argv[++i];
    else if (!strcmp(argv[i],"-o")) ofn=argv[++i];
    else if (!strcmp(argv[i],"-dim")) dim=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-dihed")) dihed=1;
  }

  if (dim<1) {
    fprintf(stderr,"ERROR: must specify CV-space dimensionality with -dim\n");
    exit(-1);
  }
  if (!ifn || !ofn) {
    fprintf(stderr,"Usage: FEPFromForces -i <image forces> -o <output> -dim <number of CVs> [-dihed]\n");
    exit(-1);
  }
  if (!(ifp=fopen(ifn,"r"))) {
    fprintf(stderr,"ERROR: cannot open %s\n",ifn);
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
  if (ni<2) {
    fprintf(stderr,"ERROR: %s has %d image(s); a string needs at least 2\n",ifn,ni);
    exit(-1);
  }
  if (!(ofp=fopen(ofn,"w"))) {
    fprintf(stderr,"ERROR: cannot open %s\n",ofn);
    exit(-1);
  }

  fprintf(stdout,"INFO: %d-dimensional data read from %d images\n",dim,ni);


  compute_tangents(D);
  write_dat(D,ofp);  

  fclose(ofp);
}
