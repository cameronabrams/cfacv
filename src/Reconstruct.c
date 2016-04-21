/* Reconstruct.c -- solves the linear algebraic system 
 * 
 * B a = c
 *
 * where a is the vector of coefficients used in the
 * radial basis function decomposition of free energy; i.e.,
 *
 * A(z) = \sum_k a_k \phi_sigma(|z-z_k|);
 *
 * B is the matrix of gradient cross-products; i.e.,
 *
 * B_{k,k^\prime} = \sum_{k^{\prime\prime}}\grad_z\phi_\sigma(|z_k-z_{k^{\prime\prime}}|)\grad_z\phi_\sigma(|z_{k^{\prime\prime}}-z_{k^\prime}|);
 *
 * and c is
 * 
 * c_k = -\sum_{k^\prime}\grad_z\phi_\sigma(|z_k-z_{k^{\prime}})f_{k^\prime};
 * 
 * where f_k is the vector force at point k measured from restrained MD.
 *
 * Note that the GNU scientific library must be installed 
 *
 * Reference: Maragliano and Vanden-Eijnden, JCP 128:184110 (2008).
 *
 * (c) 2010-2016
 * Cameron F Abrams
 * Drexel University
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include "basisrep.h"

#ifndef NDATA
#define NDATA 1000
#endif
#ifndef EPSILON
#define EPSILON 1.0e-6
#endif

char * next_word ( char * p ) {
  if (p&&*p) {
    while (*p&&!isspace(*(p++)));
    while (*p&&isspace(*(p++)));
    if (*p) return p;
  }
  return NULL;
}

int sscanf_stringList ( char * str, char ** slist, int n ) {
  if (str) {
    char * p=str;
    int i=0;
    sscanf(p,"%s",slist[i++]);
    while (p=next_word(p)) sscanf(p,"%s",slist[i++]);
    return (i==n);
  }
}

/* EXPERIMENTAL:  Monte Carlo optimization of quadrature points and sigma */
typedef struct MC {
  int nc;
  unsigned long int Seed;
  double dR;
  double dS;
  char flag[5];
} mct;

int _modified_single_sweep_ = 0;

/* read_samples
 *
 * reads in data for centers from a formatted ASCII input file.
 * Lines at beginning of file that begin w/ # are comments.
 * The first non-comment line is of the form
 *
 * NCV NDATA
 *
 * where NCV is the integer dimensionality of CV space and NDATA 
 * is the number of centers.
 *
 * The next line has NCV keywords indicating the type of CV component.
 * Currently, the only keyword that matters is DIHED; other values
 * are ignored.
 * 
 * The next line is blank. Very important.
 *
 * Following the blank lines is a set of NDATA lines, each of the format:
 *
 * r_0  r_1 ... r_(NCV-1)  f_0  f_1  ... f_(NCV-1)
 *
 * where r_i is the i'th component of the CV and f_i is the i'th
 * component of the force.
 * 
 */
typedef struct FORCE_SAMPLE {
  int ncv;
  double * r; // location (dim)
  double * f; // force (kcal/mol/dim)
  int i;       // unique index
  int u;       // use-flag; 1 means point is used, 0 means not used
  int z;       // zeroed for begin below threshold
} fSamp;

fSamp * read_samples ( int * n, char * fn, int D, int quiet ) {
  FILE * fp = stdin;
  char scrln[2048];
  int i=0;
  fSamp * data, * d;
  char * c;
  int j;
  int ndata, ncv;
  char ** keyword;
  double * datalist;

  if (fn) {
    fp = fopen(fn,"r");
  }

  /* read past comments */
  while (fgets(scrln,20,fp) && scrln[0]=='#');
  /* read the required line and the blank line */
  sscanf(scrln,"%i %i",&ncv,&ndata);
  if (ncv!=D) {
    fprintf(stderr,"ERROR: mismatch of CV space dimensionality\n");
    exit(-1);
  }
  keyword=(char**)malloc(ncv*sizeof(char*));
  for (j=0;j<ncv;j++) keyword[i]=(char*)malloc(100*sizeof(char));
  /* read and parse the line that provides a keyword description for each CV */
  fgets(scrln,1000,fp);
  if (sscanf_stringList(scrln,keyword,ncv)) {
    fprintf(stderr,"ERROR: could not read keyword list from %s.\n",fn);
    exit(-1);
  }
  /* allocate the data array */
  data = (fSamp*)malloc(ndata*sizeof(fSamp));
  datalist = (double*)malloc(2*ncv*sizeof(double));

  for (i=0;i<ndata;i++) {
    d=&data[i];
    d->ncv=ncv;
    d->r=(double*)malloc(ncv*sizeof(double));
    d->f=(double*)malloc(ncv*sizeof(double));
    /* read next line from file */
    fgets(scrln,2048,fp);
    /* parse line into center and force data arrays */
    if (sscanf_doubleList(scrln,datalist,2*ncv)) {
      fprintf(stderr,"ERROR: could not read data on line %i from %s.\n",i,fn);
      exit(-1);
    }
    if (j=0;j<ncv;j++) {
      d->r[j]=datalist[j];
      d->f[j]=datalist[ncv+j];
    }
  }
  *n=ndata;
  fclose(fp);
  free(datalist);
  return data;
}

/* writes out center/forces samples in the same format as the input */

int write_samples ( fSamp  * data, int n, char * fn, int quiet, char ** cv_keywords ) {
  FILE * fp = stdout;
  int i=0,j;
  fSamp * c;
  int ncv;

  fp=fopen(fn,"w");
  if (!fp) return -1;
  fprintf(fp,"# forces.dat\n");
  fprintf(fp,"# generated by Reconstruct::write_samples()\n");
  fprintf(fp,"#\n");
  fprintf(fp,"%i %i\n",ncv=data[0].ncv,n); 
  for (j=0;j<ncv;j++) fprintf(fp,"%s",cv_keywords[j]);
  for (i=0;i<n;i++) {
    c=&data[i];
    for (j=0;j<ncv;j++) {
      fprintf(fp,"%.8lf ",c->r[j]);
    }
    for (j=0;j<ncv;j++) {
      fprintf(fp,"%.8lf ",c->f[j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  return 0;
}

/* Sets up the set of basis function centers, where each center is located
 * on a distinct sample site (a "sample-congruent" basis set).
 */
bCntr * SetupSampleCongruentBasisCenters ( fSamp * data, int nD, double sig, int * nC ) {
  int i;
  bCntr * bsc = (bCntr*)malloc(nD*sizeof(bCntr));
  fSamp * dp;
  bCntr * bp;

  for (i=0;i<nD;i++) {
    dp=&data[i];
    bp=&bsc[i];
    memcpy(bp->r,dp->r,dp->ncv*sizeof(double));
    //fprintf(stderr,"# basis center %i at %.5lf %.5lf %.5lf\n",i,bp->r[0],bp->r[1],bp->r[2]);
    bp->s=sig;
    bp->a=0.0;
    bp->i=i;
  }
  *nC=nD;
  return bsc;
}

bCntr * ReadBasisCenters ( char * fn, double sig, int * nC ) {
  bCntr * bsc = NULL;
  bCntr * bp;
  char scrln[255];
  FILE * fp = fopen(fn,"r");
  int i=0;

  *nC=0;
  if (fp) {
    while (fgets(scrln,255,fp)) {
      if (scrln[0]!='#') (*nC)++;
    }
    rewind(fp);
    
    bsc = (bCntr*)malloc((*nC)*sizeof(bCntr));
    i=0;
    while (fgets(scrln,255,fp)) {
      if (scrln[0]!='#') {
	bp=&bsc[i];
	sscanf(scrln,"%lf %lf %lf",bp->r[0],bp->r[1],bp->r[2]);
	bp->s=sig;
	bp->a=0.0;
	bp->i=i;
	i++;
      }
    }
    fclose(fp);
  }
  fprintf(stderr,"# Read %i basis function centers from %s.\n",
	  *nC,fn);
  return bsc;
}


/* creates, populates, and returns the B-matrix:
 *
 *
 * B_{ik} = B_{ki} = \sum_{j=1}{n_D}\nabla_{\vec{z}_j}\phi_i({\vec{z}_j})\dot\nabla_{\vec{z}_j}\phi_k({\vec{z}_j})
 *
 */
double * B_matrix ( fSamp * data, int nD, bCntr * bsc, int nC ) {
  int i,j,k,ind;
  double * B = (double*)malloc(nC*nC*sizeof(double));
  double Zjk[data[0].ncv], Zji[data[0].ncv];
  double zjk, zji;
  double phi_jk, phi_ji;
  double gphi_jk, gphi_ji, dot;
  double ff=1.0;
  fSamp * jp;
  bCntr * ip, * kp;

  // populating element Bik
  for (i=0;i<nC;i++) {
    ip=&bsc[i];
    for (k=0;k<nC;k++) {
      kp=&bsc[k];
      ind=i*nC+k;
      B[ind]=0.0; // B[i][k]
      for (j=0;j<nD;j++) {
	jp=&data[j];
	vecdiff(Zjk,jp->r,kp->r,jp->ncv,jp->per,jp->doms);
	vecdiff(Zji,jp->r,ip->r,jp->ncv,jp->per,jp->doms);
	zjk=norm(Zjk,jp->ncv);
	zji=norm(Zji,jp->ncv);
	dot=(zjk>0&&zji>0)?vecdot(Zjk,Zji,jp->ncv)/zjk/zji:0.0;
	kernel(zjk,kp->s,&phi_jk,&gphi_jk);
	kernel(zji,ip->s,&phi_ji,&gphi_ji);
	//#if MODIFIED_SINGLE_SWEEP
	if (_modified_single_sweep_) {
	  ff=sqrt(norm(jp->f,jp->ncv))+EPSILON;
	//#else
	} else {
	  ff=1.0;
	}
	//#endif
	B[ind]+=dot*gphi_jk*gphi_ji/ff;
      }
    }
  }
  
  //  fprintf(stderr,"# constructed %i x %i B-matrix.\n",n,n);
  return B;
}

/* creates, populates, and returns the c-vector. */
double * c_vector ( fSamp * data, int nD, bCntr * bsc, int nC ) {
  int k,j;
  double * c = (double*)malloc(nC*sizeof(double));
  double Zjk[SDIM], zjk, dot;
  double nm;
  double phi;
  double gphi;
  double ff=1.0;
  fSamp * jp;
  bCntr * kp;

  for (k=0;k<nC;k++) {
    kp=&bsc[k];
    c[k]=0.0;
    for (j=0;j<nD;j++) {
      jp=&data[j];
      // to do -- use minimum image convention
      vecdiff(Zjk,jp->r,kp->r,SDIM);
      zjk=norm(Zjk,SDIM);
      dot=(zjk>0)?vecdot(jp->f,Zjk,SDIM)/zjk:0.0;
      kernel(zjk,kp->s,&phi,&gphi);
      //#if MODIFIED_SINGLE_SWEEP
      if (_modified_single_sweep_) {
	ff=sqrt(norm(jp->f,SDIM))+EPSILON;
      } else {
	//#else
	ff=1.0;
      }
      //#endif
      c[k]+=-dot*gphi/ff;
    }
  }
  return c;

}

/* a silly converter to produce a simple vector from a gsl-vector */
double * a_vector ( gsl_vector * A, int n ) {

  int i;
  double * a = malloc(n*sizeof(double));

  for (i=0;i<n;i++) {
    a[i]=gsl_vector_get(A,i);
  }
  return a;

}

/* Must reconstruct free energy at each sample site in order to compute
   the error in the forces. */
lRecon * ReconstructAtSamples ( fSamp * data, int nD, bCntr * bsc, int nC, double * relErr ) {
  int i,j;
  fSamp * ip;
  lRecon * recon, * irp;
  double f2s=0.0;
  double * B_data, * c_data, * a_data;
  gsl_matrix_view A, B;
  gsl_vector_view c;
  gsl_vector * a, * residual;

  B_data=B_matrix(data,nD,bsc,nC);
  A=gsl_matrix_view_array(B_data,nC,nC);
  B=gsl_matrix_view_array(B_data,nC,nC);
  c_data=c_vector(data,nD,bsc,nC);
  c=gsl_vector_view_array(c_data,nC);
  a=gsl_vector_alloc(nC);
  residual=gsl_vector_alloc(nC);

  /* Maragliano uses QR */
  gsl_linalg_QR_decomp(&B.matrix,residual);
  gsl_linalg_QR_solve(&B.matrix,residual,&c.vector,a);

  a_data=a_vector(a,nC);
  for (i=0;i<nC;i++) {
    bsc[i].a=a_data[i];
  }

  gsl_vector_free(a);

  gsl_vector_free(residual);
  free(B_data);
  free(c_data);
  free(a_data);

  recon=(lRecon*)malloc(nC*sizeof(lRecon));
  *relErr=0.0;
  for (i=0;i<nC;i++) {
    ip=&data[i];
    irp=&recon[i];
    irp->i=i;
    memcpy(irp->r,ip->r,SDIM*sizeof(double));
    localReconstruct(irp,bsc,nC);
    for (j=0;j<SDIM;j++) {
      *relErr+=(ip->f[j]-irp->f[j])*(ip->f[j]-irp->f[j]);
      f2s+=ip->f[j]*ip->f[j];
    }
  }

  *relErr/=f2s;

  return recon;

}

int generate_dxmap ( bCntr * bsc, int nC, double S, char * fn ) {
  int i,j,k,m,ii;
  FILE * fp;
  double thisz[SDIM],d[SDIM],nm;
  double mx[3]={-1.e99,-1.e99,-1.e99},mn[3]={1.e99,1.e99,1.e99};
  double sp[SDIM];
  int isp[SDIM];
  lRecon gridpnt, *gp=&gridpnt;
  double thisA,*A,Amin;
  int counter;
  int nGrid=0;
  bCntr * mp;
  double delta = 1.0;
  double sig;

  sig=0.0;
  for (i=0;i<nC;i++) {
    mp=&bsc[i];
    for (j=0;j<SDIM;j++) {
      if (mp->r[j]>mx[j]) mx[j]=mp->r[j];
      if (mp->r[j]<mn[j]) mn[j]=mp->r[j];
      if (mp->s>sig) sig=mp->s;
    }
  }

  mn[0]-=S*sig;
  mn[1]-=S*sig;
  mn[2]-=S*sig;
  mx[0]+=S*sig;
  mx[1]+=S*sig;
  mx[2]+=S*sig;
  sp[0]=mx[0]-mn[0];
  sp[1]=mx[1]-mn[1];
  sp[2]=mx[2]-mn[2];
  isp[0]=(int)(sp[0]/delta)+1;
  isp[1]=(int)(sp[1]/delta)+1;
  isp[2]=(int)(sp[2]/delta)+1;
  nGrid=isp[0]*isp[1]*isp[2];

  fprintf(stdout,"# map generation on a %i x %i x %i grid (%u bytes)\n",isp[0],isp[1],isp[2],nGrid*sizeof(double));

  A=malloc(nGrid*sizeof(double));

  /* data output here */
  ii=0;
  Amin=1.e99;
  for (i=0;i<isp[0];i++) {
    gp->r[0]=(double)(i+mn[0]);
    for (j=0;j<isp[1];j++) {
      gp->r[1]=(double)(j+mn[1]);
      for (k=0;k<isp[2];k++) {
	gp->r[2]=(double)(k+mn[2]);
	localReconstruct(gp,bsc,nC);
	A[ii++]=gp->e;
	if (thisA<Amin) Amin=thisA;
      }
    }
  }

  fp=fopen(fn,"w");

  /* *.dx file format for VMD */  
  fprintf(fp,"# Data calculated by Reconstruct.c (cfa) 2010\n");
  fprintf(fp,"object 1 class gridpositions counts %i %i %i\n",isp[0],isp[1],isp[2]);
  fprintf(fp,"origin %.1lf %.1lf %.1lf\n",mn[0],mn[1],mn[2]);
  fprintf(fp,"delta %.1lf 0 0\n",delta);
  fprintf(fp,"delta 0 %.1lf 0\n",delta);
  fprintf(fp,"delta 0 0 %.1lf\n",delta);
  fprintf(fp,"object 2 class gridconnections counts %i %i %i\n",isp[0],isp[1],isp[2]);
  fprintf(fp,"object 3 class array type double rank 0 items %i data follows\n",nGrid);

  counter=0;
  for (i=0;i<nGrid;i++) {
    //    fprintf(fp,"%.5lf ",A[i]-Amin);
    fprintf(fp,"%.5lf ",A[i]);
    counter++;
    if (counter==3) {
      fprintf(fp,"\n");
      counter=0;
    }
  }

  fprintf(fp,"\n");
  fprintf(fp,"object \"free energy reconstruction (kcal/mol)\" class field\n");
  close(fp);
  free(A);


  return 0;

}

bCntr * DisplaceBasisCenters_All ( bCntr * bsc, int nC, double dR, double dS, gsl_rng * r ) {
  bCntr * b = malloc(nC*sizeof(bCntr));
  int i;
  double dx;
  double dy;
  double dz;
  double ds;

  memcpy(b,bsc,nC*sizeof(bCntr));
  
  for (i=0;i<nC;i++) {
    dx = dR*(0.5-gsl_rng_uniform(r));
    dy = dR*(0.5-gsl_rng_uniform(r));
    dz = dR*(0.5-gsl_rng_uniform(r));
    ds = dS*(0.5-gsl_rng_uniform(r));
/*     fprintf(stderr,"#MC: displacing center %i by (% .3lf,% .3lf,% .3lf)\n", */
/* 	    i,dx,dy,dz); */
/*     fflush(stderr); */
    b[i].r[0]+=dx;
    b[i].r[1]+=dy;
    b[i].r[2]+=dz;
    b[i].s+=ds;
  }

  return b;
}

bCntr * DisplaceBasisCenters_Random ( bCntr * bsc, int nC, double dR, double dS, gsl_rng * r ) {
  bCntr * b = malloc(nC*sizeof(bCntr));

  int i = (int)(nC*(double)gsl_rng_uniform(r));
  double dx = dR*(0.5-gsl_rng_uniform(r));
  double dy = dR*(0.5-gsl_rng_uniform(r));
  double dz = dR*(0.5-gsl_rng_uniform(r));
  double ds = dS*(0.5-gsl_rng_uniform(r));

  memcpy(b,bsc,nC*sizeof(bCntr));

  fflush(stderr);
  b[i].r[0]+=dx;
  b[i].r[1]+=dy;
  b[i].r[2]+=dz;
  b[i].s+=ds;
  fprintf(stderr,"#MC: displacing center %i by (% .3lf,% .3lf,% .3lf; s % .4lf)\n",
	  i,dx,dy,dz,b[i].s);
  return b;
}

int main ( int argc, char * argv[] ) {
  fSamp * data; // force samples from input file
  lRecon * recon; // reconstructed forces
  bCntr * bsc; // basis centers
  int D; // number of CV's, aka, dimensionality of CV space
  int nD; // number of force samples
  int nC; // number of basis centers (need not equal nD, but usually does)
  int nR; // number of reconstructed force samples
  int i;
  double sig=2.5;
  double sig1,dsig=0.0;

  double Error;

  int s;
  int map=0;
  double S=2.0;
  int quiet=0;
  int D = -1;
  double fmag_lowthresh=0.0;
  double fmag_hithresh=0.0;

  char * cvinp_fn = "cv.inp";
  char * fsinp_fn = "forces.dat";
  char * bsout_fn = "basis_set.out";
  char * rcout_fn = "reconstruction.out";
  char * dxout_fn = "A_map.dx";

  /* EXPERIMENTAL:  Monte-carlo optimization */
  mct * mc = NULL;

  /* Functionality for input of basis centers that may be 
     different from the centers at which forces were measured. */
  char * bsin_fn = NULL;
  int inputBasisCenters = 0;

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-sig")) sig=atof(argv[++i]);
    else if (!strcmp(argv[i],"-sigscan")) 
      sscanf(argv[++i],"%lf,%lf,%lf",&sig,&sig1,&dsig);
    else if (!strcmp(argv[i],"-map")) map=1;
    else if (!strcmp(argv[i],"-q")) quiet=1;
    else if (!strcmp(argv[i],"-f_lowthresh")) fmag_lowthresh=atof(argv[++i]);
    else if (!strcmp(argv[i],"-f_hithresh")) fmag_hithresh=atof(argv[++i]);
    else if (!strcmp(argv[i],"-map_margin_factor")) S=atof(argv[++i]);
    else if (!strcmp(argv[i],"-cv")) cvinp_fn=argv[++i];
    else if (!strcmp(argv[i],"-fr")) fsinp_fn=argv[++i];
    else if (!strcmp(argv[i],"-bs")) bsout_fn=argv[++i];
    else if (!strcmp(argv[i],"-rc")) rcout_fn=argv[++i];
    else if (!strcmp(argv[i],"-dx")) dxout_fn=argv[++i];
    else if (!strcmp(argv[i],"-mod")) _modified_single_sweep_=1;
    else if (!strcmp(argv[i],"-i")) {
      inputBasisCenters=1;
      bsin_fn=argv[++i];
    }
    /* EXPERIMENTAL:  Monte-carlo optimization */
    else if (!strcmp(argv[i],"-mc")) {
      mc = malloc(sizeof(mct));
      sscanf(argv[++i],"%i,%lu,%lf,%lf,%s",&mc->nc,&mc->Seed,&mc->dR,&mc->dS,mc->flag);
    }
  }

  fprintf(stdout,"# Reconstruct.c begins.\n");

  D=get_dimensionality(cvinp_fn);
  fprintf(stdout,"# Dimensionality of CV space from %s: %i\n",cvinp_fn,D);

  /* Read data from the formatted forces.dat file from gather.tcsh */
  if (data=read_samples(&nD,fsinp_fn,D,quiet)) {
    fprintf(stderr,"ERROR reading force samples from %s.\n",fn);
    exit(-1);
  }
  
  if (inputBasisCenters) {
    bsc=readBasisSet(&nC,bsin_fn);
  } else {
    bsc=SetupSampleCongruentBasisCenters(data,nD,sig,&nC);
  }
  recon=ReconstructAtSamples(data,nD,bsc,nC,&Error);
  nR=nD;

  if (!quiet) fprintf(stdout,"# sigma  rel-resid-error (kcal/mol)\n");
  fprintf(stdout,"%.8lf %.8lf\n",sig,sqrt(Error));

  if (dsig>0) {
    sig+=dsig;
    for (;sig<=sig1;sig+=dsig) {
      free(bsc);
      bsc=SetupSampleCongruentBasisCenters(data,nD,sig,&nC);
      free(recon);
      recon=ReconstructAtSamples(data,nD,bsc,nC,&Error);
      fprintf(stdout,"%.8lf %.8lf\n",sig,sqrt(Error));fflush(stdout);
    }
  }

  if (mc) {
    int c;
    lRecon * oldrecon, * newrecon;
    bCntr * oldbsc, * newbsc;
    double newError;
    gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

    gsl_rng_set(r,mc->Seed);

    fprintf(stderr,"Entering MC loop; %i cycles; seed %u; dR %.5lf; dS %.5lf\n",mc->nc,mc->Seed,mc->dR,mc->dS);
    fflush(stderr);
    for (c=0;c<mc->nc;c++) {
      if (!strcmp(mc->flag,"all")) {
	newbsc=DisplaceBasisCenters_All(bsc,nC,mc->dR,mc->dS,r);
      } else {
	newbsc=DisplaceBasisCenters_Random(bsc,nC,mc->dR,mc->dS,r);
      }
      newrecon=ReconstructAtSamples(data,nD,newbsc,nC,&newError);
      if (newError<Error) {
	oldrecon=recon;
	recon=newrecon;
	free(recon);
	oldbsc=bsc;
	bsc=newbsc;
	free(oldbsc);
	Error=newError;
      } else {
	free(newrecon);
	free(newbsc);
      }
      fprintf(stderr,"#MC %i: %.5lf (%.5lf)\n",c,sqrt(Error),sqrt(newError));
    }
  }

  writeBasisSet(bsc,nC,bso_fn);
  writeReconstruction(recon,nR,rcn_fn);

  if (map) {
    generate_dxmap(bsc,nC,s,map_fn);
  }

  if (!quiet) fprintf(stdout,"# program ends.\n");
}
