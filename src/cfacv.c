#include "cfacv.h"

//#define _PARANOIA_ 1

void cfacvBanner ( void ) {
  printf("CFACV/C) Version %.2f\n",CFACVC_VERSION);
}

FILE * my_fopen ( char * name, char * code ) {
  if (strcmp(name,"stdout")) return fopen(name,code);
  else return stdout;
}

FILE * my_binfopen ( char * name, char * code, unsigned int outputLevel, DataSpace * ds ) {
  FILE * fp=fopen(name,code);
  fwrite(&outputLevel,sizeof(unsigned int),1,fp);
  fwrite(&(ds->iK),sizeof(int),1,fp);
  return fp;
}

char * CVSTRINGS[NULL_CV] = {"BOND", "ANGLE", "DIHED", "CARTESIAN_X", "CARTESIAN_Y", "CARTESIAN_Z"};
int cv_getityp ( char * typ ) {
  int i;
  for (i=0;i<NULL_CV&&strcmp(typ,CVSTRINGS[i]);i++);
  if (i<NULL_CV) return i;
  else return -1;
}

char * cv_getstyp ( int ityp ) {
  if (ityp<NULL_CV) return CVSTRINGS[ityp];
  else return "NOT_FOUND";
}

char * RFSTRINGS[NULL_RF] = {"HARMONIC", "PERIODIC", "VORONOI-CENTER", "SPHERE-CENTER"};
int rf_getityp ( char * typ ) {
  int i;
  for (i=0;i<NULL_RF&&strcmp(typ,RFSTRINGS[i]);i++);
  if (i<NULL_RF) return i;
  else return -1;
}

char * BADRESSTRINGMESSAGE = "NOT FOUND";
char * rf_getstyp ( int ityp ) {
  if (ityp>-1&&ityp<NULL_RF) return RFSTRINGS[ityp];
  else return BADRESSTRINGMESSAGE;
}

cvStruct * New_cvStruct ( int typ, int nC, int * ind ) {
  cvStruct * newcv=malloc(sizeof(cvStruct));
  int i;
  newcv->typ=typ;
  newcv->nC=nC;
  newcv->val=0.0;
  newcv->ind=calloc(nC,sizeof(int));
  if (ind) for (i=0;i<nC;i++) newcv->ind[i]=ind[i];
  newcv->gr=(double**)malloc(nC*sizeof(double*));
  for (i=0;i<nC;i++) newcv->gr[i]=(double*)malloc(3*sizeof(double));
  return newcv;
}

metricTensorStruct * New_metricTensorStruct ( int M ) {
  int i;
  metricTensorStruct * mt = (metricTensorStruct*)malloc(sizeof(metricTensorStruct));
  mt->MM=(double*)malloc(M*M*sizeof(double*));
  mt->n=0;
  mt->m=M;
  /* sparse handler */
  mt->nn=0;
  mt->cva=mt->cvb=NULL;
  mt->nca=NULL;
  mt->ca=NULL;
  return mt;
}

int DataSpace_metricTensor_Reset ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    int m=ds->mt->m;
    for (i=0;i<m;i++) for (j=0;j<m;j++)
	ds->mt->MM[i*m+j]=0.0;
    ds->mt->n=0;
  }
}

int DataSpace_forceAccumulators_Reset ( DataSpace * ds ) {
  if (ds) {
    int i;
    restrStruct * r;
    forceAccumStruct * facc;
    for (i=0;i<ds->K;i++) {
      r=ds->restr[i];
      facc=r->facc;
      facc->f=0.0;
      facc->n=0;
    }
  }
}

void DataSpace_metricTensor_fprintf ( DataSpace * ds, FILE * fp ) {
  if (ds) {
    int k,l;
    metricTensorStruct * mt = ds->mt;
    fprintf(fp,"CFACV/C) Metric tensor: %i non-zero elements (%i off-diagonals)\n",mt->nn,mt->nn-ds->M);
    fprintf(fp,"         [linear-index]  [cva  cvb]  (value (1/amu)) (#-common-centers) : list-of-center-pair-indices [global-ctr-index:(cva-i,cvb-i)]\n");
    for (k=0;k<mt->nn;k++) {
      fprintf(fp,"    [%i]  [%i %i] : (%.5le) : (%i) : ",k,mt->cva[k],mt->cvb[k],mt->nca[k],mt->MM[mt->cva[k]*mt->m+mt->cvb[k]]);
      for (l=0;l<mt->nca[k];l++) fprintf(fp,"[%i:(%i,%i)]",mt->ca[k][l].a,mt->ca[k][l].i,mt->ca[k][l].j);
      fprintf(fp,"\n");
    }
    fflush(fp);
  }
}

int DataSpace_metricTensor_Setup ( DataSpace * ds ) {
  if (ds) {
    // visit each pair of CV's. If any pair of CV's has one or more
    // element in common in the ind[] array, then that pair of CV's
    // may have a non-zero metric tensor element.
    cvStruct * cva, * cvb;
    int i,j,hit,ii,jj,k,l;
    metricTensorStruct * mt = ds->mt;
    
    if (mt) {
      /* count the number of elements that are not identically zero.
	 These are (1) diagonal elements, and (2) elements
	 corresponding to pairs of CV's that share dependence on one or
	 more configurational variables.  Condition 2 is met 
	 if two CV's have common centers in
	 their lists of centers AND, if both are Cartesian, the same direction */
      fprintf(stdout,"CFACV/C) Setting up metric tensor structure...\n");fflush(stdout);
      DataSpace_metricTensor_Reset(ds);
      mt->nn=0;
      for (i=0;i<ds->M;i++) {
	cva=ds->cv[i];
	mt->nn++; // diagonal elements assumed to be non-zero
	for (j=i+1;j<ds->M;j++) {
	  cvb=ds->cv[j];
	  /* for this pair of CV's, search through both list in center
	     indices for common elements */
	  hit=0; // assume this pair of CV's does not have a common center
	  for (ii=0;ii<cva->nC&&!hit;ii++) {
	    for (jj=0;jj<cvb->nC&&!hit;jj++) {
	      if (cva->ind[ii]==cvb->ind[jj]) {
		// these two CV's share a common center, but this does not mean
		// they necessarily share dependence on one or more configurational
		// variable.  E.g., if both are CARTESIAN, they cannot.
		if (!((cva->typ==CARTESIAN_X||cva->typ==CARTESIAN_Y||cva->typ==CARTESIAN_Z)&&(cvb->typ==CARTESIAN_X||cvb->typ==CARTESIAN_Y||cvb->typ==CARTESIAN_Z))) {
		  hit=1;
		}
	      }
	    }
	  }
	  if (hit) mt->nn++; // pair-ij of CV's is a hit, so increment the counter
	}
      }
      fprintf(stdout,"CFACV/C) mt has %i non-zero elements\n",mt->nn);fflush(stdout);
      // allocate
      mt->cva=(int*)malloc(mt->nn*sizeof(int));
      mt->cvb=(int*)malloc(mt->nn*sizeof(int));
      mt->nca=(int*)malloc(mt->nn*sizeof(int));
      mt->ca=(mtmappair**)malloc(mt->nn*sizeof(mtmappair*));
      // populate
      k=0;
      for (i=0;i<ds->M;i++) {
	cva=ds->cv[i];
	/* diagonal element */
	mt->cva[k]=i;
	mt->cvb[k]=i;
	mt->nca[k]=cva->nC;
	mt->ca[k]=(mtmappair*)malloc(cva->nC*sizeof(mtmappair));
	for (l=0;l<cva->nC;l++) {
	  mt->ca[k][l].a=cva->ind[l];
	  mt->ca[k][l].i=l;
	  mt->ca[k][l].j=l;
	}
	k++;
	/* search (again) for partner CV's */
	for (j=i+1;j<ds->M;j++) {
	  hit=0;
	  cvb=ds->cv[j];
	  for (ii=0;ii<cva->nC;ii++) {
	    for (jj=0;jj<cvb->nC;jj++) {
	      if (cva->ind[ii]==cvb->ind[jj]) {
		if (!((cva->typ==CARTESIAN_X||cva->typ==CARTESIAN_Y||cva->typ==CARTESIAN_Z)&&(cvb->typ==CARTESIAN_X||cvb->typ==CARTESIAN_Y||cvb->typ==CARTESIAN_Z))) {
		  hit++;
		}
	      }
	    }
	  }
	  /* 'hit' now counts the number of common centers between
	     these two CV's */
	  if (hit) {
	    mt->cva[k]=i;
	    mt->cvb[k]=j;
	    mt->nca[k]=hit;
	    mt->ca[k]=(mtmappair*)malloc(cva->nC*sizeof(mtmappair));
	    l=0;
	    for (ii=0;ii<cva->nC;ii++) {
	      for (jj=0;jj<cvb->nC;jj++) {
		if (cva->ind[ii]==cvb->ind[jj]) {
		  if (!((cva->typ==CARTESIAN_X||cva->typ==CARTESIAN_Y||cva->typ==CARTESIAN_Z)&&(cvb->typ==CARTESIAN_X||cvb->typ==CARTESIAN_Y||cvb->typ==CARTESIAN_Z))) { 
		    mt->ca[k][l].a=cva->ind[ii];
		    mt->ca[k][l].i=ii;
		    mt->ca[k][l].j=jj;
		    l++;
		  }
		}
	      }
	    }
	    mt->ca[k]=realloc(mt->ca[k],l*sizeof(mtmappair));
	    k++; // increment the index
	  }
	}
      }
    }
    fprintf(stdout,"CFACV/C) Metric tensor structure is set up:\n");
    DataSpace_metricTensor_fprintf(ds,stdout);
    return 0;
  }
  return -1;
}

double rf_Harmonic ( double k, double v, double z, double half_domain ) {
  return -k*(v-z);
}

double re_Harmonic ( double k, double v, double z, double half_domain ) {
  return 0.5*k*(v-z)*(v-z);
}

double rf_Periodic ( double k, double v, double z, double half_domain ) {
  double dz=v-z;
#ifdef _PARANOIA_
  if (_PARANOIA_) {
    if ((dz!=dz)||(k!=k)||(half_domain!=half_domain)) {
      fprintf(stderr,"CFACV/C/PARANOIA) Tripped in rf_Periodic\n");
      fprintf(stderr,"CFACV/C/PARANOIA) k %.5lf v %.5lf z %.5lf dz %.5lf half_domain %.5lf\n",
	      k,v,z,dz,half_domain);
      fprintf(stderr,"Program exits\n");
      fflush(stderr);
      exit(-1);
    }
  }
#endif
  if (dz<-half_domain) dz+=2*half_domain;
  else if (dz>half_domain) dz-=2*half_domain;
  
  return -k*dz;
}

double re_Periodic ( double k, double v, double z, double half_domain ) {
  double dz=v-z;
  if (dz<-half_domain) dz+=2*half_domain;
  else if (dz>half_domain) dz-=2*half_domain;
  return 0.5*k*dz*dz;
}

restrStruct * New_restrStruct ( double k, double z, int nCV, double * cvc, char * rftypstr, double zmin, double zmax ) {
  restrStruct * newr=malloc(sizeof(restrStruct));
  int i;
  newr->rfityp=rf_getityp(rftypstr);
  newr->forceFunc=NULL;
  newr->energyFunc=NULL;
  newr->k=k;
  newr->z=z;
  newr->u=newr->f=0.0;
  newr->nCV=nCV;
  newr->cvc=(double*)malloc(nCV*sizeof(double));
  if (cvc) for (i=0;i<nCV;i++) newr->cvc[i]=cvc[i];
  else for (i=0;i<nCV;i++) newr->cvc[i]=0;
  newr->tamdOpt=NULL;
  newr->tamd_noise=0.0;
  newr->tamd_restraint=0.0;
  newr->tmdOpt=NULL;
  newr->min=zmin;
  newr->max=zmax;
  newr->half_domain=0.5*(zmax-zmin);
  newr->sbc=NULL;
  newr->sbr=0.0;
  if (newr->rfityp==HARMONIC) {
    newr->forceFunc = rf_Harmonic;
    newr->energyFunc = re_Harmonic;
  }
  else if (newr->rfityp==PERIODIC) {
    newr->forceFunc = rf_Periodic;
    newr->energyFunc = re_Periodic;
  }

  newr->facc=New_forceAccumStruct();
  
  return newr;
}

forceAccumStruct * New_forceAccumStruct ( void ) {
  forceAccumStruct * facc=malloc(sizeof(forceAccumStruct));
  facc->f=0.0;
  facc->n=0;
  return facc;
}

tamdOptStruct * New_tamdOptStruct ( double g, double kt, double dt, int riftyp, double half_domain ) {
  tamdOptStruct * tamd=malloc(sizeof(tamdOptStruct));
  tamd->kT=kt;
  tamd->gamma=g;
  tamd->dt=dt;
  if (g>0.0) tamd->ginv=1.0/g;
  else tamd->ginv=0.0;
  tamd->noise = sqrt(2.0*kt*dt*tamd->ginv);
  tamd->periodic=(riftyp==PERIODIC);
  tamd->half_domain=half_domain;
  return tamd;
}

void forceAccumulate ( double f, forceAccumStruct * facc ) {
  if (facc) {
    facc->f+=f;
    facc->n++;
  }
}

double forceAvg ( forceAccumStruct * facc ) {
  if (facc) {
    return facc->n?(facc->f/facc->n):0.0;
  } else return 0.0;
}

double forceStartAccum ( double f, forceAccumStruct * facc ) {
  if (facc) {
    facc->f=f;
    facc->n=1;
  }
}

/* Diffusive dynamics! */
double tamdUpdate ( DataSpace * ds, double z, double f, tamdOptStruct * tamd, double * noise, double * res ) {
  double dd = tamd->ginv*tamd->dt*f;
  double rd = tamd->noise*my_whitenoise(ds->Xi);
  *noise=rd;
  *res=dd;
  if (tamd->periodic) {
    double newrawz=z+dd+rd;
    if (newrawz>tamd->half_domain) newrawz-=2*tamd->half_domain;
    if (newrawz<-tamd->half_domain) newrawz+=2*tamd->half_domain;
    return newrawz;
  }
  return z+dd+rd; 
}

double tmdUpdate ( double z, double increment, int OK ) {
  return z+OK*increment;
}

double tmdUpdate_Periodic ( double z, double increment, int OK ) {
  double newz=z+increment;
  if (newz>M_PI) newz-=2*OK*M_PI;
  else if (newz<-M_PI) newz+=2*OK*M_PI;
  return newz;
}

tmdOptStruct * New_tmdOptStruct ( double target, int t0, int t1, int periodic ) {
  tmdOptStruct * tmd=malloc(sizeof(tmdOptStruct));
  tmd->t0=t0;
  tmd->t1=t1;
  tmd->invinterval = 1.0/(t1-t0);
  tmd->target=target;
  if (periodic) tmd->update=tmdUpdate_Periodic;
  else tmd->update=tmdUpdate;
  return tmd;
}

int tmdOptInit ( tmdOptStruct * tmd, double initval, int periodic ) {
  double dist=0.0;
  if (tmd) {
    tmd->initval=initval;
    dist = tmd->target-tmd->initval;
    if (periodic) {
      if (dist<-M_PI) {
	dist+=2*M_PI;
      }
      else if (dist > M_PI) {
	dist-=2*M_PI;
      } 
    }
    tmd->increment=dist*tmd->invinterval;
    printf("CFACV/C) tmd increment is %.5lf based on target %.5lf and initval %.5lf\n",
	   tmd->increment,tmd->target,tmd->initval);
    return 0;
  }
  return -1;
}

DataSpace * NewDataSpace ( int N, int M, int K, long int seed ) {
  DataSpace * ds = malloc(sizeof(DataSpace));
  int i;
  ds->N=N;
  ds->M=M;
  ds->K=K;
  ds->iN=ds->iM=ds->iK=0;

  ds->Xi=(unsigned short *)malloc(3*sizeof(unsigned short));
  ds->Xi[0]=(seed & 0xffff00000000) >> 32;
  ds->Xi[1]=(seed & 0x0000ffff0000) >> 16;
  ds->Xi[2]= 0x330e;

  ds->R=(double**)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) ds->R[i]=calloc(3,sizeof(double));

  ds->ac = (atomCenterStruct**)malloc(N*sizeof(atomCenterStruct*));

  ds->cv=(cvStruct**)malloc(M*sizeof(cvStruct*));
  ds->restr=(restrStruct**)malloc(K*sizeof(restrStruct*));
  ds->z=(double*)malloc(K*sizeof(double));
  ds->oldz=(double*)malloc(K*sizeof(double));
  ds->f=(double*)malloc(K*sizeof(double));
  ds->MM=(double*)malloc(K*K*sizeof(double));

  /* metric tensor -- only use of K == M (one restraint per CV) */
  ds->mt=NULL;
  if (M>0&&M==K) {
    ds->mt=New_metricTensorStruct(M);
  }

  return ds;
}

int DataSpace_getN ( DataSpace * ds ) {
  if (ds) return ds->N;
  else return -1;
}

double * DataSpace_centerPos ( DataSpace * ds, int i ) {
  if (ds) {
    if (i<ds->N) {
      return ds->R[i];
    }
    return NULL;
  }
  return NULL;
}

double * DataSpace_z ( DataSpace * ds ) {
  if (ds) {
    return ds->z;
  }
  return NULL;
}

double * DataSpace_g ( DataSpace * ds ) {
  if (ds) {
    return ds->f; // "f" means force, but in the sense that this is on the atoms -- this is really the gradient dF/dz
  }
  return NULL;
}

double * DataSpace_MM ( DataSpace * ds ) {
  if (ds) {
    return ds->MM;
  }
  return NULL;
} 

int DataSpace_nz ( DataSpace * ds ) {
  if (ds) {
    return ds->K;
  }
  return 0;
}

int DataSpace_checkdata ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    for (i=0;i<ds->N;i++) {
      for (j=0;j<3;j++) {
	if (ds->R[i][j]!=ds->R[i][j]) {
	  return 1;
	}
      }
    }
    return 0;
  }
  return 0;
}

int DataSpace_dump ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    for (i=0;i<ds->N;i++) {
      fprintf(stderr,"%i ",i);
      for (j=0;j<3;j++)
	fprintf(stderr,"%.5lf ",ds->R[i][j]);
      fprintf(stderr,"\n");
      fflush(stderr);
    }
  }
  return 0;
}

int DataSpace_AddCenter_MassOnly ( DataSpace * ds, double m ) {
  if (ds) {
    if (ds->iN<ds->N) {
      int i,j;
      ds->ac[ds->iN++]=New_atomCenterStruct(1);
      i=ds->iN-1;
      ds->ac[i]->ind[0]=-1;
      ds->ac[i]->m[0]=m;
      ds->ac[i]->M=m;
    }
  }
}

int DataSpace_AddAtomCenter ( DataSpace * ds, int n, int * ind, double * m ) {
  if (ds) {
    if (ds->iN<ds->N) {
      int i,j;
      ds->ac[ds->iN++]=New_atomCenterStruct(n);
      i=ds->iN-1;
      for (j=0;j<n;j++) {
	ds->ac[i]->ind[j]=ind[j];
	ds->ac[i]->m[j]=m[j];
	ds->ac[i]->M+=m[j];
      }
      return (ds->iN-1);
    }
    return -1;
  }
  return -1;
}

int DataSpace_AddCV ( DataSpace * ds, char * typ, int nind, int * ind ) {
  if (ds) {
    int ityp=cv_getityp(typ);
    if (ds->iM<ds->M) {
      ds->cv[ds->iM++]=New_cvStruct(ityp,nind,ind);
      return (ds->iM-1);
    }
    return -1;
  }
  return -1;
}

int DataSpace_AddRestr ( DataSpace * ds, double k, double z, int nCV, double * cvc, char * rftypstr, 
			 double zmin, double zmax ) {
  if (ds) {
    if (ds->iK<ds->K) {
      ds->restr[ds->iK++]=New_restrStruct(k,z,nCV,cvc,rftypstr,zmin,zmax);
      return (ds->iK-1);
    }
    return -1;
  }
  return -1;
}


int DataSpace_updateRestraintValue ( DataSpace * ds, int i, double rv ) {
  if (ds) {
    ds->restr[i]->z=rv;
    ds->z[i]=rv;
    ds->restr[i]->facc->f=0.0;
    ds->restr[i]->facc->n=0;
    return 0;
  }
  return -1;
}


int DataSpace_AddSphericalBoundaryToRestr ( DataSpace * ds, int ir, double * ctr, double rad ) {
  if (ds) {
    restrStruct * r = ds->restr[ir];
    int i;
    if (r&&!r->sbc) {
      r->sbc=(double*)malloc(r->nCV*sizeof(double));
      for (i=0;i<r->nCV;i++) {
	r->sbc[i]=ctr[i];
      }
      r->sbr=rad;
    }
  }
  return -1;
}

int DataSpace_AddForceAccum ( DataSpace * ds, int ir ) {
  if (ds) {
    if (ir < ds->iK ) {
      ds->restr[ir]->facc=New_forceAccumStruct();
      return ir;
    }
    return -1;    
  }
  return -1;
}

int DataSpace_AddTamdOpt ( DataSpace * ds, int ir, double g, double kt, double dt ) {
  if (ds) {
/*     printf("DB: adding tamd options (% 7.3f% 7.3f% 7.3f) to r %i out of %i...\n", */
/* 	   g,kt,dt,ir,ds->iK); */
    if (ir < ds->iK ) {
      ds->restr[ir]->tamdOpt=New_tamdOptStruct(g,kt,dt,ds->restr[ir]->rfityp,ds->restr[ir]->half_domain);
      return ir;
    }
    return -1;
  }
  return -1;
}

int DataSpace_AddTmdOpt  ( DataSpace * ds, int ir, double target, int t0, int t1 ) {
  if (ds) {
    if (ir < ds->iK ) {
      ds->restr[ir]->tmdOpt=New_tmdOptStruct(target,t0,t1,(int)(ds->restr[ir]->forceFunc==rf_Periodic));
      return ir;
    }
    return -1;
  }
  return -1;
}

int DataSpace_SetRestraints ( DataSpace * ds, double * rval ) {
  int i;
  restrStruct * r;
  for (i=0;i<ds->iK;i++) {
    r=ds->restr[i];
    r->z=rval[i];
    // printf("CFACV/C) Reinitialize Z[%i] = %.5lf\n",i,r->z);
    ds->z[i]=r->z;
  }
  return 0;
}

/* compute instantaneous values of CV's from center positions */
int DataSpace_ComputeCVs ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    cvStruct * cvi, * cvj;

    for (i=0;i<ds->iM;i++) {
      cvi=ds->cv[i];
      if (cvi->typ==BOND) {
	cvi->val=my_getbond(ds->R[cvi->ind[0]],ds->R[cvi->ind[1]],cvi->gr[0],cvi->gr[1]);
      } 
      else if (cvi->typ==ANGLE) {
	cvi->val=my_getangle(ds->R[cvi->ind[0]],ds->R[cvi->ind[1]],ds->R[cvi->ind[2]],
			     cvi->gr[0],        cvi->gr[1],        cvi->gr[2]);
      }
      else if (cvi->typ==DIHED) {
	cvi->val=my_getdihed(ds->R[cvi->ind[0]],ds->R[cvi->ind[1]],ds->R[cvi->ind[2]],ds->R[cvi->ind[3]],
			     cvi->gr[0],        cvi->gr[1],        cvi->gr[2],        cvi->gr[3]);
#ifdef _PARANOIA_
	if (_PARANOIA_) {
	  if (cvi->val!=cvi->val) {
	    fprintf(stderr,"CFACV/C/PARANOIA) Tripped at dihed cvi->val %.5lf\n",cvi->val);
	    fprintf(stderr,"Program exits.\n");
	    fflush(stderr);
	    exit(-1);
	  }
	}
#endif
      }
      else if (cvi->typ==CARTESIAN_X) {
	cvi->val=ds->R[cvi->ind[0]][0];
	cvi->gr[0][0]=1.0;
	cvi->gr[0][1]=0.0;
	cvi->gr[0][2]=0.0;
      }
      else if (cvi->typ==CARTESIAN_Y) {
	cvi->val=ds->R[cvi->ind[0]][1];
	cvi->gr[0][0]=0.0;
	cvi->gr[0][1]=1.0;
	cvi->gr[0][2]=0.0;
      }
      else if (cvi->typ==CARTESIAN_Z) {
	cvi->val=ds->R[cvi->ind[0]][2];
	cvi->gr[0][0]=0.0;
	cvi->gr[0][1]=0.0;
	cvi->gr[0][2]=1.0;
      }
    }

    // update metric tensor tally
    if (ds->mt) {
//      fprintf(stdout,"CFACV/C) metric tensor [%x] update\n",ds->mt);fflush(stdout);
      metricTensorStruct * mt = ds->mt;
      int k,d;
      int ii,jj;
      double incr, thismte;
/* #ifdef SHORTMT */
      for (i=0;i<mt->nn;i++) {
	cvi=ds->cv[mt->cva[i]];
	cvj=ds->cv[mt->cvb[i]];
	// for each common center, tally 1/mass * ( grad(cva,rx)*grad(cvb,rx) + grad(cva,ry)*grad(cvb,ry) + grad(cva,rz)*grad(cvb,rz) )
        thismte=0.0;
	for (k=0;k<mt->nca[i];k++) {
          incr=0.0; 
	  for (d=0;d<3;d++) {
            incr+=1.0/ds->ac[mt->ca[i][k].a]->M * ( cvi->gr[mt->ca[i][k].i][d] * cvj->gr[mt->ca[i][k].j][d] );
          }
          thismte+=incr;
	  mt->MM[mt->cva[i]*mt->m+mt->cvb[i]] += incr;
	  if (mt->cvb[i]!=mt->cva[i]) mt->MM[mt->cvb[i]*mt->m+mt->cva[i]] += incr;
#ifdef PARANOIA
          fprintf(stdout,"CFACV/C/DEBUG: mt n.z.e. %i common-atom %i (%i in cv-%i, %i in cv-%i) out of %i\n",
                   i,k, mt->ca[i][k].i, mt->cva[i], mt->ca[i][k].j, mt->cvb[i], mt->nca[i]);fflush(stdout);
          fprintf(stdout,"              (%i,%i): mass %.5lf\n",
                   mt->cva[i],mt->cvb[i],ds->ac[mt->ca[i][k].a]->M);fflush(stdout);
          fprintf(stdout,"              grad-%i %.5lf %.5lf %.5lf grad-%i %.5lf %.5lf %.5lf\n",
                   mt->cva[i],cvi->gr[mt->ca[i][k].i][0],cvi->gr[mt->ca[i][k].i][1],cvi->gr[mt->ca[i][k].i][2],
                   mt->cvb[i],cvj->gr[mt->ca[i][k].j][0],cvj->gr[mt->ca[i][k].j][1],cvj->gr[mt->ca[i][k].j][2]);fflush(stdout);
          fprintf(stdout,"              metric tensor increment %.5lf\n",
                   1.0/ds->ac[mt->ca[i][k].a]->M*(cvi->gr[mt->ca[i][k].i][0]*cvj->gr[mt->ca[i][k].j][0]+
                   cvi->gr[mt->ca[i][k].i][1]*cvj->gr[mt->ca[i][k].j][1]+
                   cvi->gr[mt->ca[i][k].i][2]*cvj->gr[mt->ca[i][k].j][2]));fflush(stdout);
#endif
        }
#ifdef PARANOIA
        fprintf(stdout,"CFACV/C/DEBUG: mt n.z.e. %i thiselementvalue %.5lf\n",i,thismte);fflush(stdout);
#endif
      }
/* #else */
//       for (i=0;i<ds->M;i++) {
// 	for (j=0;j<ds->M;j++) {
// 	  for (ii=0;ii<ds->cv[i]->nC;ii++) {
// 	    for (jj=0;jj<ds->cv[j]->nC;jj++) {
// 	      if (ds->cv[i]->ind[ii]==ds->cv[j]->ind[jj]) {
// 		for (d=0;d<3;d++) {
// 		  mt->MM[i*mt->m+j]+=1.0/ds->ac[ds->cv[i]->ind[ii]]->M * (ds->cv[i]->gr[ii][d]*ds->cv[j]->gr[jj][d]);
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//      }
/* #endif */
      // this is one full time-step of update, so..
      mt->n++;
      //DataSpace_metricTensor_fprintf(ds,stdout);
    }

    return 0;
  }
  return -1;
}

int DataSpace_Tally ( DataSpace * ds ) {
//  fprintf(stderr,"CFACV/C) SM entered DataSpace_Tally()...\n");fflush(stderr);
  if (ds) {
  //  fprintf(stderr,"CFACV/C) SM tallying gradients and MT dim %d match %d...\n",ds->K,ds->mt->m);fflush(stderr);
    int i,j;
    for (i=0;i<ds->K;i++) {
      ds->f[i]=ds->restr[i]->facc->f/ds->restr[i]->facc->n;
    //  fprintf(stderr,"  -> f[%d] %.6lf\n",i,ds->f[i]);
      if (ds->mt) for (j=0;j<ds->K;j++) {
        ds->MM[i*ds->mt->m+j] = ds->mt->MM[i*ds->mt->m+j]/ds->mt->n;
      //  fprintf(stderr,"     -> M[%d][%d] %.6lf\n",i,j,ds->MM[i*ds->mt->m+j]);fflush(stderr);
      }
    }  
  }
}

int DataSpace_AssignZsFromCVs ( DataSpace * ds ) {
  if (ds) {
    int i,j;
    restrStruct * r;
    cvStruct * cv;
    double * cvc;
    /* Loop over all restraints */
    for (i=0;i<ds->iK;i++) {
      r=ds->restr[i];
      cvc=r->cvc;
      r->val=0.0;
      for (j=0;j<r->nCV;j++) {
	cv=ds->cv[j];
	r->val+=cvc[j]*cv->val;
      }
      r->z=r->val;
      ds->z[i]=r->z;
    }
  }
  return 0;
}

int DataSpace_RestrainingForces ( DataSpace * ds, int first, int timestep ) {
  if (ds) {
    int i,j,jj,k;
    double v, u, f, kr, z, dz;
    cvStruct * cv;
    restrStruct * r;
    double * cvc;

    //   printf("CFACV/C) RestrainingForces begins.\n");

    /* clear out the position arrays to hold forces to be communicated
       back to tclForces */
    for (i=0;i<ds->N;i++) {
      ds->R[i][0]=0.0;
      ds->R[i][1]=0.0;
      ds->R[i][2]=0.0;
    }


    /* Loop over all restraints */
    for (i=0;i<ds->iK;i++) {
      r=ds->restr[i];
      cvc=r->cvc;
      r->val=0.0;
      for (j=0;j<r->nCV;j++) {
	cv=ds->cv[j];
	r->val+=cvc[j]*cv->val;
#ifdef _PARANOIA_
	if (_PARANOIA_) {
	  if (r->val!=r->val) {
	    fprintf(stderr,"CFACV/C/PARANOIA) Tripped at r->val %.5lf\n",r->val);
	    fprintf(stderr,"CFACV/C/PARANOIA) cvc[%i] %.5lf cv->val %.5lf\n",j,cv->val);
	    fprintf(stderr,"Program exits.\n");
	    fflush(stderr);
	    exit(-1);
	  }
	}
#endif
      }
      
      if (first) {
	if (r->tamdOpt) {
	  r->z=r->val;
	}
	if (r->tmdOpt) {
	  r->z=r->val;
	  tmdOptInit(r->tmdOpt,r->val,(int)(r->forceFunc==rf_Periodic));	  
	}
      }
      //printf("CFACV/C) %i restraint %i val %.4lf targ %.4lf\n",timestep,i,r->val,r->z);
      r->f=r->forceFunc(r->k,r->val,r->z,r->half_domain);
      forceAccumulate(r->f,r->facc);
      
      r->u=r->energyFunc(r->k,r->val,r->z,r->half_domain);

      //printf("CFACV/C) restraint %i: f %.5lf u %.5lf\n",i,r->f,r->u);


      /* accumulate force increments */
      /* cv->gr[jj][kk] = \partial CV / \partial (kk-coord of center jj of CV) */
      for (j=0;j<r->nCV;j++) {
	cv=ds->cv[j];
	/* if the j'th cv referenced by this restraint contributes to this restraint: */
	if (r->cvc[j]) {
	  /* ... for each center in this CV ... */
	  for (jj=0;jj<cv->nC;jj++) {
	    /* increment the forces of the jj'th cartesian center on the j'th cv of this restraint */
	    ds->R[cv->ind[jj]][0]+=r->cvc[j]*cv->gr[jj][0]*r->f;
	    ds->R[cv->ind[jj]][1]+=r->cvc[j]*cv->gr[jj][1]*r->f;
	    ds->R[cv->ind[jj]][2]+=r->cvc[j]*cv->gr[jj][2]*r->f;
#ifdef _PARANOIA_
	    if (_PARANOIA_) {
	      if ((ds->R[cv->ind[jj]][0]!=ds->R[cv->ind[jj]][0])
		  ||(ds->R[cv->ind[jj]][1]!=ds->R[cv->ind[jj]][1])
		  ||(ds->R[cv->ind[jj]][2]!=ds->R[cv->ind[jj]][2])) {
		fprintf(stderr,"CFACV/C/PARANOIA) Tripped.\n");
		fprintf(stderr,"CFACV/C/PARANOIA) r->cvc[%i] %.5lf\n",j,r->cvc[j]);
		fprintf(stderr,"CFACV/C/PARANOIA) cv->gr[%i][0] %.5lf\n",jj,cv->gr[jj][0]);
		fprintf(stderr,"CFACV/C/PARANOIA) cv->gr[%i][1] %.5lf\n",jj,cv->gr[jj][1]);
		fprintf(stderr,"CFACV/C/PARANOIA) cv->gr[%i][2] %.5lf\n",jj,cv->gr[jj][2]);
		fprintf(stderr,"CFACV/C/PARANOIA) r->f %.5lf\n",r->f);
		fprintf(stderr,"Program exits\n");
		fflush(stderr);
		exit(-1);
	      }
	    }
#endif
	  }
	}
      }
    }

    

    // update of z's
    for (i=0;i<ds->iK;i++) {
      r=ds->restr[i];
      if (r->tamdOpt) {
	r->z=tamdUpdate(ds,r->z,-r->f,r->tamdOpt,&r->tamd_noise,&r->tamd_restraint);
      }
      if (r->tmdOpt) {
	r->z=r->tmdOpt->update(r->z,r->tmdOpt->increment,
			       (r->tmdOpt->t0<=timestep)&&(r->tmdOpt->t1>=timestep));
      }
      ds->z[i]=r->z;
    }
    return 0;
  }
  return -1;
}


double DataSpace_RestraintEnergy ( DataSpace * ds ) {
  if (ds) {
    double u=0.0;
    int i;
    restrStruct * r;
    for (i=0;i<ds->iK;i++) {
      r=ds->restr[i];
      u+=r->u;
    }
    return u;
  }
  return 0.0;    
}

int DataSpace_VoronoiCellCheck ( DataSpace * ds, int home_cell ) {
  if (ds) {
    int i,j,jj,k;
    double v, u, f, kr, z, dz;
    double cvpt[ds->M];
    restrStruct * r;
    cvStruct * cv;
    double * this_vc;
    double min_d2=1.e69;
    int min_i;
    double this_d, this_d2;
    //printf("CFACV/C) whichCell begins.\n");

    // printf("CFACV/C) cv-pt at ");
    for (i=0;i<ds->M;i++) {
      cv=ds->cv[i];
      cvpt[i]=cv->val;
      // printf("%.5lf ",cvpt[i]);
    }
 //  printf("\n");

    /* Loop over all "restraints" */
    for (i=0;i<ds->iK;i++) {
      r=ds->restr[i];
      this_vc=r->cvc;
      this_d2=0.0;
      for (j=0;j<ds->M;j++) {
	this_d=cvpt[j]-this_vc[j];
	this_d2+=this_d*this_d;
      }
      if (this_d2 < min_d2) {
	min_d2=this_d2;
	min_i=i;
      }
      //printf("CFACV/C/MIL) distance to center of cell %i is %.5lf\n",
      //     i,sqrt(this_d2));
    }
    //printf("CFACV/C) whichCell ends: CV-pt is in cell %i\n",min_i);
    /* if we are nearest to the home-cell center, then we need to check
       if we have left the spherical boundary of that cell, if it has one. */
    if (min_i==home_cell) {
      r=ds->restr[min_i];
      this_vc=r->cvc;

      // if this cell has a sphere boundary defined
       if (r->sbc) {
	 double sbr2=r->sbr*r->sbr;
	 //fprintf(stderr,"CFACV/MIL/DEBUG) checking to see if CV is outside sphere of radius %.5lf centered in cell %i at ",min_i,r->sbr);
	this_d2=0.0;
	for (j=0;j<ds->M;j++) {
	  this_d=cvpt[j]-r->sbc[j];
	  this_d2+=this_d*this_d;
	  //fprintf(stderr,"%.5lf ",r->sbc[j]);
	}
	//fprintf(stderr,"\nCFACV/MIL/DEBUG) distance = %.5lf\n",sqrt(this_d2));fflush(stderr);
	if (this_d2>sbr2) {
	  /* we are outside the sphere */
	  min_i=-1;
	}
      }
    }

    return min_i;
  }
}

void DataSpace_WriteCVtoArray ( DataSpace * ds, int * active, double * res ) {
  if (ds) {
    int i;
    for (i=0;i<ds->iM;i++) {
      if (active[i]) res[i]=ds->cv[i]->val;
      else res[i]=0.0;
    }
  }
}

void DataSpace_ReportAll ( DataSpace * ds ) {
  if (ds) {
    int i,j,k;
    printf("CFACV/C) DataSpace: %i Centers, %i CV's, %i Restraints\n",
	   ds->N, ds->iM, ds->iK);
    for (i=0;i<ds->N;i++) {
      printf("CFACV/C)  Center %i : %.5lf %.5lf %.5lf\n",
	     i,ds->R[i][0],ds->R[i][1],ds->R[i][2]);
    }
    for (i=0;i<ds->iM;i++) {
      printf("CFACV/C)      CV %i : typ %s val %.5lf ind ",i,cv_getstyp(ds->cv[i]->typ),
	     ds->cv[i]->val);
      for (j=0;j<ds->cv[i]->nC;j++) printf("%i ",ds->cv[i]->ind[j]);
      printf(" grads: ");
      for (j=0;j<ds->cv[i]->nC;j++) 
	for (k=0;k<3;k++) printf("%.5lf ",ds->cv[i]->gr[j][k]);
      printf("\n");
    }
    for (i=0;i<ds->iK;i++) {
      printf("CFACV/C)       R %i : k %.3lf z %.3lf cvc ",i,ds->restr[i]->k,ds->restr[i]->z);
      for (j=0;j<ds->restr[i]->nCV;j++) printf("%.5lf ",ds->restr[i]->cvc[j]);
      if (ds->restr[i]->tamdOpt) {
	printf("TAMD(kT=%.2f,gamma=%.2f,dt=%.3f) ",
	       ds->restr[i]->tamdOpt->kT,
	       ds->restr[i]->tamdOpt->gamma,
	       ds->restr[i]->tamdOpt->dt);
      }
      printf(" typ %s min %.5lf max %.5lf\n",rf_getstyp(ds->restr[i]->rfityp),ds->restr[i]->min,ds->restr[i]->max);
      printf("\n");
    }
    
  }
}

void DataSpace_ReportCVs ( DataSpace * ds, int step ) {
  int i;
  if (ds) {
    cvStruct * cv;
    fprintf(stdout,"CFACV/MIL/C) %i ",step);
    for (i=0;i<ds->M;i++) {
      cv=ds->cv[i];
      fprintf(stdout,"%i %.5lf ",i,cv->val);
    }
    fprintf(stdout,"\n");fflush(stdout);
  }
}

void DataSpace_Report_Z_and_F ( DataSpace * ds, int step, FILE * fp ) {
  int i;
  fprintf(fp,"CFACV/C) Z  % 10i ",step);
  for (i=0;i<ds->K;i++) {
    fprintf(fp,"% 11.5lf",ds->z[i]);
  }
  fprintf(fp,"\n");
  fprintf(fp,"CFACV/C) F  % 10i ",step);
  for (i=0;i<ds->K;i++) {
    fprintf(fp,"% 11.5lf",ds->f[i]);
  }
  fprintf(fp,"\n");
}

void DataSpace_ReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp ) {
  int i;
  
  if (outputlevel & 1) {
    fprintf(fp,"CFACV/C) Z  % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5lf",ds->restr[i]->z);
    }
    fprintf(fp,"\n");
  }
  if (outputlevel & 2) {
    fprintf(fp,"CFACV/C) Th % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5lf",ds->restr[i]->val);
    }
    fprintf(fp,"\n");
  }
  if (outputlevel & 4) {
    fprintf(fp,"CFACV/C) FD % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5lf",ds->restr[i]->f);//tamd_restraint);
    }
    fprintf(fp,"\n");
  }
  if (outputlevel & 8) {
    fprintf(fp,"CFACV/C) ND % 10i ",step);
    for (i=0;i<ds->iK;i++) {
      fprintf(fp,"% 11.5lf",ds->restr[i]->tamd_noise);
    }
    fprintf(fp,"\n");
  }
  fflush(fp);
/*   printf("CFACV/C) INFO %i ",step); */
/*   for (i=0;i<ds->iK;i++) { */
/*     printf(" | typ %s min %.5lf max %.5lf ",rf_getstyp(ds->restr[i]->rfityp),ds->restr[i]->min,ds->restr[i]->max); */
/*   } */
/*   printf("\n"); */
}

void DataSpace_BinaryReportRestraints ( DataSpace * ds, int step, int outputlevel, FILE * fp ) {
  int i;
  
  if (outputlevel & 1) {
    for (i=0;i<ds->iK;i++) {
      fwrite(&(ds->restr[i]->z),sizeof(double),1,fp);
    }
  }
  if (outputlevel & 2) {
    for (i=0;i<ds->iK;i++) {
      fwrite(&(ds->restr[i]->val),sizeof(double),1,fp);
    }
  }
  if (outputlevel & 4) {
    for (i=0;i<ds->iK;i++) {
      fwrite(&(ds->restr[i]->tamd_restraint),sizeof(double),1,fp);
    }
  }
  if (outputlevel & 8) {
    for (i=0;i<ds->iK;i++) {
      fwrite(&(ds->restr[i]->tamd_noise),sizeof(double),1,fp);
    }
  }
  if (outputlevel & 16) {
    double force;
    for (i=0;i<ds->iK;i++) {
      force=ds->restr[i]->k*(ds->restr[i]->z-ds->restr[i]->val);
      fwrite(&(force),sizeof(double),1,fp);
    }
  }
  fflush(fp);
}

/* Update image 'img' based on gradients and metric tensors */
int SMDataSpace_RawImageUpdate ( SMDataSpace * sm, double stepsize, int img ) {
  if (sm) {
    int ni=sm->ni;
//    int tst;
    if (sm->dual) ni/=2;
//    tst=img + ((sm->dual&&(img>=ni))?(-ni):0);
//    if (tst<0 || tst>sm->ni) {
//      fprintf(stderr,"ERROR: img-%i (of %i) index error (%i=%i+%i) for sm->g[] [%s] tst\n",
//         img,ni,tst,img,(sm->dual&&(img>=ni))?(-ni):0,
//         sm->dual?"DUAL":"SINGLE");fflush(stderr);exit(-1);
//    }
    if (img>=0 && img<ni) {
      double * z=sm->z[img];
      double * oldz=sm->oldz[img];
      /* if we are using dual-image mode, then if this image's index is less than ni,
         use the metric tensor from it's partner, shifted by ni;  if this image's
         index is greater than or equal to ni, use the gradient from the (-ni)-shifted
         partner */
      double * g=sm->g[ img + ((sm->dual&&(img>=ni))?(-ni):0)];
      double * M=sm->MM[img + ((sm->dual&&(img< ni))?  ni :0)];
      // put in abililty to see full update
      int n = sm->nz;
      int i,j;
      double tmp=0.0;
      for (i=0;i<n;i++) {
        tmp=0.0;
        for (j=0;j<n;j++) tmp+=g[j]*M[i*n+j];
        oldz[i]=z[i];
	z[i]-=tmp*stepsize;
	//z[i]-=g[i]*stepsize;
      }
    }
  }
}

int SMDataSpace_MoveString ( SMDataSpace * sm, double stepsize ) {
  if (sm) {
    int img;
    int ni=sm->ni;

    //fprintf(stderr,"CFACV/C) SM moving [%s]-string with stepsize %.6lf\n",sm->dual?"DUAL":"SINGLE",stepsize);fflush(stderr);

    // perform raw update of each z-position
    if (sm->evolve_ends) SMDataSpace_RawImageUpdate(sm,stepsize,0);
    if (sm->evolve_ends) SMDataSpace_RawImageUpdate(sm,stepsize,ni-1);
    for (img=1;img<(ni-1);img++) SMDataSpace_RawImageUpdate(sm,stepsize,img);

    // reparameterize
    SMDataSpace_reparameterize(sm);

    // done
  } 
}

SMDataSpace * New_stringMethod_Dataspace ( int ni, int nz, int outputlevel, double nu, int evolve_ends, int dual ) {
  SMDataSpace * sm = (SMDataSpace*)malloc(sizeof(SMDataSpace));
  int i,j;

  fprintf(stderr,"CFACV/C) Allocating new string method dataspace (%i[%s],%i)\n",ni,dual?"DUAL":"SINGLE",nz);
  fflush(stderr);

  if (!nz) {
    fprintf(stderr,"ERROR/SM: cannot allocate a 0-dimensional CV space!\n");
    exit(-1);
  }
 
  sm->ni=ni;
  sm->nz=nz;
  sm->outputlevel=outputlevel;
  sm->evolve_ends=evolve_ends;

  sm->dual=dual;

  sm->z=(double**)malloc(ni*sizeof(double*));
  for (i=0;i<ni;i++) sm->z[i]=(double*)malloc(nz*sizeof(double));
  sm->oldz=(double**)malloc(ni*sizeof(double*));
  for (i=0;i<ni;i++) sm->oldz[i]=(double*)malloc(nz*sizeof(double));
  sm->zn=(double**)malloc(ni*sizeof(double*));
  for (i=0;i<ni;i++) sm->zn[i]=(double*)malloc(nz*sizeof(double));

  sm->g=(double**)malloc(ni*sizeof(double));
  for (i=0;i<ni;i++) sm->g[i]=(double*)malloc(nz*sizeof(double));
  sm->MM=(double**)malloc(ni*sizeof(double));
  for (i=0;i<ni;i++) sm->MM[i]=(double*)malloc(nz*nz*sizeof(double));
 
  sm->L=(double*)malloc(ni*sizeof(double));

  sm->s=(double*)malloc(ni*sizeof(double));

  sm->ztyp=(int*)malloc(nz*sizeof(int));

  sm->nu = nu;

  fprintf(stderr,"CFACV/C) Finished string-method dataspace allocation.\n");
  fflush(stderr);

  return sm;
}

double * SMDataSpace_image_z ( SMDataSpace * sm, int i ) {
  if (sm) {
    if (i>=0 && i<sm->ni) 
      return sm->z[i];
    else return NULL;
  } else return NULL;
}

double * SMDataSpace_image_oldz ( SMDataSpace * sm, int i ) {
  if (sm) {
    if (i>=0 && i<sm->ni) 
      return sm->oldz[i];
    else return NULL;
  } else return NULL;
}

double * SMDataSpace_image_g ( SMDataSpace * sm, int i ) {
  if (sm) {
    if (i>=0 && i<sm->ni) return sm->g[i];
    else return NULL;
  } else return NULL;
}

double * SMDataSpace_image_M ( SMDataSpace * sm, int i ) {
  if (sm) {
    if (i>=0 && i<sm->ni) return sm->MM[i];
    else return NULL;
  } else return NULL;
}

double sm_euc ( double * a, double * b, int D ) {
  int i;
  double d=0.0;
  for (i=0;i<D;i++) d+=(b[i]-a[i])*(b[i]-a[i]);
  return sqrt(d);
}

void sm_sub ( double * d, double * a, double * b, int D ) {
  int i;
  for (i=0;i<D;i++) d[i]=a[i]-b[i];
}

void sm_scl ( double * d, double x, int D ) {
  int i;
  for (i=0;i<D;i++) d[i]*=x;
}

void sm_add ( double * d, double * a, double * b, int D ) {
  int i;
  for (i=0;i<D;i++) d[i]=a[i]+b[i];
}

void sm_cop ( double * d, double * a, int D ) {
  int i;
  for (i=0;i<D;i++) d[i]=a[i];
}

double sm_mag ( double * a, int D ) {
  double m=0.0;
  int i;
  for (i=0;i<D;i++) m+=a[i]*a[i];
  return sqrt(m);
}

double sm_dot ( double * a, double * b, int D ) {
  double m=0.0;
  int i;
  for (i=0;i<D;i++) m+=a[i]*b[i];
  return m;
}

int SMDataSpace_climb ( SMDataSpace * sm ) {
  if (sm) {
    int ni=sm->ni;
    int nz=sm->nz;
    double * dzt=(double*)malloc(nz*sizeof(double));
    double * dzi=(double*)malloc(nz*sizeof(double));
    double mdzi, fac;
    
    // compute the instantaneous local tangent vector prior to update
    sm_sub(dzi,sm->oldz[ni-1],sm->oldz[ni-2],nz);
    mdzi=sm_mag(dzi,nz);
    sm_scl(dzi,1.0/mdzi,nz);
    fprintf(stderr,"CFACV/C/SM) tangent at climbing end: %.5lf %.5lf\n",dzi[0],dzi[1]);
    fflush(stderr);
    // compute the last update to the end's z; note that this
    // reflects a raw update and therefore one application of hDMDG
    sm_sub(dzt,sm->z[ni-1],sm->oldz[ni-1],nz);
    fprintf(stderr,"CFACV/C/SM) dz at climbing end: %.5lf %.5lf\n",dzt[0],dzt[1]);
    fflush(stderr);
    // dot this into tangent to find the component of the update that was tangent
    fac=sm_dot(dzt,dzi,nz);
    sm_scl(dzi,-fac*sm->nu,nz);
    // climb!
    sm_add(sm->z[ni-1],sm->z[ni-1],dzi,nz);
    fprintf(stderr,"CFACV/C/SM) climbing force: %.5le %.5le\n",dzi[0],dzi[1]);
    fflush(stderr);
    free(dzi);
    free(dzt);
    return 0;
  }
  return -1;
}

int SMDataSpace_set_reparam_tol ( SMDataSpace * sm, double reparam_tol, int maxiter ) {
  if (sm) {
    sm->reparam_tol=reparam_tol;
    sm->reparam_maxiter=maxiter;
  }
  return 0;
}

int SMDataSpace_reparameterize ( SMDataSpace * sm ) {
  if (sm) {
    int ni=sm->ni;
    int nz=sm->nz;
    double * L=sm->L;
    double * s=sm->s;
    int i,k,tk;
    double del=0.0;
    int iter=0;
    double dL;
    double * tmp1=(double*)malloc(nz*sizeof(double));
    double * tmp2=(double*)malloc(nz*sizeof(double));
    double err=1000.0,sum, sum2, d;
    int no_iter=0;

    /* if we are using dual-images, then the first half of the images are the reference string; and the
       second half are just slaved congruently */
    if (sm->dual) {
      ni/=2;
      fprintf(stdout,"CFACV/C) reparam of dual-image string; using %i images out of %i as reference\n",
        ni,sm->ni);
    }

    /* preserve the ends; by definition the ends don't move as a result of reparameterization! */
    sm_cop(sm->zn[0],sm->z[0],nz);
    sm_cop(sm->zn[ni-1],sm->z[ni-1],nz);

    if (sm->reparam_tol==0.0) no_iter=1;
    while (no_iter || (iter < sm->reparam_maxiter && err > sm->reparam_tol)) {
    /* 1. compute actual running arc length of the un-reparameterized string */
    L[0]=0.0;
    for (i=1;i<ni;i++) L[i]=L[i-1]+sm_euc(sm->z[i],sm->z[i-1],nz);
    /* 2. compute desired (equidistant) running arc length */
    s[0]=0.0;
    dL=L[ni-1]/(ni-1);
    for (i=1;i<ni;i++) s[i]=i*dL;
 
    if (sm->outputlevel>0) {
      fprintf(stdout,"CFACV/C) string pre-reparam interimage stepsizes: |");
      for (i=1;i<ni;i++) {
        fprintf(stdout,"%.5lf|",sm_euc(sm->z[i-1],sm->z[i],nz));
      }
      fprintf(stdout,"\n");
    }

    /* 3. for each point on raw string, excluding ends... */
    for (i=1;i<(ni-1);i++) {
      tk=-1;
      for (k=1;tk==-1&&k<ni;k++) {
        if ((L[k-1])<s[i] && s[i]<=(L[k]+del)) {
          tk=k;
        }
      }
      /* tk is now the index of the point on the un-reparameterized string
         for which the actual arc length is the ceiling of the desired arc
         length for the i'th point */
      if (tk!=-1) {
        /* compute vector displacement between tk and tk-1 */
        sm_sub(tmp1,sm->z[tk],sm->z[tk-1],nz);
        /* scale result */
        sm_scl(tmp1,(s[i]-L[tk-1])/(L[tk]-L[tk-1]),nz);
        /* add to previous position */
        sm_add(tmp2,tmp1,sm->z[tk-1],nz);
        sm_cop(sm->zn[i],tmp2,nz);
      } else {
        fprintf(stderr,"ERROR: could not reparameterize at image %i\n",i);
        exit(1);
      }
    }

    for (i=0;i<ni;i++) sm->z[i]=sm->zn[i];
    if (sm->dual) for (i=ni;i<sm->ni;i++) sm->z[i]=sm->zn[i]=sm->z[i-ni];
 
    /* check the convergence of the reparameterization */
    err=0.0;
    sum=0.0;
    sum2=0.0;
    for (i=1;i<ni;i++) {
      d=sm_euc(sm->z[i-1],sm->z[i],nz);
      sum+=d;
      sum2+=d*d;
    }
    err=sqrt((sum2-sum*sum/(ni-1))/(ni-1));
    if (sm->outputlevel>0) 
      fprintf(stdout,"CFACV/C) reparam iter %i err %.5le tol %.5le\n",iter,err,sm->reparam_tol);
      if (no_iter) no_iter=0;
      iter++;
    }
    if (iter==sm->reparam_maxiter) {
      fprintf(stderr,"ERROR: string reparameterization did not converge after %i iterations (%.5le > %.5le)\n",iter,err,sm->reparam_tol);
      exit(1);
    }

    if (sm->outputlevel>0) {
      fprintf(stdout,"CFACV/C) string post-reparam interimage stepsizes after %i iterations (%.5le): |",iter,err);
      for (i=1;i<ni;i++) {
	fprintf(stdout,"%.5lf|",sm_euc(sm->z[i-1],sm->z[i],nz));
      }
      fprintf(stdout,"\n");
      if (sm->dual) {
        fprintf(stdout,"CFACV/C)    plus partner images in dual-image configuration |");
        for (i=ni+1;i<sm->ni;i++) {
          fprintf(stdout,"%.5lf|",sm_euc(sm->z[i-1],sm->z[i],nz));
        }
        fprintf(stdout,"\n");
      }
      fflush(stdout);
    }
    free(tmp1);
    free(tmp2);
    return 0;
  }
  return -1;
}

