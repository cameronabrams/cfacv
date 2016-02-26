#include "centers.h"
/* bin[i] contains the cluster id of the i'th atom
   x[i] is x-coordinate, y[i] is y-coordinate, z[i] is z-coordinate,
   N is number of atoms, nc is number of cycles

   This routine will reorder bin[] such that the centers of mass are 
   as close as equidistant as possible.
*/
atomCenterStruct * New_atomCenterStruct ( int n ) {
  atomCenterStruct * newac=malloc(sizeof(atomCenterStruct));
  newac->n=n;
  newac->ind=calloc(n,sizeof(int));
  newac->m=calloc(n,sizeof(double));
  newac->M=0.0;
  return newac;
}

int bin_sort ( int * bin, double * x, double * y, double * z, int nAtoms, int nCenters, int nCycles, 
	       unsigned int Seed ) {
  int i,j,nsucc;
  double FUZZ=0.0;

  int * cnt=malloc(nCenters*sizeof(int)), c, ii, jj, ib, jb;
  double ** cm=malloc(nCenters*sizeof(double*));
  double * rg=malloc(nCenters*sizeof(double*)), srg, srg0;
  for (i=0;i<nCenters;i++) {
    cm[i]=calloc(3,sizeof(double));
  }

  srand(Seed);

  for (i=0;i<nCenters;i++) {
    cm[i][0]=cm[i][1]=cm[i][2]=0.0;
    rg[i]=0.0;
    cnt[i]=0;
  }
  for (i=0;i<nAtoms;i++) {
    j=bin[i]-1;
    cm[j][0]+=x[i];
    cm[j][1]+=y[i];
    cm[j][2]+=z[i];
    cnt[j]++;
  }
  for (i=0;i<nCenters;i++) {
    for (j=0;j<3;j++) cm[i][j]/=cnt[i];
  }
  for (i=0;i<nAtoms;i++) {
    rg[bin[i]-1]+=(cm[bin[i]-1][0]-x[i])*(cm[bin[i]-1][0]-x[i]);
    rg[bin[i]-1]+=(cm[bin[i]-1][1]-y[i])*(cm[bin[i]-1][1]-y[i]);
    rg[bin[i]-1]+=(cm[bin[i]-1][2]-z[i])*(cm[bin[i]-1][2]-z[i]);
  }
  srg=0;
  for (i=0;i<nCenters;i++) {
    rg[i]=(cnt[i]>0.0)?(sqrt(rg[i]/cnt[i])):0.0;
    srg+=rg[i];
  }
  srg0=srg;

  fprintf(stdout,"CFACV/C) bin_sort initial center assignments: ");
  for (i=0;i<nAtoms;i++) fprintf(stdout,"%i ",bin[i]);
  fprintf(stdout,"\n");
  fflush(stdout);
  fprintf(stdout,"CFACV/C) bin_sort initial center rg's: ");
  for (i=0;i<nCenters;i++) fprintf(stdout,"%.3lf ",rg[i]);
  fprintf(stdout," average %.3lf\n",srg/nCenters);
  fflush(stdout);


  nsucc=0;

  for (c=0;c<nCycles;c++) {
    ii = (int)(nAtoms * (rand()/(RAND_MAX+1.0)));
    jj = (int)(nAtoms * (rand()/(RAND_MAX+1.0)));

    while (bin[ii]==bin[jj]) {
      ii = (int)(nAtoms * (rand()/(RAND_MAX+1.0)));
      jj = (int)(nAtoms * (rand()/(RAND_MAX+1.0)));
    }
    ib=bin[ii];
    jb=bin[jj];
    bin[jj]=ib;
    bin[ii]=jb;

/*     fprintf(stderr,"CFACV/C) cycle %i swap %i(%i) %i(%i)...\n", */
/* 	    c,ii,ib,jj,jb); fflush(stderr); */

    for (i=0;i<nCenters;i++) {
      cm[i][0]=cm[i][1]=cm[i][2]=0.0;
      rg[i]=0.0;
      cnt[i]=0;
    }

    for (i=0;i<nAtoms;i++) {
      cm[bin[i]-1][0]+=x[i];
      cm[bin[i]-1][1]+=y[i];
      cm[bin[i]-1][2]+=z[i];
      cnt[bin[i]-1]++;
    }

    for (i=0;i<nCenters;i++) {
      for (j=0;j<3;j++) cm[i][j]/=cnt[i];
    }
    for (i=0;i<nAtoms;i++) {
      rg[bin[i]-1]+=(cm[bin[i]-1][0]-x[i])*(cm[bin[i]-1][0]-x[i]);
      rg[bin[i]-1]+=(cm[bin[i]-1][1]-y[i])*(cm[bin[i]-1][1]-y[i]);
      rg[bin[i]-1]+=(cm[bin[i]-1][2]-z[i])*(cm[bin[i]-1][2]-z[i]);
    }

    srg=0.0;
    for (i=0;i<nCenters;i++) {
      rg[i]=sqrt(rg[i]/cnt[i]);
      srg+=rg[i];
    }

    if ((srg-srg0)<FUZZ) {
      srg0=srg;
      fprintf(stdout,"CFACV/C) cycle %i <<rg>> %.5lf\n",c,srg/nCenters);
      nsucc++;
    }
    else {
      bin[ii]=ib;
      bin[jj]=jb;
    }
  }
  
  for (i=0;i<nCenters;i++) {
    cm[i][0]=cm[i][1]=cm[i][2]=0.0;
    rg[i]=0.0;
    cnt[i]=0;
  }
  for (i=0;i<nAtoms;i++) {
    cm[bin[i]-1][0]+=x[i];
    cm[bin[i]-1][1]+=y[i];
    cm[bin[i]-1][2]+=z[i];
    cnt[bin[i]-1]++;
  }
  for (i=0;i<nCenters;i++) {
    for (j=0;j<3;j++) cm[i][j]/=cnt[i];
  }
  for (i=0;i<nAtoms;i++) {
    rg[bin[i]-1]+=(cm[bin[i]-1][0]-x[i])*(cm[bin[i]-1][0]-x[i]);
    rg[bin[i]-1]+=(cm[bin[i]-1][1]-y[i])*(cm[bin[i]-1][1]-y[i]);
    rg[bin[i]-1]+=(cm[bin[i]-1][2]-z[i])*(cm[bin[i]-1][2]-z[i]);
  }
  srg=0;
  for (i=0;i<nCenters;i++) {
    rg[i]=sqrt(rg[i]/cnt[i]);
    srg+=rg[i];
  }
  fprintf(stdout,"CFACV/C) bin_sort final center assignments: ");
  for (i=0;i<nAtoms;i++) fprintf(stdout,"%i ",bin[i]);
  fprintf(stdout,"\n");
  fflush(stdout);
  fprintf(stdout,"CFACV/C) bin_sort final center rg's: ");
  for (i=0;i<nCenters;i++) fprintf(stdout,"%.3lf ",rg[i]);
  fprintf(stdout," average %.3lf\n",srg/nCenters);
  fflush(stdout);
  
  for (i=0;i<nCenters;i++) {
    fprintf(stdout,"CFACV/C) cm %i (%i): %.5lf %.5lf %.5lf : rg %.5lf\n",
	    i+1,cnt[i],cm[i][0],cm[i][1],cm[i][2], rg[i]);
  }
  fprintf(stdout,"CFACV/C) successful swaps: %i out of %i cycles\n",nsucc,nCycles);
  fflush(stdout);
  free(cm);
  free(cnt);
}


/** The CenterStruct Tree -- useful for sorting residues into 
    centers */

centerStruct * New_centerStruct ( int id, int maxN ) {
  centerStruct * c = malloc(sizeof(centerStruct));
  c->id=id;
  c->maxN=maxN;
  c->iN=0;
  c->rg=0.0;
  c->cm[0]=c->cm[1]=c->cm[2]=0.0;
  c->mList=(int*)calloc(maxN,sizeof(int));
  c->left=NULL;
  c->right=NULL;
  return c;
}

void centerStruct_addMember ( centerStruct * c, int i ) {
  if (c) {
    c->mList[c->iN++]=i;
  }
}

void centerStruct_rg ( centerStruct * c, double * x, double * y, double * z ) {
  if (c) {
    int i;
    c->rg=0.0;
    c->cm[0]=c->cm[1]=c->cm[2]=0.0;

    for (i=0;i<c->iN;i++) {
      c->cm[0]+=x[c->mList[i]];
      c->cm[1]+=y[c->mList[i]];
      c->cm[2]+=z[c->mList[i]];
    }
    c->cm[0]/=c->iN;    c->cm[1]/=c->iN;    c->cm[2]/=c->iN;
    for (i=0;i<c->iN;i++) {
      c->rg+=(x[c->mList[i]]-c->cm[0])*(x[c->mList[i]]-c->cm[0])+
	(y[c->mList[i]]-c->cm[1])*(y[c->mList[i]]-c->cm[1])+
	(z[c->mList[i]]-c->cm[2])*(z[c->mList[i]]-c->cm[2]);
    }
    c->rg=sqrt(c->rg/c->iN);
  }
}

void centerStruct_randomSwapMember (centerStruct * c, centerStruct * d, int * i, int * j ) {
  if (c&&d) {
    int tmp;
    *i=(int)(c->iN * (rand()/(RAND_MAX+1.0)));
    *j=(int)(d->iN * (rand()/(RAND_MAX+1.0)));
    tmp=c->mList[*i];
    c->mList[*i]=d->mList[*j];
    d->mList[*j]=tmp;
  }
}

void centerStruct_SwapMembers ( centerStruct * c, centerStruct * d, int i, int j ) {
  if (c&&d) {
    int tmp;
    tmp=c->mList[i];
    c->mList[i]=d->mList[j];
    d->mList[j]=tmp;
  }
}

void centerStruct_split ( centerStruct * c, double * x, double * y, double * z, int nAtom ) {
  if (c) {
    int nacc=0;
    int nCycles = 100000;
    double rgsum,rgsum0;
    int i;
    int left, right;
    c->left=New_centerStruct(0,c->iN/2+1);
    c->right=New_centerStruct(0,c->iN/2+1);
    for (i=0;i<(c->iN/2);i++) centerStruct_addMember(c->left,c->mList[i]);
    for (; i<c->iN;i++)  centerStruct_addMember(c->right,c->mList[i]);  
    centerStruct_rg(c->left,x,y,z);
    centerStruct_rg(c->right,x,y,z);
    rgsum0=c->left->rg+c->right->rg;
    printf("CFACV/C) split; before cycles %.5lf %.5lf = %.5lf\n",c->left->rg,c->right->rg,rgsum0);
    
    for (i=0;i<nCycles;i++) {
      centerStruct_randomSwapMember(c->left,c->right,&left,&right);
      centerStruct_rg(c->left,x,y,z);
      centerStruct_rg(c->right,x,y,z);
      rgsum=c->left->rg+c->right->rg;
      //     printf("CFACV/C) swapped %i and %i : %.5lf %.5lf = %.5lf\n",left,right,c->left->rg,c->right->rg,rgsum);
      if (rgsum<rgsum0) {
	rgsum0=rgsum;
	nacc++;
      }
      else {
	centerStruct_SwapMembers(c->left,c->right,left,right);
	centerStruct_rg(c->left,x,y,z);
	centerStruct_rg(c->right,x,y,z);
	//	printf("CFACV/C) unswapped %i and %i : %.5lf %.5lf = %.5lf\n",left,right,c->left->rg,c->right->rg,rgsum);
      }
    }
    printf("CFACV/C) split; after cycles (%i/%i): %.5lf %.5lf = %.5lf\n",nacc,nCycles,c->left->rg,c->right->rg,rgsum);
  }
}

void centerStruct_binnify ( centerStruct * c, int * bin, int * id, double * rg ) {
  if (c) {
    //    printf("CFACV/C) traversal left %x right %x\n",c->left,c->right);
    if (!c->left && !c->right) {
      int i;
      // printf("CFACV/C) at a leaf id %i\n",*id);
      for (i=0;i<c->iN;i++) {
	//printf("CFACV/C) assigning bin id of atom %i to %i\n",c->mList[i],*id);
	bin[c->mList[i]]=*id+1;
	rg[*(id)]=c->rg;
      }
      (*id)++;
    }
    else {
      centerStruct_binnify(c->left,bin,id,rg);
      centerStruct_binnify(c->right,bin,id,rg);
    }
  }
}

int rgyr_sort (centerStruct * c, int * bin, double * x, double * y, double * z, int nAtom, int minAtom, double * rg, unsigned int Seed  ) {
  if (c) {
    if ((c->iN)/2 < minAtom) {
      return 0;
    }
    else {
      centerStruct_split(c,x,y,z,nAtom);
      rgyr_sort(c->left,bin,x,y,z,nAtom,minAtom,rg,Seed);
      rgyr_sort(c->right,bin,x,y,z,nAtom,minAtom,rg,Seed);
      return 0;
    }
  }
  else {
    int i,j;
    srand(Seed);
    c=New_centerStruct(0,nAtom);
    for (i=0;i<nAtom;i++) {
      centerStruct_addMember(c,i);
      bin[i]=0;
    }
    centerStruct_rg(c,x,y,z);
    printf("CFACV/C) init rg %.3lf\n",c->rg);
    rgyr_sort(c,bin,x,y,z,nAtom,minAtom,rg,Seed);
    i=0;
    centerStruct_binnify(c,bin,&i,rg);
    //    printf("CFACV/C) bins %i:\n",i);
    // for (j=0;j<nAtom;j++) printf("CFACV/C) %i %i\n",j,bin[j]);
    return i;
  }
}

centerStruct * Null_centerStruct ( void ) {
  return NULL;
}
