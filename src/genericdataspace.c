#include "genericdataspace.h"

GenDataSpace * NewGenDataSpace ( int N, int M ) {
  GenDataSpace * gds = malloc(sizeof(GenDataSpace));
  int i,j;
  gds->n=N;
  gds->m=M;
  
  gds->A=(double**)malloc(N*sizeof(double*));
  for (i=0;i<N;i++) gds->A[i]=(double*)malloc(M*sizeof(double));
  for (i=0;i<N;i++) for (j=0;j<M;j++) gds->A[i][j]=0.0; 
  return gds;
}

void FreeGenDataSpace ( GenDataSpace * gds ) {
  int i;
  for (i=0;i<gds->n;i++) free(gds->A[i]);
  free(gds->A);
  free(gds);
}

void GenDataSpace_Set ( GenDataSpace * gds, int i, int j, double val ) {
  gds->A[i][j]=val;
}

double GenDataSpace_Get ( GenDataSpace * gds, int i, int j ) {
  return gds->A[i][j];
}

void GenDataSpace_Scale ( GenDataSpace * gds, double x ) {
  int i,j;
  for (i=0;i<gds->n;i++) {
    for (j=0;j<gds->m;j++) {
      gds->A[i][j]*=x;
    }
  }
}

void GenDataSpace_Accumulate ( GenDataSpace * gds, int i, int j, double val ) {
  gds->A[i][j]+=val;
}

void GenDataSpace_Diff_RunningAverage ( GenDataSpace * gds ) {
  int i,j;
  double * sums=calloc(gds->m,sizeof(double));
  for (j=0;j<gds->m;j++) {
    sums[j]=0.0;
  }
  for (i=1;i<gds->n;i++) {
    for (j=0;j<gds->m;j++) {
      sums[j]+=gds->A[i][j]-gds->A[0][j];
      gds->A[i][j]=sums[j]/i;
      //      printf("C: i %i j %i th %.5lf z %.5lf\n",i,j,gds->A[i][j],gds->A[0][j]);
    }
  }
  for (j=0;j<gds->m;j++) {
    gds->A[0][j]=sums[j]/(gds->n-1);
  }
  free(sums);
}

IGenDataSpace * NewIGenDataSpace ( int N, int M ) {
  IGenDataSpace * igds = malloc(sizeof(GenDataSpace));
  int i;
  igds->n=N;
  igds->m=M;
  
  igds->I=(int**)malloc(N*sizeof(int*));
  for (i=0;i<N;i++) igds->I[i]=(int*)malloc(M*sizeof(int));
  return igds;
}

double * GenDataSpace_getAddr ( GenDataSpace * gds, int i ) {
  if (gds && i < gds->n) {
    return gds->A[i];
  }
  else return NULL;
}

int * IGenDataSpace_getAddr ( IGenDataSpace * igds, int i ) {
  if (igds && i < igds->n) {
    return igds->I[i];
  }
  else return NULL;
}

void GenDataSpace_WriteToFile ( GenDataSpace * gds, char * filename ) {
  int i,j;
  FILE * fp = fopen(filename,"w");
  fprintf(fp,"# gendataspace output (%i x %i)\n",gds->n,gds->m);
  for (i=0;i<gds->n;i++) {
    for (j=0;j<gds->m;j++) {
      fprintf(fp,"% 10.5lf",gds->A[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

void GenDataSpace_WriteToFile_Gnuplot ( GenDataSpace * gds, char * filename, int outerpad, char * directive ) {
  int i,j,j0;
  FILE * fp = fopen(filename,"w");
  fprintf(fp,"# gendataspace output (%i x %i)\n",gds->n,gds->m);
  for (i=0;i<gds->n;i++) {
    if (directive) {
      if (!strcmp(directive,"full")) j0=0;
      else j0=i;
    }
    for (j=j0;j<gds->m;j++) {
      fprintf(fp,"% 5i% 5i% 10.5lf\n",i,j,gds->A[i][j]);
    }
    if (outerpad) fprintf(fp,"% 5i% 5i% 10.5lf p\n",i,j,gds->A[i][j-1]);
    fprintf(fp,"\n");
  }
/*   if (outerpad) { */
/*     for (j=0;j<gds->m;j++) { */
/*       fprintf(fp,"% 5i% 5i% 10.5lf p\n",i,j,gds->A[i-1][j]); */
/*     } */
/*     fprintf(fp,"% 5i% 5i% 10.5lf p\n",i,j,gds->A[i-1][j-1]); */
/*     fprintf(fp,"\n"); */
/*   } */
  fclose(fp);
}
void GenDataSpace_WriteToFile_Gnuplot_IntIndices ( GenDataSpace * gds, int * I, int * J, char * filename, int outerpad, char * directive ){
  int i,j,j0;
  FILE * fp = fopen(filename,"w");
  fprintf(fp,"# gendataspace output (%i x %i)\n",gds->n,gds->m);
  for (i=0;i<gds->n;i++) {
    if (directive) {
      if (!strcmp(directive,"full")) j0=0;
      else j0=i;
    }
    for (j=j0;j<gds->m;j++) {
      fprintf(fp,"% 5i% 5i% 10.5lf% 5i% 5i\n",I[i],J[j],gds->A[i][j],i,j);
    }
    if (outerpad) fprintf(fp,"% 5i% 5i% 10.5lf% 5i% 5i p\n",I[i],J[j-1]+1,gds->A[i][j-1],i,j);
    fprintf(fp,"\n");
  }
}

void GenDataSpace_DistMap ( GenDataSpace * map, GenDataSpace * r1, GenDataSpace * r2 ) {
  /* map is n x m -- assume r1 is n x 3 and r2 is m x 3 */
  int i,j;
  double dx,dy,dz;
  
  if (map->n!= r1->n || map->m != r2->n) {
    fprintf(stderr,"CFACV) ERROR: DistMap coordinates don't match size of map.\n");
    exit(-1);
  }

/*   fprintf(stdout,"CFACV) GenDataSpace_DistMap: coordset 1: %ix%i 2: %ix%i map %ix%i\n", */
/* 	  r1->n,r1->m,r2->n,r2->m,map->n,map->m); */
/*   fflush(stdout); */

  for (i=0;i<r1->n;i++) {
    for (j=0;j<r2->n;j++) {
      dx=r1->A[i][0]-r2->A[j][0];
      dy=r1->A[i][1]-r2->A[j][1];
      dz=r1->A[i][2]-r2->A[j][2];
      map->A[i][j]=sqrt(dx*dx+dy*dy+dz*dz);
    }
  }
}

void GenDataSpace_DiffMap ( GenDataSpace * map1, GenDataSpace * map0) {
  int i,j;
  for (i=0;i<map1->n;i++) {
    for (j=0;j<map1->m;j++) {
      map1->A[i][j]-=map0->A[i][j];
    }
  }
}

void GenDataSpace_AddToMap ( GenDataSpace * map1, GenDataSpace * map0) {
  int i,j;
  for (i=0;i<map1->n;i++) {
    for (j=0;j<map1->m;j++) {
      map1->A[i][j]+=map0->A[i][j];
    }
  }
}

/* computes correlation coefficients */
void GenDataSpace_CorrCoeff (GenDataSpace * timedat, GenDataSpace * corrdat) {
  int i,j,k;
  double * ssum = malloc(sizeof(double)*timedat->n);
  double sum;
  
  for (i=0;i<timedat->n;i++) {
    sum=ssum[i]=0.0;
    for (j=0;j<timedat->m;j++) {
      sum+=timedat->A[i][j];
      ssum[i]+=timedat->A[i][j]*timedat->A[i][j];
    }
  }

  for (i=0;i<timedat->n;i++) {
    for (j=0;j<=i;j++) {
      sum=0.0;
      for (k=0;k<timedat->m;k++) {
	sum+=timedat->A[i][k]*timedat->A[j][k];
      }
      corrdat->A[i][j]=sum/sqrt((ssum[i]*ssum[j]));
    }
  }
  for (i=0;i<timedat->n;i++) {
    for (j=i+1;j<timedat->n;j++) {
      corrdat->A[i][j]=corrdat->A[j][i];
    }
  }

  free(ssum);
}
