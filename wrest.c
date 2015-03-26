//#################################################################
//
// wrest.c
//
// subroutine for store restart file
//
// Written by Yong Su Jung
//##################################################################
#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

//====================================================================
//
//===================================================================
void wrest(GRID *g,SOLN *s, int n, int nn)
{
  int i;
  double fdummy;
  FILE *fp;
  char fname[80];
  fdummy=0;

  sprintf(fname,"./QuadData/urest%d.tc1",nn);
  fp=fopen(fname,"w");

  for(i=0;i<g->ncells;i++) 
  {
    fprintf(fp,"%lf %lf %lf %lf\n",
    s->q[NVAR*i],s->q[NVAR*i+1],s->q[NVAR*i+2],s->q[NVAR*i+3]);
  }
  fprintf(fp,"%d\n",n+1); //total number of iteration

  close(fp);
  printf("Writing restart files.......\n");
}
