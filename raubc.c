// ##################################################################
//
// raubc.c
//
// ##################################################################
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <string.h>
#include <stdio.h>


void raubc(GRID *g, SOLN *s)
{
  int i,k,ibtype;
  int rightCell,leftCell,node1,node2;
  double x1,y1,x2,y2;

  ibtype = 1; // characteristic inflow/outflow

  if(ibtype==1) 
  {
    
    for (i=0;i<g->nfaces;i++)
    {
       rightCell = g->faces[6*i+4];
       leftCell  = g->faces[6*i+2];

       if(rightCell == -1) //far-boundary
       {
         for(k=0;k<NVAR;k++)
         {
           s->qb[k] = s->q[NVAR*leftCell+k]; 
         }
         node1 = g->faces[6*i];
         node2 = g->faces[6*i];

         x1 = g->x[2*node1];
         y1 = g->x[2*node1+1];
         x2 = g->x[2*node2];
         y2 = g->x[2*node2+1];
         
         s->qb[4] =  (y2-y1);
         s->qb[5] = -(x2-x1);

         //calculate riemann invariant
         rie1d(g,s);

         // set again
         for(k=0;k<NVAR;k++)
         {
           s->q[NVAR*leftCell+k] = s->qb[k]; 
         }
      
       
       
       }      
    }
  }  
  else // other boundary condition.. 
  {
    printf("Need to be implemented\n");
  }




}
// ####################################
// END OF FILE
// ###################################

