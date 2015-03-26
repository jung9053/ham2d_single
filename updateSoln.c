// ##################################################################
//
// updateSoln.c
//
// Update the solution scheme (1st order EE)
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
void updateSoln(double *qsrc,  // "source/initial" states
                double *qdest, // "destination/updated" states
                double *res,   // residual
                double *sigma, // spectral radius
	              double dt,     // timestep
                double CFL,    // CFL number
                double coef,   // scaling factor
                int    ncells) // total number of cells
{
  int    i,j,m;
  double dtfac;
  double dtfac_min;
  //
  // 1st order euler explicit for now
  //
  m = 0;
  //
   dtfac_min = 1.e9;
  for(i=0;i<ncells;i++)
  {
    for(j=0;j<NVAR;j++)
	  {
      dtfac    = coef*CFL/sigma[i];
      
      qdest[m] = qsrc[m]+dtfac*res[m];
      m++;
	  }
  }

}
// ##################################################################
// END OF FILE
// ##################################################################
