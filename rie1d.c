// ##################################################################
//
// rie1d.c
//
// Determine far-field boundary data using quasi 1-d 
// characteristic relations.
// Written by Yong Su Jung
// ##################################################################

#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void rie1d(GRID *g, SOLN *s)
{
  int i,j;
  double gm1,rho0,uout,vout,aout,s0;
  double xgm1,xnorm,ynorm,dnorm;
  double rhoi,ui,vi,pi,unormi,unormo;
  double u2,a2,ai,rplus,rminus,unorm,a;
  double u,v,ss,rho,p;
  
  gm1 = gamm-1;
  rho0 = rinf;
  uout = s->uinf;
  vout = s->vinf;
  aout = 1.0;
  s0   = aout*aout/(gamm*pow(rho0,gm1));

  xgm1 = 1.0/gm1; 
  xnorm = s->qb[4];
  ynorm = s->qb[5];
  
  dnorm = sqrt(xnorm*xnorm+ynorm*ynorm);
  xnorm = xnorm/dnorm;
  ynorm = ynorm/dnorm;

  rhoi   = s->qb[0];
  ui     = s->qb[1]/rhoi;
  vi     = s->qb[2]/rhoi;
  pi     = (gamm-1)*(s->qb[3]-0.5*(s->qb[1]*s->qb[1]+s->qb[2]*s->qb[2])/rhoi);
  unormi = ui*xnorm + vi*ynorm;
  unormo = uout*xnorm + vout*ynorm;
  u2     = ui*ui + vi*vi;
  a2     = gamm*pi/rhoi;
  a      = sqrt(a2);
  ai     = a;
  rplus  = unormi + 2.0*a/gm1;
  rminus = unormo - 2.0*aout/gm1;
  unorm  = 0.5*(rplus + rminus);
  a      = 0.25*gm1*(rplus - rminus);

  // if unorm > 0 this is outflow: take variables from inside
  // if unorm > 0 this is outflow: take variables from inside
  
  if(unorm>0.0) 
  {
    u = ui + xnorm*(unorm - unormi);
    v = vi + ynorm*(unorm - unormi);
    ss = ai*ai/(gamm*pow(rhoi,gm1));
  }
  else
  {
    u = uout + xnorm*(unorm - unormo);
    v = vout + ynorm*(unorm - unormo);
    ss = s0;
  }


  rho = pow((a*a/(gamm*ss)),xgm1);
  p   = rho*a*a/gamm;

  s->qb[0] = rho;
  s->qb[1] = u*rho;
  s->qb[2] = v*rho;
  s->qb[3] = p/gm1 + 0.5*rho*(u*u+v*v);
}
