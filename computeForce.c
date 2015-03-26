#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"

void computeForce(GRID *g,SOLN *s)
{
  int i,n;
  double x0,x1,x2,x3,x4,y0,y1,y2,y3,y4;
  int iface,iface2,node1,node2,node3,node4,icell;
  double rho,rhou,rhov,e,p,temp;
  double fac,facmu;
  double cs,ss;
  double c2b=0.3678, rgas=1./1.4;
  double mu,tau,dis,vel;
  double c1=4.0/3.0, c2=2.0/3.0;
  FILE   *fp;
  s->cx=s->cy=s->cl=s->cd=0;

  //fp=fopen("output/skinfriction.dat","w");

  for (i=0;i<g->nbfaces;i++)
  {
    iface=g->bfaces[i];
    node1=g->faces[6*iface];
    node2=g->faces[6*iface+1];
    x1=g->x[2*node1];
    y1=g->x[2*node1+1];
    x2=g->x[2*node2];
    y2=g->x[2*node2+1];
    icell=g->faces[6*iface+2];
    //
    rho=s->q[NVAR*icell];
    rhou=s->q[NVAR*icell+1];
    rhov=s->q[NVAR*icell+2];
    e=s->q[NVAR*icell+3];
    //
    p=(gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov)/rho);
    //
    s->cy-=p*(x2-x1);
    s->cx+=p*(y2-y1);

    //skin friction
    if(g->visc)
    {
       iface2 = g->c2f[4*icell+2];
       // calculate wall normal distance
       node3=g->faces[6*iface2];
       node4=g->faces[6*iface2+1];

       x3=g->x[2*node3];
       y3=g->x[2*node3+1];
       x4=g->x[2*node4];
       y4=g->x[2*node4+1];

       x0 = (x1+x2+x3+x4)*0.25;
       y0 = (y1+y2+y3+y4)*0.25;

       dis = fabs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1));
       dis = dis/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
       
       // sutherland's law for laminar viscosity
       temp=p/rho/rgas;
       mu=(c2b+1.)*temp*sqrt(temp)/(c2b+temp);       
       vel = sqrt((rhou/rho)*(rhou/rho)+(rhov/rho)*(rhov/rho));

       tau = 1./(s->rey)*mu*vel/(dis/g->vol[icell]);
       
       s->cy += tau*(y2-y1);
       s->cx += tau*(x2-x1);

       //printf("tau:%lf\n",tau);
       //fprintf(fp,"%lf %lf\n",x0,tau/(0.5*s->mach*s->mach));

    }
  }
 
  //fclose(fp);
  fac=0.5*rinf*s->mach*s->mach;
  s->cy/=fac;
  s->cx/=fac;
  cs=cos(s->alpha*deg2rad);
  ss=sin(s->alpha*deg2rad);
  s->cl=s->cy*cs-s->cx*ss;
  s->cd=s->cy*ss+s->cx*cs;      
}


