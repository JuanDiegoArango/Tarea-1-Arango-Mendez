#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "riemann.h"
#define min(a, b) (((a) < (b)) ? (a) : (b))
/*
  Translated from a C++ Riemann solver written by Richard
  J. Gonsalves,  which in turn is rewrite of a Riemann solver
  found in Laney's textbook Computational Gasdynamics (Cambridge
  University Press.)     
 */ 

inline double fg(double x) {
    const double gamma = 1.4;
    const double g2 = (gamma + 1) / (2 * gamma);
    return (x-1) / sqrt(g2 * (x - 1) + 1);
}


void Riemann(double *U1, double *U4, double *F) {
    int i;
    float rhoRl;
    float uRl;
    float hRl;
    float aRl;
    float lambda[3];
    float deltaV[3];
    float r[3][3];
    
  const double gamma = 1.4;
  const double g1 = (gamma - 1) / (2 * gamma);
  const double g2 = (gamma + 1) / (2 * gamma);
  const double g3 = (gamma + 1) / (gamma - 1);
  const double tol = 1.0*1e-10;
  
    // compute primitive variables
  float rho1 = U1[0];
  float u1 = U1[1] / rho1;
  float p1 = (U1[2] - rho1 * u1 * u1 / 2.0) * (gamma - 1.0);
  float h1=(U1[0]+p1)/rho1;


    
    
  float rho4 = U4[0];
  float u4 = U4[1] / rho4;
  float p4 = (U4[2] - rho4 * u4 * u4 / 2.0) * (gamma - 1.0);
  float h4=(U4[0]+p4)/rho4;
    
   printf("%f,%f,%f,%f,%f,%f,%f,%f\n",rho1,rho4,u1,u4,p1,p4,h1,h4);

    

  // switch states if necessary so high pressure is on left
  int revflag = FALSE;
  if (p4 < p1) {
    double swap = p1; p1 = p4; p4 = swap;
    swap = u1; u1 = -u4; u4 = -swap;
    swap = rho1; rho1 = rho4; rho4 = swap;
    revflag = TRUE;
  }
    
    float delta_rho=rho1-rho4;
    float delta_p=p1-p4;
    float delta_u=u1-u4;


    rhoRl =sqrt(rho1*rho4);
    uRl=(sqrt(rho1)*u1+sqrt(rho4)*u4)/(sqrt(rho1)+sqrt(rho4));
    hRl=(sqrt(rho1)*h1+sqrt(rho4)*h4)/(sqrt(rho1)+sqrt(rho4));
    aRl=sqrt((gamma - 1)*(hRl-0.5*uRl*uRl));
    
    lambda[0]=uRl;
    lambda[1]=uRl-aRl;
    lambda[2]=uRl+aRl;
    
    deltaV[0]=delta_rho-delta_p/(aRl*aRl);
    deltaV[1]=delta_u+delta_p/(aRl*rhoRl);
    deltaV[2]=delta_u-delta_p/(aRl*rhoRl);
    
    double a=rhoRl/(2*aRl);
    
    
    r[0][0]=1.0;
    r[0][1]=uRl;
    r[0][2]=0.5*uRl*uRl;
    
    r[1][0]=a;
    r[1][1]=a*(uRl+aRl);
    r[1][2]=a*(hRl+aRl*uRl);
    
    r[2][0]=-a;
    r[2][1]=-a*(uRl-aRl);
    r[2][2]=-a*(hRl-aRl*uRl);
    

    F[0]=r[0][0]*min(0,lambda[0])*deltaV[0]+r[1][0]*min(0,lambda[1])*deltaV[1]+r[2][0]*min(0,lambda[2])*deltaV[2]+(rho1*u1);
    F[1]=r[0][1]*min(0,lambda[0])*deltaV[0]+r[1][1]*min(0,lambda[1])*deltaV[1]+r[2][1]*min(0,lambda[2])*deltaV[2]+(rho1*u1*u1+p1)*0.5;
    F[2]=r[0][2]*min(0,lambda[0])*deltaV[0]+r[1][2]*min(0,lambda[1])*deltaV[1]+r[2][2]*min(0,lambda[2])*deltaV[2]+(rho1*h1*u1);

//printf("%f,%f,%f,%f\n",F[0],uRl );


}
