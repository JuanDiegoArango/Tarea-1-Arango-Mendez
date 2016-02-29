/* roe_fds.c *****************************************\
*                                                     * 
*  First-order scheme for 1D Conservation Laws of     *
*  Gas Dynamics with Roe flux difference splitting    *
*  (FDS)                                              *
*  author: Andrei Chernousov                          *
*          E-mail: <andrei99@iname.com>               *
*          http://www.geocities.com/andrei_chernousov *
*  November 04, 2001                                  *
*                                                     *   
\*****************************************************/

/*

 References:

 1. Roe P. L.
        Approximate Riemann solvers, parameter vectors and difference scheme.
	// J. Comp. Phys. 1981. v. 43. p. 357-472.

 2. Roe P. L. Characteristic-based schemes for the Euler equations.
	// Ann. Rev. Fluid. Mech. 1986. v. 8. p. 337-369.

 3. file "roe_fds.ps" in the current directory or at url: 
        http://www.geocities.com/andrei_chernousov/CODES/roe_fds.tar.gz
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* combinations of the ratio of specific heats */
#define K     1.40
#define K_1   0.40
#define K_K_1 3.50

/* alias names for fluxes */
#define F1 _U1
#define F2 _U2
#define F3 _U3

float
  /* arrays */
  *U1,  *U2,  *U3,
  *_U1, *_U2, *_U3,
  *U1_, *U2_, *U3_,
  /* spacing and time step */
  deltaX, deltaT, deltaT_X;

unsigned
  LEN,           /* number of cells */
  LENN,          /* LEN+1 */
  i, i_, _i,     /* indices */
  step, numstep; /* current time step and overall number of time steps */

/* for flux-difference splitting procedure of Roe */
float
  /* left state */
  rL, uL, pL,
  /* right state */
  rR, uR, pR,
  /* coefficients */
  D, C,
  /* enthalpies */
  HL, HR,
  /* Roe-averaged quantities */
  r, u, H, a, a2,
  /* right and left fluxes */
  FR1, FR2, FR3,   FL1, FL2, FL3,
  /* differences of primitive variables */
  dr, dp, du,
  /* differences of fluxes */
  dF11, dF12, dF13,   dF41, dF42, dF43,   dF51, dF52, dF53;


/************
*  Array1D  *   Allocates memory for array
************/
float *Array1D(unsigned size)
{
  float *x;

  if ((x = (float *)malloc(size * sizeof(float))) == NULL) {
    fprintf(stderr, "Cannot allocate memory\n");
    exit(-1);
  }
  return x;
   
} /* end Array1D() */


/*************
* Initialize *
*************/
void Initialize(void)
{
  float P, R, U;
  char str[80];
  FILE *pF;

  /* loading data from "roe_fds.ini" file */
  if ((pF = fopen("roe_fds.ini", "r")) == NULL) {
    fprintf(stderr, "Cannot open file \"roe_fds.ini\"\n");
    exit(-1);
  }
  i = 0;
  do {
    if (fgets(str, 80, pF) != NULL) {
      switch (i) {
      case 0:  LEN     = atoi(str); break;
      case 1:  deltaX  = atof(str); break;
      case 2:  deltaT  = atof(str); break;
      case 3:  numstep = atoi(str); break;
      default: break;
      }
      i++;
    }
    else {
      fprintf(stderr,"Error while reading file \"roe_fds.ini\"\n");
      exit(-1);
    }
  } while (i < 4);
  fclose(pF);

  LENN = LEN+1;
  deltaT_X = deltaT / deltaX;

  /* memory allocation */
  U1  = Array1D(LEN + 2);  U2 = Array1D(LEN + 2);  U3 = Array1D(LEN + 2);
  _U1 = Array1D(LEN + 1); _U2 = Array1D(LEN + 1); _U3 = Array1D(LEN + 1);
  U1_ = Array1D(LEN + 1); U2_ = Array1D(LEN + 1); U3_ = Array1D(LEN + 1);

  /* initial conditions */
  U = 0.; R = 2.; P = 200000.;
  for (i = 1; i <= LEN/2; i++) {
    U1[i] = R;
    U2[i] = R * U;
    U3[i] = P / K_1 + 0.5 * R * U * U;
  }
  R = 1.; P = 100000.;
  for (; i < LENN; i++) {
    U1[i] = R;
    U2[i] = R * U;
    U3[i] = P / K_1 + 0.5 * R * U * U;
  }

} /* end Initialize() */


/***********************
* BounCondInGhostCells *   boundary conditions in ghost cells
***********************/
void BounCondInGhostCells(void)
{
  /* reflective BC */
  U1[0] =   U1[1];
  U2[0] = - U2[1];
  U3[0] =   U3[1];
  U1[LENN] =   U1[LEN];
  U2[LENN] = - U2[LEN];
  U3[LENN] =   U3[LEN];
} /* end BounCondInGhostCells() */


/***********
*  Fluxes  *
***********/
void Fluxes(void)
{
  /* flux-difference splitting procedure of Roe */
  for (i = 0, i_ = 1; i < LENN; i++, i_++) {
    /* left state */
    rL =  U1[i];
    uL =  U2[i] / rL;
    pL = (U3[i] - 0.5 * rL * uL * uL) * K_1;
    /* right state */
    rR =  U1[i_];
    uR =  U2[i_] / rR;
    pR = (U3[i_] - 0.5 * rR * uR * uR) * K_1;
    /* total enthalpies */
    HL = K_K_1 * pL / rL + 0.5 * uL * uL;
    HR = K_K_1 * pR / rR + 0.5 * uR * uR;
    /* Roe averaging */
    r = rL * (D = sqrt(rR / rL));
    u = (uL + uR * D) / (1 + D);
    H = (HL + HR * D) / (1 + D);
    a = sqrt(a2 = K_1 * (H - 0.5 * u * u));
    /* right fluxes */
    FR1 = rR * uR;
    FR2 = rR * uR * uR + pR;
    FR3 = rR * uR * HR;
    /* left fluxes */
    FL1 = rL * uL;
    FL2 = rL * uL * uL + pL;
    FL3 = rL * uL * HL;
    /* differences */
    /* of primitive variables */
    dr = rR - rL;
    dp = pR - pL;
    du = uR - uL;
    /* 1 */
    dF11 = C = fabs(u) * (dr - dp / a2);
    dF12 = C * u;
    dF13 = C * 0.5 * u * u;
    /* 4 */
    dF41 = C = fabs(u + a) * ((dp + r * a * du) / (a2 + a2));
    dF42 = C * (u + a);
    dF43 = C * (H + u * a);
    /* 5 */
    dF51 = C = fabs(u - a) * ((dp - r * a * du) / (a2 + a2));
    dF52 = C * (u - a);
    dF53 = C * (H - u * a);
    /* fluxes */
    F1[i] = 0.5 * (FL1 + FR1 - dF11 - dF41 - dF51);
    F2[i] = 0.5 * (FL2 + FR2 - dF12 - dF42 - dF52);
    F3[i] = 0.5 * (FL3 + FR3 - dF13 - dF43 - dF53);
  }
} /* end Fluxes( ) */


/**************
*  Evolution  *
**************/
void Evolution(void)
{
  for (i = 1, _i = 0; i <= LEN; i++, _i++) {
    U1[i] += deltaT_X * (F1[_i] - F1[i]);
    U2[i] += deltaT_X * (F2[_i] - F2[i]);
    U3[i] += deltaT_X * (F3[_i] - F3[i]);
  }
} /* end Evolution() */


/***********
*  OUTPUT  *
***********/
void Output(void)
{
  float P, R, U;
  FILE *pF;

  /* output to "input.dat" file */
  if ((pF = fopen("input.dat", "w")) == NULL) {
    fprintf(stderr, "Cannot open file \"input.dat\"\n");
    exit(-1);
  }
//  fprintf(pF, "%d\n%d\n", LEN + 2, 3);
//  fprintf(pF, "0. 90000. 0.9 0. \n");
  for (i = 1; i <= LEN; i++) {
    U =  U2[i] / (R = U1[i]);
    P = (U3[i] - 0.5 * R * U * U) * K_1;
    fprintf(pF,"%3.3f %6.0f %1.4f %3.1f \n", (i - 0.5) * deltaX, P, R, U);
  }
//  fprintf(pF,"%d 210000. 2.1 100. \n", LEN);
//  fclose(pF );

} /* end Output() */


/*********
*  MAIN  *
**********/
int main(void)
{
  /*--- initialization ---*/
  Initialize();
   
  /*---*/
  while (step < numstep) {
    /*--- increment of step count ---*/
    //printf("%d\n", ++step);
    step++;
    /*--- BC in ghost cells ---*/
    BounCondInGhostCells();
    /*--- fluxes at cell interfaces ---*/
    Fluxes();
    /*--- evolution stage ---*/
    Evolution();
  }
   
  /*--- output of the flowfield ---*/
  Output();
   
  return 0;
   
}  /* end main() */

/*--- end roe_fds.c ---*/
