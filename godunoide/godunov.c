/* godunov.c *************************************************\
*                                                             *
*  Classical first order Godunov scheme for the system        *
*  of conservation laws of 1D gasdynamics,                    *
*  uses exact Newton-type iterative Riemann solver of Godunov *
*                                                             *
*  written by  Andrei A. Chernousov                           *
*  E-mail: <andrei99@iname.com>                               *
*  http://www.geocities.com/andrei_chernousov                 *
*  October 11, 2002                                          * 
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double
  *U1, *U2, *U3, /* conservative variables */
  *F1, *F2, *F3, /* fluxes */

  R, U, P, /* primitive variables */
  R_1, U_1, P_1,
  R_2, U_2, P_2,

  deltaX, deltaT, deltaT_X;

unsigned LEN, LENN, i, i_, step, numstep;

FILE *pF;

#include "nondim.c" // for nondimensionalizing

/*--- exact Newton-type iterative Riemann solver of Godunov ---*/
#include "solver.c"
//#include "godunoff.c"

/************
*  Array1D  *
************/
double *Array1D(unsigned size)
{
double *x;

  if ((x = (double *)malloc(size * sizeof(double))) == NULL) {
    fprintf(stderr, "Can't allocate memory\n");
    exit(-1);
  }
  return x;
   
} /* end Array1D() */


/*************
* Initialize *
*************/
void Initialize(void)
{
  char str[80];
  double p1, r1, p4, r4;

  /* loading data from file "input" */
  if ((pF = fopen("input", "r")) == NULL) {
    fprintf(stderr, "Can't open file \"input\"\n\n");
    exit(-1);
  }
  i = 0;
  do{
    if (fgets(str, 80, pF) != NULL) {
      switch (i) {
      case 0:  LEN    = atoi(str); break;
      case 1:  deltaX = atof(str); break;
      case 2:  deltaT = atof(str); break;
      case 3:  numstep= atoi(str); break;
      case 4:  p4     = atof(str); break;
      case 5:  r4     = atof(str); break;
      case 6:  p1     = atof(str); break;
      case 7:  r1     = atof(str); break;
      default: break;
      } /* end switch() */
      i++;
    } /* end if() */
    else {
      fprintf(stderr, "Error reading file \"input\"\n\n");
      exit(-1);
    }
  } while (i < 8);
  fclose(pF);

  deltaT_X = deltaT / deltaX;
  LENN = LEN + 1;

  /* memory allocation */
  U1 = Array1D(LEN + 2); U2 = Array1D(LEN + 2); U3 = Array1D(LEN + 2);
  F1 = Array1D(LEN + 1); F2 = Array1D(LEN + 1); F3 = Array1D(LEN + 1);

  /* initial conditions */
  U = 0.; 
  R = r4, P = p4;
  for (i = 1; i <= LEN/2; i++) {
    U1[i] = R;
    U2[i] = R * U;
    U3[i] = P / K_1 + 0.5 * R * U * U;
  } /* end for() */
  R = r1, P = p1;
  for(; i <= LEN; i++) {
    U1[i] = R;
    U2[i] = R * U;
    U3[i] = P / K_1 + 0.5 * R * U * U;
  } /* end for() */

} /* end Initialize() */


/***********************
* BounCondInGhostCells *
***********************/
void BounCondInGhostCells(void)
{
  /* solid walls */
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
  for (i = 0, i_ = 1; i < LENN; i++, i_++) {
    /* leftmost parameters */
    U_1 =  U2[i ] / (R_1 = U1[i ]);
    P_1 = (U3[i ] - 0.5 * R_1 * U_1 * U_1) * K_1;
    /* rightmost parameters */
    U_2 =  U2[i_] / (R_2 = U1[i_]);
    P_2 = (U3[i_] - 0.5 * R_2 * U_2 * U_2) * K_1;
    /* exact Riemann solver of Godunov, see "solver.c" */
    Godunov(R_1, U_1, P_1,  R_2, U_2, P_2,  &R,  &U,  &P);
    //Godunoff(R_1, U_1, P_1,  R_2, U_2, P_2,  &R,  &U,  &P);
    /* fluxes */
    F1[i] = R * U;
    F2[i] = R * U * U + P;
    F3[i] = U * (K_K_1 * P + 0.5 * R * U * U);
  }
} /* end Fluxes() */


/**************
*  Evolution  *
**************/
void Evolution(void)
{
  for (i = 1; i < LENN; i++) {
    U1[i] += deltaT_X * (F1[i-1] - F1[i]);
    U2[i] += deltaT_X * (F2[i-1] - F2[i]);
    U3[i] += deltaT_X * (F3[i-1] - F3[i]);
  }
} /* end Evolution() */


/***********
*  Output  *
***********/
void Output(void)
{
  /* output to "godunov.dat" file */
  if ((pF = fopen("godunov.dat", "w")) == NULL) {
    fprintf(stderr, "Can't open \"godunov.dat\"\n\n");
    exit(-1);
  }
  for (i = 1; i <= LEN; i++) {
    U =  U2[i] / (R = U1[i]);
    P = (U3[i] - 0.5 * R * U * U) * K_1;
    fprintf(pF,"%g %g %g %g\n", (i-0.5) * deltaX, 
            NDimP(P), NDimR(R), NDimU(U));
  }
  fclose(pF);
   
} /* end Output() */


/*********
*  MAIN  *
**********/
int main(int argc, char *argv[])
{
  /*=== initialization ===*/
  Initialize();
   
   /*=== time integration ===*/
  while (step++ < numstep) {
    /*---  ---*/
//  printf("%d\n", step);
    /*--- BC in ghost cells ---*/
    BounCondInGhostCells();
    /*--- fluxes at cell interfaces ---*/
    Fluxes();
    /*--- evolution ---*/
    Evolution();
  } /* end while */
   
  /*--- output of the flowfield ---*/
  Output();

  /**/ 
  if (argc > 0)
    printf("%s: OK!\n\n", argv[0]);   
  return 0;
   
} /* end main() */

/*--- end of "godunov.c" ---*/
