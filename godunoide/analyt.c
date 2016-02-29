#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define k 1.4
#define R 287.1

#include "nondim.c" // for nondimensionalizing the output


// gasdynamical functions
double alfa__   (double M ) { return     1. / (1. + 0.5 * (k - 1.) * M); }
double epsilon__(double M ) { return pow(1. / (1. + 0.5 * (k - 1.) * M), 2.     / (k - 1.)); }
double pi__     (double M ) { return pow(1. / (1. + 0.5 * (k - 1.) * M), 2. * k / (k - 1.)); }


double
  t1, p1, r1, a1,
  t4, p4, r4, a4,

  x0, x = 50., t = 0.064;

double 
  p2, M3, pi, u2, u3, a3;


/********
*  FUN  *
********/
double Fun(double My)
{
  u2 = (2. * a1 / (k + 1.)) * (My - 1/My);
  p2 = p1 * (1. + (2. * k / (k + 1.)) * (My * My - 1.));
  pi = p2 / p4;
  M3 = (2./ (k - 1.)) * (1. / pow(pi, (k - 1.) / (2. * k)) - 1.);
  a3 = a4 / (1. + ((k - 1.) / 2.) * M3);
  u3 = a3 * M3;
  return u2 - u3;
} // end Fun()


/**********
*  DIXOT  * 
**********/
double Dixot(double a,
	     double b,
	     double epsilon)
{
  double c, fb;

  if ((Fun(a)) * (fb = Fun(b)) > 0.) {
    printf("wrong [a, b] range in Dixot()\n");
    exit(-1);
  }
  while ((b - a) > epsilon) {
    c = 0.5 * (b + a);
    if (Fun(c) * fb >= 0.) 
      b = c;
    else 
      a = c;
  }
  return c;

} // end Dixot()


/*********
*  MAIN  *
*********/
int main(int argc, char *argv[])
{
  int N = 15, i;
  double dM, M, My, r, a, u, c, w, p;
  FILE *pF;
  char str[80];
  // extra parameters read from "input" file
  double deltaX, deltaT;
  int numstep, LEN;

    
  // loading initial data from file "input"
  if ((pF = fopen("input", "r")) == NULL) {
    fprintf(stderr, "Cannot open file \"input\"\n\n");
    exit(-1);
  }
  i = 0;
  do {
    if (fgets(str, 80, pF) != NULL) {
      switch (i) {
      case 0: LEN    = atoi(str); break;
      case 1: deltaX = atof(str); break;
      case 2: deltaT = atof(str); break;
      case 3: numstep= atoi(str); break;
	// left
      case 4: p4     = atof(str); break;
      case 5: r4     = atof(str); break;
	// right 
      case 6: p1     = atof(str); break;
      case 7: r1     = atof(str); break;
      default: break;
      }
      i++;
    } // end if()
    else {
      fprintf(stderr, "Error reading \"input\" file \n");
      exit(-1);
    }
  } while (i < 8);
  fclose(pF);

  // domain
  x0 = LEN * deltaX / 2.;
  t = deltaT * numstep;  

  // other initial parameters
  t4 = p4 / (R * r4);
  a4 = sqrt(k * R * t4);
  t1 = p1 / (R * r1);
  a1 = sqrt(k * R * t1);

  // find the My for the exact ("analytical") solution
  My = Dixot(1.02, 4.0, 0.00000000001);

  // writing profiles of parameters from the exact solution to file
  pF = fopen("analyt.dat", "w");
  //
  fprintf(pF, "0.000001 %f %f %f \n", NDimP(p4), NDimR(r4), NDimU(0.) );
  for (M = 0, dM = M3/N; M <= M3; M += dM) {
    a = a4 * alfa__   (M);
    r = r4 * epsilon__(M);
    p = p4 * pi__     (M);
    u = M * a;
    c = u - a;
    x = x0 + c * t;
    fprintf(pF, "%f %f %f %f \n", x, NDimP(p), NDimR(r), NDimU(u));
  }
//printf("My = %.10f, M3 = %f, u = %f \n", My, u / (c = sqrt(k * p /r)), u);
  //
  fprintf(pF, "%f %f %f %f \n", x0 + u3 * t, NDimP(p), NDimR(r), NDimU(u3));
  w = a1 * My;
  r = - r1 * w / (u3 - w);
  fprintf(pF, "%f %f %f %f \n", x0 + u3 * t + 0.000001, NDimP(p),  NDimR(r),   NDimU(u3));
  fprintf(pF, "%f %f %f %f \n", x0 + w  * t,            NDimP(p),  NDimR(r),   NDimU(u3));
  fprintf(pF, "%f %f %f %f \n", x0 + w  * t + 0.000001, NDimP(p1), NDimR(r1),  NDimU(0.));
  fprintf(pF, "%f %f %f %f \n", 2 * x0 - 0.000001,      NDimP(p1), NDimR(r1),  NDimU(0.));
  fclose(pF);

  // success!	
  if (argc > 0)
    printf("%s: OK!\n\n", argv[0]); 
  return 0;

} // end main()

// end of "analyt.c"
