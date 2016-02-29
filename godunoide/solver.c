/* solver.c **********************************************\
*                                                         *
*  Exact Newton-type iterative Riemann solver of Godunov  *
*                                                         *
*  written by  Andrei A. Chernousov                       *
*  E-mail: <andrei99@iname.com>                           *
*  http://www.geocities.com/andrei_chernousov             *
*  October 11, 2002                                      *
**********************************************************/

   /* complexes with the ratio of specific heats */
#define K        1.4
#define K_1_2    0.2
#define K_1      0.4
#define K_K      0.142857142857
#define _2K_K_1  7.0
#define K__1_2   1.2
#define K__1_2K  0.857142857143
#define K_1_2K   0.142857142857
#define K_1_K__1 0.166666666667
#define K_K_1    3.5
#define _2_K__1  0.833333333333
#define K__1     2.4
#define _3K_1    3.2
#define _2_K_1   5.0
#define _4K      5.6

/************
*  Godunov  *   Exact Newton-type Riemann solver of Godunov
************/
void Godunov(
  double r1,  double u1, double p1,  /* parameters in left cell */
  double r2,  double u2, double p2,  /* parameters in left cell */
  double *R,  double *U, double *P) /* parameters on the interface */
{
  static double
    c1, c2,    /* sound speeds in left and right cells */
    a1, a2,    /* mass velocities */
    c1_, c2_,  /* sound speeds in rarefaction waves */
    P_,        /* iterated pressure at the C(ontact)D(iscontinuity) */
    D1, D2,    /* velocties of the fronts of S(hock)W(ave)s and R(arefaction)Ws */
    D1_, D2_,  /* velocties of the tails of RWs */
    R1, R2,    /* densities to leftand right from CD */
    f, ff,     /* function and its derivative */
    Regim;     /* auxiliary variable */

  /* speeds of sound for left and right states */
  c1 = sqrt(K * p1 / r1);
  c2 = sqrt(K * p2 / r2);

  /* corrected parameters */
  if (p1 > p2) {
    Regim = p2; p2 = p1; p1 = Regim;
    Regim = r2; r2 = r1; r1 = Regim;
    Regim = u2; u2 = u1; u1 = Regim;
    Regim = c2; c2 = c1; c1 = Regim;
    Regim = - 1.;
  }
  else Regim = 1.;

  /* vacuum solution */
  if ((u1 - u2) <= - _2_K_1 * (c1 + c2) * 0.99999) {

    /* backward correction of parameters */
    if (Regim == -1.) {
      Regim = p2; p2 = p1; p1 = Regim;
      Regim = r2; r2 = r1; r1 = Regim;
      Regim = u2; u2 = u1; u1 = Regim;
      Regim = c2; c2 = c1; c1 = Regim;
    }

    /* velocities of vacuum characteristics */
    D1_ = _2_K_1 * c1 + u1;
    D2_ = u2 - _2_K_1 * c2;

    /* parameters from vacuum zone */
    if ((D1_ < 0.) && (D2_ > 0.)) {
      *P = *R = *U = 0.;
      return;
    }

    /* velocities of characteristics of disturbed zone */
    D1 = u1 - c1;
    D2 = u2 + c2;

    /* parameters from left RW */
    if ((D1 < 0.) && (D1_ >= 0.)) {
      *U = (c1_ =  _2_K__1 * c1 + K_1_K__1 * u1);
      *P = p1 * pow(c1_ / c1, _2K_K_1);
      *R = K * *P / (c1_ * c1_);
      return;
    }

    /* parameters from right RW */
    if ((D2_ <= 0.) && (D2 > 0.)) {
      *U = - (c2_ = _2_K__1 * c2 - K_1_K__1 * u2);
      *P = p2 * pow(c2_ / c2, _2K_K_1);
      *R = K * *P / (c2_ * c2_);
      return;
    }

    if (D1 >= 0.) { /* parameters from left stste */
      *P = p1, *R = r1, *U = u1; 
      return;
    }
    else {           /* parameters from right state */
      *P = p2, *R = r2, *U = u2; 
      return;
    }

  } /* end of vacuum solution */

  /* acoustic initial guess */
  P_ = (p1 * r2 * c2 + p2 * r1 * c1 + ( u1 - u2 ) * r1 * c1 * r2 * c2) / (r1 * c1 + r2 * c2);

  /* ? */
  if (P_ < 0.) 
    P_ = 10.;

  /* iterative Newton-type solution procedure */
  /* F(Pi-1) ( 13.16 ), F'(Pi-1) ( 13.17 ) */
  do {

    *P = P_;

    /* left wave */
    if (*P > p1) { /* SW */
      f   = (*P - p1) / (r1 * c1 * sqrt(K__1_2K * (*P / p1) + K_1_2K));
      ff  = (K__1 * (*P / p1) + _3K_1) / 
               (_4K * r1 * c1 * pow(K__1_2K * (*P / p1) + K_1_2K, 1.5));
    }
    else {         /* RW */
      f   = _2_K_1 * c1 * (pow(*P / p1, K_1_2K) - 1.);
      ff  = c1 * pow(*P / p1, K_1_2K) / (*P * K);
    }

    /* right wave */
    if (*P > p2) { /* SW */
      f  += (*P - p2) / (r2 * c2 * sqrt(K__1_2K * (*P / p2) + K_1_2K));
      ff += (K__1 * (*P / p2) + _3K_1) / 
               (_4K * r2 * c2 * pow(K__1_2K * (*P / p2) + K_1_2K, 1.5));
    }
    else {          /* RW */
      f  += _2_K_1 * c2 * (pow(*P / p2, K_1_2K) - 1.);
      ff += c2 * pow(*P / p2, K_1_2K) / (*P * K);
    }

    f -= (u1 - u2) * Regim;

    /* new P */
    P_ -= f / ff;

  } while (fabs(*P - P_) > 2.0); /* check convergence */

  /* final pressure */
  *P = P_;

  /* backward correction of parameters */
  if (Regim == -1.) {
    Regim = p2; p2 = p1; p1 = Regim;
    Regim = r2; r2 = r1; r1 = Regim;
    Regim = u2; u2 = u1; u1 = Regim;
    Regim = c2; c2 = c1; c1 = Regim;
  }

  /* computation of mass velocities */
    /* left wave is SW */
  if (*P > p1)
    a1 = sqrt(r1 * (K__1_2 * (*P) + K_1_2 * p1));
    /* left wave is RW */
  else if ((p1 - *P) > 300.) /* 300 Pa ! */
    a1 = K_1_2K * r1 * c1 * (1. - *P / p1) / (1.- pow(*P / p1, K_K));
  else 
    a1 = r1 * c1;
  /* rigt wave is SW */
  if (*P > p2)
    a2 = sqrt(r2 * (K__1_2 * (*P) + K_1_2 * p2));
  /* rigt wave is RW */
  else if ((p2 - *P) > 300.) /* 300 Pa ! */
    a2 = K_1_2K * r2 * c2 * (1. - *P / p2) / (1.- pow(*P / p2, K_K));
  else 
    a2 = r2 * c2;

  /* velocity of CD */
  *U = (a1 * u1 + a2 * u2 + p1 - p2) / (a1 + a2);

  /* characteristic velocities */
    /* left wave is SW */
  if (*P > p1) {
    D1 = u1 - a1 / r1;
    D1_ = D1;
    R1 = r1 * a1 / (a1 - r1 * (u1 - *U));
  }
   /* left wave is RW */
  else {          
    D1 = u1 - c1;
    c1_ = c1 + K_1_2 * (u1 - *U);
    D1_ = *U - c1_;
    R1 = K * *P / (c1_ * c1_);
  }
    /* right wave is SW */
  if (*P > p2) {  
    D2 = u2 + a2 / r2;
    D2_ = D2;
    R2 = r2 * a2 / (a2 + r2 * (u2 - *U));
  }
    /* right wave is RW*/
  else {
    D2 = u2 + c2;
    c2_ = c2 - K_1_2 * (u2 - *U);
    D2_ = *U + c2_;
    R2 = K * *P / (c2_ * c2_);
  }


  // parameters from zone 3 or 4
  if ((D1_ < 0.) && (D2_ > 0.)) {
    if (*U >= 0.)
      *R = R1; // 3 
    else
      *R = R2; // 4
    return;
  }

  // parameters from zone 2
  if ((D1 < 0.) && (D1_ >= 0.)) {
    *U = (c1_ = _2_K__1 * c1 + K_1_K__1 * u1);
    *P = p1 * pow(c1_/c1, _2K_K_1);
    *R = K * *P / (c1_ * c1_);
    return;
  }

  // parameters from zone 5
  if ((D2_ <= 0.) && (D2 > 0.)) {
    *U = - (c2_ = _2_K__1 * c2 - K_1_K__1 * u2);
    *P = p2 * pow(c2_/c2, _2K_K_1);
    *R = K * *P / (c2_ * c2_);
    return;
  }

  // parameters from the left stste - zone 1
  if (D1 >= 0.) { 
    *U = u1; 
    *P = p1;
    *R = r1; 
    return;
  }

  // parameters from the right state - 6
  *U = u2;
  *P = p2;
  *R = r2;
  return; 

} /* end Godunov() */

/*--- end solver.c ---*/
