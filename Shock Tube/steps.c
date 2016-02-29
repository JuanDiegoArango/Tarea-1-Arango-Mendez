#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "steps.h"
#include "riemann.h"

/*
 C adaptation from C++ code written by Richard J. Gonsalves.
 */


double L = 4.0;                    // length of shock tube
double gama = 1.4;               // ratio of specific heats
int N = 10000;                     // number of grid points

double CFL = 0.4;                // Courant-Friedrichs-Lewy number

double **U = NULL;                      // solution with 3 components
double **F = NULL;                      // flux with 3 components

double h;                        // lattice spacing
double tau;                      // time step

int step;


void allocate() {
    int j;
    
    U = malloc(N * sizeof(double *));
    F = malloc(N * sizeof(double *));
    for (j = 0; j < N; j++) {
        U[j] = malloc(3 * sizeof(double));
        F[j] = malloc(3 * sizeof(double));
    }
}

void initialize() {
    int j;
    double rho,p,u,e;
    allocate();
    h = 1.0 * L / (N - 1);
    for (j = 0; j < N; j++) {
        u = 0;
        
        if (j <= N / 2)
        {
            rho = 1.0;
            p = 1.0;
        }
        
        else if (j > N / 2)
        {
            rho = 0.125;
            p = 0.1;
        }
        
        
        e = p / (gama - 1) + (rho * u * u) /2.0;
        U[j][0] = rho;
        U[j][1] = rho * u;
        U[j][2] = e;
    }
    tau = CFL * h / cMax();
    step = 0;
}


double cMax() {
    double uMax = 0;
    double rho, u, p, c;
    int i;
    for (i = 0; i < N; i++) {
        if (U[i][0] == 0)
            continue;
        
        rho = U[i][0];
        u = U[i][1] / rho;
        p = (U[i][2] - rho * u * u / 2) * (gama - 1);
        c = sqrt(gama * fabs(p) / rho);
        
        if (uMax < (c + fabs(u)))
            uMax = c + fabs(u);
    }
    return uMax;
}


void boundaryConditions(double **U) {
    
    // reflection boundary conditions at the tube ends
    U[0][0] = U[1][0];
    U[0][1] = -U[1][1];
    U[0][2] = U[1][2];
    U[N - 1][0] = U[N - 2][0];
    U[N - 1][1] = -U[N - 2][1];
    U[N - 1][2] = U[N - 2][2];
}

void upwindGodunovStep() {
    int i, j;
    // find fluxes using Riemann solver
    for (j = 0; j < N - 1; j++){
        Riemann(U[j], U[j + 1], F[j]);
    }
    // update U
    
    
    for (j = 1; j < N - 1; j++){
    
        U[j][0] =   U[j][0]-tau*(F[j][0] - F[j - 1][0])/h;
        U[j][1] =   U[j][1]-tau*(F[j][1] - F[j - 1][1])/h;
        U[j][2] =   U[j][2]-tau*(F[j][2] - F[j - 1][2])/h;
        
        
        
    }
    
}



