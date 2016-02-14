//
//  Problema de los tres cuerpos.c
//  
//
//  Created by Juan Diego Arango on 2/13/16.
//
//

#include <stdio.h>
#include <math.h>

#define C1 1/(2(2-2^(1/3)))
#define C2 (1-2^(1/3))/(2(2-2^(1/3)))
#define C3 (1-2^(1/3))/(2(2-2^(1/3)))
#define C4 1/(2(2-2^(1/3)))
#define D1 1/(2-2^(1/3))
#define D2 2^(1/3)/(2-2^(1/3))
#define D3 1/(2-2^(1/3))
#define D4 0



typedef double (*funcion)(double t, double y);

float ecuacion_diferencial_1 ( float p1,float p2,float p3, float q1,float q2,float q3,float epsilon);

float ecuacion_diferencial_3 ( float p1,float p2,float p3, float q1,float q2,float q3,float epsilon);

void main()


{
    double dt;
    double tiempo_total;
    double num_iteraciones=tiempo_total/dt;
    
    float *posicion1;
    float *posicion2;
    float *posicion3;
    
    float *momento1;
    float *momento2;
    float *momento3;

}



float ecuacion_diferencial_1 ( float p1,float p2,float p3, float q1,float q2,float q3,float epsilon )
{
    return -2.0*q1/(4q1^2+epsilon^2)^(1.5);

}


float ecuacion_diferencial_3( float p1,float p2,float p3, float q1,float q2,float q3,float epsilon
 )

{
    return (q1-q3)/((q1-q3)^2+epsilon^2/4)^(1.5)+(q1+q3)/((q1+q3)^2+epsilon^2/4)^(1.5);
    
}



void integrado_simplectico (float p,float q, funcion f)


{



}

