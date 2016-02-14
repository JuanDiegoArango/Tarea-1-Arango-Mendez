//
//  Problema de los tres cuerpos.c
//  
//
//  Created by Juan Diego Arango and Juan Felipe Mendez
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>





typedef double (*funcion_aceleracion)(double q);
typedef double (*funcion_posicion_1)(double q1,double epsilon);
typedef double (*funcion_posicion_3)(double q3,double q1,double epsilon );


double ecuacion_diferencial_q1 (double p1,double p2,double p3, double q1,double q2,double q3,double epsilon);
double ecuacion_diferencial_q3 (double p1,double p2,double p3, double q1,double q2,double q3,double epsilon);

double ecuacion_diferencial_p1 (double p1,double p2,double p3, double q1,double q2,double q3,double epsilon);
double ecuacion_diferencial_p3 (double p1,double p2,double p3, double q1,double q2,double q3,double epsilon);


void integrado_simplectico_1 (double dt,double *p1,double *q1, double epsilon, funcion_posicion_1 f1, funcion_aceleracion f2);
void integrado_simplectico_3 (double dt,double *p3,double *q3,double *q1,double epsilon, funcion_posicion_3 f1, funcion_aceleracion f2);
void imprimir (double *p1,double *q1, double num_iteraciones);



int main()


{
    const double  C1= 1.0/(2.0*(2.0-pow(2.0,(1.0/3.0))));
    const double  C2 =(1.0-pow(2.0,(1.0/3.0))) /( 2.0*(2.0-pow(2.0,(1.0/3.0))));
    const double  C3 =(1.0-pow(2.0,(1.0/3.0)))/ (2.0*(2.0-pow(2.0,(1.0/3.0))));
    const double  C4 =1.0/(2.0*(2.0-pow(2.0,(1.0/3.0))));
    const double  D1 =1.0/(2.0-pow(2.0,(1.0/3.0)));
    const double  D2=pow(2.0,(1.0/3.0))/(2-pow(2.0,(1.0/3.0)));
    const double  D3=1.0/(2.0-pow(2.0,(1.0/3.0)));
    const double  D4=0.0;

    double dt=0.006;
    double tiempo_total=3000.0;
    double num_iteraciones=tiempo_total/dt;
    
    double *posicion1;
    double *posicion3;
    double *momento1;
    double *momento3;
    
     posicion1=malloc((int)num_iteraciones * sizeof(double));
     posicion3=malloc((int)num_iteraciones * sizeof(double));
     momento1=malloc((int)num_iteraciones * sizeof(double));
     momento3=malloc((int)num_iteraciones * sizeof(double));
    
    posicion1[0]=0.003;
    posicion3[0]=0.003;
    momento1[0]=0.003;
    momento3[0]=0.003;


    
    for (int i=1; i< num_iteraciones; i=i+1)
    {
        integrado_simplectico_1(dt,posicion1[i],momento1[i], epsilon, ecuacion_diferencial_p1,  ecuacion_diferencial_q1);
        integrado_simplectico_3(dt,posicion3[i],momento3[i],posicion3[i], epsilon,  ecuacion_diferencial_p3,  ecuacion_diferencial_q3);
        
    }
    
    
    
    imprimir (posicion1,posicion3,  num_iteraciones);

    return 0;
}



double ecuacion_diferencial_p1 (double q1,double epsilon )
{
    return -2.0*q1/pow((4*pow(q1,2)+pow(epsilon,2)),(1.5));

}


double ecuacion_diferencial_p3(double q3,double q1,double epsilon )
{
    return (q1-q3)/pow(((q1-q3)*(q1-q3)+epsilon*epsilon/4),(1.5))+ (q1+q3)/pow(((q1+q3)*(q1+q3)+epsilon*epsilon/4),(1.5));
    
}


double ecuacion_diferencial_q1 (double q1)
{
    return q1;
    
}


double ecuacion_diferencial_q3(double q3)

{
    return q3;
    
}




void integrado_simplectico_1 (double dt,double *p1,double *q1, double epsilon, funcion_posicion_1 f1, funcion_aceleracion f2)


{
    q1=q1+C1*f1(q1,epsilon)*dt;
    p1=p1+D1*f2(q1)*dt;
    q1=q1+C2*f1(q1,epsilon)*dt;
    p1=p1+D2*f2(q1)*dt;
    q1=q1+C3*f1(q1,epsilon)*dt;
    p1=p1+D3*f2(q1)*dt;
    q1=q1+C4*f1(q1,epsilon)*dt;
    p1=p1+D4*f2(q1)*dt;


    


}



void integrado_simplectico_3 (double dt,double *p3,double *q3,double *q1,double epsilon, funcion_posicion_3 f1, funcion_aceleracion f2)


{
   
    
    
    q3=q3+C1*f1(q3,q1,epsilon)*dt;
    p3=p3+D1*f2(q3)*dt;
    q3=q3+C2*f1(q3)*dt;
    p3=p3+D2*f2(q3)*dt;
    q3=q3+C3*f1(q3,q1,epsilon)*dt;
    p3=p3+D3*f2(q3)*dt;
    q3=q3+C4*f1(q3,q1,epsilon)*dt;
    p3=p3+D4*f2(q3)*dt;
    
}


void imprimir (double *p1,double *q1, double num_iteraciones)

{   for (int i 0; i<num_iteraciones; i=i+1 )
    printf("%f, %f\n", p1[i],q1[i]);



}
