
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef double (*derivative)(double t, double w, double x, double y, double z, double epsilon);
double func_q1(double t, double w, double x, double y, double z, double epsilon);
double func_p1(double t, double w, double x, double y, double z, double epsilon);
double func_q2(double t, double w, double x, double y, double z, double epsilon);
double func_p2(double t, double w, double x, double y, double z, double epsilon);
void RK4_step(double step, double t, double *w1, double *x1, double *y1, double *z1, double epsilon, derivative func_q1, derivative func_q2, derivative func_p1, derivative func_p2);
void simplectic_step(double step, double t, double *w2, double *x2, double *y2, double *z2, double epsilon, derivative func_q1, derivative func_q2,derivative func_p1, derivative func_p2);


int main(void){
	double q1_RK4, q2_RK4, p1_RK4, p2_RK4, epsilon;
	double q1_S, q2_S, p1_S, p2_S;
	double a=1/(2*sqrt(2));
	double t;
	int i,j;
	
	double T=2800.0;
	double step=0.006;
	int n_step = (int)(T/step);
	int n_runs=83;
	
	FILE *in0;
	char filename1[100];
	sprintf(filename1, "orbitasSimp.dat");
	in0 = fopen(filename1,"w");
	
	if(!in0){
		printf("problems opening the file %s\n", filename1);
		exit(1);
	}
	
	FILE *in1;
	char filename2[100];
	sprintf(filename2, "orbitasRK4.dat");
	in1 = fopen(filename2,"w");
	
	if(!in1){
		printf("problems opening the file %s\n", filename2);
		exit(1);
	}
	

/*******************/
/*******************/


t=0.0;
q1_RK4=a;
q2_RK4=-0.031;
p1_RK4=0.0;
p2_RK4=0.01;
q1_S=q1_RK4;
q2_S=q2_RK4;
p1_S=p1_RK4;
p2_S=p2_RK4;
epsilon=1.0;

for(i=0;i<n_step;i++){
	if (i%100==0) {
		fprintf(in1,"%f %.15e %.15e %.15e %.15e\n", t, q1_RK4, p1_RK4, q2_RK4, p2_RK4);
		fprintf(in0,"%f %.15e %.15e %.15e %.15e\n", t, q1_S, p1_S, q2_S, p2_S);
	}
	
	RK4_step( step,  t, &q1_RK4, &q2_RK4, &p1_RK4, &p2_RK4, epsilon, func_q1, func_q2, func_p1, func_p2);
	simplectic_step( step,  t, &q1_S, &q2_S, &p1_S, &p2_S, epsilon, func_q1, func_q2, func_p1, func_p2);
	t += step;
}
/*******************/
/*******************/



return 0;
}



double func_q1(double t, double w, double x, double y, double z, double epsilon){
	double v;
	v = y;
	return v;
}

double func_p1(double t, double w, double x, double y, double z, double epsilon){
	double v;
	v = -2*w/pow(4*w*w + epsilon*epsilon, 3.0/2.0);
	return v;
}

double func_q2(double t, double w, double x, double y, double z, double epsilon){
	double v;
	v = z;
	return v;
}

double func_p2(double t, double w, double x, double y, double z, double epsilon){
	double v;
	v = (w-x)/pow((w-x)*(w-x) + epsilon*epsilon/4.0, 3.0/2.0)-(w+x)/pow((w+x)*(w+x) + epsilon*epsilon/4.0, 3.0/2.0);
	return v;
}



/*runge kutta 4*/



void RK4_step(double step, double t, double *w1, double *x1, double *y1, double *z1, double epsilon, derivative func_q1, derivative func_q2, derivative func_p1, derivative func_p2){
	
	double k1, k2, k3, k4;
	double l1, l2, l3, l4;
	double m1, m2, m3, m4;
	double n1, n2, n3, n4;
	double q1_in;
	double p1_in;
	double q2_in;
	double p2_in;
	q1_in = *w1;
	q2_in = *x1;
	p1_in = *y1;
	p2_in = *z1;
	
	
	k1 = func_q1( t, q1_in, q2_in, p1_in, p2_in, epsilon);
	l1 = func_q2( t, q1_in, q2_in, p1_in, p2_in, epsilon);
	m1 = func_p1( t, q1_in, q2_in, p1_in, p2_in, epsilon);
	n1 = func_p2( t, q1_in, q2_in, p1_in, p2_in, epsilon);
	
	k2 = func_q1( t+ step*0.5, q1_in+ k1*step*0.5, q2_in+ l1*step*0.5, p1_in+ m1*step*0.5, p2_in+ n1*step*0.5, epsilon);
	l2 = func_q2( t+ step*0.5, q1_in+ k1*step*0.5, q2_in+ l1*step*0.5, p1_in+ m1*step*0.5, p2_in+ n1*step*0.5, epsilon);
	m2 = func_p1( t+ step*0.5, q1_in+ k1*step*0.5, q2_in+ l1*step*0.5, p1_in+ m1*step*0.5, p2_in+ n1*step*0.5, epsilon);
	n2 = func_p2( t+ step*0.5, q1_in+ k1*step*0.5, q2_in+ l1*step*0.5, p1_in+ m1*step*0.5, p2_in+ n1*step*0.5, epsilon);
	
	k3 = func_q1( t+ step*0.5, q1_in+ k2*step*0.5, q2_in+ l2*step*0.5, p1_in+ m2*step*0.5, p2_in+ n2*step*0.5, epsilon);
	l3 = func_q2( t+ step*0.5, q1_in+ k2*step*0.5, q2_in+ l2*step*0.5, p1_in+ m2*step*0.5, p2_in+ n2*step*0.5, epsilon);
	m3 = func_p1( t+ step*0.5, q1_in+ k2*step*0.5, q2_in+ l2*step*0.5, p1_in+ m2*step*0.5, p2_in+ n2*step*0.5, epsilon);
	n3 = func_p2( t+ step*0.5, q1_in+ k2*step*0.5, q2_in+ l2*step*0.5, p1_in+ m2*step*0.5, p2_in+ n2*step*0.5, epsilon);
	
	k4 = func_q1( t+ step, q1_in+ k3*step, q2_in+ l3*step, p1_in+ m3*step, p2_in+ n3*step, epsilon);
	l4 = func_q2( t+ step, q1_in+ k3*step, q2_in+ l3*step, p1_in+ m3*step, p2_in+ n3*step, epsilon);
	m4 = func_p1( t+ step, q1_in+ k3*step, q2_in+ l3*step, p1_in+ m3*step, p2_in+ n3*step, epsilon);
	n4 = func_p2( t+ step, q1_in+ k3*step, q2_in+ l3*step, p1_in+ m3*step, p2_in+ n3*step, epsilon);
	
	q1_in += (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*step;
	q2_in += (l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0)*step;
	p1_in += (m1/6.0 + m2/3.0 + m3/3.0 + m4/6.0)*step;
	p2_in += (n1/6.0 + n2/3.0 + n3/3.0 + n4/6.0)*step;
	
	*w1=q1_in;
	*x1=q2_in;
	*y1=p1_in;
	*z1=p2_in;
}


/*integrador simplectico*/


void simplectic_step(double step, double t, double *w2, double *x2, double *y2, double *z2, double epsilon, derivative func_q1, derivative func_q2,derivative func_p1, derivative func_p2){
	double q1_in;
	double p1_in;
	double q2_in;
	double p2_in;
	double alpha0=-pow(2,1.0/3.0)/(2.0-pow(2,1.0/3.0));
	double alpha1=1.0/(2.0-pow(2,1.0/3.0));
	q1_in = *w2;
	q2_in = *x2;
	p1_in = *y2;
	p2_in = *z2;
	
	
	/*cosa masiva*/
	
	/*1kick*/
	p1_in += 0.5 * func_p1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*1drift*/
	q1_in += 1.0 * func_q1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*1kick*/
	p1_in += 0.5 * func_p1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*2kick*/
	p1_in += 0.5 * func_p1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	/*2drift*/
	q1_in += 1.0 * func_q1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	/*2kick*/
	p1_in += 0.5 * func_p1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	/*3kick*/
	p1_in += 0.5 * func_p1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*3drift*/
	q1_in += 1.0 * func_q1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*3kick*/
	p1_in += 0.5 * func_p1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	
	
	/*cosa chiquita*/
	
	
	/*1kick*/
	p2_in += 0.5 * func_p2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*1drift*/
	q2_in += 1.0 * func_q2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*1kick*/
	p2_in += 0.5 * func_p2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*2kick*/
	p2_in += 0.5 * func_p2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	/*2drift*/
	q2_in += 1.0 * func_q2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	/*2kick*/
	p2_in += 0.5 * func_p2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	/*3kick*/
	p2_in += 0.5 * func_p2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*3drift*/
	q2_in += 1.0 * func_q2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*3kick*/
	p2_in += 0.5 * func_p2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	
	
	*w2=q1_in;
	*x2=q2_in;
	*y2=p1_in;
	*z2=p2_in;
}


