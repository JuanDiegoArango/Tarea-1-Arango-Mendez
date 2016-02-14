
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef double (*derivative)(double t, double w, double x, double y, double z, double epsilon);
double RK4_step(double step, double t, double w, double x, double y, double z, double epsilon, derivative dev);
double func_q1(double t, double w, double x, double y, double z, double epsilon);
double func_p1(double t, double w, double x, double y, double z, double epsilon);
double func_q2(double t, double w, double x, double y, double z, double epsilon);
double func_p2(double t, double w, double x, double y, double z, double epsilon);
void simplectic_step(double step, double t, double *w, double *x, double *y, double *z, double epsilon, derivative fq1, derivative fq2,derivative fp1, derivative fp2);

int main(void){
	double q1_RK4, q2_RK4, p1_RK4, p2_RK4, epsilon;
	double q1_S, q2_S, p1_S, p2_S;
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
	
	t=0.0;
	q1_RK4=0.4325;
	q2_RK4=50;
	p1_RK4=0.0;
	p2_RK4=10;
	q1_S=0.4325;
	q2_S=0.3;
	p1_S=0.0;
	p2_S=0.4;
	epsilon=1.0;
	
	for(i=0;i<n_step;i++){
		
		fprintf(in1,"%f %.15e %.15e %.15e %.15e\n", t, q1_RK4, p1_RK4, q2_RK4, p2_RK4);
		fprintf(in0,"%f %.15e %.15e %.15e %.15e\n", t, q1_S, p1_S, q2_S, p2_S);
		q1_RK4 += RK4_step(step, t, q1_RK4, q2_RK4, p1_RK4, p2_RK4, epsilon, func_q1);
		q2_RK4 += RK4_step(step, t, q1_RK4, q2_RK4, p1_RK4, p2_RK4, epsilon, func_q2);
		p1_RK4 += RK4_step(step, t, q1_RK4, q2_RK4, p1_RK4, p2_RK4, epsilon, func_p1);
		p2_RK4 += RK4_step(step, t, q1_RK4, q2_RK4, p1_RK4, p2_RK4, epsilon, func_p2);
		simplectic_step( step,  t, &q1_S, &q2_S, &p1_S, &p2_S, epsilon, func_q1, func_q2, func_p1, func_p2);
		t += step;
	}
	
	
	
	return 0;
}


double RK4_step(double step, double t, double w, double x, double y, double z, double epsilon, derivative dev){
	double k1, k2, k3, k4;
	double y_step;
	k1 = dev( t, w, x, y, z, epsilon);
	k2 = dev( t+ step*0.5, w+ step*0.5, x+ step*0.5, y+ step*0.5, z+ step*0.5, epsilon);
	k3 = dev( t+ step*0.5, w+ step*0.5, x+ step*0.5, y+ step*0.5, z+ step*0.5,epsilon);
	k4 = dev( t+ k3*step, w+ k3*step, x+ k3*step, y+ k3*step, z+ k3*step,epsilon);
	y_step = (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*step;
	return y_step;
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

void simplectic_step(double step, double t, double *w, double *x, double *y, double *z, double epsilon, derivative fq1, derivative fq2,derivative fp1, derivative fp2){
	double q1_in;
	double p1_in;
	double q2_in;
	double p2_in;
	double alpha1=-pow(2,1.0/3.0)/(2.0-pow(2,1.0/3.0));
	double alpha0=1.0/(2.0-pow(2,1.0/3.0));
	q1_in = *w;
	q2_in = *x;
	p1_in = *y;
	p2_in = *z;
	
	/*1kick*/
	p1_in += 0.5 * fp1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	p2_in += 0.5 * fp2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*1drift*/
	q1_in += 1.0 * fq1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	q2_in += 1.0 * fq2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*1kick*/
	p1_in += 0.5 * fp1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	p2_in += 0.5 * fp2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*2kick*/
	p1_in += 0.5 * fp1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	p2_in += 0.5 * fp2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	/*2drift*/
	q1_in += 1.0 * fq1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	q2_in += 1.0 * fq2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	/*2kick*/
	p1_in += 0.5 * fp1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	p2_in += 0.5 * fp2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha0 * step;
	/*3kick*/
	p1_in += 0.5 * fp1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	p2_in += 0.5 * fp2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*3drift*/
	q1_in += 1.0 * fq1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	q2_in += 1.0 * fq2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	/*3kick*/
	p1_in += 0.5 * fp1( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	p2_in += 0.5 * fp2( t, q1_in, q2_in, p1_in, p2_in, epsilon) * alpha1 * step;
	
	
	*w=q1_in;
	*x=q2_in;
	*y=p1_in;
	*z=p2_in;
}
