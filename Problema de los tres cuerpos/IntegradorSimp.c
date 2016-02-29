
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef double (*derivative)(double t, double w, double x, double y, double z, double epsilon);
double func_q1(double t, double w, double x, double y, double z, double epsilon);
double func_p1(double t, double w, double x, double y, double z, double epsilon);
double func_q2(double t, double w, double x, double y, double z, double epsilon);
double func_p2(double t, double w, double x, double y, double z, double epsilon);
void simplectic_step(double step, double t, double *w2, double *x2, double *y2, double *z2, double epsilon, derivative func_q1, derivative func_q2,derivative func_p1, derivative func_p2);


int main(int argc, char ** argv){
	double epsilon;
	double q1_S, q2_S, p1_S, p2_S;
	double a=atof(argv[1]);
	double t;
	double dp=0.001;
	int i,j;
	
	double T=2800.0;
	double step=0.006;
	int n_step = (int)(T/step);
	int n_runs=100;
	
	FILE *in0;
	char filename1[100];
	sprintf(filename1, "Simplectico_%s.dat", argv[1]);
	in0 = fopen(filename1,"w");
	
	if(!in0){
		printf("problems opening the file %s\n", filename1);
		exit(1);
	}
	

	srand48(n_runs);
	
	/*******************/
	/***Condiciones iniciales diminutas***/
	/*******************/
	
	for (j=0; j<n_runs/8; j++) {
		t=0.0;
		q1_S=a;
		q2_S=0.1*(2*drand48()-1);
		p1_S=0.0;
		p2_S=0.1*(2*drand48()-1);
		epsilon=1.0;
		
		for(i=0;i<n_step;i++){
			
			if (fabs(p1_S)<dp) {
				fprintf(in0,"%f %.15e %.15e %.15e %.15e\n", t, q1_S, p1_S, q2_S, p2_S);
			}
			
			simplectic_step( step,  t, &q1_S, &q2_S, &p1_S, &p2_S, epsilon, func_q1, func_q2, func_p1, func_p2);
			t += step;
		}
	}
	
	/*******************/
	/*******************/
	
	/*******************/
	/***Condiciones iniciales pequeñas***/
	/*******************/
	
	for (j=0; j<n_runs/8; j++) {
		t=0.0;
		q1_S=a;
		q2_S=0.5*(2*drand48()-1);
		p1_S=0.0;
		p2_S=0.5*(2*drand48()-1);
		epsilon=1.0;
		
		for(i=0;i<n_step;i++){
			
			if (fabs(p1_S)<dp) {
				fprintf(in0,"%f %.15e %.15e %.15e %.15e\n", t, q1_S, p1_S, q2_S, p2_S);
			}
			
			simplectic_step( step,  t, &q1_S, &q2_S, &p1_S, &p2_S, epsilon, func_q1, func_q2, func_p1, func_p2);
			t += step;
		}
	}
	
	/*******************/
	/*******************/
	
	/*******************/
	/***Condiciones iniciales medio grandes***/
	/*******************/
	
		for (j=0; j<n_runs/4; j++) {
			t=0.0;
			q1_S=a;
			q2_S=(2*drand48()-1);
			p1_S=0.0;
			p2_S=(2*drand48()-1);
			epsilon=1.0;
			
			for(i=0;i<n_step;i++){
				
				if (fabs(p1_S)<dp) {
					fprintf(in0,"%f %.15e %.15e %.15e %.15e\n", t, q1_S, p1_S, q2_S, p2_S);
				}
				
				simplectic_step( step,  t, &q1_S, &q2_S, &p1_S, &p2_S, epsilon, func_q1, func_q2, func_p1, func_p2);
				t += step;
			}
		}
	
	/*******************/
	/*******************/
	
	
	/*******************/
	/***Condiciones iniciales grandes***/
	/*******************/
	
			for (j=0; j<n_runs/8; j++) {
				t=0.0;
				q1_S=a;
				q2_S=5*(2*drand48()-1);
				p1_S=0.0;
				p2_S=5*(2*drand48()-1);
				epsilon=1.0;
				
				for(i=0;i<n_step;i++){
					
					if (fabs(p1_S)<dp) {
						fprintf(in0,"%f %.15e %.15e %.15e %.15e\n", t, q1_S, p1_S, q2_S, p2_S);
					}
					
					simplectic_step( step,  t, &q1_S, &q2_S, &p1_S, &p2_S, epsilon, func_q1, func_q2, func_p1, func_p2);
					t += step;
				}
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