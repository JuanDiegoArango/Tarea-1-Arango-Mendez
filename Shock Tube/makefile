all: comparacion.pdf UpwindGodunov_step_0.dat UpwindGodunov_step_1.dat UpwindGodunov_step_2.dat UpwindGodunov_step_3.dat UpwindGodunov_step_4.dat UpwindGodunov_step_5.dat

comparacion.pdf: grafica.py a.out
	python grafica.py
	./a.out

UpwindGodunov_step_0.dat:a.out
	./a.out

UpwindGodunov_step_1.dat:a.out
	./a.out

UpwindGodunov_step_2.dat:a.out
	./a.out

UpwindGodunov_step_3.dat:a.out
	./a.out

UpwindGodunov_step_4.dat:a.out
	./a.out

UpwindGodunov_step_5.dat:a.out
	./a.out

a.out: shocktube.c riemann.c steps.c
	cc shocktube.c riemann.c steps.c
