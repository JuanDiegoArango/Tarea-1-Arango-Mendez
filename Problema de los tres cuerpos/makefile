all:IntegradorSimp.x IntegradorSimpKutta.x RK4_0.3535.dat energy.dat Simplectico_0.3535.dat  Simplectico_0.4325.dat  Simplectico_0.425.dat Simplectico_0.46.dat Simplectico_0.3535_RK4_0.3535.pdf Simplectico_0.4325.pdf  Simplectico_0.425.pdf Simplectico_0.46.pdf Simplectico_0.3535_RK4_0.3535zoom1.pdf Simplectico_0.4325zoom1.pdf Simplectico_0.3535_RK4_0.3535zoom2.pdf Simplectico_0.4325zoom2.pdf ENE.pdf

IntegradorSimp.x: IntegradorSimp.c
	gcc IntegradorSimp.c -lm -o IntegradorSimp.x
IntegradorSimpKutta.x: IntegradorSimpKutta.c
	gcc IntegradorSimpKutta.c -lm -o IntegradorSimpKutta.x
Simplectico_0.3535.dat:IntegradorSimpKutta.x
	./IntegradorSimpKutta.x 0.3535
RK4_0.3535.dat:IntegradorSimpKutta.x
	./IntegradorSimpKutta.x 0.3535
energy.dat:IntegradorSimpKutta.x
	./IntegradorSimpKutta.x 0.3535
Simplectico_0.4325.dat:IntegradorSimp.x
	./IntegradorSimp.x 0.4325
Simplectico_0.425.dat:IntegradorSimp.x
	./IntegradorSimp.x 0.425
Simplectico_0.46.dat:IntegradorSimp.x
	./IntegradorSimp.x 0.46
Simplectico_0.3535_RK4_0.3535.pdf:plot2.py
	python plot2.py Simplectico_0.3535.dat RK4_0.3535.dat
Simplectico_0.4325.pdf:plot1.py
	python plot1.py Simplectico_0.4325.dat
Simplectico_0.425.pdf:plot1.py
	python plot1.py Simplectico_0.425.dat
Simplectico_0.46.pdf:plot1.py
	python plot1.py Simplectico_0.46.dat
Simplectico_0.3535_RK4_0.3535zoom1.pdf:plot2zoom1.py
	python plot2zoom1.py Simplectico_0.3535.dat RK4_0.3535.dat
Simplectico_0.4325zoom1.pdf:plot1zoom1.py
	python plot1zoom1.py Simplectico_0.4325.dat
Simplectico_0.3535_RK4_0.3535zoom2.pdf:plot2zoom2.py
	python plot2zoom2.py Simplectico_0.3535.dat RK4_0.3535.dat
Simplectico_0.4325zoom2.pdf:plot1zoom1.py
	python plot1zoom1.py Simplectico_0.4325.dat
ENE.pdf:plotEne.py
	python plotEne.py energy.dat

clean:
	rm -f IntegradorSimp.x IntegradorSimpKutta.x RK4_0.3535.dat energy.dat Simplectico_0.3535.dat  Simplectico_0.4325.dat  Simplectico_0.425.dat Simplectico_0.46.dat Simplectico_0.3535_RK4_0.3535.pdf Simplectico_0.4325.pdf  Simplectico_0.425.pdf Simplectico_0.46.pdf Simplectico_0.3535_RK4_0.3535zoom1.pdf Simplectico_0.4325zoom1.pdf Simplectico_0.3535_RK4_0.3535zoom2.pdf Simplectico_0.4325zoom2.pdf ENE.pdf
