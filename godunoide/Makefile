#
# Makefile for project 'godunov' 
#

CC = gcc
EXECS = analyt godunov
FLAGS = -lm -O2

all: $(EXECS)

analyt: analyt.c nondim.c
	$(CC) $(FLAGS) -o analyt analyt.c

godunov: godunov.c nondim.c solver.c
	$(CC) $(FLAGS) -o godunov godunov.c

.PHONY:

clean:
	rm *~ $(EXECS) *.dat  *.r *.u *.p

arch:
	tar -cvf ../godunov.tar ../godunov
	gzip ../godunov.tar

### end of Makefile ###