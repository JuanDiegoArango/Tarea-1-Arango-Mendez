# compilation:
roe_fds: roe_fds.c
	gcc -lm -O2 -o roe_fds roe_fds.c

# plotting:
plot:
	gnuplot roe_fds.plot

# making "roe_fds.dvi"
dvi:
	latex roe_fds

# making "roe_fds.ps"
ps:
	dvips roe_fds -o roe_fds.ps

# cleaning up:
clean:
	rm -f *.dvi *.aux *.log *~