
import numpy as np
import matplotlib.pyplot as plt
import sys

import matplotlib.mlab as mlab

if len(sys.argv) != 3:
	print "Este programa necesita un (2) argumentos para funcionar"
	sys.exit()

t=sys.argv[2]
datosRK4 = np.loadtxt(""+t+"")
a=sys.argv[1]
datosS = np.loadtxt(""+a+"")


q3RK = datosRK4[:,3]
p3RK = datosRK4[:,4]
q3S = datosS[:,3]
p3S = datosS[:,4]

figura = plt.figure(figsize=(12.0, 17.0))


ax = figura.add_subplot(2, 1, 1)
plt.ylim([-0.4,0.4])
plt.xlim([-0.2, 0.7])
plt.scatter(q3S, p3S,   s=1, facecolor='0.5', lw = 0)
ax.set_title("$\mathrm{Mapa}$"+" "+"$\mathrm{de}$"+" "+"$\mathrm{Poincare}$"+" "+"$\mathrm{Simplectico}$")
ax.set_xlabel("$q3$")
ax.set_ylabel("$p3$")

ax = figura.add_subplot(2, 1, 2)
plt.ylim([-0.4,0.4])
plt.xlim([-0.2, 0.7])
plt.scatter(q3RK, p3RK,   s=1, facecolor='0.5', lw = 0)
ax.set_title("$\mathrm{Mapa}$"+" "+"$\mathrm{de}$"+" "+"$\mathrm{Poincare}$"+" "+"$\mathrm{RK4}$")
ax.set_xlabel("$q3$")
ax.set_ylabel("$p3$")

t = t[:-4]
a = a[:-4]

plt.savefig(t+"_"+a+".pdf")



