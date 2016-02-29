
import numpy as np
import matplotlib.pyplot as plt
import sys

import matplotlib.mlab as mlab

if len(sys.argv) != 2:
	print "Este programa necesita un (1) argumento para funcionar"
	sys.exit()

a=sys.argv[1]
datosS = np.loadtxt(""+a+"")


q3S = datosS[:,3]
p3S = datosS[:,4]

figura = plt.figure()

ax = figura.add_subplot(1, 1, 1)
plt.ylim([-0.1,0.1])
plt.xlim([-0.08, 0.08])
plt.scatter(q3S, p3S,   s=1, facecolor='0.5', lw = 0)
ax.set_title("$\mathrm{Mapa}$"+" "+"$\mathrm{de}$"+" "+"$\mathrm{Poincare}$"+" "+"$\mathrm{Simplectico}$")
ax.set_xlabel("$q3$")
ax.set_ylabel("$p3$")

a = a[:-4]

plt.savefig(a+".pdf")