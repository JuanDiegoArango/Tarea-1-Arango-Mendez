
import numpy as np
import matplotlib.pyplot as plt
import sys

import matplotlib.mlab as mlab

if len(sys.argv) != 2:
	print "Este programa necesita un (1) argumento para funcionar"
	sys.exit()

a=sys.argv[1]
datosS = np.loadtxt(""+a+"")


q3S = datosS[:,0]
p3S = datosS[:,1]

figura = plt.figure()

ax = figura.add_subplot(1, 1, 1)
plt.ylim([-3,3])
plt.xlim([-2.5, 2.5])
plt.scatter(q3S, p3S,   s=1, facecolor='0.5', lw = 0)
ax.set_title("$\mathrm{Energy}$")
ax.set_xlabel("$E$")
ax.set_ylabel("$T$")

a = a[:-4]

plt.savefig("ENE.pdf")
