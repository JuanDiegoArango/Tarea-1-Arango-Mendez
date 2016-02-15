
import numpy as np
import matplotlib.pyplot as plt
import sys

import matplotlib.mlab as mlab


a=sys.argv[1]
data = np.loadtxt(""+a+"")
plt.xlabel("p")
plt.ylabel("q")
plt.scatter(data[:,3][0:10000], data[:,4][0:10000])


plt.savefig("grafica.pdf")





