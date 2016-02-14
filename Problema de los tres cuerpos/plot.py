
import numpy as np
import matplotlib.pyplot as plt
import sys

import matplotlib.mlab as mlab


a=sys.argv[1]
data = np.loadtxt(""+a+"")
plt.xlabel("p")
plt.ylabel("q")
plt.scatter(data[:,0], data[:,1])


plt.savefig("grafica.pdf")





