import numpy as np
import matplotlib.pyplot as plt

v = np.loadtxt("G.dat", usecols = 0)
G = np.loadtxt("G.dat", usecols = 1)


plt.figure(dpi=400)
plt.grid()
plt.xlabel("$v$")
plt.ylabel("$G(v)$")
plt.title("Plot of $G(v)$")
plt.plot(v,G,'b', lw = 0.5)
plt.savefig("zeros.pdf")
plt.show()



