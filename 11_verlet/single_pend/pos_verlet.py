import numpy as np
import matplotlib.pyplot as plt

N = 1

t = np.loadtxt("pos_verlet.dat", usecols = 0)
x = np.loadtxt("pos_verlet.dat", usecols = 1) 


def log(x):
	return np.log10(x)




plt.figure(dpi=400)
plt.grid()
plt.xlabel("$\Delta t$")
plt.ylabel("$\Delta x$")
plt.title("residue for different values of $\Delta t$")
plt.loglog(t,x,'b:', lw = 0.5)
plt.savefig("pendulum_verlet.pdf")
plt.show()




