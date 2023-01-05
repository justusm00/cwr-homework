import numpy as np
import matplotlib.pyplot as plt

N = 30
t = np.loadtxt("time.dat")
E = np.loadtxt("energy_verlet.dat")

delta_t = np.loadtxt("delta.dat")


plt.figure(dpi = 400)

plt.grid()
plt.ylabel("$E$")
plt.xlabel("$t$")
plt.plot(t, E, 'k' ,lw = 0.5)
plt.title("simulation using verlet method and $\Delta t =%g$ s" %delta_t)
plt.savefig("energy_verlet.pdf")
plt.show()


