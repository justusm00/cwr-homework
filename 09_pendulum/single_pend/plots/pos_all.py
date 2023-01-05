import numpy as np
import matplotlib.pyplot as plt

N = 1

t = np.loadtxt("pos_rk2.dat", usecols = 0)
x2 = np.loadtxt("pos_rk2.dat", usecols = 1) ##rk2
x4 = np.loadtxt("pos_rk4.dat", usecols = 1) ##rk4
xe = np.loadtxt("pos_euler.dat", usecols = 1) ##euler



def log(x):
	return np.log10(x)




plt.figure(dpi=400)
plt.grid()
plt.xlabel("$\Delta t$")
plt.ylabel("$\Delta x$")
plt.title("residue for different values of $\Delta t$")
plt.loglog(t,x2,'b:', lw = 0.5,label = "rk2")
plt.loglog(t,x4,'r:', lw = 0.5,label = "rk4")
#plt.loglog(t,xe,'g:', lw = 0.5,label = "Euler")
plt.legend(loc = "upper right")
plt.savefig("pendulum_all.png")
plt.show()




