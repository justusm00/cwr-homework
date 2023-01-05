import numpy as np
import matplotlib.pyplot as plt

N = 30
t = np.loadtxt("time.dat")
E_eu = np.loadtxt("euler_energy.dat")
E_rk2 = np.loadtxt("rk2_energy.dat")
E_rk4 = np.loadtxt("rk4_energy.dat")
delta_t = np.loadtxt("delta.dat")



fig, axs = plt.subplots(3, figsize=(10,12), dpi = 400)
plt.subplots_adjust(hspace = 0.5)


axs[0].grid()
axs[0].set_ylabel("$x$")
axs[0].set_xlabel("$t$")
axs[0].plot(t, E_eu, 'k' ,lw = 0.5)
axs[0].set_title("Euler")

axs[1].grid()
axs[1].set_ylabel("$x$")
axs[1].set_xlabel("$t$")
axs[1].plot(t, E_rk2, 'k',lw = 0.5)
axs[1].set_title("RK2")

axs[2].grid()
axs[2].set_ylabel("$x$")
axs[2].set_xlabel("$t$")
axs[2].plot(t, E_rk4,'k', lw = 0.5)
axs[2].set_title("RK4")

fig.suptitle("energy for $\Delta t=%g$ s" %delta_t)
fig.savefig("mult_pend_energy.pdf")



