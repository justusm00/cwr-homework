import numpy as np
import matplotlib.pyplot as plt

ev1 = np.loadtxt("ev.dat",usecols = 0)
ev2 = np.loadtxt("ev.dat",usecols = 1)
ev3 = np.loadtxt("ev.dat",usecols = 2)

x = np.linspace(0,1,1000)


fig, axs = plt.subplots(3, figsize=(12,14), dpi = 200)
plt.subplots_adjust(hspace = 0.5)


axs[0].grid()
axs[0].set_ylabel("$\Psi$")
axs[0].set_xlabel("x")
axs[0].plot(x, ev1, lw = 1)
axs[0].set_title("plot of $\Psi_1$")

axs[1].grid()
axs[1].set_ylabel("$\Psi$")
axs[1].set_xlabel("x")
axs[1].plot(x, ev2, lw = 1)
axs[1].set_title("plot of $\Psi_2$")

axs[2].grid()
axs[2].set_ylabel("$\Psi$")
axs[2].set_xlabel("x")
axs[2].plot(x, ev3, lw = 1)
axs[2].set_title("plot of $\Psi_3$")

fig.savefig("ev.pdf")

