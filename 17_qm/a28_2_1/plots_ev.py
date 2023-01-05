import numpy as np
import matplotlib.pyplot as plt

ev = np.loadtxt("ev.dat")
ev = np.transpose(ev)
L = 5
x = np.linspace(-L,L,500)


fig, axs = plt.subplots(3, figsize=(12,14), dpi = 200)
plt.subplots_adjust(hspace = 0.5)


for i in range(3):
	axs[i].grid()
	axs[i].set_ylabel("$|\Psi_%d|^2$"%i)
	axs[i].set_xlabel("$x$")
	axs[i].plot(x, ev[i], lw = 1)




fig.savefig("ev_anharm.pdf")

