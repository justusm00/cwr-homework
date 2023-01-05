import numpy as np
import matplotlib.pyplot as plt

dx = np.loadtxt("dx.dat")
ew = np.loadtxt("ew.dat")
ew = np.transpose(ew)



fig, axs = plt.subplots(4, figsize=(12,14), dpi = 200)
plt.subplots_adjust(hspace = 0.5)

for i in range(4):
	axs[i].grid()
	axs[i].set_ylabel("$\Delta E_{rel}(\Delta x_i)$")
	axs[i].set_xlabel("$\Delta x_i$")
	axs[i].plot(dx,ew[i], 'b', lw = 0.5)


fig.savefig("ew_err.pdf")
plt.show()

