import numpy as np
import matplotlib.pyplot as plt

dx = np.loadtxt("dx.dat")
ew = np.loadtxt("ew.dat")
ew = np.transpose(ew)
print(ew[0])
def w_anal(n):
	return 0.5 * n**2 * np.pi**2

fig, axs = plt.subplots(4, figsize=(12,14), dpi = 200)
plt.subplots_adjust(hspace = 0.5)

for i in range(4):
	axs[i].grid()
	axs[i].set_ylabel("Error $E_{numerical}/E_{analytical}$")
	axs[i].set_xlabel("$\Delta x$")
	axs[i].plot(dx,ew[i] , lw = 1)
	axs[i].set_title("Plot for eigenvalue %g"%(w_anal(i+1)))


fig.savefig("ev.pdf")

