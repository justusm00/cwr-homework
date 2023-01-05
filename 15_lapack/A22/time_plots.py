import numpy as np
import matplotlib.pyplot as plt

n = np.loadtxt("cpu_time.dat",usecols = 0)
t = np.loadtxt("cpu_time.dat",usecols = 1)

plt.figure(dpi = 400)
plt.grid()
plt.title("cpu computation time for different dimensions $n$")
plt.xlabel("$n$")
plt.ylabel("$t$ in s")
plt.plot(n,t,'b.', lw = 0.5)
plt.legend(loc = "upper right")
plt.savefig("cpu_time.pdf")
plt.show()
