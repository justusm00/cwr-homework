import numpy as np
import matplotlib.pyplot as plt

k = np.loadtxt("modulus1.dat",usecols = 0)
mod_jac = np.loadtxt("modulus1.dat",usecols = 1)
mod_gs = np.loadtxt("modulus1.dat",usecols = 2)

plt.figure(dpi = 400)
plt.grid()
plt.title("modulus of $x_k$ for each iteration $k$")
plt.xlabel("$k$")
plt.ylabel("modulus")
plt.plot(k,mod_jac, 'b',lw = 0.5,label = "jacobi")
plt.plot(k,mod_gs, 'g', lw = 0.5,label = "gauss seidel")
plt.legend(loc = "upper right")
plt.savefig("modulus1.pdf")
plt.show()
