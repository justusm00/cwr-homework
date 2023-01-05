import numpy as np
import matplotlib.pyplot as plt


ew = np.loadtxt("ew.dat")
ew = np.transpose(ew)
n = np.zeros(20)

for i in range(20):
	n[i] = i
	

plt.figure(dpi = 400)
plt.grid()
plt.xlabel("$n$")
plt.ylabel("Eigenvalue $E_n$")
plt.plot(n,ew[0],'bo', label = "numerical", markersize = 1)
plt.plot(n,ew[1],'ro',  label = "analytical", markersize = 1)
plt.legend(loc = "upper left")
plt.savefig("ew_harm.pdf")
plt.show()
