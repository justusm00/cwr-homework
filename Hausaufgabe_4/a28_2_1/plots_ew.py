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
plt.plot(n,ew,'bo', markersize = 1)
plt.savefig("ew_anharm.pdf")
plt.show()
