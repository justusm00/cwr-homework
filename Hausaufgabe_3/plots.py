import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("pos.dat")
data = np.transpose(data)


w = 1
t = 0.5
n = 10
def anal(a):
	return w + 2 * t * np.cos(a * np.pi /(n+1))

plt.figure(dpi = 400)
plt.grid()
plt.title("characteristic polynomial")
plt.xlabel("$\lambda$")
plt.ylabel("$P(\lambda)$")
plt.plot(data[0], data[1], lw = 0.5)
for i in range(10):
	plt.vlines(anal(i+1),np.min(data[1]),np.max(data[1]),'k', lw = 0.5)
plt.savefig("char.pdf")
plt.show()


