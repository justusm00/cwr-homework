import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("x1000.csv")
mu = data[:,0]
x = data[:,1:]


plt.figure(dpi=400)
plt.grid()
plt.xlabel("$\mu$")
plt.ylabel("$x_{990,..,1000}$")
plt.title("some values of $x$ for different values of $\mu$")
plt.plot(mu,x,'r',linestyle = ':', lw = 0.05)
	
	

plt.savefig("x_1000.png")
plt.show()


