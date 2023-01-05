import numpy as np
import matplotlib.pyplot as plt

x0 = 0.1
mu = np.array([0.4,0.74,0.77])

t1 = np.loadtxt("pop1.dat", usecols = 0) #i wert
x1 = np.loadtxt("pop1.dat", usecols = 1) #x wert
t2 = np.loadtxt("pop2.dat", usecols = 0) #i wert
x2 = np.loadtxt("pop2.dat", usecols = 1) #x wert
t3 = np.loadtxt("pop3.dat", usecols = 0) #i wert
x3 = np.loadtxt("pop3.dat", usecols = 1) #x wert

##analytical solution

def anal(x0, mu, i):
	return x0 * np.exp(i *(4 * mu - 1)) / (1 + 4 * mu / (4 * mu - 1) * x0 * (np.exp(i *(4 * mu - 1)) - 1))


plt.figure(dpi=400)
plt.grid()
plt.xlabel("i")
plt.ylabel("x")
plt.title("numerical solution for different values of $\mu$ and $x_0=0.1$")
plt.plot(t1,x1, 'k', lw = 0.5, label = "0.4 numerical")
plt.plot(t1,anal(x0,mu[0],t1), 'r' ,lw = 0.5, label = "0.4 analytical")
plt.plot(t2,x2, 'g', lw = 0.5, label = "0.74 numerical")
plt.plot(t2,anal(x0,mu[1],t2), 'b' ,lw = 0.5, label = "0.74 analytical")
plt.plot(t3,x3, 'c', lw = 0.5, label = "0.77 numerical")
plt.plot(t3,anal(x0,mu[2],t3), 'y' ,lw = 0.5, label = "0.77 analytical")

plt.legend(loc = "upper right")
plt.savefig("pop.png")
plt.show()


