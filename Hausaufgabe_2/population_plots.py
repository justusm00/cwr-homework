import numpy as np
import matplotlib.pyplot as plt



t1 = np.loadtxt("pop1.dat", usecols = 0) #i wert
x1 = np.loadtxt("pop1.dat", usecols = 1) #x wert
t2 = np.loadtxt("pop2.dat", usecols = 0) #i wert
x2 = np.loadtxt("pop2.dat", usecols = 1) #x wert
t3 = np.loadtxt("pop3.dat", usecols = 0) #i wert
x3 = np.loadtxt("pop3.dat", usecols = 1) #x wert




plt.figure(dpi=400)
plt.grid()
plt.xlabel("i")
plt.ylabel("x")
plt.title("numerical solution for different values of $\mu$ and $x_0=0.1$")
plt.plot(t1,x1, lw = 0.5, label = "0.4")
plt.plot(t2,x2, lw = 0.5, label = "0.74")
plt.plot(t3,x3, lw = 0.5, label = "0.77")
plt.legend(loc = "upper right")
plt.savefig("pop.png")
plt.show()


