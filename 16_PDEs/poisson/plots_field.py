import numpy as np
import matplotlib.pyplot as plt

ex = np.loadtxt("ex.dat", delimiter=",")
ey = np.loadtxt("ey.dat", delimiter=",")

grid = 100
x = np.linspace(-0.005,0.005,grid)
y = np.linspace(-0.005,0.005,grid)
rx,ry = np.meshgrid(x,y)


plt.figure(dpi=200)
plt.title("$E$-Feld in xy-Ebene")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.streamplot(rx,ry,ex,ey)
plt.savefig("field.pdf")
plt.show()
