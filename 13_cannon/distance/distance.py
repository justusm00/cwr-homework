import numpy as np
import matplotlib.pyplot as plt


d = np.loadtxt("distance.dat",usecols = 1)
angle = np.loadtxt("distance.dat",usecols = 0) * 180 / 3.14159265358979323846;



plt.figure(dpi = 400)
plt.grid()
plt.title("cannon ball distance for different initial angles")
plt.ylabel("distance $d$ in m")
plt.xlabel("angle $\Theta$ in °")
plt.plot(angle,d,lw = 0.5)
plt.savefig("distance.pdf")
plt.show()
