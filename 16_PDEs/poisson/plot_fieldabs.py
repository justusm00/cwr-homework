import numpy as np
import matplotlib.pyplot as plt

E = np.loadtxt("e.dat", delimiter=",")

plt.figure(dpi=200)
plt.title("$E$-Feld in xy-Ebene")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.imshow(E, aspect='auto',cmap = 'jet')
plt.savefig("fieldabs.pdf")
plt.show()
