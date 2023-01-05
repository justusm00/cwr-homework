import numpy as np
import matplotlib.pyplot as plt

temp = np.loadtxt("heat.dat", skiprows = 1, delimiter = ",")

plt.figure(dpi = 200)
plt.imshow(temp,aspect= 'auto')
plt.savefig("bla.pdf")
plt.show()
