import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mco

data = np.loadtxt("pot.dat", delimiter=",")
plt.imshow(data, aspect='auto',cmap = 'jet')
plt.savefig("potential.pdf")
plt.show()
