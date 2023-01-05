import numpy as np
import matplotlib.pyplot as plt


x = np.loadtxt("pos.dat")
x = np.transpose(x)
print(x)
alpha = np.loadtxt("angles.dat")





plt.figure(dpi = 400)
plt.grid()
plt.title("Trajectories for different initial angles")
plt.xlabel("$x$ in m")
plt.xlabel("$y$ in m")
for i in range(10):
	plt.plot(x[2*i],x[2*i+1],lw = 0.5, label = "%g °"%(alpha[i] * 180 / np.pi))
plt.legend(loc = "upper right")
plt.savefig("traj.pdf")
plt.show()



