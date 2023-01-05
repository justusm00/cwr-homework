import numpy as np
import matplotlib.pyplot as plt


N = np.loadtxt("N.dat")
ew = np.loadtxt("ew.dat")
ew = np.transpose(ew)





plt.figure(dpi = 200)
plt.grid()
plt.ylabel("$E^*_0$")
plt.xlabel("$N$")
plt.plot(N,ew[0], 'b.', markersize = 3)
plt.savefig("ew_0.pdf")


plt.figure(dpi = 200)
plt.grid()
plt.ylabel("$E^*_1$")
plt.xlabel("$N$")
plt.plot(N,ew[1], 'b.', markersize = 3)
plt.savefig("ew_1.pdf")



plt.figure(dpi = 200)
plt.grid()
plt.ylabel("$E^*_2$")
plt.xlabel("$N$")
plt.plot(N,ew[2], 'b.', markersize = 3)
plt.savefig("ew_2.pdf")



plt.figure(dpi = 200)
plt.grid()
plt.ylabel("$E^*_3$")
plt.xlabel("$N$")
plt.plot(N,ew[3], 'b.', markersize = 3)
plt.savefig("ew_3.pdf")


