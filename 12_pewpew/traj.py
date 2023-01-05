import numpy as np
import matplotlib.pyplot as plt


t = np.loadtxt("traj1.dat",usecols = 0)
r1 = np.loadtxt("traj1.dat", usecols = 1)
r2 = np.loadtxt("traj2.dat", usecols = 1)
r3 = np.loadtxt("traj3.dat", usecols = 1)



fig, axs = plt.subplots(3, figsize=(10,12), dpi = 400)
plt.subplots_adjust(hspace = 0.5)


axs[0].grid()
axs[0].set_ylabel("$r$")
axs[0].set_xlabel("$t$")
axs[0].plot(t, r1, 'k' ,lw = 0.5)
axs[0].set_title("$v_0=-5$m/s")

axs[1].grid()
axs[1].set_ylabel("$r$")
axs[1].set_xlabel("$t$")
axs[1].plot(t, r2, 'k' ,lw = 0.5)
axs[1].set_title("$v_0=-4$m/s")

axs[2].grid()
axs[2].set_ylabel("$r$")
axs[2].set_xlabel("$t$")
axs[2].plot(t, r3, 'k' ,lw = 0.5)
axs[2].set_title("numerically computed value for $v_0$")



fig.suptitle("simulations for different initial velocities")
fig.savefig("traj.pdf")
plt.show()



