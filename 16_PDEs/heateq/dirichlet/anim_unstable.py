from matplotlib import animation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mco


L = 1;
t = np.loadtxt("time2.dat")
temp = np.loadtxt("heat2.dat", skiprows = 1, delimiter = ",")
x = np.linspace(0,1,len(temp[0]))

fig = plt.figure(dpi = 200)

ax = plt.axes()
line, = ax.plot([], [], lw=2)


# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    
    
    ax.clear()
    ax.set_title(f'time ={t[i]} s')
    ax.set_xlim(0,1)
    ax.set_xlabel("position $x$ in m")
    ax.set_ylabel("temperatur $T$ in K")
    ax.plot(x,temp[i],'r',lw = 1)
    
    
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(t)-1, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('heateq_unstable.mp4', fps=20, extra_args=['-vcodec', 'libx264'])

plt.show()
