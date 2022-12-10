import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pyarma as pa

# Import data
p0 = pa.cube()
p0.load('../build/data/modulus_nslit0.bin')  # noslit
p0 = np.array(p0)

p1 = pa.cube()
p1.load('../build/data/modulus_nslit1.bin')  # single-slit
p1 = np.array(p1)

p2 = pa.cube()
p2.load('../build/data/modulus_nslit2.bin')  # double-slit
p2 = np.array(p2)

p3 = pa.cube()
p3.load('../build/data/modulus_nslit3.bin')  # triple-slit
p3 = np.array(p3)

dt = 0.000025
t_points = np.arange(0, 1+dt, dt)
t_min = t_points[0]
x_min, x_max = 0, 1
y_min, y_max = 0, 1

# Some settings
fontsize = 12
w = 10
h = 10
show_colorbar = False
animation_extention = ".gif"

nslit = 0
plist = [p0,p1,p2,p3]
for p in plist:
    # Create figure
    fig = plt.figure(figsize=(w,h))
    ax = plt.gca()
    
    # Create a colour scale normalization according to the max value in the frame
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(p[0]))

    # Plot the first frame
    img = ax.imshow(p[0].T, extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("inferno"), norm=norm, interpolation='Gaussian')

    # Axis labels
    plt.axis('off')

    # Add a colorbar
    if show_colorbar:
        cbar = fig.colorbar(img, ax=ax)
        cbar.set_label("probability", fontsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

    # Function that takes care of updating the data and other things for each frame
    def animation(i):
        # Normalize the colour scale to the current frame?
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(p[i]))
        img.set_norm(norm)

        img.set_data(p[i].T)

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(p), 2), repeat=False, blit=0)

    # # Save the animation
    anim.save(f'../build/plots/p{nslit}'+animation_extention, writer="ffmpeg", bitrate=-1, fps=30)
    nslit+=1
    plt.cla()