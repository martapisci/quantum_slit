import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyarma as pa

plt.rcParams.update({'figure.max_open_warning': 0}) #not sure why but needs to be here

# ---------------- IMPORT DATA -----------------

# 1 slit

p1 = pa.cube()
p1.load('../build/data/modulus_nslit1.bin')
p1 = np.array(p1)

re1 = pa.cube()
re1.load('../build/data/re_nslit1.bin')
re1 = np.array(re1)

im1 = pa.cube()
im1.load('../build/data/im_nslit1.bin')
im1 = np.array(im1)

# 2 slit

p2 = pa.cube()
p2.load('../build/data/modulus_nslit2.bin')
p2 = np.array(p2)

re2 = pa.cube()
re2.load('../build/data/re_nslit2.bin')
re2 = np.array(re2)

im2 = pa.cube()
im2.load('../build/data/im_nslit2.bin')
im2 = np.array(im2)

# 3 slit

p3 = pa.cube()
p3.load('../build/data/modulus_nslit3.bin')
p3 = np.array(p3)

re3 = pa.cube()
re3.load('../build/data/re_nslit3.bin')
re3 = np.array(re3)

im3 = pa.cube()
im3.load('../build/data/im_nslit3.bin')
im3 = np.array(im3)

# ---------------- PLOT OPTIONS ----------------

w = 5
h = 4
fs = 12
save_fig = True

def make_colormap(data, filename, label, w, h, fs, save_fig):
    fig = plt.figure(figsize=(w, h))
    ax = plt.gca()
    plt.xlabel("x", fontsize=fs)
    plt.ylabel("y", fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
    norm = mpl.cm.colors.Normalize(vmin=0.0, vmax=np.max(data))
    img = ax.imshow(data.T, extent=[0,1,0,1], cmap=plt.get_cmap("inferno"), norm=norm)

    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(label=label, fontsize=fs)
    cbar.ax.tick_params(labelsize=fs)
    plt.savefig(filename) if save_fig else plt.show()

# -------------------- PLOT ---------------------

counter = 1
for (p, re, im) in zip([p1, p2, p3], [re1,  re2, re3], [im1, im2, im3]):
    for t in range(3):
        i = int(0.001*t/0.000025)
        make_colormap(p[i], f"../build/plots/cmap_p{counter}_t{t*0.001}.pdf", "p(x, y)", w, h, fs, save_fig)
        make_colormap(re[i], f"../build/plots/cmap_re{counter}_t{t*0.001}.pdf", "Re[Ψ](x, y)", w, h, fs, save_fig)
        make_colormap(im[i], f"../build/plots/cmap_im{counter}_t{t*0.001}.pdf", "Im[Ψ](x, y)", w, h, fs, save_fig)
    counter += 1
    