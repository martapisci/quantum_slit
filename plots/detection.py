import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

plt.rcParams.update({'figure.max_open_warning': 0}) #not sure why but needs to be here

# ---------------- IMPORT DATA -----------------

# 1 slit

p1 = pa.cube()
p1.load('../build/data/modulus_nslit1.bin')
p1 = np.array(p1)

# 2 slit

p2 = pa.cube()
p2.load('../build/data/modulus_nslit2.bin')
p2 = np.array(p2)

# 3 slit

p3 = pa.cube()
p3.load('../build/data/modulus_nslit3.bin')
p3 = np.array(p3)

# ---------------- PLOT OPTIONS ----------------

w = 8
h = 6
fs = 18
save_fig = True

def detect_onscreen(data, x_screen, filename, w, h, fs, save_fig):
    plt.figure(figsize=(w, h))
    plt.xlabel("y", fontsize=fs)
    plt.ylabel("p(y|x=0.8; t=0.002)", fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
    y = np.linspace(0,1, 200)
    prob = data[int(x_screen*200),:]
    prob /= np.linalg.norm(prob)
    
    plt.plot(y, prob, 'k')
    plt.savefig(filename) if save_fig else plt.show()

# -------------------- PLOT ---------------------
plist = [p1,p2,p3]
for i in range(3):
    p = plist[i]
    index = int(0.002/0.000025)
    detect_onscreen(p[index], 0.8, f"../build/plots/detection{i+1}.pdf", w, h, fs, save_fig)

    