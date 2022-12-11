import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

#plt.rcParams.update({'figure.max_open_warning': 0}) #not sure why but needs to be here

# ---------------- IMPORT DATA -----------------

p0 = pa.cube()
p0.load('../build/data/modulus_nslit0.bin')
p0 = np.array(p0)

# ---------------- PLOT OPTIONS ----------------
M = 200
w = 8
h = 6
fs = 18
save_fig = True

for t in range(3):
    plt.figure(figsize=(w, h))
    plt.xlabel("x", fontsize=fs)
    plt.ylabel(f"p(x|y=0.5; t={0.001*t})", fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)

    y = np.linspace(0,1, M)
    prob = p0[int(t*0.001/2.5e-5),:,M//2]
    prob /= np.linalg.norm(prob)

    plt.plot(y[1:], prob[1:], 'k')
    #plt.vlines(0.49, 0, 1, 'grey')
    #plt.vlines(0.51, 0, 1, 'grey')
    plt.axvspan(0.49, 0.51, alpha=0.5, color='grey')
    plt.savefig(f"../build/plots/tunnel_t{t}.pdf") if save_fig else plt.show()    