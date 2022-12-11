import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

#plt.rcParams.update({'figure.max_open_warning': 0}) #not sure why but needs to be here

# ---------------- IMPORT DATA -----------------

p0 = pa.cube()
p0.load('../build/data/modulus_nslit0.bin')
p0 = np.array(p0)

# ---------------- PLOT OPTIONS ----------------
Nt, M, _ = p0.shape
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

    p_tunnel = max(prob[M//2+20:])
    x_tunnel = np.where(prob == p_tunnel)[0][0]/200.

    print(x_tunnel)
    plt.plot(y[1:], prob[1:], 'k')
    plt.axvspan(0.49, 0.51, alpha=0.5, color='grey')
    
    plt.text(0.8, 0.05, f"p(x={x_tunnel:.2g}) = {p_tunnel:.2g}", color="black", horizontalalignment="center", verticalalignment="center", fontsize=fs)
    
    
    plt.savefig(f"../build/plots/tunnel_t{t}.pdf") if save_fig else plt.show()    