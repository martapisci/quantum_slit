import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

psi0 = pd.read_csv("../build/data/probability_M200_nslit0.txt", names= ['t', 'p'], sep=" ")
#psi1 = pd.read_csv("../build/data/probability_M200_nslit1.txt", names= ['t', 'p'], sep=" ")
psi2 = pd.read_csv("../build/data/probability_M200_nslit2.txt", names= ['t', 'p'], sep=" ")
#psi3 = pd.read_csv("../build/data/probability_M200_nslit3.txt", names= ['t', 'p'], sep=" ")

w = 5
h = 4
fs = 12
save_fig = True

# epsilon(T)
plt.figure(figsize=(w,h))
plt.plot(psi0['t'], np.abs(psi0['p']), '.', markersize=2, label='no slit')
#plt.plot(psi1['t'], np.abs(psi1['p']), '.', markersize=2, label='1 slit')
plt.plot(psi2['t'], np.abs(psi2['p']), '.', markersize=2, label='2 slits')
#plt.plot(psi3['t'], np.abs(psi3['p']), '.', markersize=2, label='3 slits')

plt.xscale('log')
plt.legend(loc='best', fontsize = fs)
plt.xlabel("$t$")
plt.ylabel("$\Delta p(t)$")
plt.savefig("../build/plots/probability_conservation.pdf") if save_fig else plt.show()