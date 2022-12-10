import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns
#state0 = pd.read_csv("../build/data/probability_M5_no_slit.txt", names= ['t', 'p'], sep=" ")
#state = pd.read_csv("../build/data/probability_M10_no_slit.txt", names= ['t', 'p'], sep=" ")
#state2 = pd.read_csv("../build/data/probability_M20_no_slit.txt", names= ['t', 'p'], sep=" ")
#state25 = pd.read_csv("../build/data/probability_M25.txt", names= ['t', 'p'], sep=" ")
#state3 = pd.read_csv("../build/data/probability_M30.txt", names= ['t', 'p'], sep=" ")
#state35 = pd.read_csv("../build/data/probability_M35.txt", names= ['t', 'p'], sep=" ")
#state4 = pd.read_csv("../build/data/probability_M40.txt", names= ['t', 'p'], sep=" ")
#state5 = pd.read_csv("../build/data/probability_M50.txt", names= ['t', 'p'], sep=" ")
#state100_no = pd.read_csv("../build/data/probability_M100_no_slit.txt", names= ['t', 'p'], sep=" ")
#state200 = pd.read_csv("../build/data/probability_M200_slit.txt", names= ['t', 'p'], sep=" ")
#state200_no = pd.read_csv("../build/data/probability_M200_no_slit.txt", names= ['t', 'p'], sep=" ")
psi0 = pd.read_csv("../build/data/probability_M200_nslit0.txt", names= ['t', 'p'], sep=" ")
psi1 = pd.read_csv("../build/data/probability_M200_nslit1.txt", names= ['t', 'p'], sep=" ")
psi2 = pd.read_csv("../build/data/probability_M200_nslit2.txt", names= ['t', 'p'], sep=" ")
psi3 = pd.read_csv("../build/data/probability_M200_nslit3.txt", names= ['t', 'p'], sep=" ")

#state = state.astype(float)
w = 5
h = 4
fs = 12
cmap = plt.get_cmap('tab20')
save_fig = True

# epsilon(T)
plt.figure(figsize=(w,h))
#plt.plot(state0['t'], state0['p'], label='M=5')
#plt.plot(state['t'], state['p'], label='M=10')
#plt.plot(state2['t'], state2['p'], label='M=20')
#plt.plot(state25['t'], state25['p'], label='M=25')
#plt.plot(state3['t'], state3['p'], label='M=30')
#plt.plot(state35['t'], state35['p'], label='M=35')
#plt.plot(state4['t'], state4['p'], label='M=40')
#plt.plot(state5['t'], state5['p'], label='M=50')
#plt.plot(state100_no['t'], state100_no['p'], label='M=100')
#plt.plot(state200['t'], state200['p'], label='slit')
#plt.plot(state200_no['t'], np.abs(state200_no['p']), '.', markersize=2, label='M=200, no slit')
#plt.plot(state200['t'], np.abs(state200['p']), '.', markersize=2, label='M=200, double slit')
plt.plot(psi0['t'], np.abs(psi0['p']), '.', markersize=2, label='no slit')
plt.plot(psi1['t'], np.abs(psi1['p']), '.', markersize=2, label='1 slit')
plt.plot(psi2['t'], np.abs(psi2['p']), '.', markersize=2, label='2 slits')
plt.plot(psi3['t'], np.abs(psi3['p']), '.', markersize=2, label='3 slits')

plt.xscale('log')
#plt.yscale('log')
plt.legend(loc='best', fontsize = fs)
plt.xlabel("$t$")
plt.ylabel("$\Delta p(t)$")
plt.savefig("../build/plots/p(t).pdf") if save_fig else plt.show()


plt.cla()
