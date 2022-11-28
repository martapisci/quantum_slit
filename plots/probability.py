import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
#import seaborn as sns

#state = pd.read_csv("../build/data/probability_M10.txt", names= ['t', 'p'], sep=" ")
#state2 = pd.read_csv("../build/data/probability_M20.txt", names= ['t', 'p'], sep=" ")
#state25 = pd.read_csv("../build/data/probability_M25.txt", names= ['t', 'p'], sep=" ")
#state3 = pd.read_csv("../build/data/probability_M30.txt", names= ['t', 'p'], sep=" ")
#state35 = pd.read_csv("../build/data/probability_M35.txt", names= ['t', 'p'], sep=" ")
#state4 = pd.read_csv("../build/data/probability_M40.txt", names= ['t', 'p'], sep=" ")
#state5 = pd.read_csv("../build/data/probability_M50.txt", names= ['t', 'p'], sep=" ")
#state100 = pd.read_csv("../build/data/probability_M100.txt", names= ['t', 'p'], sep=" ")
state200 = pd.read_csv("../build/data/probability_M200_slit.txt", names= ['t', 'p'], sep=" ")
state200_no = pd.read_csv("../build/data/probability_M200_no_slit.txt", names= ['t', 'p'], sep=" ")

#state = state.astype(float)
w = 10
h = 10
fs = 12
cmap = plt.get_cmap('tab20')
save_fig = True

# epsilon(T)
plt.figure(figsize=(w,h))
#plt.plot(state['t'], state['p'], label='M=10')
#plt.plot(state2['t'], state2['p'], label='M=20')
#plt.plot(state25['t'], state25['p'], label='M=25')
#plt.plot(state3['t'], state3['p'], label='M=30')
#plt.plot(state35['t'], state35['p'], label='M=35')
#plt.plot(state4['t'], state4['p'], label='M=40')
#plt.plot(state5['t'], state5['p'], label='M=50')
#plt.plot(state100['t'], state100['p'], label='M=100')
plt.plot(state200['t'], state200['p'], label='slit')
plt.plot(state200_no['t'], state200_no['p'], label='no slit')

plt.xscale('log')
#plt.yscale('log')
plt.legend(loc='best', fontsize = fs)
plt.xlabel("$t$")
plt.ylabel("$\Delta p(t)$")
plt.savefig("../build/plots/p(t).pdf") if save_fig else plt.show()


plt.cla()
