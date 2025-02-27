import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import sys
from matplotlib.font_manager import FontProperties


# My style preferences
plt.style.use("./style.mplstyle")

plt.figure(1, (6, 3))

ax = plt.gca()

df = pd.read_csv("../fermion/in_files/edge_modes/aw2.dat")
plt.plot(df.omega, df.S, c = 'k', zorder= 1, lw = 2)#, label = 'Artificial Spectrum')

# Change to location of output spectrum
folder = "../fermion/edge/out_files/t2_double_symm/Nw80/Ac_0.000/p_0.500"

spec = 'd'

n = 0

df= pd.read_csv(f"{folder}/{spec}w{n:03}_1.dat")#, label = '')
plt.plot(df.omega, df.S, c='b')
plt.plot(-df.omega, df.S, c='b')

plt.text(0.98, 0.95, r'$A(-\omega) = A(\omega)$', ha = 'right', va = 'top', size =12, transform = ax.transAxes)


plt.rcParams['font.family'] = 'monospace'
plt.text(0.05, 0.85, 'double_edge_symm', ha = 'left', va = 'top', size =15, transform = ax.transAxes)

plt.ylim(0, 8)
plt.xlim(-10, 10)

plt.ylabel(r"$A(\omega)$")
plt.xlabel(r"$\omega$")

plt.savefig("figs/t2.jpg", bbox_inches='tight')

plt.show()