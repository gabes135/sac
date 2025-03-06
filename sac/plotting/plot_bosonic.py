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

fname = "../free/out_files/hchain_beta2048/sw000.csv"
A_df = pd.read_csv(fname)
plt.plot(A_df.omega, A_df.S, c='r', lw=1, label = 'Free')


# fname = "../free/out_files/hchain_beta16/sw000.csv"
# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='b', lw=1)


fname = "../edge/out_files/hchain_beta2048_single/Nw80/Ac_0.000/p_0.500/dw000_1.dat"
A_df = pd.read_csv(fname)
plt.plot(A_df.omega, A_df.S, c='b', lw=1)



plt.rcParams['font.family'] = 'monospace'
plt.text(0.05, 0.85, 'single_edge \& free', ha = 'left', va = 'top', size =15, transform = ax.transAxes)


plt.xlabel(r"$\omega$")
plt.ylabel(r"$A(\omega)$")


plt.xlim(0)
plt.ylim(0, 20)



plt.savefig("figs/bosonic_edge.jpg", bbox_inches='tight')

plt.show()