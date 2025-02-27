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

df = pd.read_csv("/Users/gabeschumm/BU/Research/research_code/sac_repo/sac/process_G/synthetic/double_edge_out/aw.dat")
plt.plot(df.omega, df.S, c = 'k', zorder= 1, lw = 2)#, label = 'Artificial Spectrum')

# Change to location of output spectrum
folder = "../fermion/edge/out_files/t4_double_out/Nw80/Ac_0.000/p_0.500/Ar_0.700"

spec = 'd'
n = 0

df_R= pd.read_csv(f"{folder}/{spec}w{n:03}_1.dat")
df_L= pd.read_csv(f"{folder}/{spec}w{n:03}_2.dat")


plt.plot(df_R.omega, df_R.S, c='b')
plt.plot(-df_L.omega, df_L.S, c='b')



plt.rcParams['font.family'] = 'monospace'
plt.text(0.05, 0.85, 'double_edge_out', ha = 'left', va = 'top', size =15, transform = ax.transAxes)

plt.ylim(0, 8)
plt.xlim(-10, 10)

plt.ylabel(r"$A(\omega)$")
plt.xlabel(r"$\omega$")

# plt.savefig("figs/t4.pdf", bbox_inches='tight')

plt.show()