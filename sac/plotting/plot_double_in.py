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

# df = pd.read_csv("../fermion/in_files/edge_modes/aw3.dat")
# plt.plot(df.omega, df.S, c = 'k', zorder= 1, lw = 2)#, label = 'Artificial Spectrum')

# Change to location of output spectrum
# folder = "../edge/out_files/t3_double_in/Nw80/Ac_0.000/p_0.500/Ar_0.500"

folder = "../edge/out_files/1D_tJ/04_pi2_double_in/Nw80/Ac_0.000/p_0.500/Ar_0.500"
folder = "../edge/out_files/1D_tJ/04_pi8_double_in/Nw_200/Ac_0.000/p_0.500/Ar_0.600"

spec = 'd'
n = 0

df_R= pd.read_csv(f"{folder}/{spec}w{n:03}_1.dat")
df_L= pd.read_csv(f"{folder}/{spec}w{n:03}_2.dat")

x_interp = np.linspace(df_R.omega.min(), -df_L.omega.min(), 1000)
R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])

# plt.plot(x_interp, R_interp, alpha = 0.5, c='r')
# plt.plot(x_interp, L_interp, alpha = 0.5, c='g')

plt.plot(x_interp, R_interp+L_interp, c= 'b')


plt.rcParams['font.family'] = 'monospace'
# plt.text(0.5, 0.85, 'double_edge_in', ha = 'center', va = 'top', size =15, transform = ax.transAxes)

plt.ylim(0, 8)
plt.xlim(-2, 3)

plt.ylabel(r"$A(\omega)$")
plt.xlabel(r"$\omega$")

plt.savefig("figs/1D_tJ_04_pi8_Ar6.pdf", dpi = 300, bbox_inches='tight')

plt.show()