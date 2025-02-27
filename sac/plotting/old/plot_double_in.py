import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import sys

# My style preferences
plt.style.use("./style.mplstyle")


# Change to location of output spectra
folder = "../fermion/double_edge/out_files/synth_double_edge/double_edge_in/Nw80/Ac_0.000/p_0.500/Ar_0.500" 
folder = "../fermion/double_edge/out_files/tJ_04_pi8/double_edge_in/Nw80/Ac_0.000/p_0.500/Ar_0.500"

spec = 'd'
n = 0

df_R= pd.read_csv(f"{folder}/{spec}w{n:03}_1.dat")
df_L= pd.read_csv(f"{folder}/{spec}w{n:03}_2.dat")

x_interp = np.linspace(df_R.omega.min(), -df_L.omega.min(), 1000)
R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])

plt.plot(x_interp, R_interp, alpha = 0.5)
plt.plot(x_interp, L_interp, alpha = 0.5)

plt.plot(x_interp, R_interp+L_interp, c= 'k')


plt.ylim(0, 8)
plt.xlim(-10, 10)

plt.ylabel(r"$S(\omega)$")

plt.xlabel(r"$\omega$")

plt.show()