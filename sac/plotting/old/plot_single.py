import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import sys

# My style preferences
plt.style.use("./style.mplstyle")


# df = pd.read_csv("../fermion/double_edge/in_files/aw.dat")
# plt.plot(df.omega, df.S, c = 'k', zorder= 1)

# Change to location of output spectra
# folder = "../fermion/double_edge/out_files/synth_symm/double_edge_symm/Nw80/Ac_0.000/p_0.500"
# folder = "../fermion/double_edge/out_files/synth_single_edge/single_edge/Nw80/Ac_0.000/p_0.500"
folder = "../fermion/double_edge/out_files/tJ_04_pi8/single_edge/Nw80/Ac_0.600/p_0.500"


spec = 'd'
n = 0

df= pd.read_csv(f"{folder}/{spec}w{n:03}_1.dat")
plt.plot(df.omega, df.S, c = 'b')


plt.ylim(0, 8)
plt.xlim(-10, 10)

plt.ylabel(r"$S(\omega)$")

plt.xlabel(r"$\omega$")

plt.show()