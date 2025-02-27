import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import sys

# My style preferences
plt.style.use("./style.mplstyle")


# Change to location of output spectra
folder = "../fermion/double_edge/out_files/synth_symm/tJ_04_pi8/Nw80/Ac_0.000/p_0.500/Ar_0.500"


spec = 'd'
n = 0

df= pd.read_csv(f"{folder}/{spec}w{n:03}.dat")


plt.plot(df.omega, df.S, c= 'k')


plt.ylim(0, 8)
plt.xlim(-10, 10)

plt.ylabel(r"$S(\omega)$")

plt.xlabel(r"$\omega$")

plt.show()