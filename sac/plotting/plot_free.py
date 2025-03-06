import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("./style.mplstyle")

plt.figure(1, (6, 3))

fname = "../free/out_files/hchain_beta2048/sw000.csv"
A_df = pd.read_csv(fname)
plt.plot(A_df.omega, A_df.S, c='k', lw=1)


fname = "../free/out_files/hchain_beta16/sw000.csv"
A_df = pd.read_csv(fname)
plt.plot(A_df.omega, A_df.S, c='b', lw=1)


plt.xlabel(r"$\omega$")
plt.ylabel(r"$A(\omega)$")

# plt.xlim(0)

# plt.savefig("A_omega.pdf")

plt.show()