import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True

A_df = pd.read_csv("../output/sw000.csv")

plt.figure(1, (6, 4), 300)

plt.plot(A_df.omega, A_df.S, c='k', lw=1)
plt.plot(-A_df.omega, A_df.S, c='k', lw=1)

plt.xlabel("$\omega$")
plt.ylabel("$A(\omega)$")

plt.savefig("A_omega.pdf")