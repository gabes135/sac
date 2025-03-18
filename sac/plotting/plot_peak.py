import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


plt.style.use("./style.mplstyle")

plt.figure(1, (6, 3))
ax = plt.gca()

# Hubbard symmetric
# fname = "../peak/out_files/hubbard_0pi_symm/Np_01/A0_0.700/sw000.csv"
# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='k', lw=1)

# fname = "../peak/out_files/hubbard_0pi/Np_01/A0_0.700/sw000.csv"
# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='b', lw=1)

# # Hubbard asymmetric
# fname = '../in_files/hubbard/aw_03pi4.dat'
# A_df = pd.read_csv(fname)
# # plt.plot(A_df.omega, A_df.S/256, c='k', lw=1)

# fname = "../peak/out_files/hubbard/03pi4/Np_01/A0_0.350/sw000.csv"
# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='b', lw=1)


# fname = "../free/out_files/hubbard/03pi4/sw000.csv"
# A_df = pd.read_csv(fname)
# # plt.plot(A_df.omega, A_df.S, c='k', lw=1)


# fname = "../peak/out_files/hubbard/03pi4/Np_01/A0_0.560_v1/sw000.csv"
# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='r', lw=1)



# Bosonic
# fname = "../process_G/synthetic/bosonic_delta_peak/aw.dat"

# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='k', lw=2)

# fname = "../peak/out_files/bosonic/Np_01/A0_0.700/sw000.csv"
# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='b', lw=1)



# # Fermionic 1
# fname = "../in_files/peak/aw_f1.dat"
# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='k', lw=1, zorder = 3)


# fname = "../peak/out_files/fermionic1_symm/Np_01/A0_0.700/sw000.csv"
# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='b', lw=1)

# fname = "../peak/out_files/fermionic1/Np_01/A0_0.700/sw000.csv"
# A_df = pd.read_csv(fname)
# plt.plot(A_df.omega, A_df.S, c='r', lw=1)


# Fermionic 2
fname = "../in_files/peak/aw_f2.dat"
A_df = pd.read_csv(fname)
plt.plot(A_df.omega, A_df.S, c='k', lw=2, zorder = 1)



fname = "../peak/out_files/fermionic2/Np_01/A0_0.700/sw000_v1.csv"
A_df = pd.read_csv(fname)
plt.plot(A_df.omega, A_df.S, c='b', lw=1)






# plt.text(0.98, 0.95, r'$A(-\omega) = A(\omega)$', ha = 'right', va = 'top', size =12, transform = ax.transAxes)


plt.xlabel(r"$\omega$")
plt.ylabel(r"$A(\omega)$")

plt.ylim(0, 1)

plt.xlim(-3, 3)



plt.savefig("figs/peak_fermionic2.png")

plt.show()