#!/usr/bin/env python
# coding: utf-8

# In[208]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sc



plt.style.use("~/style.mplstyle")
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 18


# # Real Data

# # $\pi/2$

# In[2]:


plt.figure(1, (8, 4), 200)

folder = "tJ_1D_J04/pi2"
Ar = 0.5

df = pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/free/sw000.dat")
# plt.plot(df.omega, df.S, label = 'Unconstrained', alpha = 1)


# df = pd.read_csv(f"sw000.dat", sep = r'\s+', header=None)
# plt.plot(df[0], df[1], label = 'Sibin 0000', alpha = 1)

# df = pd.read_csv(f"sw001.dat", sep = r'\s+', header=None)
# plt.plot(df[0], df[1], label = 'Sibin 00001', alpha = 1)


Nw = 160
n = 0
spec = 'd'
df_R= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_1.dat")
# plt.plot(df_R.omega, df_R.S, alpha = 0.5, label = r"Rightward $\delta$'s")

df_L= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_2.dat")
# plt.plot(-df_L.omega, df_L.S, alpha = 0.5, label = r"Leftward $\delta$'s")

x_interp = np.linspace(-10, 10, 1000)
R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])

plt.plot(x_interp, R_interp+L_interp, label = fr"Double Edge, $a_R = {Ar:.2f}$", c= 'b')

# plt.plo
plt.legend()
plt.ylim(0, 8)
plt.xlim(-4, 4)

plt.ylabel(r"$S(\omega)$")

plt.xlabel(r"$\omega$")
# plt.title("1D $t-J$ model, $J/t = 0.4$")


# ### Ar Scan

# In[3]:


Nw = 160
folder = "tJ_1D_J04/pi2"
N_tau = 80
A_range =  np.arange(.01, .9, .01)#np.concatenate([[.1, .2, .3], np.arange(.4, .901, .02)])

N_anneal = 200
scan = np.zeros((N_anneal, len(A_range), 2))
scan[:, :, :] = np.nan
for a, A in enumerate(A_range):
    
    fname = f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/p_0.500/Ar_{A:.3f}/anneal.csv"
    df = pd.read_csv(fname).iloc[1:, :]

    i_max = df.index[df.theta.diff() > 0][0]-1
  
        
    df = df.loc[:i_max, :]
    scan[:df.shape[0], a, :] = df[['chi2_avg', 'chi2_min']].values
    # except:
    #     pass

chi2_min = np.nanmin(scan[:, :, 1])
a = 0.5
chi2_target = chi2_min + a * np.sqrt(2*chi2_min/N_tau)
print(chi2_min, chi2_target)

i_target = np.zeros_like(A_range, dtype = int)
for a, A in enumerate(A_range):
    i_target[a] = np.arange(N_anneal)[scan[:, a, 0] > chi2_target][-1]+1
min(i_target)


# In[4]:


plt.figure(1, (8, 4))
i_plot = min(i_target)
plt.plot(A_range, scan[i_plot // 2, :, 0], marker='o')
plt.plot(A_range, scan[i_plot, :, 0], marker='o')

plt.axhline(chi2_target, c= 'k', ls ='--')
plt.ylim(chi2_target*.95, chi2_target*1.8);


# # $\pi/8$

# In[263]:


plt.figure(1, (6, 4), 200)

folder = "tJ_1D_J04/pi8"

# df = pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/free/sw000.dat")
# plt.plot(df.omega, df.S, label = 'Unconstrained', alpha = 1)


# df = pd.read_csv(f"sw000.dat", sep = r'\s+', header=None)
# plt.plot(df[0], df[1], label = 'Sibin 0000', alpha = 1)



Nw = 200
n = 0


Ac = 0.7
if Ac > 0:
    Ac_folder = f"Ac_{Ac:.3f}/"
else:
    Ac_folder = ""
Ar = 0.5


spec = 'd'

fig, ax = plt.subplots(2, 1, figsize = (8, 8), dpi = 100)

df = pd.read_csv(f"sw001.dat", sep = r'\s+', header=None)
ax[0].plot(df[0], df[1], label = 'Unconstrained', alpha = 1)
# ax[1].plot(df[0], df[1], label = 'Unconstrained', alpha = 1)



df_R= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/{Ac_folder}p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_1.dat")
ax[1].plot(df_R.omega, df_R.S_edge, alpha = 1, label = r"Left Edges")
ax[1].plot(df_R.omega, df_R.S_cont, alpha = 1, label = r"Continuum")


df_L= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/{Ac_folder}p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_2.dat")
ax[1].plot(-df_L.omega, df_L.S_edge, alpha = 1, label = r"Right Edge")

x_interp = np.linspace(-5, 5, 10000)
R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])

ax[0].plot(x_interp, R_interp+L_interp, label = r"$A_{\mathrm{left}} = " +rf"{0.5:.2f}$," + fr"$A_c = {Ac:.2f}$")
ax[0].plot(x_interp, spec_opt, label = r"$A_{\mathrm{left}} = " +rf"{0.6:.2f}$," + fr" $A_c = 0$" )

ax[1].plot(x_interp, R_interp+L_interp, label = "Combined", alpha = 1)


# # plt.plo
for a in ax:
    
    a.legend()
    a.set_ylim(0, 5)
    a.set_xlim(-2, 2)

    a.set_ylabel(r"$S(\omega)$")

    a.set_xlabel(r"$\omega$")

plt.savefig(f"figs/tJ_1D/pi8_Ac_{Ac:.2f}.pdf")


# ## Ar Scan

# In[247]:


Nw = 200
N_tau = 55
Ar_range =  np.arange(.01, .9, .01)#np.concatenate([[.1, .2, .3], np.arange(.4, .901, .02)])

N_anneal = 200
scan = np.zeros((N_anneal, len(Ar_range), 3))
scan[:, :, :] = np.nan
for a, A in enumerate(Ar_range):
    
    fname = f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/p_0.500/Ar_{A:.3f}/anneal.csv"
    df = pd.read_csv(fname).iloc[1:, :]

    try:
        i_max = df.index[df.theta.diff() > 0][0]-1
    except:
        i_max = None
        
    df = df.loc[:i_max, :]
    scan[:df.shape[0], a, :] = df[['chi2_avg', 'chi2_min', 'chi2_sigma']].values


chi2_min = np.nanmin(scan[:, :, 1])
a = 0.5
chi2_target = chi2_min + a * np.sqrt(2*chi2_min/N_tau)
print(chi2_min, chi2_target)

i_target = np.zeros_like(Ar_range, dtype = int)
for a, A in enumerate(Ar_range):
    i_target[a] = np.arange(N_anneal)[scan[:, a, 0] > chi2_target][-1]+1
print(min(i_target))

##################################################



plt.figure(1, (8, 4))
i_plot = min(i_target)
# plt.plot(A_range, scan[i_plot // 4, :, 0], marker='o')
plt.errorbar(Ar_range, scan[i_plot//2, :, 0],scan[i_plot//2, :, 2], capsize = 3, fmt='-o')

plt.errorbar(Ar_range, scan[i_plot, :, 0],scan[i_plot, :, 2], capsize = 3, fmt='-o')

print(A_range[np.nanargmin( scan[i_plot, :, 0])])

plt.axhline(chi2_target, c= 'k', ls ='--')
plt.ylim(chi2_target*.95, chi2_target*1.8);


plt.xlabel(r"$A_{\mathrm{left}}$")
plt.ylabel(r"$\chi^2/N_{\tau}$")

plt.savefig("figs/tJ_1D/pi8_Aleft_scan.pdf")


# In[248]:


Nw = 200
folder = "tJ_1D_J04/pi8"
N_tau = 55
Ac_range =  np.arange(.01, .901, .01)#np.concatenate([[.1, .2, .3], np.arange(.4, .901, .02)])
Ar = .5
N_anneal = 200
scan = np.zeros((N_anneal, len(Ac_range), 3))
scan[:, :, :] = np.nan
for a, A in enumerate(Ac_range):
    
    fname = f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/Ac_{A:.3f}/p_0.500/Ar_{Ar:.3f}/anneal.csv"
    df = pd.read_csv(fname).iloc[1:, :]
    try:
        i_max = df.index[df.theta.diff() > 0][0]-1
    except:
        i_max = None
        
    df = df.loc[:i_max, :]
    scan[:df.shape[0], a, :] = df[['chi2_avg', 'chi2_min', 'chi2_sigma']].values
    # except:
    #     pass

chi2_min = np.nanmin(scan[:, :, 1])
a = 0.5
chi2_target = chi2_min + a * np.sqrt(2*chi2_min/N_tau)
print(chi2_min, chi2_target)

i_target = np.zeros_like(Ac_range, dtype = int)
for a, A in enumerate(Ac_range):
    i_target[a] = np.arange(N_anneal)[scan[:, a, 0] > chi2_target][-1]+1
min(i_target)


# In[249]:


fig, ax = plt.subplots(1, 1, figsize = (8, 4), dpi = 100)

i_plot = min(i_target)

i_plot_high = 40

# plt.errorbar(Ac_range, scan[i_plot//4, :, 0], scan[i_plot//4, :, 2],  capsize = 3,fmt='-o')
plt.errorbar(Ac_range, scan[i_plot_high, :, 0],scan[i_plot_high, :, 2],  capsize = 3,fmt='-o')
plt.errorbar(Ac_range, scan[i_plot, :, 0],scan[i_plot, :, 2], capsize = 3, fmt='-o')

plt.axhline(chi2_target, c= 'k', ls ='--')
# plt.ylim(chi2_target*.95, chi2_target*1.5);
plt.ylim(min(scan[i_plot, :, 0])*.95, min(scan[i_plot, :, 0])*1.5);

plt.xlabel("$A_c$")
plt.ylabel(r"$\chi^2/N_{\tau}$")


plt.savefig("figs/tJ_1D/pi8_Ac_scan.pdf")


# In[251]:


plt.figure(1, (6, 4), 200)

folder = "tJ_1D_J04/pi8"


Nw = 200
n = 0


Ac = 0
if Ac > 0:
    Ac_folder = f"Ac_{Ac:.3f}/"
else:
    Ac_folder = ""
Ar = 0.5


spec = 'd'

fig, ax = plt.subplots(1, 1, figsize = (8, 4), dpi = 100)

df = pd.read_csv(f"sw001.dat", sep = r'\s+', header=None)
ax.plot(df[0], df[1], label = 'Unconstrained', alpha = 1)
# ax[1].plot(df[0], df[1], label = 'Unconstrained', alpha = 1)


x_interp = np.linspace(-5, 5, 10000)

Ar = 0.5
df_R= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/{Ac_folder}p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_1.dat")
df_L= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/{Ac_folder}p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_2.dat")

R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])
ax.plot(x_interp, R_interp+L_interp, label = r"$A_{\mathrm{left}} = " +rf"{Ar:.2f}$")

Ar = 0.6
df_R= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/{Ac_folder}p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_1.dat")
df_L= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/{Ac_folder}p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_2.dat")

R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])
spec_opt = R_interp+L_interp

ax.plot(x_interp, R_interp+L_interp, label = r"$A_{\mathrm{left}} = " +rf"{Ar:.2f}$")


Ar = 0.7
df_R= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/{Ac_folder}p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_1.dat")
df_L= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/{Ac_folder}p_0.500/Ar_{Ar:.03f}/{spec}w{n:03}_2.dat")

R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])
ax.plot(x_interp, R_interp+L_interp, label = r"$A_{\mathrm{left}} = " +rf"{Ar:.2f}$")


ax.legend()
ax.set_ylim(0, 5)
ax.set_xlim(-2, 2)

ax.set_ylabel(r"$S(\omega)$")

ax.set_xlabel(r"$\omega$")
# plt.title("1D $t-J$ model, $J/t = 0.4$")

plt.savefig("figs/tJ_1D/pi8_no_cont.pdf")


# # Synth

# ## Double Edge Cont

# In[268]:


fig, ax = plt.subplots(2, 1, figsize = (8, 6), dpi = 200, sharex=True)

df = pd.read_csv("../in_files/tJ_double_edge_cont/aw.dat")
df.loc[0, 'S'] = 0
for a in ax:
    a.plot(df.omega, df.S, c= 'k', alpha = 1, lw = 2)

    a.axvline(0, c= 'k', lw=0.5)
    a.axhline(0, c= 'k', lw=0.5)
    a.set_ylim(0, 10)
    a.set_ylim(0, 10)
    a.set_xlim(-2, 3)


# df= pd.read_csv("/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_double_edge_cont/free/sw000.dat")
# # ax[0].plot(df.omega, df.S, lw = 1, label = 'Unconstrained')


Nw = 200
n = 0
spec = 's'

Ac = 0.4
# if Ac > 0:
Ac_folder = f"Ac_{Ac:.3f}/"
# else:
#     Ac_folder = ""

Ar = .5
df_R= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_double_edge_cont/edge/Nw_{Nw}/{Ac_folder}p_0.500/Ar_{Ar:.3f}/{spec}w{n:03}_1.dat")
ax[1].plot(df_R.omega, df_R.S_edge, alpha = 0.5, label = r"Rightward $\delta$'s")

ax[1].plot(df_R.omega, df_R.S_cont, alpha = 0.5, label = r"Rightward $\delta$'s")


df_L=  pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_double_edge_cont/edge/Nw_{Nw}/{Ac_folder}/p_0.500/Ar_{Ar:.3f}/{spec}w{n:03}_2.dat")
ax[1].plot(-df_L.omega, df_L.S, alpha = 0.5, label = r"Leftward $\delta$'s")

x_interp = np.linspace(min(df_R.omega), max(-df_L.omega), 1000)
R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])


ax[0].plot(x_interp, R_interp+L_interp, label = "Combined", c= 'y')


# In[265]:


Nw = 200
folder = "tJ_double_edge_cont"
N_tau = 80
Ac_range =  np.arange(.00, .901, .01)#np.concatenate([[.1, .2, .3], np.arange(.4, .901, .02)])

A_r = 0.5

N_anneal = 200
scan = np.zeros((N_anneal, len(Ac_range), 3))
scan[:, :, :] = np.nan
for a, A in enumerate(Ac_range):
    # print(A)
    # if A == .86:
    #     continue
    
    fname = f"/projectnb/qspin/gschumm/sac/fermion/out_files/{folder}/edge/Nw_{Nw}/Ac_{A:.3f}/p_0.500/Ar_{A_r:.3f}/anneal.csv"
    df = pd.read_csv(fname).iloc[1:, :]
    try:
        i_max = df.index[df.theta.diff() > 0][0]-1
    except:
        i_max = None
    df = df.loc[:i_max, :]
    scan[:df.shape[0], a, :] = df[['chi2_avg', 'chi2_min', 'chi2_sigma']].values
    # except:
    #     pass


# In[266]:


chi2_min = np.nanmin(scan[:, :, 1])
a = 0.5
chi2_target = chi2_min + a * np.sqrt(2*chi2_min/N_tau)
print(chi2_min, chi2_target)

i_target = np.zeros_like(Ac_range, dtype = int)
for a, A in enumerate(Ac_range):
    try:
        i_target[a] = np.arange(N_anneal)[scan[:, a, 0] > chi2_target][-1]+1
    except:
        i_target[a] = 1e10
        
print(np.nanmin(i_target))

plt.figure(1, (8, 4))
i_plot = min(i_target)


plt.errorbar(Ac_range, scan[i_plot//4, :, 0], scan[i_plot//4, :, 2],  capsize = 3,fmt='o')
plt.errorbar(Ac_range, scan[i_plot//2, :, 0],scan[i_plot//2, :, 2],  capsize = 3,fmt='o')
plt.errorbar(Ac_range, scan[i_plot, :, 0],scan[i_plot, :, 2], capsize = 3, fmt='o')

plt.axhline(chi2_target, c= 'k', ls ='--')
# plt.ylim(chi2_target*.95, chi2_target*1.5);
plt.ylim(min(scan[i_plot, :, 0])*.95, min(scan[i_plot, :, 0])*2);


# ## Double Edge Pos

# In[28]:


fig, ax = plt.subplots(2, 1, figsize = (8, 6), dpi = 200, sharex=True)

df = pd.read_csv("../in_files/tJ_double_edge_pos/aw.dat")
df.loc[0, 'S'] = 0
for a in ax:
    a.plot(df.omega, df.S, c= 'k', alpha = 1, lw = 2)

    a.axvline(0, c= 'k', lw=0.5)
    a.axhline(0, c= 'k', lw=0.5)

   


df= pd.read_csv("/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_double_edge_pos/free/sw000.dat")
ax[0].plot(df.omega, df.S, lw = 1, label = 'Unconstrained')



Nw = 80
n = 0
spec = 'd'
df_R= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_double_edge_pos/edge/Nw_{Nw}/p_0.500/{spec}w{n:03}_1.dat")
ax[1].plot(df_R.omega, df_R.S, alpha = 0.5, label = r"Rightward $\delta$'s")

df_L= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_double_edge_pos/edge/Nw_{Nw}/p_0.500/{spec}w{n:03}_2.dat")
ax[1].plot(-df_L.omega, df_L.S, alpha = 0.5, label = r"Leftward $\delta$'s")

x_interp = np.linspace(-3, 3, 1000)
R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])

ax[1].plot(x_interp, R_interp+L_interp, label = "Combined", c= 'y')

# plt.plot(-df_L.omega.values[::-1], df_L.S.values[::-1])

for a in ax:
    a.set_xlim(0, 4)
    a.set_ylim(0, 12)
    a.legend()
    a.set_ylabel(r"$S(\omega)$")

ax[1].set_xlabel(r"$\omega$")



# ## Double Edge

# In[21]:


fig, ax = plt.subplots(2, 1, figsize = (8, 6), dpi = 200, sharex=True)

df = pd.read_csv("../in_files/tJ_double_edge/aw.dat")
for a in ax:
    a.plot(df.omega, df.S, c= 'k', alpha = 1, lw = 2)

    a.axvline(0, c= 'k', lw=0.5)
    a.axhline(0, c= 'k', lw=0.5)

   


df= pd.read_csv("/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_double_edge/free/sw000.dat")
ax[0].plot(df.omega, df.S, lw = 1, label = 'Unconstrained')



Nw = 80
n = 0
spec = 'd'
df_R= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_double_edge/edge/Nw_{Nw}/p_0.500/{spec}w{n:03}_1.dat")
ax[1].plot(df_R.omega, df_R.S, alpha = 0.5, label = r"Rightward $\delta$'s")

df_L= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_double_edge/edge/Nw_{Nw}/p_0.500/{spec}w{n:03}_2.dat")
ax[1].plot(-df_L.omega, df_L.S, alpha = 0.5, label = r"Leftward $\delta$'s")

x_interp = np.linspace(-3, 3, 1000)
R_interp = np.interp(x_interp, df_R.omega.values, df_R.S.values)
L_interp = np.interp(x_interp, -df_L.omega.values[::-1], df_L.S.values[::-1])

ax[1].plot(x_interp, R_interp+L_interp, label = "Combined", c= 'y')

# plt.plot(-df_L.omega.values[::-1], df_L.S.values[::-1])

for a in ax:
    a.set_xlim(-2, 3)
    a.set_ylim(0, 12)
    a.legend()
    a.set_ylabel(r"$S(\omega)$")

ax[1].set_xlabel(r"$\omega$")



# ## Single Edge

# In[ ]:


plt.figure(1, (8, 4), 200)

df = pd.read_csv("../in_files/tJ_edge/aw.dat")

plt.plot(df.omega, df.S, c= 'k')
plt.axvline(0, c= 'k', lw=0.5)
plt.axhline(0, c= 'k', lw=0.5)

df= pd.read_csv("/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_edge/free/sw000.dat")
plt.plot(df.omega, df.S)

Nw = 80
spec = 'd'
df= pd.read_csv(f"/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_edge/edge/Nw_{Nw}/p_0.500/{spec}w000_1.dat")
plt.plot(df.omega, df.S)

# df= pd.read_csv("/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_edge/edge_fixed/Nw_120/p_0.500/dw000_1.dat")
# plt.plot(df.omega, df.S)
# # plt.figure()

# df = pd.read_csv("/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_edge/G.dat", header=None, sep =r'\s+')

# plt.plot(df[0], df[1])
# plt.yscale('log')

plt.ylim(0, 20)
plt.xlim(-3, 3)


# In[ ]:


plt.figure(1, (8, 4), 200)

df = pd.read_csv("../in_files/tJ_edge/aw.dat")

plt.plot(df.omega, df.S, c= 'k')
plt.axvline(0, c= 'k', lw=0.5)
plt.axhline(0, c= 'k', lw=0.5)

df= pd.read_csv("/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_edge/edge/Nw_20/p_0.500/sw000_1.dat")
# plt.plot(df.omega, df.S)

df= pd.read_csv("/projectnb/qspin/gschumm/sac/fermion/out_files/tJ_edge/edge/Nw_20/p_0.500/dw000_1.dat")
plt.plot(df.omega, df.S)
# plt.figure()

# df = pd.read_csv("/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_edge/G.dat", header=None, sep =r'\s+')

# plt.plot(df[0], df[1])
# plt.yscale('log')

plt.ylim(0, 20)


# # t.in

# In[26]:


# f_synth = '/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_edge/t.in'
# N_tau = 80
# df = pd.read_csv(f_synth, skiprows = 1, header=None, sep = r"\s+", nrows =N_tau)
# plt.plot(df[0], df[1])
# # plt.plot(df[0], df[21/df[1])
# plt.yscale("log")

f_sib = '/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_1D_J2/t.in'
N_tau = 80
df = pd.read_csv(f_sib, skiprows = 1, header=None, sep = r"\s+", nrows =N_tau)
# plt.plot(df[0], df[1])
plt.plot(df[0], df[1])


f_synth = '/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_double_edge_pos/t.in'
N_tau = 80
df = pd.read_csv(f_synth, skiprows = 1, header=None, sep = r"\s+", nrows =N_tau)
plt.plot(df[0], df[1])
# # plt.plot(df[0], df[2]/df[1])
plt.yscale("log")


# In[27]:


[(np.log(df[1].values[-1]) - np.log(df[1].values[-2]))/.1,
 (np.log(df[1].values[1]) - np.log(df[1].values[0]))/.1]


# In[28]:


f_sib = '/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_1D_J2/t.in'
N_tau = 80
df = pd.read_csv(f_sib, skiprows = 1, header=None, sep = r"\s+", nrows =N_tau)
# plt.plot(df[0], df[1])
plt.plot(df[0], df[2]/df[1])


f_synth = '/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_double_edge/t.in'
N_tau = 80
df = pd.read_csv(f_synth, skiprows = 1, header=None, sep = r"\s+", nrows =N_tau)
# plt.plot(df[0], df[1])
plt.plot(df[0], df[2]/df[1])

plt.yscale("log")

# plt.ylim(1e-1, 1e1)


# In[29]:


f_sib = '/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_1D_J2/t.in'
N_tau = 80
df = pd.read_csv(f_sib, skiprows = 1, header=None, sep = r"\s+", nrows =N_tau)
# plt.plot(df[0], df[1])
plt.plot(df[0], df[2])


f_synth = '/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_double_edge_pos/t.in'
N_tau = 80
df = pd.read_csv(f_synth, skiprows = 1, header=None, sep = r"\s+", nrows =N_tau)
# plt.plot(df[0], df[1])
plt.plot(df[0], df[2])

plt.yscale("log")

# plt.ylim(1e-1, 1e1)


# In[30]:


f_sib = '/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_1D_J2/t.in'
N_tau = 80
df = pd.read_csv(f_sib, skiprows = 1, header=None, sep = r"\s+", nrows =N_tau)
# plt.plot(df[0], df[1])
plt.plot(df[0], df[3])


f_synth = '/projectnb/qspin/gschumm/sac/fermion/in_files/tJ_double_edge_pos/t.in'
N_tau = 80
df = pd.read_csv(f_synth, skiprows = 1, header=None, sep = r"\s+", nrows =N_tau)
# plt.plot(df[0], df[1])
plt.plot(df[0], df[3])

plt.yscale("log")

# plt.ylim(1e-1, 1e1)


# In[ ]:




