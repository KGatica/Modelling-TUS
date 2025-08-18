# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:07:12 2024

@author: mgatica
"""

import numpy as np
from scipy import signal
import scipy.io as sio
import matplotlib.pyplot as plt
import Hopf_model_multi as HM

# Number of subjects
Nsubj = 19
# Number of brain areas 
modules = 84

# Range of global coupling values (G)
G_vals =  np.arange(0,0.51,0.05)#np.arange(0,0.51,0.01)
# Kuramoto order parameter per subject and G
kop_sim = np.zeros((len(G_vals),Nsubj))
# Correlation between empirical and simulated FC
corr_FC = np.zeros((len(G_vals),Nsubj))

# Load empirical data
triu_emp_all = np.load("../datasets/triu_emp_subset.npy")# Empirical FC (upper triangular)
kop_emp = np.load("../datasets/kop_emp_subset.npy")# Empirical KOP


for n in range(Nsubj):
    print(n)
    triu_emp = triu_emp_all[:,n]
    for (idG,G) in enumerate(G_vals):
        BOLD_filt_sim = sio.loadmat('./homogeneous_CN/subj%i/Bold_G%.2f.mat'%(n+1,G))["BOLD_original"]
        resolution = 0.1 * 2
        Fmin, Fmax = 0.01, 0.1 #Allowed frequencies
        a0,b0 = signal.bessel(3,[2 * resolution * Fmin, 2 * resolution * Fmax],btype='bandpass')# filteres signals
        BOLD_filt_sim = signal.filtfilt(a0,b0,BOLD_filt_sim,axis=0)
        BOLD_filt_sim = BOLD_filt_sim[500:-500,:]
        phases_signal = np.angle(signal.hilbert(BOLD_filt_sim, axis = 0)) # Phases
        kop_sim[idG,n] = np.mean(HM.simple_order_parameter(phases_signal)) # KOP for simulations
        
        FC_sim = np.corrcoef(BOLD_filt_sim.T) # Simulated FC
        triu_sim = FC_sim[np.triu_indices(modules, k = 1)] # Upper triangle of the simulated FC matrix
        corr_FC[idG,n] = np.corrcoef(triu_sim,triu_emp)[0,1] # Correlation between empirical and simulated FC

        break

#%% Visualizing
diff_kop = abs(kop_emp - kop_sim) # difference between simulated and empirical KOP
plt.figure()

plt.plot(G_vals, np.mean(corr_FC,1))
plt.fill_between(G_vals, np.mean(corr_FC,1)-np.std(corr_FC,1), np.mean(corr_FC,1)+np.std(corr_FC,1),alpha=0.2)

plt.plot(G_vals,np.mean(diff_kop,1))
plt.fill_between(G_vals, np.mean(diff_kop,1)-np.std(diff_kop,1), np.mean(diff_kop,1)+np.std(diff_kop,1),alpha=0.2)

# plt.ylabel("ErrKOP FCfit")
plt.xlabel("global coupling (G)")
plt.legend(["FCfit","","ErrKOP",""])

plt.ylim(0,0.5)
#%% This is the input for the CN heterogeneous model
mean_ErrKOP = np.mean(diff_kop,1)
min_ErrKOP =np.min(mean_ErrKOP)
G_vals[mean_ErrKOP == min_ErrKOP]
