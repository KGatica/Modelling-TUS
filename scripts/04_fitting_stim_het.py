# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:07:12 2024

@author: mgatica
"""

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


#print(sys.argv)
#Select the model: communicability (CMY) or distance
stim = "CMY" #str(sys.argv[1])
simpath = "heterogeneous_CN_%s"%(stim)

# Number of subjects and brain areas 
Nsubj = 19
# Number of brain areas 
modules = 84


# Testing some bias and scale values, the original ones are commented
bias_vals = np.round(np.arange(-1.1,0.1,0.1), 2)#np.round(np.arange(-1.1,0.1,0.03), 2)
scale_vals =  np.arange(0,1.1,0.1)#np.arange(0,1.1,0.03)

# Initialize result matrices (scale x bias x subject) to store:
# kop_sim: simulated KOP
# corr_FC: correlation between simulated and empirical FC
# kop_diff: difference between simulated and empirical KOP
kop_sim = np.zeros((len(scale_vals),len(bias_vals),Nsubj))
corr_FC = np.zeros((len(scale_vals),len(bias_vals),Nsubj))
kop_diff = np.zeros((len(scale_vals),len(bias_vals),Nsubj))

# Load empirical data
triu_emp_all = np.load("../datasets/triu_emp_subset.npy")#Empirical FC (upper triangular part), per subject
kop_emp = np.load("../datasets/kop_emp_subset.npy")#Empirical KOP, per subject

for n in range(Nsubj):
    print(n)
    triu_emp = triu_emp_all[:,n]
    #Computing the FC and KOP per simulation
    for (idscale,scale) in enumerate(scale_vals):
        for (idbias,bias) in enumerate(bias_vals):
            kop_sim[idscale,idbias,n] = np.mean(sio.loadmat('./%s/subj%i/Bold_bias%.2f_scale%.2f.mat'%(simpath,n+1,bias,scale))["kop"]) 
            BOLD_filt_sim = sio.loadmat('./%s/subj%i/Bold_bias%.2f_scale%.2f.mat'%(simpath,n+1,bias,scale))["BOLD"]
            FC_sim = np.corrcoef(BOLD_filt_sim.T)  # Simulated FC
            triu_sim = FC_sim[np.triu_indices(modules, k = 1)] # Upper triangle of the simulated FC matrix
            corr_FC[idscale,idbias,n] = np.corrcoef(triu_sim,triu_emp)[0,1] # Correlation between empirical and simulated FC
            kop_diff[idscale,idbias,n] = abs(kop_sim[idscale,idbias,n] - kop_emp[0,n]) # difference between simulated and empirical KOP
            break
        break
#%% To visualize and find the fitting

# Compute optimization score combining FC fit and KOP error (lower KOP error and higher FC correlation are better)
opt_ErrKOP_FC = (1- np.mean(kop_diff,axis=2))* np.mean(corr_FC,axis=2)
plt.figure(figsize=(15,4))
plt.subplot(133)
plt.imshow(opt_ErrKOP_FC);
plt.colorbar()
plt.xticks(np.linspace(0,len(bias_vals)-1,len(bias_vals[0::5])),bias_vals[0::5])
plt.yticks(np.linspace(0,len(scale_vals)-1,len(scale_vals[0::5])),scale_vals[0::5])

# Identify and mark optimal (scale, bias) combination
dic_opt = {"opt_ErrKOP_FC":opt_ErrKOP_FC}

i, j = np.unravel_index(opt_ErrKOP_FC.argmax(), opt_ErrKOP_FC.shape)
plt.text(j, i, "*", fontsize=14,
         horizontalalignment='center',
         verticalalignment='center_baseline') # mark optimal point

opt_bias = bias_vals[j]
opt_scale = scale_vals[i]
plt.ylabel("scale")
plt.xlabel("bias")

# Subplot: Mean KOP error across subjects
plt.subplot(131)
mean_ErrKOP = np.mean(kop_diff,axis=2)

plt.imshow(np.mean(kop_diff,axis=2),cmap="Reds")
plt.text(j, i, "*", fontsize=14,
         horizontalalignment='center',
         verticalalignment='center_baseline') # mark optimal point
plt.colorbar()
plt.title("ErrKOP")
plt.xticks(np.linspace(0,len(bias_vals)-1,len(bias_vals[0::5])),bias_vals[0::5])
plt.yticks(np.linspace(0,len(scale_vals)-1,len(scale_vals[0::5])),scale_vals[0::5])

# Subplot: Mean FC correlation across subjects
plt.subplot(132)
plt.imshow(np.mean(corr_FC,axis=2),cmap="Blues")
plt.colorbar()
plt.title("FCfit")
plt.xticks(np.linspace(0,len(bias_vals)-1,len(bias_vals[0::5])),bias_vals[0::5])
plt.yticks(np.linspace(0,len(scale_vals)-1,len(scale_vals[0::5])),scale_vals[0::5])

plt.text(j, i, "*", fontsize=14,
         horizontalalignment='center',
         verticalalignment='center_baseline')#coords (x,y)
