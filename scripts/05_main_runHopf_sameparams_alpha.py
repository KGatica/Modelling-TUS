# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 14:57:25 2018

Script for running a fMRI-BOLD whole-brain simulation using the Hopf model [1,2].

[1] Deco, G., Cruzat, J., Cabral, J., Tagliazucchi, E., Laufs, H., Logothetis, N. K., 
& Kringelbach, M. L. (2019). Awakening: Predicting external stimulation to force 
transitions between different brain states. Proceedings of the National Academy 
of Sciences, 116(36), 18088-18097.

[2] Escrichs, A., Perl, Y. S., Uribe, C., Camara, E., Türker, B., 
Pyatigorskaya, N., ... & Deco, G. (2022). Unifying turbulent dynamics 
framework distinguishes different brain states. Communications Biology,
 5(1), 638.

@author: Carlos Coronel
@modified by: Marilyn Gatica

"""

import numpy as np
from scipy import signal
import time
import Hopf_model_multi as HM
import matplotlib.pyplot as plt
import importlib
importlib.reload(HM)
import scipy.io as sio
import os
from multiprocessing import Pool
import sys

#%%
# plt.close("all")
#Simulation parameters
HM.dt = 1E-1 #Integration step (seconds)
HM.teq = 100 #Equilibrium time (seconds)
HM.tmax = 1100 #Signals' length (seconds)
HM.downsamp = 2#2 #This reduces the number of points of the signals by a factor X
HM.seed = 0 #random seed
HM.F0 = 0 * np.ones((HM.nnodes,1)) #External stimulation

#Network parameters
sc_aux = np.load("../datasets/sc_subset.npy")#empirical SC

#print(sys.argv)
ncores =  4#int(sys.argv[1])
# Model selection: communicability (CMY) or distance.
model =  "CMY"#str(sys.argv[2])
# Stimulation intensity
alpha = 0.3#float(sys.argv[4]) #In the paper: alpha \in np.arange(-0.8,0.9,0.1)
print(model)
global_het = sio.loadmat("../datasets/global_%s.mat"%model)["global_measure"].T

# target: parstriangularis (inferior frontal cortex), thalamusproper (thalamus)
stim =  "parstriangularis"#str(sys.argv[3])
local_het = sio.loadmat("../datasets/discrete_local_%s.mat"%(stim))["local_measure"].T

try:
    os.mkdir("./alpha%.3f_heterogeneous_%s_%s"%(alpha,model,stim))
except:
    print("not created")

def main_Hopf_parallel(n):
    try:
        os.mkdir("./alpha%.3f_heterogeneous_%s_%s/subj%i"%(alpha,model,stim,(n+1)))
    except:
        print("")
            
    
    HM.C = sc_aux[:,:,n]
    HM.C = HM.C/np.max(HM.C)#normalize
    
    HM.nnodes = len(HM.C)
    HM.norm = 1 #global normalization
    
    #Model parameters
    HM.w = 0.05 * 2 * np.pi #Oscillatory frequency
    HM.beta = 0.02#noise scaling factor (background noise) 
    
    G = 0.16#Fit in the script 02_fitting_CN.py
    HM.G = G #Global coupling

    ###We fit the parameters in the previous script (04_fitting_stim_het.py)
    if model == "CMY":
        bias = -0.08
        scale = 0.21
    elif model == "distance":
        bias = -0.17
        scale = 0.24
    
    # Local stimulation at different intensities
    # Heterogeneity is shaped by global CMY or distance model, and local stimulation modulated via 'alpha'
    HM.a = bias + scale*(global_het + alpha*local_het)#Bifurcation parameter
    
    #Update parameters that aren't inputs of model's functions (e.g. oscillatory frequency)
    HM.update()
        
    #Simulation starts here
    init = time.time()
    y, time_vector = HM.Sim(verbose = False) #True showd time, for anxious people
    BOLD_signals = y[:,:,0] #fMRI BOLD-like output of the model
    end = time.time()
    #Simulation ends here
    # print('Sim time: %.3f'%(end-init))
    
    #Filtering BOLD signals
    resolution = HM.dt * HM.downsamp
    Fmin, Fmax = 0.01, 0.1 #Allowed frequencies

    a0,b0 = signal.bessel(3,[2 * resolution * Fmin, 2 * resolution * Fmax],btype='bandpass')
    BOLD_filt = signal.filtfilt(a0,b0,BOLD_signals,axis=0) #filteres signals

        
    phases_signal = np.angle(signal.hilbert(BOLD_filt, axis = 0)) #Phases
    kop = HM.simple_order_parameter(phases_signal) #Kuramoto order parameter
    data={'BOLD':BOLD_filt[500:-500],'kop':kop[500:-500]}
    sio.savemat('./alpha%.3f_heterogeneous_%s_%s/subj%i/Bold_bias%.2f_scale%.2f.mat'%(alpha,model,stim,n+1,bias,scale),data) 


#print(sys.argv)   
#ncores =  int(sys.argv[1])
#model =  str(sys.argv[2])
#print(model)
#global_het = sio.loadmat("global_%s.mat"%model)["global_measure"].T

if __name__ == '__main__':
    with Pool(processes=ncores) as pool:
        results = pool.map(main_Hopf_parallel, range(19))
        


