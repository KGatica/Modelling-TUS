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
#import matplotlib.pyplot as plt
import importlib
importlib.reload(HM)
import scipy.io as sio
import os
from multiprocessing import Pool
#import sys

#%%
# plt.close("all")
#Simulation parameters
HM.dt = 1E-1 #Integration step (seconds)
HM.teq = 100 #Equilibrium time (seconds)
HM.tmax = 1100 #Signals' length (seconds)
HM.downsamp = 2#2 #This reduces the number of points of the signals by a factor X
HM.seed = 0 #random seed

HM.F0 = 0 * np.ones((HM.nnodes,1)) #External stimulation

HM.a = -0.02 * np.ones((HM.nnodes,1)) #Bifurcation parameter

#Network parameters
sc_aux = np.load("../datasets/sc_subset.npy")


try:
    os.mkdir("homogeneous_CN")
except:
    print("")

def main_Hopf_parallel(n):
    try:
        os.mkdir("homogeneous_CN/subj%i"%(n+1))
    except:
        print("")
            

    HM.C = sc_aux[:,:,n]
    HM.C = HM.C/np.max(HM.C)#normalize
    
    HM.nnodes = len(HM.C)
    HM.norm = 1 #global normalization
    
    #Model parameters
    HM.w = 0.05 * 2 * np.pi #Oscillatory frequency
    #Note: frequencies could be different between brain regions and also obtained using the real data
    HM.beta = 0.02#noise scaling factor (background noise) 
    
    for G in np.arange(0,0.51,0.05):#np.arange(0,0.51,0.01):
        HM.G = G #Global coupling
            
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
        # print(np.mean(kop))

        phases_original = np.angle(signal.hilbert(BOLD_signals,axis=0))
        kop_original = HM.simple_order_parameter(phases_original)
                
        #break
        data={'BOLD':BOLD_filt[500:-500],'kop':kop[500:-500],'BOLD_original':BOLD_signals,'kop_original':kop_original}
        sio.savemat('homogeneous_CN/subj%i/Bold_G%.2f.mat'%(n+1,G),data)
        break

#print(sys.argv)   
ncores =  2#int(sys.argv[1])
if __name__ == '__main__':
    with Pool(processes=ncores) as pool:
        results = pool.map(main_Hopf_parallel, range(19))
        



