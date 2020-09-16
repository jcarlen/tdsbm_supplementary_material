# Assumes you're in the mixed_model_implementation_python folder

from lib import gradient_mixed
from lib import visualize
from lib import mixmem_aux

import pandas as pd
import math
import random
import numpy as np
import scipy as scp
import scipy.stats as stats
import matplotlib.pyplot as plt

gradient_descent = gradient_mixed.gradient_descent
realistic_initial = gradient_mixed.realistic_initial
likelihood = mixmem_aux.likelihood
mu_likelihood = mixmem_aux.mu_likelihood

scatterplot_city = visualize.scatterplot_city
block_omegas = visualize.block_omegas

### 
# Variables across models
"""
    - Model Parameters:
        - numblocks: the number of blocks.
        - r: variables with r refer to the block to block time interconnectivity parameter
        - they are np.arrays of dimensions numblocks x numblocks x 24
        - r0: random initialization of omega for the gradient descent
        - Cityr_fit (e.g. LAr_fit) is computed maximum likelihood estimate omega for the network of City
        - c: variables with c  refer to the block membership parameter C
        - they are np.arrays of dimension number-of-stations x numblocks
        - Cityc_fit (e.g. Manhattanc_fit) is computed maximum likelihood estimate C for the network of City
        - c0: random initializtion of C for the gradient descent
    - Data Variables (replace City with the city of the network (e.g. LAdf, NY_stations))
        - Citydf: the pandas dataframe for the network of City
        - City_stations: the list of station IDs
        - City_lat,City_lon: the lists of latitudes and longitudes of stations
        - City_N: number of station in the City
        - City_info: a dictionary mapping station id's to (latitude,longitude) pairs.
    - Suggestions:
        - The number of iterations should suffice to find a meaningful local optimum
        - Sometimes less iterations (600 works for Manhattan) can be used to see if the optimization
        is headed towards a meaningful local optimum.
        - This can be refined by plugging in r0=Cityr_fit c0=Cityc_fit and running the optimization for more steps.
        
"""


# LA models
##############################################################################################################
# Load data in array (matrix) format
SIM_matrix = pd.read_csv("")
SIM_T = 16

num_runs=10#number of times we re-run the algorithm.
numblocklist = np.array(range(2,6)) #blocks we try is 2,3,4,5,6
# SIM TDMM-SBM for # Blocks in numblocklist
SIM_ll=-np.inf*np.ones((num_runs,np.size(numblocklist,0)))
SIMnumparameters=SIM_N*numblocklist- numblocklist + SIM_T*numblocklist**2
SIM_fit=[];
numiter=3000

for n_run in range(0,num_runs):
    SIM_fit_sample=[];
    print("n_run "+str(n_run+1))
    for blocks_i in range(0,np.size(numblocklist,0)):
        numblocks=numblocklist[blocks_i]
        print("numblocks "+str(numblocks))
        r0,c0=realistic_initial(SIM_matrix,SIM_N,numblocks)
        SIMr_fit, SIMc_fit, wr, wc = gradient_descent(SIM_matrix, r0, c0, N=numiter,sig_digs=4, N_stable=600, nblocks = numblocks, verbose_time = 100)
        SIM_fit_sample.append((SIMr_fit,SIMc_fit))
        SIM_ll[n_run,blocks_i]=likelihood(SIMc_fit,SIMr_fit,SIM_matrix)  
    SIM_fit.append(SIM_fit_sample)
SIM_maxind=np.argmax(SIM_ll,0)
SIM_llmax=np.max(SIM_ll,0)

#Save
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    pd.DataFrame(LA_fit[LA_maxind[blocks_i]][blocks_i][0].reshape([numblocks**2, SIM_T])).to_csv("../mixed_model_results/LA_"+str(numblocks)+"_omega.csv",  index = False)
    pd.DataFrame(LA_fit[LA_maxind[blocks_i]][blocks_i][1], index = LA_stations).to_csv("../mixed_model_results/LA_"+str(numblocks)+"_roles.csv", index = LA_stations)
