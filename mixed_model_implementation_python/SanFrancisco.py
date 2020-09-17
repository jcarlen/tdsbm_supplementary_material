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
        - Cityr_fit (e.g. SFr_fit) is computed maximum likelihood estimate omega for the network of City
        - c: variables with c  refer to the block membership parameter C
        - they are np.arrays of dimension number-of-stations x numblocks
        - Cityc_fit (e.g. Manhattanc_fit) is computed maximum likelihood estimate C for the network of City
        - c0: random initializtion of C for the gradient descent
    - Data Variables (replace City with the city of the network (e.g. SFdf, NY_stations))
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

# San Fransisco models
###############################################################################################################
# Load SF Data

SFdf = pd.read_csv("../data/cleaned/SF16_cleaned_final_no_weekend.csv")
SF_stations_dups = list(SFdf['start_station_id'])
SF_lat = list(SFdf['start_lat'])
SF_lon = list(SFdf['start_lon'])
SF_station_info = dict(zip(SF_stations_dups, zip(SF_lat, SF_lon)))
SF_stations = list(set(SF_stations_dups))
SF_stations.sort()
SF_N = len(SF_stations)
SF_matrix = np.zeros([SF_N, SF_N, 24])
for index, row in SFdf.iterrows():
    # print(index)
    #if index < 60000:
    i = SF_stations.index(row['start_station_id'])
    j = SF_stations.index(row['end_station_id'])
    t = int(row['Start.Date'].split(':')[0].split(' ')[1])
    #if i != j:
    SF_matrix[i, j, t] += 1
    # print(SF_matrix[i,j,t])


num_runs=10#number of times we re-run the algorithm.
numblocklist = np.array(range(2,5)) #blocks we try is 2,3,4
# SF TDMM-SBM for # Blocks in numblocklist
SF_ll=-np.inf*np.ones((num_runs,np.size(numblocklist,0)))
SFnumparameters=SF_N*numblocklist- numblocklist + 24*numblocklist**2
SF_fit=[];
numiter=3000

for n_run in range(0,num_runs):
    SF_fit_sample=[];
    print("n_run "+str(n_run+1))
    for blocks_i in range(0,np.size(numblocklist,0)):
        numblocks=numblocklist[blocks_i]
        print("numblocks "+str(numblocks))
        r0,c0=realistic_initial(SF_matrix,SF_N,numblocks)
        SFr_fit, SFc_fit, wr, wc = gradient_descent(SF_matrix, r0, c0, N=numiter,sig_digs=4, N_stable=600, nblocks = numblocks, verbose_time = 100)
        SF_fit_sample.append((SFr_fit,SFc_fit))
        SF_ll[n_run,blocks_i]=likelihood(SFc_fit,SFr_fit,SF_matrix)  
    SF_fit.append(SF_fit_sample)
SF_maxind=np.argmax(SF_ll,0)
SF_llmax=np.max(SF_ll,0)
SF_AIC=2*SFnumparameters-2*SF_llmax
fig = plt.figure(figsize=(20,20))
plt.scatter(numblocklist,SF_AIC)
plt.xlabel("number of blocks",fontsize=20)
plt.ylabel("AIC",fontsize=20)
plt.savefig("img/SF/SFAIC.png")
print("Worst log-likelihoods"+str(np.min(SF_ll,0)))
print("Best log-likelihoods"+str(np.max(SF_ll,0)))
#Visualize 
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    SFr_fit=SF_fit[SF_maxind[blocks_i]][blocks_i][0]
    SFc_fit=SF_fit[SF_maxind[blocks_i]][blocks_i][1]
    
    SF_weights=np.sum(SF_matrix,(1,2))+np.sum(SF_matrix,(0,2))
    SF_sizes=SF_weights/15;
    SF_colors=SFc_fit/np.max(SFc_fit,0,keepdims=True)
    scatterplot_city(SFc_fit,SF_N,numblocks,np.array([list(SF_station_info[x]) for x in SF_stations])
                     ,sizes=SF_sizes,filename="img/SF/SF"+str(numblocks)+"blocks"+str(numiter)+"iter")
    block_omegas(SFr_fit,numblocks,filename="img/SF/SF"+str(numblocks)+"blocks"+str(numiter)+"iter")

#Save
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    pd.DataFrame(SF_fit[SF_maxind[blocks_i]][blocks_i][0].reshape([numblocks**2, 24], order = 'F')).to_csv("../mixed_model_results/SF_"+str(numblocks)+"_omega.csv",  index = False)
    pd.DataFrame(SF_fit[SF_maxind[blocks_i]][blocks_i][1], index = SF_stations).to_csv("../mixed_model_results/SF_"+str(numblocks)+"_roles.csv", index = SF_stations)


