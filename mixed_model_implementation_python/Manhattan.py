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
        - Cityr_fit (e.g. Manhattanr_fit) is computed maximum likelihood estimate omega for the network of City
        - c: variables with c  refer to the block membership parameter C
        - they are np.arrays of dimension number-of-stations x numblocks
        - Cityc_fit (e.g. Manhattanc_fit) is computed maximum likelihood estimate C for the network of City
        - c0: random initializtion of C for the gradient descent
    - Data Variables (replace City with the city of the network (e.g. Manhattandf, NY_stations))
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


# Manhattan home subset models
###############################################################################################################
# Load Manhattan Data


Manhattandf = pd.read_csv("../data/cleaned/ny1610_hm_no_weekend.csv")
Manhattan_stations_dups = list(Manhattandf['Start.Station.ID'])
Manhattan_lat = list(Manhattandf['Start.Station.Longitude'])
Manhattan_lon = list(Manhattandf['Start.Station.Latitude'])
Manhattan_station_info = dict(zip(Manhattan_stations_dups, zip(Manhattan_lat, Manhattan_lon)))
Manhattan_stations = list(set(Manhattan_stations_dups))
Manhattan_stations.sort()
Manhattan_N = len(Manhattan_stations)
Manhattan_matrix = np.zeros([Manhattan_N, Manhattan_N, 24])

for index, row in Manhattandf.iterrows():
    #if ~row['weekend']:
    i = Manhattan_stations.index(row['Start.Station.ID'])
    j = Manhattan_stations.index(row['End.Station.ID'])
    t = int(row['Start.Time'].split(':')[0].split(' ')[1])
    Manhattan_matrix[i, j, t] += 1


# Run TDMM-SBM with 2 blocks

num_runs=10#number of times we re-run the algorithm.
numblocklist = np.array(range(2,7)) #blocks we try is 2,3,4,5,6
# Manhattan TDMM-SBM for # Blocks in numblocklist
Manhattan_ll=-np.ones((num_runs,np.size(numblocklist,0)))*np.inf
Manhattannumparameters=Manhattan_N*numblocklist- numblocklist + 24*numblocklist**2
Manhattan_fit=[];
numiter=3000

for n_run in range(0,num_runs):
    Manhattan_fit_sample=[];
    print("n_run "+str(n_run+1))
    for blocks_i in range(0,np.size(numblocklist,0)):
        numblocks=numblocklist[blocks_i]
        print("numblocks "+str(numblocks))
        r0,c0=realistic_initial(Manhattan_matrix,Manhattan_N,numblocks)
        Manhattanr_fit, Manhattanc_fit, wr, wc = gradient_descent(Manhattan_matrix, r0, c0, N=numiter,sig_digs=4, N_stable=600, nblocks = numblocks, verbose_time = 100)
        Manhattan_fit_sample.append((Manhattanr_fit,Manhattanc_fit))
        Manhattan_ll[n_run,blocks_i]=likelihood(Manhattanc_fit,Manhattanr_fit,Manhattan_matrix)  
    Manhattan_fit.append(Manhattan_fit_sample)
Manhattan_maxind=np.argmax(Manhattan_ll,0)
Manhattan_llmax=np.max(Manhattan_ll,0)
Manhattan_AIC=2*Manhattannumparameters-2*Manhattan_llmax
fig = plt.figure(figsize=(20,20))
plt.scatter(numblocklist,Manhattan_AIC)
plt.xlabel("number of blocks",fontsize=20)
plt.ylabel("AIC",fontsize=20)
plt.savefig("img/Manhattan/ManhattanAIC.png")
print("Worst log-likelihoods"+str(np.min(Manhattan_ll,0)))
print("Best log-likelihoods"+str(np.max(Manhattan_ll,0)))
#Visualize 
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    Manhattanr_fit=Manhattan_fit[Manhattan_maxind[blocks_i]][blocks_i][0]
    Manhattanc_fit=Manhattan_fit[Manhattan_maxind[blocks_i]][blocks_i][1]
    
    Manhattan_weights=np.sum(Manhattan_matrix,(1,2))+np.sum(Manhattan_matrix,(0,2))
    Manhattan_sizes=Manhattan_weights/15;
    Manhattan_colors=Manhattanc_fit/np.max(Manhattanc_fit,0,keepdims=True)
    scatterplot_city(Manhattanc_fit,Manhattan_N,numblocks,np.array([list(Manhattan_station_info[x]) for x in Manhattan_stations])
                     ,sizes=Manhattan_sizes,filename="img/Manhattan/Manhattan"+str(numblocks)+"blocks"+str(numiter)+"iter")
    block_omegas(Manhattanr_fit,numblocks,filename="img/Manhattan/Manhattan"+str(numblocks)+"blocks"+str(numiter)+"iter")

#Save
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    pd.DataFrame(Manhattan_fit[Manhattan_maxind[blocks_i]][blocks_i][0].reshape([numblocks**2, 24],order = 'F')).to_csv("../mixed_model_results/Manhattan_"+str(numblocks)+"_omega.csv",  index = False)
    pd.DataFrame(Manhattan_fit[Manhattan_maxind[blocks_i]][blocks_i][1], index = Manhattan_stations).to_csv("../mixed_model_results/Manhattan_"+str(numblocks)+"_roles.csv", index = Manhattan_stations)
