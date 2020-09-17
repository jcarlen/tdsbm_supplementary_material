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
        - Cityr_fit (e.g. NYr_fit) is computed maximum likelihood estimate omega for the network of City
        - c: variables with c  refer to the block membership parameter C
        - they are np.arrays of dimension number-of-stations x numblocks
        - Cityc_fit (e.g. Manhattanc_fit) is computed maximum likelihood estimate C for the network of City
        - c0: random initializtion of C for the gradient descent
    - Data Variables (replace City with the city of the network (e.g. NYdf, NY_stations))
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


# NY Models
###############################################################################################################
# Load NY Data

NYdf = pd.read_csv("../data/cleaned/Oct16_nyc_cleaned_final_no_weekend.csv")
NY_stations_dups = list(NYdf['Start.Station.ID'])
NY_lat = list(NYdf['Start.Station.Longitude'])
NY_lon = list(NYdf['Start.Station.Latitude'])
NY_station_info = dict(zip(NY_stations_dups, zip(NY_lat, NY_lon)))
NY_stations = list(set(NY_stations_dups))
NY_stations.sort()
NY_N = len(NY_stations)
NY_matrix = np.zeros([NY_N, NY_N, 24])
for index, row in NYdf.iterrows():
    #if ~row['weekend']:
    i = NY_stations.index(row['Start.Station.ID'])
    j = NY_stations.index(row['End.Station.ID'])
    t = int(row['Start.Time'].split(':')[0].split(' ')[1])
    NY_matrix[i, j, t] += 1

num_runs=10#number of times we re-run the algorithm.
numblocklist = np.array(range(2,5)) #blocks we try is 2,3,4
# NY TDMM-SBM for # Blocks in numblocklist
NY_ll=-np.inf*np.ones((num_runs,np.size(numblocklist,0)))
NYnumparameters=NY_N*numblocklist- numblocklist + 24*numblocklist**2
NY_fit=[];
numiter=3000

for n_run in range(0,num_runs):
    NY_fit_sample=[];
    print("n_run "+str(n_run+1))
    for blocks_i in range(0,np.size(numblocklist,0)):
        numblocks=numblocklist[blocks_i]
        print("numblocks "+str(numblocks))
        r0,c0=realistic_initial(NY_matrix,NY_N,numblocks)
        NYr_fit, NYc_fit, wr, wc = gradient_descent(NY_matrix, r0, c0, N=numiter,sig_digs=4, N_stable=600, nblocks = numblocks, verbose_time = 100)
        NY_fit_sample.append((NYr_fit,NYc_fit))
        NY_ll[n_run,blocks_i]=likelihood(NYc_fit,NYr_fit,NY_matrix)  
    NY_fit.append(NY_fit_sample)
NY_maxind=np.argmax(NY_ll,0)
NY_llmax=np.max(NY_ll,0)
NY_AIC=2*NYnumparameters-2*NY_llmax
fig = plt.figure(figsize=(20,20))
plt.scatter(numblocklist,NY_AIC)
plt.xlabel("number of blocks",fontsize=20)
plt.ylabel("AIC",fontsize=20)
plt.savefig("img/NY/NYAIC.png")
print("Worst log-likelihoods"+str(np.min(NY_ll,0)))
print("Best log-likelihoods"+str(np.max(NY_ll,0)))
#Visualize 
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    NYr_fit=NY_fit[NY_maxind[blocks_i]][blocks_i][0]
    NYc_fit=NY_fit[NY_maxind[blocks_i]][blocks_i][1]
    
    NY_weights=np.sum(NY_matrix,(1,2))+np.sum(NY_matrix,(0,2))
    NY_sizes=NY_weights/15;
    NY_colors=NYc_fit/np.max(NYc_fit,0,keepdims=True)
    scatterplot_city(NYc_fit,NY_N,numblocks,np.array([list(NY_station_info[x]) for x in NY_stations])
                     ,sizes=NY_sizes,filename="img/NY/NY"+str(numblocks)+"blocks"+str(numiter)+"iter")
    block_omegas(NYr_fit,numblocks,filename="img/NY/NY"+str(numblocks)+"blocks"+str(numiter)+"iter")

#Save
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    pd.DataFrame(NY_fit[NY_maxind[blocks_i]][blocks_i][0].reshape([numblocks**2, 24], order = 'F')).to_csv("../mixed_model_results/NY_"+str(numblocks)+"_omega.csv",  index = False)
    pd.DataFrame(NY_fit[NY_maxind[blocks_i]][blocks_i][1], index = NY_stations).to_csv("../mixed_model_results/NY_"+str(numblocks)+"_roles.csv", index = NY_stations)

