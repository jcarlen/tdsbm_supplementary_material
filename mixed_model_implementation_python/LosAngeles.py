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
# Load LA data
LAdf = pd.read_csv("../data/cleaned/LA16_cleaned_final_no_weekend.csv")
LA_stations_dups = list(LAdf['start_station_id'])
LA_lat = list(LAdf['start_lat'])
LA_lon = list(LAdf['start_lon'])
LA_station_info = dict(zip(LA_stations_dups, zip(LA_lat, LA_lon)))
LA_stations = list(set(LA_stations_dups))
LA_stations.sort()
LA_N = len(LA_stations)
LA_matrix = np.zeros([LA_N, LA_N, 24])
for index, row in LAdf.iterrows():
    i = LA_stations.index(row['start_station_id'])
    j = LA_stations.index(row['end_station_id'])
    t = int(row['start_time'].split(':')[0].split(' ')[1])
    LA_matrix[i, j, t] += 1
    
num_runs=10#number of times we re-run the algorithm.
numblocklist = np.array(range(2,11)) #blocks we try is 2,3,4,5,6,7,8,9,10
# LA TDMM-SBM for # Blocks in numblocklist
LA_ll=-np.inf*np.ones((num_runs,np.size(numblocklist,0)))
LAnumparameters=LA_N*numblocklist- numblocklist + 24*numblocklist**2
LA_fit=[];
numiter=10000

for n_run in range(0,num_runs):
    LA_fit_sample=[];
    print("n_run "+str(n_run+1))
    for blocks_i in range(0,np.size(numblocklist,0)):
        numblocks=numblocklist[blocks_i]
        print("numblocks "+str(numblocks))
        r0,c0=realistic_initial(LA_matrix,LA_N,numblocks)
        LAr_fit, LAc_fit, wr, wc = gradient_descent(LA_matrix, r0, c0, N=numiter,sig_digs=4, N_stable=600, nblocks = numblocks, verbose_time = 100)
        LA_fit_sample.append((LAr_fit,LAc_fit))
        LA_ll[n_run,blocks_i]=likelihood(LAc_fit,LAr_fit,LA_matrix)  
    LA_fit.append(LA_fit_sample)
LA_maxind=np.argmax(LA_ll,0)
LA_llmax=np.max(LA_ll,0)
LA_AIC=2*LAnumparameters-2*LA_llmax
fig = plt.figure(figsize=(20,20))
plt.scatter(numblocklist,LA_AIC)
plt.xlabel("number of blocks",fontsize=20)
plt.ylabel("AIC",fontsize=20)
plt.savefig("img/LA/LAAIC.png")
print("Worst log-likelihoods"+str(np.min(LA_ll,0)))
print("Best log-likelihoods"+str(np.max(LA_ll,0)))
#Visualize 
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    LAr_fit=LA_fit[LA_maxind[blocks_i]][blocks_i][0]
    LAc_fit=LA_fit[LA_maxind[blocks_i]][blocks_i][1]
    
    LA_weights=np.sum(LA_matrix,(1,2))+np.sum(LA_matrix,(0,2))
    LA_sizes=LA_weights/15;
    LA_colors=LAc_fit/np.max(LAc_fit,0,keepdims=True)
    scatterplot_city(LAc_fit,LA_N,numblocks,np.array([list(LA_station_info[x]) for x in LA_stations])
                     ,sizes=LA_sizes,filename="img/LA/LA"+str(numblocks)+"blocks")
    block_omegas(LAr_fit,numblocks,filename="img/LA/LA"+str(numblocks)+"blocks")

#Save
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    pd.DataFrame(LA_fit[LA_maxind[blocks_i]][blocks_i][0].reshape([numblocks**2, 24],order = 'F')).to_csv("../mixed_model_results/LA_"+str(numblocks)+"_omega.csv",  index = False)
    pd.DataFrame(LA_fit[LA_maxind[blocks_i]][blocks_i][1], index = LA_stations).to_csv("../mixed_model_results/LA_"+str(numblocks)+"_roles.csv", index = LA_stations)


