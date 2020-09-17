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
numiter=3000

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
                     ,sizes=LA_sizes,filename="img/LA/LA"+str(numblocks)+"blocks"+str(numiter)+"iter")
    block_omegas(LAr_fit,numblocks,filename="img/LA/LA"+str(numblocks)+"blocks"+str(numiter)+"iter")

#Save
for blocks_i in range(0,np.size(numblocklist,0)):
    numblocks=numblocklist[blocks_i]
    pd.DataFrame(LA_fit[LA_maxind[blocks_i]][blocks_i][0].reshape([numblocks**2, 24], order = 'F')).to_csv("../mixed_model_results/LA_"+str(numblocks)+"_omega.csv",  index = False)
    pd.DataFrame(LA_fit[LA_maxind[blocks_i]][blocks_i][1], index = LA_stations).to_csv("../mixed_model_results/LA_"+str(numblocks)+"_roles.csv", index = LA_stations)

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
    pd.DataFrame(Manhattan_fit[Manhattan_maxind[blocks_i]][blocks_i][0].reshape([numblocks**2, 24], order = 'F')).to_csv("../mixed_model_results/Manhattan_"+str(numblocks)+"_omega.csv",  index = False)
    pd.DataFrame(Manhattan_fit[Manhattan_maxind[blocks_i]][blocks_i][1], index = Manhattan_stations).to_csv("../mixed_model_results/Manhattan_"+str(numblocks)+"_roles.csv", index = Manhattan_stations)

# San Fransisco models
###############################################################################################################
# Load SF Data
    
# Assumes you're in the mixed_model_implementation_python folder

# A maintained and documented version of the lib can be found here https://github.com/jaumededios/Bike_Networks
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

