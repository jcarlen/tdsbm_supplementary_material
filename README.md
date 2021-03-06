Description of the contents of tdsbm\_supplementary\_material:

The material in this folder is hosted at: https://github.com/jcarlen/tdsbm_supplementary_material

1) tdsbm\_data\_plots.R 

- This is the main replication script of the paper.  It contains all code to clean the data as described in the paper, estimate model parameters and create all plots in the paper. (It has commented out code to call the TDMM-SBM estimation (in Python) directly from R, but this code can also be accessed directly in the mixed\_model\_implementation\_python folder.)

- At beginning of the script, set your working directory as the path to your tdsbm\_supplementary\_material folder by modifying the line:
    
    ````setwd("~/Documents/285J/tdsbm_supplementary_material")````
    
- One way to run this script in the terminal is to navigate to the tdsbm\_supplementary\_material folder and then (assuming you have R installed) enter 

````> Rscript tdsbm_data_plots.R````
   
- Output from the discrete (TDD-SBM) models run by this script is stored in the discrete\_model\_results\_folder.
    - If models are very time consuming to run (for New York City), code to run the model has been commented out and pre-run model results are loaded instead from the discrete\_model\_results folder.

2) IMG

-  Destintion folder for saved plots generated by tdsbm\_data\_plots.R. (This includes more plots than in the paper.) 

3) data 

- LA, SF, NY - city specific trip histories and (for LA and SF) station data tables.

- maps - Background maps for some plots. Code to generate the maps is in the replication script, but commented out because it requires an individual google API key.

- zoning - Shape files for LA and New York City zoning maps used as background for some figures.

- cleaned - This is a folder where cleaned versions of the data sets will be stored once created in the tdsbm\_data\_plots.R script. Cleaned data sets without weekend trips for LA, SF and our Manhattan subnetwork of New York City are given as examples, but additional data sets will be populated by running the tdsbm\_data\_plots.R script. Note the complete set of cleaned data sets is about 400mb. Also note that the python replication script for the TDMM-SBM models, tdmm\_sbm\_model\_replication.py,  uses the cleaned data (without weekends) as input, so it's recommended to run the tdsbm\_data\_plots.R scripts before tdmm\_sbm\_model\_replication.py. 


4) mixed\_model\_implementation\_python

- tdmm\_sbm\_model\_replication.py - This contains all the code to estimate parameters of the TDMM-SBM models referenced in the paper.
- For convenience the code in this file is also saved in files for each city in the paper: LosAngeles.py, Manhattan.py, NewYork.py, and SanFrancisco.py.  
    -  Output from the mixed-membership (TDMM-SBM) models run in this script is stored in the mixed\_model\_results\_folder.
    -  *.png files are produced for visualizing the models and are saved in city-specific folders of an "img" folder that is created in mixed\_model\_implementation\_python.
- The scripts are set up to be run from the mixed\_model\_implementation\_python folder.
    
- lib contains the Python package of functions to estimate parameters  of the TDMM-SBM and calculate likelihood of a given model as described in the paper.
        
5) discrete\_model\_results

- Various output from TDD-SBM models labeled by city (note that ny\_hm refers to the Manhattan subnetwork of the New York City network) and number of blocks, which is the first number in each name. The second number in each name,  3, indicates that the type of degree correction applied was the type described in the paper. Results with names ending in \_T indicated time-aggregated data (i.e., the results of a time-independent SBM).  All objects in this folder are in R data format (.RDS), which can be loaded into R or Python.  

6) mixed\_model\_results

- Output from the TDMM-SBM models run by the tdmm\_sbm\_model\_replication.py script. Each model output is listed by city (note that Manhattan refers to the Manhattan subnetwork of the New York City network) and number of blocks. Each model run stores both a "role" and "omega" object. The role object contains the estimated C parameters for the model, listed with rownames equal to the corresponding station ID. The omega object is a (K\*K) x T matrix of the time-dependent block-to-block activity parameters, where the K x K matrix of block-to-block parameters for an individual time step has been converted to a column (of length K\*K) of the matrix.

7) rstan\_mixed\_implementation

- Code in R which calls the Stan language to perform an alternate (Hamiltonian MCMC) estimation of our TDMM-SBM models. This is meant as a check on the output of our gradient-descent mixed-model estimation implementation in Python.  
    
- The tdmm\_sbm\_stan.R script implements this estimation method for several two-block TDMM-SBM models. The script uses data objects created in the main replication script, tdsbm\_data\_plots.R. 
    
- The tdmm\_sbm.stan script is called by tdmm\_sbm\_stan.R, and does not need to be accessed directly.
    
- stan\_results contains the estimated omega and C parameters for two-block TDMM-SBM of LA, SF, and the Manhattan subnetwork of New York City. (We did not complete estimation in R/Stan of the entire NYC network with an adequate sample size.)



