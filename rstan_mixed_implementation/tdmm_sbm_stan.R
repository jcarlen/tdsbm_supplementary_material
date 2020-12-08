# Code to fit the TDMM-SBM uses Stan's HMC algorithm
# Also plots the output and compares to Python gradient descent output

# Jane Carlen
# Created: 11-5-18

# Notes:
#   Data needed by this script, e.g. la_byhour, created in tdsbm_data_plots.R
#   Necessary packages are also installed (as needed) in that script. 
#
# TO DO
#     Smarter initializeation
#     See if Stan model can recover the generating values of the parameters when fit to this dataset.
# 
######################################################################################################  

#Set to the path to your tdsbm_supplementary_material folder:
setwd("~/Documents/285J/tdsbm_supplementary_material")

# O. Setup ---------------------------------------------

library(sbmt)
library(rstan)
library(MASS)
library(boot)
library(dplyr)
library(ggplot2)
library(reshape2)

Sys.setenv(USE_CXX14 = 1)
rstan_options(auto_write = T)
options(mc.cores = ceiling(parallel::detectCores()/2))

# 1. Data  -------------------------------------
# See data processed in from tdsbm_data_plots.R

#   la ####

# Convert hourly trips to hourly matrices
A_la = edgelist_to_adj(la_byhour, selfEdges = TRUE, as.array = TRUE, directed = TRUE)

# Data for Stan 
la_data <- list(
  N = dim(A_la)[1],
  K = 2,
  Time = length(la_byhour),
  A = A_la,
  selfEdges = as.numeric(TRUE)
)

#   sf ####

# Convert hourly trips to hourly matrices
A_sf = edgelist_to_adj(sf_byhour, selfEdges = TRUE, as.array = TRUE, directed = TRUE)

# Data for Stan 
sf_data <- list(
  N = dim(A_sf)[1],
  K = 2,
  Time = length(sf_byhour),
  A = A_sf,
  selfEdges = as.numeric(TRUE)
)

#   ny_hm (manhattan (home) subset) ####

# Convert hourly trips to hourly matrices
A_ny_hm = edgelist_to_adj(ny_byhour_hm, selfEdges = TRUE, as.array = TRUE, directed = TRUE)

# Data for Stan 
ny_hm_data <- list(
  N = dim(A_ny_hm)[1],
  K = 2,
  Time = length(ny_byhour_hm),
  A = A_ny_hm,
  selfEdges = as.numeric(TRUE)
)

#   ny (too big) ####
# #
# # Convert hourly trips to hourly matrices
# A_ny = edgelist_to_adj(ny_byhour, selfEdges = FALSE, as.array = TRUE, directed = TRUE)
# #
# # Data for Stan 
# ny_data <- list(
#   N = dim(A_ny)[1],
#   K = 2,
#   Time = length(ny_byhour),
#   A = A_ny,
#   selfEdges = as.numeric(TRUE)
# )

# 2. RUN -------------------------------------

#outputs of these runs are saved as fit_la_stan.RDS, fit_sf_stan.RDS, fit_ny_hm_stan.RDS, fit_ny_stan_init_run.RDS, fit_ny_stan.RDS

stanc("rstan_mixed_implementation/tdmm_sbm.stan")

fit.la = stan("rstan_mixed_implementation/tdmm_sbm.stan", model_name = "la", data = la_data, chains = 2, seed = 397031630, #had to recover seed
              iter = 10000, save_warmup = FALSE, sample_file = "rstan_mixed_implementation/lasample.csv", thin = 10)

fit.sf = stan("rstan_mixed_implementation/tdmm_sbm.stan", model_name = "sf", data = sf_data, chains = 2, seed = 1471495729,
              iter = 10000, save_warmup = FALSE, sample_file = "rstan_mixed_implementation/sfsample.csv", thin = 10)

fit.ny.hm = stan("rstan_mixed_implementation/tdmm_sbm.stan", model_name = "ny.hm", data = ny_hm_data, 
                  control = list(adapt_delta = 0.9), chains = 2, seed = 2,
                  pars = c("C","omega"), iter = 5000, save_warmup = FALSE, 
                  sample_file = "rstan_mixed_implementation/nyhmsample.csv", thin = 10) #about 3 hours with 2 cores

# Run below is commented out because it's very slow, about 24 hours with relatively  few iterations
# It was initialized with the results from a single chain with 25 samples (300 iter, 50 warmup, thin 10, seed 1783015463)
# cores1 = round(parallel::detectCores()/2) #6 on my machine
# fit.ny = stan("tdmm_sbm.stan", model_name = "ny", data = ny_data, chains = cores1, cores = cores1,
#              seed = 1783015463,
#              pars = c("C","omega"),
#              iter = 300, warmup = 100, save_warmup = FALSE, sample_file = "nysample.csv", thin = 10,
#              init = list(list(C=C, omega = ny_omega), list(C=C, omega = ny_omega), list(C=C, omega = ny_omega),
#                          list(C=C, omega = ny_omega), list(C=C, omega = ny_omega), list(C=C, omega = ny_omega)))

# 3. Output -------------------------------------
# can look bad if chains have switched labellings
# All values for R̂  should be ~1, less than 1.1.

#   la ----
#     extract results ----
la_extract = rstan::extract(fit.la)
la_summary <- summary(fit.la)
summary.function = median #found higher likelihood with median than mean

#omega
la_omega = apply(la_extract$omega, 2:4, summary.function)
par(mfrow = c(2,2))
plot(la_omega[,1,1], type = "l")
plot(la_omega[,1,2], type = "l")
plot(la_omega[,2,1], type = "l")
plot(la_omega[,2,2], type = "l")
# C
la_C = apply(la_extract$C, c(2,3), summary.function)
colnames(la_C) = rownames(A_la)

#     plot with meaningful node size and color ---- 
la.output = left_join(la.station, data.frame(id = rownames(A_la), color = la_C[1,]/colSums(la_C), size = colSums(la_C))) 

ggplot(data = la.output, aes(y = lat, x = lon, color = color, size = size)) + #size = degree
  geom_point() 

#     check results agree with gradient descent ----
cor(la.output$color, la_continuous.roles$X1/(la_continuous.roles$X0+la_continuous.roles$X1)) # 0.999737

# Extremely similar likelihood (grad. desc. teeny bit better)
K = 2
tdmm_sbm_llik(A = A_la, t(la_C), array(apply(la_omega, 1, matrix), c(K,K,24))) #stan
tmp.roles = data.frame(la_continuous.roles[,-1], row.names =  la_continuous.roles[,1]) 
tdmm_sbm_llik(A = A_la, C = tmp.roles, omega = array(unlist(la_continuous.omega))) #grad desc.


#     save ----
la_continuous_stan = list(omega = la_omega, station_params = data.frame(id = rownames(A_la), t(la_C), stringsAsFactors = FALSE))
saveRDS(la_continuous_stan, "rstan_mixed_implementation/stan_results/la_continuous.RDS")

#   sf ----

#     extract results ----
sf_extract = rstan::extract(fit.sf)
sf_summary <- summary(fit.sf)
summary.function = median

#omega
sf_omega = apply(sf_extract$omega, 2:4, summary.function)
par(mfrow = c(2,2))
plot(sf_omega[,1,1], type = "l")
plot(sf_omega[,1,2], type = "l")
plot(sf_omega[,2,1], type = "l")
plot(sf_omega[,2,2], type = "l")
# C
sf_C = apply(sf_extract$C, c(2,3), summary.function)
colnames(sf_C) = rownames(A_sf)

#     plot with meaningful node size and color ----

sf.output = left_join(sf.station, data.frame(id = as.numeric(rownames(A_sf)), color = sf_C[1,]/colSums(sf_C), size = colSums(sf_C))) 

ggplot(data = sf.output, aes(y = lat, x = lon, color = color, size = size)) + #size = degree
  geom_point()

#     check results agree with gradient descent ----
cor(sf.output$color, sf_continuous.roles$X1/(sf_continuous.roles$X0+sf_continuous.roles$X1)) #0.99997

# Extremely similar likelihood (grad. desc. teeny bit better)
K = 2
tdmm_sbm_llik(A = A_sf, t(sf_C), array(apply(sf_omega, 1, matrix), c(K,K,24))) #stan
tmp.roles = data.frame(sf_continuous.roles[,-1], row.names =  sf_continuous.roles[,1]) 
tdmm_sbm_llik(A = A_sf, tmp.roles, array(unlist(sf_continuous.omega))) #grad desc.

#     save ----

sf_continuous_stan = list(omega = sf_omega, station_params = data.frame(id = rownames(A_sf), t(sf_C), stringsAsFactors = FALSE))

saveRDS(sf_continuous_stan, "rstan_mixed_implementation/stan_results/sf_continuous.RDS")

#   manhattan subnetwork ----

#     extract results ----
ny_hm_extract = rstan::extract(fit.ny.hm)
ny_hm_summary <- summary(fit.ny.hm)
summary.function = median

#omega
ny_hm_omega = apply(ny_hm_extract$omega, 2:4, summary.function)
par(mfrow = c(2,2))
plot(ny_hm_omega[,1,1], type = "l")
plot(ny_hm_omega[,1,2], type = "l")
plot(ny_hm_omega[,2,1], type = "l")
plot(ny_hm_omega[,2,2], type = "l")
# C
ny_hm_C = apply(ny_hm_extract$C, c(2,3), summary.function)
colnames(ny_hm_C) = rownames(A_ny_hm)

#     plot with meaningful node size and color ----

ny_hm.output = left_join(ny1610_hm.station, data.frame(id = as.numeric(rownames(A_ny_hm)),
                                                       color = ny_hm_C[2,]/colSums(ny_hm_C), 
                                                       size = colSums(ny_hm_C))) 

ggplot(data = ny_hm.output, aes(y = lat, x = lon, color = color, size = size)) + #size = degree
  geom_point() 

#     check results agree with gradient descent ----
cor(ny_hm.output$color, ny_hm_continuous.roles.2$X1/(ny_hm_continuous.roles.2$X0+ny_hm_continuous.roles.2$X1))
plot(ny_hm.output$color, ny_hm_continuous.roles.2$X1/(ny_hm_continuous.roles.2$X0+ny_hm_continuous.roles.2$X1))


# almost identical likelihood to gradient descent's LOCAL max result
tdmm_sbm_llik(A = A_ny_hm, t(ny_hm_C), array(apply(ny_hm_omega, 1, matrix), c(2,2,24))) #stan
tmp.roles = data.frame(ny_hm_continuous.roles.2[,-1], row.names =  ny_hm_continuous.roles.2[,1])
tdmm_sbm_llik(A = A_ny_hm, tmp.roles, array(unlist(ny_hm_continuous.omega.2))) #grad desc.

#     save ----

ny_hm_continuous_stan = list(omega = ny_hm_omega, station_params = data.frame(id = rownames(A_ny_hm), 
                                t(ny_hm_C), stringsAsFactors = FALSE)) #switch c column order for image coloring later
saveRDS(ny_hm_continuous_stan, "rstan_mixed_implementation/stan_results/ny_hm_continuous.RDS")

# # ny (commented out -- takes too long to run and still get small samples) ####
# ny_extract = rstan::extract(fit.ny)
# ny_summary <- summary(fit.ny)
# ny_summary$c_summary[,,1]
# summary.function = median #for exact colsums
# 
# # These use all chains, but some were divergent se we'll use a subset.
# ny_omega = apply(ny_extract$omega, 2:4, summary.function)
# C = apply(ny_extract$C, c(2,3), summary.function) 
# 
# # take only the chains that seem to align
# C2 = apply(cbind(ny_summary$c_summary[,1,2], ny_summary$c_summary[,1,3]), 1, summary.function)
# 
# # Compare effective sample sizes when chains removed
# plot(apply(as.matrix(fit.ny)[c(21:60, 81:100),], 2, effectiveSize), apply(as.matrix(fit.ny)[21:100,], 2, effectiveSize))
# table(round(apply(as.matrix(fit.ny)[c(21:60, 81:100),], 2, effectiveSize)))


