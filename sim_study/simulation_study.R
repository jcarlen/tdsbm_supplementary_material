# OUTLINE
# 1. Run tdd simulation to evaluate parameter estimation
# 2. ppsbm to show impact of degree correction 
# 3. mixed-membership to show impact of potential model mis-specification
# ---------------------------------------------------------------------------------------------------------------
# TO DO
# - add generate_multilayer_array to package to facilitate simulation?
# - Note: degree correcton can also lead to a more parsimonious and interpretable model representation where there is degree heterogeneity because a unique class is not needed for each degree-activity level
#   add harder cases with smaller amplitude sin curves so hard to distinguish. if it gets harder try with more nodes.
# - More heterogeneous degree correcton case?
# - N = 120?
# ---------------------------------------------------------------------------------------------------------------
# new functions (generate_multilayer_array) ----

# generate N x N x T array (time-sliced adjacency matrices) based on TDD-SBM (type == "discrete) or TDMM-SBM (type = "mixed") model
# defaults to no degree correction (all dc_factors == 1)
# for type "discrete" roles are a length-N vector of block assignments. For type "mixed" roles are a G x N matrix of assignment weights (and dc_factor is ignored)
# (add this as a simulation function to abmt package?)
generate_multilayer_array <- function(N, Time, roles, omega, dc_factors = rep(1, N), type = "discrete") {
  # checks
  if (type == "discrete" & !identical(dc_factors, rep(1, N))) {
   if (! identical(aggregate(dc_factors ~ roles, FUN = "sum")[,2], rep(1, length(unique(roles))))) {
     warning("degree correction factors not normalized to sum 1 by group")
   }
  }
  
  edge_array = array(0, dim = c(N, N, Time))
  for (i in 1:N) {
    role_i = roles[i]
    for (j in 1:N) {
      role_j = roles[j]
      for (time in 1:Time) {
          if (type == "discrete") { ijt = rpois(dc_factors[i]*dc_factors[j]*omega[role_i, role_j, time], n= 1) }
          if (type == "mixed") { ijt = rpois(t(roles[, i])%*% omega[, , time] %*% roles[, j], n= 1) }
          edge_array[i, j, time] = ijt
      } 
    }
  }
  rownames(edge_array) = 1:N
  return(edge_array)
}


# libraries ----
# devtools::install_github("jcarlen/sbm", subdir = "sbmt") 
library(sbmt)
library(fossil) #for adj rand index

# ---------------------------------------------------------------------------------------------------------------
# 1. Run tdd simulation to evaluate parameter estimation ----
# omega curves ----

a = 10
b = 5 
Time = 16 #use a power of two for compatibility with ppsbm hist method
ymax =  max(2*a, 2*b)
x = seq(.5, Time)
  
# show two-block curves
par(mfrow = c(2, 2))
curve(b*sin(x*pi/Time) + b, 0, Time, ylim = c(0, ymax))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(b*sin(x*pi/Time)+ b, 0, Time, ylim = c(0, ymax))

# show three-block curves
par(mfrow = c(3, 3))
curve(b*sin(x*pi/Time) + b, 0, Time, ylim = c(0, ymax))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(a*sin(x*2*pi/Time-1)/2+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(b*sin(x*pi/Time)+ b, 0, Time, ylim = c(0, ymax))
curve(a*sin(x*2*pi/Time-2)/4+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time-1)/2+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time-2)/4+a, 0, Time, ylim = c(0, ymax))
curve(0*x+b, 0, Time, ylim = c(0, ymax))

# show overlapping
par(mfrow = c(1,1))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, add = TRUE, lty = 2)
curve(a*sin(x*2*pi/Time-1)/2+a, 0, Time, col = "blue", add = TRUE)
curve(-a*sin(x*2*pi/Time-1)/2+a, 0, Time, col = "blue", add = TRUE, lty = 2)
curve(a*sin(x*2*pi/Time-2)/4+a, 0, Time, col = "green", add = TRUE)
curve(-a*sin(x*2*pi/Time-2)/4+a, 0, Time, col = "green", add = TRUE, lty = 2)
curve(b*sin(x*pi/Time) + b, 0, Time, col = "red", add = TRUE)
curve(0*x+b, 0, Time, col = "red", add = TRUE)

omega_11 = b*sin(x*pi/Time) + b
omega_12 =  a*sin(x*2*pi/Time)+a
omega_21 = -a*sin(x*2*pi/Time)+a
omega_22 = omega_11
omega_13 =  a*sin(x*2*pi/Time-1)/2+a
omega_31 = -a*sin(x*2*pi/Time-1)/2+a
omega_23 =  a*sin(x*2*pi/Time-2)/4+a
omega_32 = -a*sin(x*2*pi/Time-2)/4+a
omega_33 = 0*x + b

omega_2 = array(rbind(omega_11, omega_21, omega_12, omega_22), dim = c(2, 2, Time)) #left most index moves fastest
omega_3 = array(rbind(omega_11, omega_21, omega_31,
                     omega_12, omega_22, omega_13,
                     omega_13, omega_23, omega_33), dim = c(3, 3, Time)) #left most index moves fastest

omega_list = list(omega_2, omega_3)

# parameters  ----

K= c(2,3)
Time = 16 #use a power of two for compatibility with ppsbm hist method
V = c(30, 60)
N_sim = 10
set.seed(1)

# tdd-sbm simulation function ----
#note, assumes directed = TRUE and selfEdges = TRUE
simulated_tdd <- function(N, n_roles, omega, Time, dc = 0, N_sim = 10, directed = TRUE, kl = 10) {
  
  #set roles and other params ----
  roles_discrete = rep(1:n_roles, length.out = N)
  names(roles_discrete) = 1:N
  
  # for degree corrected ----
  
  if (dc == 3) {
    
    dc_factor = seq(0,1,length.out = n_roles+1)[-1]
    #dc_factor = seq(0,1,length.out = n_roles+5)[-1]
    
    dc_factors = rep(dc_factor, each = round(N/n_roles))[1:N]
    names(dc_factors) = 1:N
    
    #apply sum to 1 constraint
    f.tmp = function(v) {v/sum(v)}
    dc_factors = as.vector(aggregate(dc_factors ~ roles_discrete, FUN = f.tmp)[,-1])
    
    # check 
    #identical(aggregate(dc_factors ~ roles_discrete, FUN = "sum")[,2], rep(1, n_roles))
  }
  
  # adjusted omega (to account for normalization in degree correction) is expected total weight from r to s instead of expected edge weight from single edge from block r to block s, 
  # useful for mixed omega also
  block_omega = omega*array(table(roles_discrete) %*% t(table(roles_discrete)), dim = c(n_roles, n_roles, Time))
  
  # ---------------------------------------------------------------------------------------------------------------
  # - fit sbmt ----
  
  role_results = 1:N_sim
  tdd_sbm_ari = 1:N_sim
  tdd_sbm_true = 1:N_sim
  tdd_sbm_sim = 1:N_sim
  
  for (s in 1:N_sim) {
    
    if (dc == 0) { # no degree correction case  ----
      discrete_edge_array = generate_multilayer_array(N, Time, roles_discrete, omega, type = "discrete")
    }
    if (dc == 3) { # degree correction case  ----
      discrete_edge_array = generate_multilayer_array(N, Time, roles_discrete, block_omega, dc_factors, type = "discrete")
    }
    discrete_edge_list = adj_to_edgelist(discrete_edge_array, directed = TRUE, selfEdges = TRUE)
    tdd_sbm = sbmt(discrete_edge_list, maxComms = n_roles, degreeCorrect = dc, directed = TRUE, klPerNetwork = kl)
    plot(tdd_sbm)
    # adjusted rand index
    tdd_sbm_ari[s] = adj.rand.index(tdd_sbm$FoundComms[order(as.numeric(names(tdd_sbm$FoundComms)))], roles_discrete)
    # compare likelihood for true vs. fit parameters for simulated data
    tdd_sbm_true[s] = 
      tdd_sbm_llik(discrete_edge_array, roles = roles_discrete-1, omega = block_omega, directed = TRUE, selfEdges = TRUE)
    tdd_sbm_sim[s] = 
      tdd_sbm_llik(discrete_edge_array, roles = tdd_sbm$FoundComms, omega = tdd_sbm$EdgeMatrix, directed = TRUE, selfEdges = TRUE)
  }
  
  # - results ----
  results = list(
    # role detection
    tdd_sbm_ari = tdd_sbm_ari,
    # compare likelihood of data under fit vs. true model
    tdd_sbm_true = tdd_sbm_true,
    tdd_sbm_sim = tdd_sbm_sim,
    tdd_true_vs_sim = tdd_sbm_true - tdd_sbm_sim
    # add omega accuracy with curve distances?
  )
  
  return(results)
}

# run tdd-sbm simulation ----

td_results = vector(mode = "list", length = length(K))

td_results = 
  lapply(1:length(K), function(i) {
  n_roles = K[i]
  omega = omega_list[[i]]
    lapply(V, function(N) {
      tdd_results_dc0 = simulated_tdd(N, n_roles, omega, Time, dc = 0, N_sim, kl = 1) 
      tdd_results_dc3 = simulated_tdd(N, n_roles, omega, Time, dc = 3, N_sim, kl = 1) 
    })
  })


# ---------------------------------------------------------------------------------------------------------------
# 2. ppsbm to show impact of degree correction ----
library(ppsbm)

# no degree correction case. their model works as expected ----
# - fit ppsbm ----

# Use the "hist" method because agrees more closely with out discrete time slices and requires little data manipulation

Nijk = sapply(adj_to_edgelist(discrete_edge_array, directed = TRUE, selfEdges = FALSE, removeZeros = FALSE), "[[", 3); dim(Nijk)
discrete_ppsbm = mainVEM(list(Nijk=Nijk,Time=Time), N, Qmin = 1, Qmax = 4, directed=TRUE, 
                         method='hist', d_part=5, n_perturb=10, n_random=0)

# - results ----

# number of blocks selected
selected_Q = modelSelection_Q(list(Nijk=Nijk,Time=Time), N, Qmin = 1, Qmax = 4, directed = TRUE, sparse = FALSE, discrete_ppsbm)$Qbest
selected_Q
selected_Q == n_roles #should equal n_roles

# role detection
apply(discrete_ppsbm[[selected_Q]]$tau, 2, which.max)
adj.rand.index(apply(discrete_ppsbm[[selected_Q]]$tau, 2, which.max), roles_discrete)

# omegas
par(mfrow = c(n_roles, n_roles))
apply(exp(discrete_ppsbm[[selected_Q]]$logintensities.ql), 1, plot, type = "l")

# degree correction case ----
# - fit ppsbm ----

dc_Nijk = sapply(adj_to_edgelist(dc_discrete_edge_array, directed = TRUE, selfEdges = FALSE, removeZeros = FALSE), "[[", 3); dim(dc_Nijk)
dc_discrete_ppsbm = mainVEM(list(Nijk=dc_Nijk,Time=Time), N, Qmin = 1, Qmax = 6, directed=TRUE, 
                            method='hist', d_part=5, n_perturb=10, n_random=0)
# - results ----

# number of blocks selected
selected_Q = modelSelection_Q(list(Nijk=dc_Nijk,Time=Time), N, Qmin = 1, Qmax = 6, directed = TRUE, sparse = FALSE, dc_discrete_ppsbm)$Qbest
selected_Q
selected_Q == n_roles

# role detection

# Wants a seperate class for each degree-correcton level
apply(dc_discrete_ppsbm[[selected_Q]]$tau, 2, which.max)
adj.rand.index(apply(dc_discrete_ppsbm[[selected_Q]]$tau, 2, which.max), roles_discrete)

# omegas
par(mfrow = c(selected_Q, selected_Q)); par(mai = rep(.5, 4))
apply(exp(dc_discrete_ppsbm[[selected_Q]]$logintensities.ql), 1, plot, type = "l", col = "blue")

# with true number of groups? gets it right
apply(dc_discrete_ppsbm[[n_roles]]$tau, 2, which.max)
adj.rand.index(apply(dc_discrete_ppsbm[[n_roles]]$tau, 2, which.max), roles_discrete)

par(mfrow = c(n_roles, n_roles)); par(mai = rep(.5, 4))
apply(exp(dc_discrete_ppsbm[[n_roles]]$logintensities.ql), 1, plot, type = "l", col = "blue")

# Bike example shows how degree correction ib model -> group statins with similar behavior across activity levels
# ---------------------------------------------------------------------------------------------------------------
# 3. mixed-membership to show impact of potential model mis-specification ----

roles_mixed = matrix(c(rep(.5, 2*N/3), rep(c(0,1), N/3), rep(c(1,0), N/3)), nrow = n_roles, ncol = N); roles_mixed

#apply sum to 1 constraint to mixed roles
roles_mixed = roles_mixed/rowSums(roles_mixed)

tdmm_sbm_role_sum_abs_error = 1:N_sim #use sum because values are normalized
tdmm_sbm_omega_mean_abs_error = 1:N_sim
tdmm_sbm_ll =1:N_sim
tdd_sbm_ll =1:N_sim

# assume starting from tdsbm_supplementary_material directory
setwd("mixed_model_implementation_python")

for (s in 1:N_sim) {
  
  mixed_edge_array = generate_multilayer_array(N, Time, roles_mixed, block_omega, type = "mixed")
  write.csv(mixed_edge_array, "../data/sim/mixed_edge_array.csv", row.names = FALSE)
  # run mixed
  system("python3 tdmm_sbm_sim_study.py")
  tdmm_sbm_roles_2 = read.csv("../mixed_model_results/SIM_2_roles.csv")
  tdmm_sbm_role_sum_abs_error[s] = min(sum(abs(roles_mixed - t(tdmm_sbm_roles_2))), sum(abs(roles_mixed - t(tdmm_sbm_roles_2)[2:1,]))) #try both orders
  tdmm_sbm_omega_2 = read.csv("../mixed_model_results/SIM_2_omega.csv")
  tdmm_sbm_omega_mean_abs_error[s] = min(mean(abs(unlist(tdmm_sbm_omega_2 - apply(block_omega, 3, as.vector)))), 
                                      mean(abs(unlist(tdmm_sbm_omega_2 - apply(block_omega, 3, function(x) {as.vector(t(x))})))))
  
  # try fitting discrete models to data generated from mixed membership
  mixed_edge_edglist = adj_to_edgelist(mixed_edge_array, directed = TRUE, selfEdges = TRUE, removeZeros = TRUE)
  tdd_sbm = sbmt(mixed_edge_edglist, maxComms = n_roles, degreeCorrect = 3, directed = TRUE, klPerNetwork = 10)
  #plot(tdd_sbm)
  # summarize error by finding likelihood of mixed-generated data fit with  discrete model vs. with mixed params.
  
  # evaluate liklihood of mixed edge array using params found by discrete model
  tdd_roles = matrix(0, N, n_roles); for (i in 1:N) {tdd_roles[i, tdd_sbm$FoundComms[i]+1] = 1}; rownames(tdd_roles) = 1:N
  tdd_roles = t(t(tdd_roles)/colSums(tdd_roles))
  tdd_omega = sapply(tdd_sbm$EdgeMatrix, as.vector)
  # discrete
  tdd_sbm_ll[s] = tdmm_sbm_llik(mixed_edge_array, C = tdd_roles, omega = tdd_omega, selfEdges = TRUE, directed = TRUE)
  # vs mixed
  tdmm_sbm_ll[s] = tdmm_sbm_llik(mixed_edge_array, C = tdmm_sbm_roles_2, omega = tdmm_sbm_omega_2, selfEdges = TRUE, directed = TRUE)
  
}

tdmm_vs_tdd_sbm_ll = tdmm_sbm_ll - tdd_sbm_ll #(scale by difference in parameters?)

setwd("..")


