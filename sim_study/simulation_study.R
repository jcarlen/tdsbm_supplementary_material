# OUTLINE
# 1. Run tdd simulation to evaluate parameter estimation
# 2. ppsbm to show impact of degree correction 
# 3. mixed-membership to show impact of potential model mis-specification
# ---------------------------------------------------------------------------------------------------------------
# TO DO
# - add generate_multilayer_array to sbmt package to facilitate simulation?
#     + implement & test for directed networks
# - Add harder cases with smaller amplitude sin curves so hard to distinguish? if it gets harder try with more nodes.
# - More heterogeneous degree correcton case?
# - Note: degree correcton can also lead to a more parsimonious and interpretable model representation where there is degree heterogeneity 
#   because a unique class is not needed for each degree-activity level.
# ---------------------------------------------------------------------------------------------------------------
# libraries ----
# devtools::install_github("jcarlen/sbm", subdir = "sbmt") 
library(sbmt)
library(fossil) #for adj rand index
library(gtools) #for permutations

# new functions (generate_multilayer_array) ----

#' generate N x N x Time array (Time-sliced adjacency matrices) based on TDD-SBM (type == "discrete) or TDMM-SBM (type = "mixed") model
#' @param roles If type is `discrete`, a length-N vector of block assignments. If type is `mixed` G x N matrix of block-assignment weights.
#' @param omega A K x K x Time array or Time-length list of K x K matrices. If the latter will be internall converted to array.
#' If no degree-correction (`dc_factors equal NULL or rep(1, N)`), omega represents mean block-to-block edge values. 
#' With degree-correction (dc_factors differ from default) omega represent the sum of block-to-block edges.
#' @param type `discrete` or `mixed` for tdd-sbm or tdmm-sbm, respectively. If mixed, dc_factors is ignored.
#' @param dc_factors 
# '   if not degree-corrected, should be NULL or rep(1, N)
# '   otherwise, dc_factors are a list of degree-correction factors that should be normalized to sum to 1 within blocks
# for type "mixed" roles are a G x N matrix of assignment weights 
# (add this as a simulation function to sbmt package?)
#' @examples 
#' generate_multilayer_array(roles = rep(1:2, 5), omega = array(1:12, dim = c(2,2,3)), dc_factors = NULL, type = "discrete") #output 10 x 10 x 3
#' generate_multilayer_array(roles = rep(1:2, 5), omega = array(1:12, dim = c(2,2,3)), dc_factors = 1:10, type = "discrete") #output 10 x 10 x 3
generate_multilayer_array <- function(roles, omega, dc_factors = NULL, type = "discrete") {
  
  #  infer N, K, Time,
  N = ifelse(is.null(dim(roles)), length(roles), nrow(roles))
  K = ifelse(class(omega) == "list", nrow(omega[[1]]), dim(omega)[1]) # number of blocks
  Time = ifelse(class(omega) == "list", length(omega[[1]]), dim(omega)[3])
  if(class(omega)!="array") {omega  = array(unlist(omega), dim = c(K,K,Time))} #make array if not already
  #  degree-corrected?
  if(is.null(dc_factors)) {dc_factors = rep(1, N)}
  dc = ifelse(identical(dc_factors, rep(1, N)), FALSE, TRUE)
         
  # checks
  if (length(unique(roles)) <= 1) stop("There should be at least two unique roles")
  if (type == "discrete" & dc) {
    if (! identical(aggregate(dc_factors ~ roles, FUN = "sum")[,2], rep(1, length(unique(roles))))) {
      warning("degree correction factors not normalized to sum to 1 by block. will normalize now.")
      #normalize?
      role_sums = aggregate(dc_factors, by = list(roles), sum)
      role_sums = setNames(role_sums$x, role_sums$Group.1)
      dc_factors = dc_factors/role_sums[roles] #group sum constraint
      cat("dc_factors normalized:", dc_factors, "\n")
      # check all(aggregate(dc_factors, by = list(roles), sum)$x == 1)
    }
  }
      
  # adjusted omega (to account for normalization in degree correction) 
  # entries are expected total weight from r to s instead of expected edge weight from single edge from block r to block s
  if(dc) {block_omega = omega*array(table(roles) %*% t(table(roles)), dim = c(K, K, Time))}
  
  edge_array = array(0, dim = c(N, N, Time))
  for (i in 1:N) {
    role_i = roles[i]
    for (j in 1:N) {
      role_j = roles[j]
      for (Time in 1:Time) {
          if (type == "discrete" & dc) { ijt = rpois(dc_factors[i]*dc_factors[j]*block_omega[role_i, role_j, Time], n= 1) }
          if (type == "discrete" & !dc) { ijt = rpois(dc_factors[i]*dc_factors[j]*omega[role_i, role_j, Time], n= 1) }
          if (type == "mixed") { ijt = rpois(t(roles[, i])%*% omega[, , Time] %*% roles[, j], n= 1) }
          edge_array[i, j, Time] = ijt
      } 
    }
  }
  rownames(edge_array) = 1:N
  return(edge_array)
}

# ---------------------------------------------------------------------------------------------------------------
# 0. Parameters ----
#   - K_set (# block options), N_set (# node options), Time, n_sim, seed ----

K_set = c(2,3)
N_set = c(30, 90)
Time = 16 #use a power of two for compatibility with ppsbm hist method
N_sim = 10
set.seed(1)
kl_per_network = 10

#   - omega ----

a = 10
b = 5 
ymax =  max(2*a, 2*b)
x = seq(.5, Time)

# show two-block curves
par(mfrow = c(2, 2))
curve(b*sin(x*pi/Time) + b, 0, Time, ylim = c(0, ymax))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(b*sin(x*pi/Time)+ b, 0, Time, ylim = c(0, ymax))

# show three-block curves
T = Time #set just for labeling, then remove
png("IMG/sim_omega_3.png", width = 600)
par(mai =c(.4,.7,.3,.2))
par(mfrow = c(3, 3))
curve(b*sin(x*pi/T) + b, 0, T, ylim = c(0, ymax), main = ("1 to 1"), xlab = "")
curve(a*sin(x*2*pi/T)+a, 0, T, ylim = c(0, ymax), main = ("1 to 2"), xlab = "")
curve(a*sin(x*2*pi/T-1)/2+a, 0, T, ylim = c(0, ymax), main = ("1 to 3"), xlab = "")
curve(-a*sin(x*2*pi/T)+a, 0, T, ylim = c(0, ymax), main = ("2 to 1"), xlab = "")
curve(b*sin(x*pi/T)+ b, 0, T, ylim = c(0, ymax), main = ("2 to 2"), xlab = "")
curve(a*sin(x*2*pi/T-2)/4+a, 0, T, ylim = c(0, ymax), main = ("2 to 3"), xlab = "")
curve(-a*sin(x*2*pi/T-1)/2+a, 0, T, ylim = c(0, ymax), main = ("3 to 1"), xlab = "")
curve(-a*sin(x*2*pi/T-2)/4+a, 0, T, ylim = c(0, ymax), main = ("3 to 2"), xlab = "")
curve(0*x+b, 0, T, ylim = c(0, ymax), main = ("3 to 3"), xlab = "")
dev.off()
rm(T)

# show overlapping
# par(mfrow = c(1,1))
# curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
# curve(-a*sin(x*2*pi/Time)+a, 0, Time, add = TRUE, lty = 2)
# curve(a*sin(x*2*pi/Time-1)/2+a, 0, Time, col = "blue", add = TRUE)
# curve(-a*sin(x*2*pi/Time-1)/2+a, 0, Time, col = "blue", add = TRUE, lty = 2)
# curve(a*sin(x*2*pi/Time-2)/4+a, 0, Time, col = "green", add = TRUE)
# curve(-a*sin(x*2*pi/Time-2)/4+a, 0, Time, col = "green", add = TRUE, lty = 2)
# curve(b*sin(x*pi/Time) + b, 0, Time, col = "red", add = TRUE)
# curve(0*x+b, 0, Time, col = "red", add = TRUE)

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

#   - omega list ----

#these versions of omega give the expected degree (without any degree correcton) of any edge from block r to block s at time t
omega_list = list(omega_2, omega_3)

# 1. Run tdd simulation to evaluate parameter estimation ----

# tdd-sbm simulation function ----

# generate discrete and mixed roles for use in simulation.
# currently assumes directed = TRUE and selfEdges = TRUE
# (custom roles can also be fed to the simulation)
# roles: (discrete case) is an integer, indicating the possible blocks (roles) a user can have
#        (mixed case) a matrix where each row is a C_i vector of block weights, indicating realized mixtures of blocks. (A small set for the purposes of simulation.)
# rel_freq: dictates relative frequency of the blocks (discrete case) or block-weight vectors (mixed cased), defaulting to as close as possible to equal groups

generate_roles <- function(N, roles, type = c("discrete", "mixed")[1],
                           rel_freq = ifelse(type == "discrete", list(rep(1, roles)), list(rep(1, nrow(roles))))[[1]]) {
  if (!type %in% c("discrete", "mixed")) {stop("type must be discrete or mixed")}
  rel_freq = rel_freq/sum(rel_freq)
  if (type == "discrete") {
    roles = rep(1:roles, times = ceiling(N*rel_freq))[1:N] #only up to first N elements in case of rounding
  } else { #type mixed
    times = ceiling(rel_freq*N)
    roles = roles[rep(1:length(times),times),][1:N,] #only up to first N elements in case of rounding
  }
  return(roles)
}

  # examples:
  #   generate_roles(30, roles = 2, type = "discrete", rel_freq = rep(1, 2))
  #   generate_roles(30, roles = 2, type = "discrete", rel_freq = c(.25,.75))
  #   generate_roles(30, roles = 3, type = "discrete", rel_freq = c(1/6,2/6,3/6))
  #   n_roles = 3
  #   mixed_roles = matrix(c(0,.25,.75)[permutations(n_roles, n_roles)],factorial(n_roles),n_roles)
  #   generate_roles(30, roles = mixed_roles, type = "mixed", rel_freq = rep(1,6))
  #   generate_roles(30, roles = mixed_roles, type = "mixed", rel_freq = 1:6)

# roles are the discrete roles of the nodes (which also tells us N)
# dc_levels is preset relative activity level we use for degree correction
# dc_weights is the fraction of users in at each level (to be rounded to integers)
generate_dc_factors <- function(roles, dc_levels, dc_weights) {
  
  dc_factors = lapply(table(roles), rep_len, x=dc_levels)
  dc_factors = lapply(dc_factors, function(x) {x/sum(x)}) #group sum constraint
  dc_factors = as.vector(do.call("c", dc_factors))
  return(dc_factors)
}

# A function to generate and fit data from the tdd-sbm
# currently for directed networks only
# roles: N-vector with element in 1:K (K is number of blocks)
# omega: K x K x T array
# dc_factors: degree correction factors. If NULL or 
# verbose: TRUE will show some intermediate plots


simulated_tdd <- function(roles, omega, dc_factors = NULL, fit_method = c("tdd-sbm-0", "tdd-sbm-3", "ppsbm")[1], N_sim = 10, directed = TRUE, kl = kl_per_network, verbose = FALSE) {
  
  if (!directed) stop("Not yet implemented for undirected networks")
  #  set degree correction parameter based on presence of dc_factors and dc_fit
  dc = ifelse(is.null(dc_factors), 0, 3)
  dc_fit = ifelse(is.null(dc_fit), dc, dc_fit)
  
  #  infer  K, Time
  K = dim(omega)[2] # number of bloks
  Time = dim(omega)[3]

  # ---------------------------------------------------------------------------------------------------------------
  # - fit sbmt ----
  
  tdd_sbm_ari = 1:N_sim
  tdd_sbm_mape = 1:N_sim
  tdd_sbm_true = 1:N_sim
  tdd_sbm_sim = 1:N_sim
  
  for (s in 1:N_sim) {
    
    # generate simulated networks
    
    if (dc == 0) { # no degree correction case  ----
      discrete_edge_array = generate_multilayer_array(roles, omega, type = "discrete")
    }
    if (dc == 3) { # degree correction case  ----
      discrete_edge_array = generate_multilayer_array(roles, block_omega, dc_factors, type = "discrete")
    }
    discrete_edge_list = adj_to_edgelist(discrete_edge_array, directed = TRUE, selfEdges = TRUE)
  
    # fit simulated networks
    tdd_sbm = sbmt(discrete_edge_list, maxComms = K, degreeCorrect = dc, directed = TRUE, klPerNetwork = kl_per_network)
    
    if (verbose) {
      #plot
      plot(tdd_sbm)
      #points(1:Time, block_omega[K,K,], col = "red", add = T, type = "l")
    }
    
    # adjusted rand index
    tdd_sbm_ari[s] = adj.rand.index(tdd_sbm$FoundComms[order(as.numeric(names(tdd_sbm$FoundComms)))], roles)
    
    # MAPE of omega
    tdd_sbm_mape[s] = mean(abs(sapply(tdd_sbm$EdgeMatrix, matrix) - apply(block_omega, 3, matrix))/apply(block_omega, 3, matrix))
    
    # compare likelihood for true vs. fit parameters for simulated data
    tdd_sbm_true[s] = 
      tdd_sbm_llik(discrete_edge_array, roles = roles - 1, omega = block_omega, directed = TRUE, selfEdges = TRUE)
    tdd_sbm_sim[s] = 
      tdd_sbm_llik(discrete_edge_array, roles = tdd_sbm$FoundComms, omega = tdd_sbm$EdgeMatrix, directed = TRUE, selfEdges = TRUE)
  }
  
  # - results ----
  results = list(
    # role detection
    tdd_sbm_ari = tdd_sbm_ari,
    # blcok-to-block activity detection
    tdd_sbm_mape = tdd_sbm_mape,
    # compare likelihood of data under fit vs. true model
    tdd_sbm_true = tdd_sbm_true,
    tdd_sbm_sim = tdd_sbm_sim,
    tdd_true_vs_sim = tdd_sbm_true - tdd_sbm_sim
  )
  
  return(results)
}

  # example:
  # no degree hetereogeneity in data generation 
  roles = generate_roles(30, roles = 2, type = "discrete", rel_freq = c(.5,.5))
  results = simulated_tdd(roles, omega = omega_2)
  # degree heterogeneity
  dc_factors = generate_dc_factors(roles, dc_levels = c(1,2,3))
  dc_results = simulated_tdd(roles, omega = omega_2, dc_factors)


# run tdd-sbm simulation ----

td_results = vector(mode = "list", length = length(K_set))

td_results = 
  lapply(1:length(K_set), function(i) {
  K = K_set[i]
  omega = omega_list[[i]]
  for(N in N_set) {
    for (dc in c(0,3)) { #no degree correction, undirected, time-independent degree correction
      roles = generate_roles(30, roles = 2, type = "discrete", rel_freq = rep(1/K, K))
      lapply(V, function(N) {
        tdd_results_dc0 = simulated_tdd(N, n_roles, omega, Time, dc = 0, N_sim, kl = kl_per_network) 
        tdd_results_dc3 = simulated_tdd(N, n_roles, omega, Time, dc = 3, N_sim, kl = kl_per_network) 
      })
    }
  })


# ---------------------------------------------------------------------------------------------------------------
# 2. ppsbm to show impact of degree correction ----
library(ppsbm)

# no degree correction case. their model works as expected ----
# - fit ppsbm ----

# Use the "hist" method because agrees more closely with out discrete Time slices and requires little data manipulation

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
adj.rand.index(apply(discrete_ppsbm[[selected_Q]]$tau, 2, which.max), roles)

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
adj.rand.index(apply(dc_discrete_ppsbm[[selected_Q]]$tau, 2, which.max), roles)

# omegas
par(mfrow = c(selected_Q, selected_Q)); par(mai = rep(.5, 4))
apply(exp(dc_discrete_ppsbm[[selected_Q]]$logintensities.ql), 1, plot, type = "l", col = "blue")

# with true number of groups? gets it right
apply(dc_discrete_ppsbm[[n_roles]]$tau, 2, which.max)
adj.rand.index(apply(dc_discrete_ppsbm[[n_roles]]$tau, 2, which.max), roles)

par(mfrow = c(n_roles, n_roles)); par(mai = rep(.5, 4))
apply(exp(dc_discrete_ppsbm[[n_roles]]$logintensities.ql), 1, plot, type = "l", col = "blue")

# Bike example shows how degree correction ib model -> group statins with similar behavior across activity levels
# ---------------------------------------------------------------------------------------------------------------
# 3. mixed-membership to show impact of potential model mis-specification ----

n_roles = K_set[1]
N = V[1]

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
  
  mixed_edge_array = generate_multilayer_array(roles_mixed, block_omega, type = "mixed")
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
  tdd_sbm = sbmt(mixed_edge_edglist, maxComms = n_roles, degreeCorrect = 3, directed = TRUE, klPerNetwork = kl_per_network)
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


