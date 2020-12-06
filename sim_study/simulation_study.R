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
library(ppsbm)
library(fossil) #for adj rand index
library(gtools) #for permutations

# generate_multilayer_array function ----

#' generate N x N x Time array (Time-sliced adjacency matrices) based on TDD-SBM (type == "discrete) or TDMM-SBM (type = "mixed") model
#' @param roles If type is `discrete`, a length-N vector of block assignments. If type is `mixed` G x N matrix of block-assignment weights.
#' @param omega A K x K x Time array or Time-length list of K x K matrices. If the latter will be internally converted to array.
#' If `type` is `discrete`, omega represents mean block-to-block edge values over time. 
#' If `type` is `mixed` omega represents total block-to-block activity over time.
#' @param type `discrete` or `mixed` for tdd-sbm or tdmm-sbm, respectively. If mixed, dc_factors is ignored.
#' @param dc_factors 
# '   if not degree-corrected, should be NULL
# '   otherwise, dc_factors are a list of degree-correction factors that should be normalized to sum to 1 within blocks
# for type "mixed" roles are a G x N matrix of assignment weights 
# (add this as a simulation function to sbmt package?)
#' @examples 
#' generate_multilayer_array(roles = rep(1:2, 5), omega = array(1:12, dim = c(2,2,3)), dc_factors = NULL, type = "discrete") #output 10 x 10 x 3
#' generate_multilayer_array(roles = rep(1:2, 5), omega = array(1:12, dim = c(2,2,3)), dc_factors = 1:10, type = "discrete") #output 10 x 10 x 3
generate_multilayer_array <- function(roles, omega, dc_factors = NULL, type = "discrete") {
  
  #  degree-corrected?
  dc = ifelse(is.null(dc_factors), FALSE, TRUE)
  
  # checks
  if (length(unique(roles)) <= 1) stop("There should be at least two unique roles")
  
  if (type == "discrete") {
    roles = as.numeric(as.factor(roles)) #in this case want to align role and index #assumes order of omega rows/columns corresponds to role order 
    if (dc) {
      #normalize dc
      if (! identical(aggregate(dc_factors ~ roles, FUN = "sum")[,2], rep(1, length(unique(roles))))) {
        warnings("degree correction factors not normalized to sum to 1 by block. will normalize now.")
        role_sums = aggregate(dc_factors, by = list(roles), sum)
        role_sums = setNames(role_sums$x, role_sums$Group.1)
        dc_factors = dc_factors/role_sums[roles] #group sum constraint
        cat("dc_factors normalized:", dc_factors, "\n")
        # check all(aggregate(dc_factors, by = list(roles), sum)$x == 1)
      }
    }
  }
  
  if (type == "mixed" && !identical(colSums(roles), rep(1, ncol(roles)))) { # normalize mixed
    if (min(colSums(roles))<=0) {warnings("At least one roles non-positive weight.")}
    roles = sweep(roles, 2, colSums(roles), FUN="/")
  }
  
  if ("data.frame" %in% class(roles)) {roles = as.matrix(roles)} #list structure of data frame would interfere later
  
  #  infer N, K, Time,
  N = ifelse(is.null(dim(roles)), length(roles), nrow(roles))
  K = ifelse("list" %in% class(omega), nrow(omega[[1]]), dim(omega)[1]) # number of blocks
  Time = ifelse("list" %in% class(omega), length(omega), dim(omega)[3])
  if (!"array" %in% class(omega)) {omega  = array(unlist(omega), dim = c(K,K,Time))} #make array if not already
  
  # adjusted omega (to account for normalization in degree correction) 
  # entries are expected total weight from r to s instead of expected edge weight from single edge from block r to block s
  if (dc) {block_omega = omega*array(table(roles) %*% t(table(roles)), dim = c(K, K, Time))}

  edge_array = array(0, dim = c(N, N, Time))
  for (i in 1:N) {
    role_i = roles[i]
    for (j in 1:N) {
      role_j = roles[j]
      for (tm in 1:Time) {
          if (type == "discrete" & !dc) { ijt = rpois(omega[role_i, role_j, tm], n= 1) }
          if (type == "discrete" & dc) { ijt = rpois(dc_factors[i]*dc_factors[j]*block_omega[role_i, role_j, tm], n= 1) }
          if (type == "mixed") { ijt = rpois(t(roles[i, ]) %*% omega[, , tm] %*% roles[j, ], n= 1) }
          edge_array[i, j, tm] = ijt
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

#   - omegas ----

a = 10
b = 5 
ymax =  max(2*a, 2*b)
x = seq(.5, Time)

# show two-block curves
par(mfrow = c(2, 2))
curve(b*sin(x*pi/Time) + b, 0, Time, ylim = c(0, ymax))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(b*sin(x*pi/Time)^3+ b, 0, Time, ylim = c(0, ymax))

# show three-block curves
T = Time #set just for labeling, then remove
png("IMG/sim_omega_3.png", width = 600)
par(mai =c(.4,.7,.3,.2))
par(mfrow = c(3, 3))
curve(b*sin(x*pi/T) + b, 0, T, ylim = c(0, ymax), main = ("1 to 1"), xlab = "")
curve(a*sin(x*2*pi/T)+a, 0, T, ylim = c(0, ymax), main = ("1 to 2"), xlab = "")
curve(a*sin(x*2*pi/T-1)/2+a, 0, T, ylim = c(0, ymax), main = ("1 to 3"), xlab = "")
curve(-a*sin(x*2*pi/T)+a, 0, T, ylim = c(0, ymax), main = ("2 to 1"), xlab = "")
curve(b*sin(x*pi/T)^3+ b, 0, T, ylim = c(0, ymax), main = ("2 to 2"), xlab = "")
curve(a*sin(x*2*pi/T-2)/4+a, 0, T, ylim = c(0, ymax), main = ("2 to 3"), xlab = "")
curve(-a*sin(x*2*pi/T-1)/2+a, 0, T, ylim = c(0, ymax), main = ("3 to 1"), xlab = "")
curve(-a*sin(x*2*pi/T-2)/4+a, 0, T, ylim = c(0, ymax), main = ("3 to 2"), xlab = "")
curve(0*x+b, 0, T, ylim = c(0, ymax), main = ("3 to 3"), xlab = "")
dev.off()
rm(T)

#show overlapping
par(mfrow = c(1,1))
curve(a*sin(x*2*pi/Time)+a, 0, Time, ylim = c(0, ymax))
curve(-a*sin(x*2*pi/Time)+a, 0, Time, add = TRUE, lty = 2)
curve(a*sin(x*2*pi/Time-1)/2+a, 0, Time, col = "blue", add = TRUE)
curve(-a*sin(x*2*pi/Time-1)/2+a, 0, Time, col = "blue", add = TRUE, lty = 2)
curve(a*sin(x*2*pi/Time-2)/4+a, 0, Time, col = "green", add = TRUE)
curve(-a*sin(x*2*pi/Time-2)/4+a, 0, Time, col = "green", add = TRUE, lty = 2)
curve(b*sin(x*pi/Time) + b, 0, Time, col = "red", add = TRUE)
curve(b*sin(x*pi/Time)^3+ b, 0, Time, col = "red", add = TRUE, lty = 2)
curve(0*x+b, 0, Time, col = "red", add = TRUE, lty = 3)

omega_11 = b*sin(x*pi/Time) + b
omega_12 =  a*sin(x*2*pi/Time)+a
omega_21 = -a*sin(x*2*pi/Time)+a
omega_22 = b*sin(x*pi/Time)^3+ b
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
omega_list = list(`2` = omega_2, `3` = omega_3)

# 1. tdd simulation ----

#   tdd-sbm simulation helper functions ----

# generate discrete and mixed roles for use in simulation.
# currently assumes directed = TRUE and selfEdges = TRUE
# (custom roles can also be fed to the simulation)
# role_types: (discrete case) is an integer, indicating the possible blocks (roles) a user can have
#        (mixed case) a matrix where each row is a C_i vector of block weights, indicating realized mixtures of blocks. (A small set for the purposes of simulation. Could also sample from a distribution with k-vector range)
# rel_freq: dictates relative frequency of the blocks (discrete case) or block-weight vectors (mixed cased), defaulting to as close as possible to equal groups. Will be normalized (to sum 1) if not already.

generate_roles <- function(N, role_types, type = c("discrete", "mixed")[1], rel_freq = 1) {
  
  if (!type %in% c("discrete", "mixed")) {stop("type must be discrete or mixed")}
  
  # if rel_freq == 1 (default for equal frequency) set the right rel_freq vector depending on "type"
  if(identical(rel_freq, 1)) {
    if (type=="discrete") {rel_freq = rep(1, role_types)}
    if (type=="mixed") {rel_freq = rep(1, nrow(role_types))}
  }
  rel_freq = rel_freq/sum(rel_freq)
  if (type == "discrete") {
    roles = rep(1:role_types, times = ceiling(N*rel_freq))[1:N] #only up to first N elements in case of rounding
  } else { #type mixed
    times = ceiling(rel_freq*N)
    roles = role_types[rep(1:length(times),times),][1:N,] #only up to first N elements in case of rounding
    roles = sweep(roles, 2, colSums(roles), FUN="/")
  }
  return(roles)
}

  # examples:
    # generate_roles(30, role_types = 2, type = "discrete", rel_freq = c(.5,.5))
    # generate_roles(30, role_types = 3, type = "discrete", rel_freq = c(1,2,3))
    # mixed_role_options = matrix(c(0,.25,.75)[permutations(3, 3)],factorial(3),3)
    # generate_roles(30, role_types = mixed_role_options, type = "mixed", rel_freq = rep(1,nrow(mixed_role_options)))
    # generate_roles(30, role_types = mixed_role_options, type = "mixed", rel_freq = 1:nrow(mixed_role_options))

# roles are the discrete roles of the nodes (which also tells us N)
# dc_levels is preset relative activity level we use for degree correction. They will be cycled through within each role.
generate_dc_factors <- function(roles, dc_levels) {
  
  dc_factors = lapply(table(roles), rep_len, x=dc_levels)
  dc_factors = lapply(dc_factors, function(x) {x/sum(x)}) #group sum constraint
  dc_factors = as.vector(do.call("c", dc_factors))
  return(dc_factors)
}

  # examples
  roles = generate_roles(30, role_types = 2, type = "discrete", rel_freq = c(.5,.5))
  generate_dc_factors(roles, dc_levels = c(1,2,3))

#   tdd-sbm simulation function ----
  
# A function to generate and fit data from the tdd-sbm
# currently for directed networks only
# roles: N-vector with element in 1:K (K is number of blocks) used for simulating data.
# omega: K x K x T array used for simulating data.
# dc_factors: degree correction factors for simulating data. If "tdd-sbm-0" this is not used.
# sim_method: method used to simulate data. Last numeric of tdd-sbm- is 0 for no degree correction, 3 for time independent undirected degree correction.
# fit_method: method used to fit the simulated data
# verbose: TRUE will show model fit progress and some intermediate plots
simulate_tdd <- function(roles, omega,  dc_factors = NULL, 
                         sim_method = c("tdd-sbm-0", "tdd-sbm-3")[1], 
                         fit_method = c("tdd-sbm-0", "tdd-sbm-3", "ppsbm")[1], 
                         N_sim = 10, directed = TRUE, kl = kl_per_network, verbose = FALSE) {
  
  #checks 
  if (!directed) stop("Not yet implemented for undirected networks")
  if (sim_method=="tdd-sbm-3" & is.null(dc_factors)) stop("Sim method tdd-sbm-3 requires dc_factors")
  #make omega an array if not already
  if(!"array" %in% class(omega)) {omega  = array(unlist(omega), dim = c(nrow(omega[[1]]),nrow(omega[[1]]),length(omega)))} 
  
  #  infer  K, N, Time, degree correction parameter
  K = dim(omega)[2] # number of bloks
  N = length(roles)
  Time = dim(omega)[3]
  dc_sim = as.numeric(gsub(sim_method, pattern = "tdd-sbm-", replacement=""))
  dc_fit = as.numeric(gsub(fit_method, pattern = "tdd-sbm-|ppsbm", replacement=""))
  block_omega = omega*array(table(roles) %*% t(table(roles)), dim = c(K, K, Time))
  
  # ---------------------------------------------------------------------------------------------------------------
  # - fit sbmt ----
  
  tdd_sbm_ari = 1:N_sim
  tdd_sbm_mape = 1:N_sim
  tdd_sbm_sim = 1:N_sim
  tdd_sbm_fit = 1:N_sim
  tdd_sbm_est_K = rep(NA, N_sim) #only filled in for ppsbm
    
  for (s in 1:N_sim) {
    
    # generate simulated networks
    if (sim_method == "tdd-sbm-0") { # no degree correction (when generating data) case  ----
      discrete_edge_array = generate_multilayer_array(roles, omega, type = "discrete")
    }
    if (sim_method == "tdd-sbm-3") { # degree correction case  ----
      discrete_edge_array = generate_multilayer_array(roles, omega, dc_factors, type = "discrete")
    }
    discrete_edge_list = adj_to_edgelist(discrete_edge_array, directed = TRUE, selfEdges = TRUE)
  
    # fit simulated networks
    if (grepl(fit_method, pattern = "tdd-sbm")) {
      
      # print each fit's progress?
      if (verbose) {
        tdd_sbm = sbmt(discrete_edge_list, maxComms = K, degreeCorrect = dc_fit, directed = TRUE, klPerNetwork = kl_per_network)
      } else {
        log = capture.output({tdd_sbm = sbmt(discrete_edge_list, maxComms = K, degreeCorrect = dc_fit, directed = TRUE, klPerNetwork = kl_per_network)})
      }
    
      # optional plotting
      if (verbose) {
        #plot
        plot(tdd_sbm, show_reverse = FALSE)
        #points(1:Time, block_omega[K,K,], col = "red", add = T, type = "l")
      }

      # adjusted rand index
      tdd_sbm_ari[s] = adj.rand.index(tdd_sbm$FoundComms[order(as.numeric(names(tdd_sbm$FoundComms)))], roles)
    
      # MAPE of omega
      tdd_sbm_mape[s] = mean(abs(sapply(tdd_sbm$EdgeMatrix, matrix) - apply(block_omega, 3, matrix))/apply(block_omega, 3, matrix))
    
      # compare likelihood for true vs. fit parameters for simulated data
      tdd_sbm_sim[s] = tdd_sbm_llik(discrete_edge_array, roles = roles - 1, omega = block_omega, degreeCorrect = dc_sim, directed = TRUE, selfEdges = TRUE)
      tdd_sbm_fit[s] = tdd_sbm_llik(discrete_edge_array, roles = tdd_sbm$FoundComms, omega = tdd_sbm$EdgeMatrix, degreeCorrect = dc_fit, directed = TRUE, selfEdges = TRUE)
    }
    if (fit_method == "ppsbm") {
      # - fit ppsbm ----
      
      # Use the "hist" method because agrees more closely with out discrete Time slices and requires less manipulation
      #reformat data for ppsbm fit
      Nijk = sapply(adj_to_edgelist(discrete_edge_array, directed = TRUE, selfEdges = FALSE, removeZeros = FALSE), "[[", 3)
      #dim(Nijk) #should have ncol == Time
      
      #ppsbm_fit
      # number of blocks to try -- for this study look between 1 and 4
      min_K = 1 # -- can index later by K (true) and selected_K as long as min_K is 1
      max_K = 4
      
      # print each fit's progress?
      if (verbose) {
        ppsbm = mainVEM(list(Nijk=Nijk,Time=Time), N, Qmin = min_K, Qmax = max_K, directed=TRUE,
                        method='hist', d_part=4, n_perturb=10, n_random=0)
      } else {
        log = capture.output({ppsbm = mainVEM(list(Nijk=Nijk,Time=Time), N, Qmin = min_K, Qmax = max_K, directed=TRUE,
                                              method='hist', d_part=4, n_perturb=10, n_random=0)})
      }
      
      if (verbose) {
        selected_K = modelSelection_Q(list(Nijk=Nijk,Time=Time), N, Qmin = min_K, Qmax = max_K, directed = TRUE, sparse = FALSE, ppsbm)$Qbest
      } else {
        log = capture.output({selected_K = modelSelection_Q(list(Nijk=Nijk,Time=Time), N, Qmin = min_K, Qmax = max_K, directed = TRUE, sparse = FALSE, ppsbm)$Qbest})
      }
      # selected K
      # ppsbm_roles = apply(ppsbm[[selected_K]]$tau, 2, which.max) #est roles
      # par(mfrow = c(selected_K, selected_K));  apply(exp(ppsbm[[selected_K]]$logintensities.ql), 1, plot, type = "l")
      tdd_sbm_est_K[s] = selected_K
      
      # true K
      # role detection
      ppsbm_roles_K = apply(ppsbm[[K]]$tau, 2, which.max)
      tdd_sbm_ari[s] = adj.rand.index(apply(ppsbm[[K]]$tau, 2, which.max), roles)
      # omegas 
      ppsbm_omega = exp(array(ppsbm[[K]]$"logintensities.ql", dim = c(K, K, Time)))
      ppsbm_block_omega = ppsbm_omega*array(table(ppsbm_roles_K) %*% t(table(ppsbm_roles_K)), dim = c(K, K, Time))

      tdd_sbm_mape[s] = mean(abs(ppsbm_block_omega - block_omega)/(block_omega))
      # use selfEdges FALSE since ppsbm fits without them? (but then adjust tdd_sbm block omega?)
      tdd_sbm_sim[s] = tdd_sbm_llik(discrete_edge_array, roles = roles - 1, omega = block_omega, degreeCorrect = dc_sim, directed = TRUE, selfEdges = TRUE) #roles should be 0-indexed
      tdd_sbm_fit[s] = tdd_sbm_llik(discrete_edge_array, roles = ppsbm_roles_K - 1, omega = ppsbm_block_omega, degreeCorrect = 0, directed = TRUE, selfEdges = TRUE)
    }
  
  }
  
  # - results ----
  results = list(
    # role detection
    tdd_sbm_ari = tdd_sbm_ari,
    # blcok-to-block activity detection
    tdd_sbm_mape = tdd_sbm_mape,
    # compare likelihood of data under fit vs. true model used to simulate the data
    tdd_sbm_sim = tdd_sbm_sim,
    tdd_sbm_fit = tdd_sbm_fit,
    tdd_sim_vs_fit_method = tdd_sbm_fit - tdd_sbm_sim,
    tdd_sbm_est_K = tdd_sbm_est_K
  )
  
  return(results)
}

#   - examples: ----
  # no degree hetereogeneity in data generation 
  # roles = generate_roles(30, role_types = 2, type = "discrete", rel_freq = c(.5,.5))
  # results00 = simulate_tdd(roles, omega = omega_2, sim_method = "tdd-sbm-0", fit_method = "tdd-sbm-0")
  # results03 = simulate_tdd(roles, omega = omega_2, sim_method = "tdd-sbm-0", fit_method = "tdd-sbm-3")
  # # degree heterogeneity
  # dc_fctrs = generate_dc_factors(roles, dc_levels = c(1:5))
  # results30 = simulate_tdd(roles, omega = omega_2, dc_factors = dc_fctrs, sim_method = "tdd-sbm-3", fit_method = "tdd-sbm-0")
  # results33 = simulate_tdd(roles, omega = omega_2, dc_factors = dc_fctrs, sim_method = "tdd-sbm-3", fit_method = "tdd-sbm-3")


#   run tdd-sbm simulation ----

tdd_results_30 = data.frame(expand.grid(K = K_set, N = N_set[1], sim_method = c("tdd-sbm-0", "tdd-sbm-3"), 
                                        fit_method = c("tdd-sbm-0", "tdd-sbm-3", "ppsbm")))
tdd_results_90 = data.frame(expand.grid(K = K_set, N = N_set[2], sim_method = c("tdd-sbm-0", "tdd-sbm-3"), 
                                        fit_method = c("tdd-sbm-0", "tdd-sbm-3", "ppsbm")))

tdd_results = apply(tdd_results_30, 1, function(x) {
  cat("running for parameters:", x,"\n")
  K = as.numeric(x['K'])
  roles = generate_roles(N = as.numeric(x['N']), role_types = K, type = "discrete", rel_freq = rep(1, K)) #<- note: can alter relative role frequency by changing this specification. E.g. 1:K would give increasing frequency to higher-numbered roles.
  omega = omega_list[as.character(K)][[1]]
  dc_fctrs = generate_dc_factors(roles, dc_levels = c(1:5)) #<- note: can alter degree correction levels by changing this specification. E.g. 1:2 would give only two levels of degree heterogeneity
  results = simulate_tdd(roles, omega, dc_fctrs, sim_method = x[['sim_method']], fit_method = x[['fit_method']],
                          N_sim = 10, verbose = FALSE)
   return(results)
})

#saveRDS(tdd_results, "sim_study/tdd_results.RDS")
tdd_results = readRDS("sim_study/tdd_results.RDS")

# Add results to summary table

tdd_results_30$sim_method = gsub(tdd_results_30$sim_method, pattern = "-0", replacement = "")
tdd_results_30$sim_method = gsub(tdd_results_30$sim_method, pattern = "-3", replacement = " DC")
tdd_results_30$fit_method = gsub(tdd_results_30$fit_method, pattern = "-0", replacement = "")
tdd_results_30$fit_method = gsub(tdd_results_30$fit_method, pattern = "-3", replacement = " DC")

tdd_results_mean = lapply(tdd_results, function(x) {sapply(x, mean)})
tdd_results_sd = lapply(tdd_results, function(x) {sapply(x, sd)})

tdd_results_30$ARI = paste0(round(sapply(tdd_results_mean, "[[", "tdd_sbm_ari"), 2), " (", round(sapply(tdd_results_sd, "[[", "tdd_sbm_ari"), 2), ")")
tdd_results_30$MAPE = paste0(round(sapply(tdd_results_mean, "[[", "tdd_sbm_mape"), 2), " (", round(sapply(tdd_results_sd, "[[", "tdd_sbm_mape"), 2), ")")
tdd_results_30$LLIK_sim = paste0(round(sapply(tdd_results_mean, "[[", "tdd_sbm_sim")), " (", round(sapply(tdd_results_sd, "[[", "tdd_sbm_sim"), 1), ")")
tdd_results_30$LLIK_diff = paste0(round(sapply(tdd_results_mean, "[[", "tdd_sim_vs_fit_method")), " (", round(sapply(tdd_results_sd, "[[", "tdd_sim_vs_fit_method"), 1), ")")
tdd_results_30$`Est k (ppsbm)` = paste0(round(sapply(tdd_results_mean, "[[", "tdd_sbm_est_K")))

tdd_results_30 = tdd_results_30[,!names(tdd_results_30) %in% c("sim_method","fit_method")]


# ppsbm, no degree correction case. their model works as expected ----
# Wants a seperate class for each degree-correcton level (but note the LLIK_sim and LLIK_diff results use the true K)
# Bike example (seperate script) shows how degree correction ib model -> group stations with similar behavior across activity levels

# ---------------------------------------------------------------------------------------------------------------
# 3. mixed-membership to show impact of potential model mis-specification ----

#mixed_role_options_2 = matrix(c(1,0)[permutations(2, 2)],factorial(2),2) #with data on the boundary it works
mixed_role_options_2 = matrix(c(.25,.75, 1, 0)[permutations(4, 2)], nrow(permutations(4,2)), 2) #with data on the boundary it works
mixed_role_options_2 = matrix(c(.25,.75)[permutations(2, 2)], nrow(permutations(2,2)), 2) #with data on the boundary it works
mixed_role_options_3 = matrix(c(0,.25,.75)[permutations(3, 3)], nrow(permutations(3,2)), 3)
#mixed_role_options_2 = matrix(rep(1:5, 2), ncol = 2) #identifiability issues if all evenly split
#mixed_role_options_3 = matrix(rep(1:5, 3), ncol = 3) 
mixed_role_options_list = list(mixed_role_options_2, mixed_role_options_3)


i = 2
N = N_set[i]
K = K_set[i]
omega = omega_list[[i]]
block_omega = omega*N^2/K^2 #these weights should align it with degree corrected example in this case (equally frequent roles)
mixed_role_options = mixed_role_options_list[[i]]

tdmm_sbm_role_err = 1:N_sim #sum of abs error averaged over number of blocks (since columns are sum-1 normalized)
#tdmm_sbm_omega_mean_abs_error = 1:N_sim
#tdmm_sbm_ari = 1:N_sim #divide by number of roles?
tdmm_sbm_mape = 1:N_sim
tdmm_sbm_sim =1:N_sim
tdmm_sbm_fit =1:N_sim
tdmm_discrete_fit =1:N_sim

setwd("mixed_model_implementation_python") # assume starting from tdsbm_supplementary_material directory

for (s in 1:N_sim) {
  
  #roles_mixed = matrix(sample(1:(N*K), replace = TRUE), ncol = K); roles_mixed = sweep(roles_mixed, 2, colSums(roles_mixed), FUN="/")
  #image(roles_mixed %*% block_omega[, , 1] %*% t(roles_mixed))
  roles_mixed = generate_roles(N, role_types = mixed_role_options, type = "mixed", rel_freq = rep(1,nrow(mixed_role_options))) #can also try 1:nrow(mixed_role_options)
  roles_mixed = la
  block_omega = la_omega
  mixed_edge_array = generate_multilayer_array(roles_mixed, block_omega, type = "mixed")

  #without randomness: mixed_edge_array = array(unlist(lapply(1:Time, function(i) {roles_mixed %*% block_omega[,,i] %*% t(roles_mixed)})), dim = c(N, N, Time))
  write.csv(mixed_edge_array, paste0("../sim_study/mixed_edge_array.csv"), row.names = FALSE)
  params = data.frame(K = K)
  write.csv(params, paste0("../sim_study/mixed_params.csv"), row.names = FALSE)
  # run mixed
  system("python3 tdmm_sbm_sim_study.py")
  tdmm_sbm_roles = read.csv("../sim_study/SIM_roles.csv")
  
  block_order = permutations(K, K)[which.min(apply(permutations(K, K), 1, function(x) { sum(abs(roles_mixed[,x] - tdmm_sbm_roles)) })),]
  
  tdmm_sbm_role_err[s] = sum(abs(roles_mixed[,block_order] - tdmm_sbm_roles))/K #try all orders
  block_omega_ordered = array(sapply(1:Time, function(i) {block_omega[block_order,block_order,i]}), dim = c(K, K, Time))
  tdmm_sbm_omega = read.csv("../sim_study/SIM_omega.csv")
  tdmm_sbm_omega_array = array(unlist(tdmm_sbm_omega), dim = c(K,K,Time))
  tdmm_sbm_mape[s] = mean(abs(tdmm_sbm_omega_array - block_omega_ordered)/block_omega_ordered+1)

  #compare data matrix reconstruction
  # mean(mixed_edge_array[,,16] - as.matrix(tdmm_sbm_roles) %*% tdmm_sbm_omega_array[,,16] %*% t(as.matrix(tdmm_sbm_roles)))
  # mean(mixed_edge_array[,,16] - as.matrix(roles_mixed) %*% block_omega[,,16] %*% as.matrix(t(roles_mixed)))
  
  if (verbose) {
    par(mfrow = c(K+1,K)); par(mai = rep(.6,4))
    # mixed role comparison
    axes = 1:2
    plot((roles_mixed[,block_order])[,axes], xlim = c(0,max(cbind(roles_mixed, tdmm_sbm_roles))), ylim = c(0,max(cbind(roles_mixed, tdmm_sbm_roles))))
    points(tdmm_sbm_roles[,axes], col = "red") 
    # overall activity pattern comparison
    plot(colSums(apply(block_omega, 3, unlist)), type = "l"); points(colSums(tdmm_sbm_omega), col = "red", type = "l")
    for (i in 1:K^2) { #omega comparison
      plot(as.numeric(tdmm_sbm_omega[i,]), type = "l", ylim = c(0,max(c(max(tdmm_sbm_omega), max(block_omega_ordered)))))
      points(apply(block_omega_ordered, 3, dplyr::nth, i), col = "red", type = "l")
    }
  }
  
  # vs self  (mixed)
  tdmm_sbm_sim[s] = tdmm_sbm_llik(A = mixed_edge_array, C = roles_mixed, omega = block_omega, selfEdges = TRUE, directed = TRUE)
  # vs fit
  tdmm_sbm_fit[s] = tdmm_sbm_llik(A = mixed_edge_array, C = tdmm_sbm_roles, omega = tdmm_sbm_omega, selfEdges = TRUE, directed = TRUE)
  
  # try fitting discrete models to data generated from mixed membership
  mixed_edge_edglist = adj_to_edgelist(mixed_edge_array, directed = TRUE, selfEdges = TRUE, removeZeros = TRUE)
  tdd_sbm = sbmt(mixed_edge_edglist, maxComms = K, degreeCorrect = 3, directed = TRUE, klPerNetwork = kl_per_network)
  tdd_roles = matrix(0, N, K); for (i in 1:N) {tdd_roles[i, tdd_sbm$FoundComms[i]+1] = 1}; rownames(tdd_roles) = 1:N
  # plot(tdd_sbm)
  # plot(colSums(tdmm_sbm_omega), type = "l)
  # evaluate likelihood of mixed edge array using params found by discrete model
  tdd_roles = sweep(tdd_roles, 2, colSums(tdd_roles), FUN="/") # column normalize
  tmp = sapply(tdd_sbm$EdgeMatrix, as.vector)
  tdmm_discrete_fit[s] = tdmm_sbm_llik(mixed_edge_array, C = tdd_roles, omega = tdd_sbm$EdgeMatrix, selfEdges = TRUE, directed = TRUE)
}

tdmm_vs_tdd_sbm_ll = tdmm_sbm_ll - tdd_sbm_ll #(scale by difference in parameters?)

setwd("..")

#potential for identifiability issues when omega 1->1 and 2->2 curves were identical

