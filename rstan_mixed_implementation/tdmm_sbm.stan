data {
  int <lower=2> N; //number of stations
  int <lower=2> K; //number of blocks
  int <lower=1> Time; //number of time periods
  int A[N,N,Time]; //A three-dimensional array of integers = station-to-station trips by time slice (N x N x T)
  int <lower=0, upper=1> selfEdges; 
}

parameters {
  matrix <lower=0>[K,K] omega[Time]; //array of KxK matrices of block to block parameters over time
  simplex[N] C[K]; //mixed membership parameters for each node, do not change over time. columns sum to 1. 
}

transformed parameters {
  matrix[N,N] mu[Time]; //expected number of trips from i to j
  matrix[N,K] C2;
  for (k in 1:K) {
    C2[,k] = C[k];
  }
  for (t in 1:Time) {
    mu[t] = C2*omega[t]*C2';
  }
}

// deal with diagonals ***
model {
  // real prior on omega? log normal?
  /*for (t in 1:Time) {
    for (k in 1:K) {
      omega[t][K,] ~ lognormal(0, 10); 
    }
  }*/
  for (t in 1:Time) {
    for (i in 1:N) {
      A[,i,t] ~ poisson(mu[t,,i]);   
      if (selfEdges == 0) {
        target += -poisson_lpmf(A[i,i,t] | mu[t,i,i]);
      } 
    }
  }
}


