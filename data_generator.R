######################################################################
# Code used to generate SW-CRT data
######################################################################

## Nested exchangeable + Random intervention model
gendata <- function(
    nclus, #number of cluster
    nperiod, #number of period
    wpicc, #wpicc in control
    wpicc_int, #wpicc in int
    cac, #cac
    sigma, #individual error
    nsubject, #number of indiviudals per period
    theta #standized treatment effect
){
  nseq = nperiod - 1
  nclus_seq = nclus/nseq
  X <- matrix(data = 0, nrow = nseq, ncol = nperiod)
  for(i in 1:nseq){
    X[i,(i+1):nperiod] <- 1
  }
  X <- X %x% rep(1, nclus_seq)
  
  
  
  # implied between-cluster variance
  alpha_sq <- cac * sigma^2 * wpicc / (1 - wpicc)
  
  # implied between-cluster-period variance
  gamma_sq <- sigma^2 * wpicc / (1 - wpicc) - alpha_sq
  
  # implied between-period intracluster correlation
  #BPICC <- wpicc * cac
  
  
  eta_sq <- (wpicc_int*(alpha_sq + gamma_sq + sigma^2) - alpha_sq - gamma_sq)/(1-wpicc_int)
  
  C <- rnorm(n = nclus, mean = 0, sd = sqrt(alpha_sq))
  CP <- rnorm(n = nclus * nperiod, mean = 0, sd = sqrt(gamma_sq))
  RT <- rnorm(n = nclus, mean = 0, sd = sqrt(eta_sq))
  
  Yijk <- c()
  e <- rnorm(n = nclus * nperiod * nsubject, mean = 0, sd = sigma)
  for (i in 1:nclus) {
    for (j in 1:nperiod) {
      q <- (i - 1) * nperiod + j
      for (k in 1:nsubject) {
        # l is the index into the data
        l <- (i - 1) * nperiod * nsubject + (j - 1) * nsubject + k
        Yijk[l] <- j + X[i,j] * (theta + RT[i]) + C[i] + CP[q] + e[l]
      }
    }
  }
  
  Xvec <- as.vector(kronecker(as.vector(t(X)), matrix(rep(1, nsubject), ncol=1)))
  Cvec <- as.factor(rep(1:nclus, each=nperiod*nsubject))
  CPvec <- as.factor(rep(1:(nclus*nperiod), each=nsubject))
  Pvec <- as.factor(rep(1:nperiod, each=nsubject, times=nclus))
  
  dat <- data.frame(
    response.var = Yijk,
    tx.var = Xvec,
    cluster.var = Cvec,
    time.var = Pvec,
    cp.var = CPvec
  )
  return(dat)
}
### Example
#data <- gendata(nclus = 32, #number of cluster
#                nperiod = 9, #number of period
#                wpicc = 0.05, #cluster-period var
#                wpicc_int = 0.1, #cluster-level
#                cac = 0.8,
#                sigma = 1, #individual error
#                nsubject = 100, #number of individuals per period
#                theta = 0 #standardized treatment effect
#                )





## Discrete-time Decay + Random intervention model 
library(MASS)
gendata_dtd <- function(
    nclus, #number of cluster
    nperiod, #number of period
    cac, #cluster autocorrelation (decay per period)
    wpicc, #within-period ICC
    wpicc_int,
    sigma, #individual error
    nsubject, #number of indiviudals per period
    theta #standized treatment effect
){
  nseq = nperiod - 1
  nclus_seq = nclus/nseq
  X <- matrix(data = 0, nrow=nseq, ncol=nperiod)
  for(i in 1:nseq){
    X[i,(i+1):nperiod] <- 1
  }
  X <- X %x% rep(1, nclus_seq)
  
  
  
  # implied between-cluster variance
  sigma_alpha_sq <- sigma^2 * wpicc / (1 - wpicc)
  
  # implied between-cluster-period variance
  # sigma_gamma_sq <- sigma^2 * WPICC / (1 - WPICC) - sigma_alpha_sq
  
  eta_sq <- (wpicc_int*(sigma_alpha_sq + sigma^2) - sigma_alpha_sq)/ (1 - wpicc_int)
  # implied between-period intracluster correlation
  # BPICC <- WPICC * CAC
  
  
  #C <- rnorm(n=nclus, mean=0, sd = sqrt(sigma_alpha_sq))
  #CP <- rnorm(n=nclus * nperiod, mean=0, sd = sqrt(sigma_gamma_sq))
  #RT <- rnorm(n=nclus, mean = theta, sd = eta)
  
  sigma_ed <- matrix(NA, ncol = nperiod, nrow = nperiod)
  for (j in 1:nseq)
    for (k in (j+1):nperiod) {
      sigma_ed[j,k] = sigma_alpha_sq * cac^abs(j-k); 
      sigma_ed[k,j] = sigma_ed[j, k];
    }
  
  for (i in 1:nperiod){
    sigma_ed[i,i] = sigma_alpha_sq
  }
  
  ran_dtd <- mvrnorm(n = nclus, mu = rep(0,nperiod), sigma_ed)
  RT <- rnorm(n = nclus, mean = 0, sd = sqrt(eta_sq))
  
  
  Yijk <- c()
  e <- rnorm(n = nclus * nperiod * nsubject, mean = 0, sd = sigma)
  for (i in 1:nclus) {
    for (j in 1:nperiod) {
      # q <- (i - 1) * nperiod + j
      for (k in 1:nsubject) {
        # l is the index into the data
        l <- (i - 1) * nperiod * nsubject + (j - 1) * nsubject + k
        Yijk[l] <- j + X[i,j] * (theta + RT[i]) + ran_dtd[i, j] + e[l]
      }
    }
  }
  
  Xvec <- as.vector(kronecker(as.vector(t(X)), matrix(rep(1, nsubject), ncol=1)))
  Cvec <- as.factor(rep(1:nclus, each=nperiod*nsubject))
  CPvec <- as.factor(rep(1:(nclus*nperiod), each=nsubject))
  Pvec <- as.factor(rep(1:nperiod, each=nsubject, times=nclus))
  
  dat <- data.frame(
    response.var = Yijk,
    tx.var = Xvec,
    cluster.var = Cvec,
    time.var = Pvec,
    cp.var = CPvec
  )
  return(dat)
}

### Example
#data <- gendata_dtd(nclus = 8, #number of cluster
#                    nperiod = 5, #number of period
#                    cac = 0.8, #cluster autocorrelation
#                    wpicc = 0.01, #within-period ICC in control group
#                    wpicc_int = 0.05,
#                    sigma = 1, #individual error
#                    nsubject = 100, #number of indiviudals per period
#                    theta = 0 #standized treatment effect
#                    )
