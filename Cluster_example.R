library(lmerTest)
library(clubSandwich)
library(foreach)
library(doParallel)
library(MASS)
library(glmmTMB)

# load data_generator
gendata_dtd <- readRDS(file = "gendata_dtd.rds")

# load function to run the true model (discrete time decay + RI or nested exchangeable + RI)
fit_AR <- readRDS(file = "fit_AR.rds")
fit_HR <- readRDS(file = "fit_HR.rds")

# load function to calculate confidence intervals
cal_confint <- readRDS(file = "cal_confint.rds")

# load the function used to summarize the simulation results
get_summary <- readRDS(file = "get_summary.rds")


set.seed(2022)
data <- replicate(2000, gendata_dtd(nclus = 8, #number of cluster
                                   nperiod = 5, #number of period
                                   cac = 0.8, #cluster autocorrelation
                                   wpicc = 0.01, #within-period ICC in control group
                                   wpicc_int = 0.05,
                                   sigma = 1, #individual error
                                   nsubject = 10, #number of indiviudals per period
                                   theta = 0 #standized treatment effect
), simplify = FALSE)


# caluclate robust variance estimators
get_rve <- function(m1, m2){
  rob_hh_cr0 <- sqrt(diag(clubSandwich::vcovCR(m1, type = "CR0")))[2]
  rob_hg_cr0 <- sqrt(diag(clubSandwich::vcovCR(m2, type = "CR0")))[2]
  rob_hh_cr1 <- sqrt(diag(clubSandwich::vcovCR(m1, type = "CR1")))[2]
  rob_hg_cr1 <- sqrt(diag(clubSandwich::vcovCR(m2, type = "CR1")))[2]
  rob_hh_cr1p <- sqrt(diag(clubSandwich::vcovCR(m1, type = "CR1p")))[2]
  rob_hg_cr1p <- sqrt(diag(clubSandwich::vcovCR(m2, type = "CR1p")))[2]
  rob_hh_cr1S <- sqrt(diag(clubSandwich::vcovCR(m1, type = "CR1S")))[2]
  rob_hg_cr1S <- sqrt(diag(clubSandwich::vcovCR(m2, type = "CR1S")))[2]
  rob_hh_cr2 <- sqrt(diag(clubSandwich::vcovCR(m1, type = "CR2")))[2]
  rob_hg_cr2 <- sqrt(diag(clubSandwich::vcovCR(m2, type = "CR2")))[2]
  rob_hh_cr3 <- sqrt(diag(clubSandwich::vcovCR(m1, type = "CR3")))[2]
  rob_hg_cr3 <- sqrt(diag(clubSandwich::vcovCR(m2, type = "CR3")))[2]
  return(c(rob_hh_cr0, rob_hh_cr1, rob_hh_cr1p, rob_hh_cr1S, rob_hh_cr2, rob_hh_cr3,
           rob_hg_cr0, rob_hg_cr1, rob_hg_cr1p, rob_hg_cr1S, rob_hg_cr2, rob_hg_cr3)
  )
}



fit <- function(dat){
  # fit exchangeable model
  fit_HH <- lmerTest::lmer(response.var ~  tx.var + time.var + (1|cluster.var), data = dat)
  # fit nested exchangeable model
  fit_HG <- lmerTest::lmer(response.var ~  tx.var + time.var + (1|cluster.var) + (1|time.var:cluster.var), data = dat) 
  s1 <- summary(fit_HH)
  s2 <- summary(fit_HG)
  rve <- get_rve(m1=fit_HH, m2=fit_HG)
  result <- c(as.numeric(s1$coefficients["tx.var",]), rve[1:6], as.numeric(s2$coefficients["tx.var",]), rve[7:12])
  return(result)
}


ncores <- Sys.getenv("SLURM_CPUS_PER_TASK")
ncores <- as.numeric(ncores) 
ncores
registerDoParallel(ncores)

results <- foreach(i = 1:length(data), .combine = rbind, .packages=c("lmerTest","clubSandwich")) %dopar% {tryCatch({fit(data[[i]])}, error = function(e) print(e))} 
sumdtd <- round(get_summary(results),4)

saveRDS(sumdtd , file = "sumdtd_1.rds")
saveRDS(results, file = "resultdtd_1.rds")
saveRDS(data, file = "datadtd_1.rds")


ture_model <- fit_AR(data)
saveRDS(ture_model, file = "true_1.rds")

