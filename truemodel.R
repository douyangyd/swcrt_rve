# fit DTD+RI model
fit_AR <- function(data){
  foreach(i = 1:length(data), .combine = rbind, .packages=c("glmmTMB","lmerTest")) %dopar%  {
    tryCatch({
      fit_AR <- glmmTMB::glmmTMB(response.var ~  tx.var + time.var + ar1(time.var + 0 | cluster.var) + (0 + tx.var|cluster.var), REML=T, 
                                 data = data[[i]], family = gaussian)
      V_AR <- VarCorr(fit_AR)
      res1 <- summary(fit_AR)
      wpicc_ar <- round(V_AR$cond$cluster.var[1,1]/sum(V_AR$cond$cluster.var[1,1],(attr(V_AR$cond,"sc"))^2),8)
      ar1 <- round(attr(V_AR$cond$cluster.var,"correlation")[1,2],8)
      res <- c(res1$coefficients$cond["tx.var",], wpicc_ar, ar1)
      #res2 <- c(res1$coefficients$cond["tx.var",])
      #res <- rbind(res,res2)
    }, warning=function(e) NULL)
  }
}


# fit NE+RI model
fit_HG <- function(data){
  foreach(i = 1:length(data), .combine = rbind, .packages=c("glmmTMB","lmerTest")) %dopar%  {
    tryCatch({
      fit_RT <- lmer(response.var ~  tx.var + time.var + (1|cluster.var) + (0 + tx.var|cluster.var) + (1|time.var:cluster.var), data = data[[i]])
      res1 <- summary(fit_RT)
      V_RT <- as.data.frame(VarCorr(fit_RT))
      res <- c(res1$coefficients["tx.var",], V_RT[,"vcov"])
    }, warning=function(e) NULL)
  }
}
  
#saveRDS(fit_AR, file = "fit_AR.rds")
#saveRDS(fit_HG, file = "fit_HG.rds")
