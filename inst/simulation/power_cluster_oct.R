#install.packages("Matrix")
library(Matrix)
#install.packages("mvtnorm")
library(mvtnorm)
#install.packages("HIMA")
library(HIMA)
#install.packages("freebird")
library(freebird)
library(parallel)
library(ncvreg)
source("../HP2016.R")
source("../sample.splitting.dim.reduction.R")
source("../PS5_mediation.R")
source("../multispit_pvalue.R")
source("../PS5_multi_split.R")

#################################################
### simulation function
power.sim <- function(
    N, 
    p,
    cov, 
    Alpha.effect,
    num.ts,
    replication = 100,
    multi.split = 100,
    X.type = "conti"
){
  # power setting
  Beta.effect <- 1
  Alpha.c <- rep(0.0, p)								                            # C -> M (the association between covariate and mediators)
  Beta.x <- 0.3 									                                  # X -> Y (the association between exposure and outcome [conditional on mediators])
  Beta.c <- 0.1                                                     # C -> Y (the association between covariate and outcome)
  Alpha.x <- rep(0, p)
  Alpha.x[1:num.ts] <- Alpha.effect
  Beta.m <- c(rep(Beta.effect, 50), rep(0, p-50))
  if(cov == 0.9){
    # create moderate correlation structure
    exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
    v <- 0^exponent
    # ultra-high correlation among pair of mediators (e.g. M1 and M2, M3 and M4)
    for(i in 1:floor(num.ts/2)){
      v[(2*i-1),2*i] <- v[2*i,(2*i-1)] <- cov
    }
  }else{
    exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p - 1))
    v <- cov^exponent                                                 # covariance matrix of residual (AR1 model)
  }
  # three measures
  All.PVAL <- data.frame(matrix(NA, ncol = 4, nrow = replication, 
                                dimnames=list(NULL, c("PS5", "HP2016", "HIMA", "HILMA"))))
  global.me <- data.frame(matrix(NA, ncol = 4, nrow = replication, 
                                 dimnames=list(NULL, c("PS5", "HP2016", "HIMA", "HILMA"))))
  HIMA.med.contri.pval <- data.frame(matrix(1, ncol = p, nrow = replication, 
                                            dimnames=list(NULL, paste0("M", c(1:p)))))
  PS5.med.contri.pval <- data.frame(matrix(1, ncol = p, nrow = replication, 
                                           dimnames=list(NULL, paste0("M", c(1:p)))))
  message("True global ME: ", sum(Alpha.x*Beta.m))
  true.me <- sum(Alpha.x*Beta.m)
  for(rep.time in 1:replication){
    print(rep.time)
    set.seed(rep.time*1357)
    epsilon <- rmvnorm(n=N, mean=rep(0, p), as.matrix(v))
    if(X.type == "conti"){
      exposure <- rnorm(N, mean = 0, sd = 1)
    }else{
      exposure <- c(rep(0:1, each=N/3), rep(2, N/3))
      exposure <- c(exposure, rep(0, N-length(exposure)))
    }
    conf <- c(rep(0, ceiling(N/6)), rep(1, floor(N/3)), rep(0, floor(N/3)), rep(1, ceiling(N/6)))
    conf2<-runif(N, 0, 1)
    mediator <- 0.3 + (conf%o%Alpha.c+exposure%o%Alpha.x+epsilon)
    colnames(mediator) <- paste0("M", c(1:p))
    error <- rnorm(N, 0, 1)
    diz <- 0.2 + conf*Beta.c + exposure*Beta.x + mediator%*%Beta.m + error  # no interaction
    conf.mat <- matrix(cbind(conf,conf2), ncol = 2)
    exposure.mat <- matrix(exposure, ncol = 1)
    # multi sample split (L2)
    PS5.result <- PS5.multi.split(M = mediator, 
                                  X = exposure.mat, 
                                  Y = diz, 
                                  C = conf.mat,
                                  n.draw = 10000,
                                  dim.reduction = T,
                                  multi.num = multi.split,
                                  cores = 7)
    selected.mediators <- sapply(row.names(PS5.result$mediation.contri), function(x) which(x == paste0("M", c(1:p))))
    PS5.med.contri.pval[rep.time,selected.mediators] <- PS5.result$mediation.contri$pval.Bonf
    
    # Huang and Pan 2016
    adaptive.ratio <- 0.8
    HP2016.result <- HP2016(M=mediator,
                            X=exposure,
                            Y=diz,
                            C=conf.mat,
                            n.draw=1000,
                            adaptive=adaptive.ratio)
    
    # HIMA method
    HIMA.result <- hima(X = exposure,
                        Y = diz,
                        M = mediator,
                        COV.XM = conf.mat,
                        COV.MY = conf.mat,
                        Y.family = "gaussian",
                        M.family = "gaussian",
                        penalty = "MCP")
    if(is.null(HIMA.result)){
      HIMA.pval <- 1
    }else{
      HIMA.pval <- min(HIMA.result$Bonferroni.p)
    }
    selected.mediators <- sapply(row.names(HIMA.result), function(x) which(x == paste0("M", c(1:p))))
    HIMA.med.contri.pval[rep.time,selected.mediators] <- HIMA.result$Bonferroni.p
    #HILMA
    hilma.result <-  hilma(Y = diz, 
                           G = mediator, 
                           S = exposure.mat)
    ######################################################################
    ############################### power ################################ 
    ######################################################################
    All.PVAL[rep.time,] <- c(PS5.result$detail.global.test[multi.split], 
                             HP2016.result$HP2016,
                             HIMA.pval,
                             hilma.result$pvalue_beta_hat)
    global.me[rep.time,] <- c(median(PS5.result$global.me),
                              HP2016.result$global.me,
                              sum(HIMA.result[,4]),
                              hilma.result$beta_hat)
  }
  power <- apply(All.PVAL<0.05, 2, function(x) mean(x, na.rm=T))
  ME.bias <- apply(global.me, 2, function(x) mean(abs(x-true.me)/true.me, na.rm=T))
  return(list(power = power, 
              ME.bias = ME.bias,
              PS5.contri.pval = PS5.med.contri.pval,
              HIMA.contri.pval = HIMA.med.contri.pval))
}

### simulation start
cov.setting <- 0       # c(0,0.5,0.9)
num.ts.setting <- 5    # c(5,10,30)
effect.a.setting <- seq(0, 0.4, 0.01)
effect.a.setting <- effect.a.setting[-1]
N.setting <- 500
p.setting <- 1000
# fix beta and change alpha
count <- 0
power.result <- NULL
for(effect.a in effect.a.setting){
  count <- count + 1
  print(paste0('alpha = ', effect.a, '; cor = ', cov.setting))
  power.result[[count]] <- power.sim(N = N.setting,
                                     p = p.setting,
                                     cov = cov.setting,
                                     Alpha.effect = effect.a,
                                     num.ts = num.ts.setting,
                                     replication = 100,
                                     multi.split = 100)
  save(power.result, file = paste0("power_N", N.setting, "p", p.setting,"_cor", cov.setting, "ts", num.ts.setting, ".RData"))
}
