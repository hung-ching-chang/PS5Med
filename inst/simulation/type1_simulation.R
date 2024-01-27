library(Matrix)
library(mvtnorm)
library(HIMA)
library(freebird)
library(parallel)
library(ncvreg)
library(PS5Med)
source("HP2016.R")

########################################
type1.sim <- function(
    N,
    p,
    cov,
    replication,
    scenario = c("I", "II", "III", "IV"),
    X.type = "conti"   # "conti" or "discrete"
){
  ### setting
  exponent <- abs(matrix(1:p - 1,
                         nrow = p,
                         ncol = p,
                         byrow = TRUE) - (1:p - 1))
  v <- cov^exponent
  Alpha.c <- rep(0.0, p)	                                          # C -> M (the association between covariate and mediators)
  Beta.x <- 0.3 					                                          # X -> Y (the association between exposure and outcome [conditional on mediators])
  Beta.c <- 0.1                                                     # C -> Y (the association between covariate and outcome)

  # result
  All.PVAL <- data.frame(matrix(NA, ncol = 4, nrow = replication,
                                dimnames = list(NULL, c("PS5", "HP2016", "HIMA", "HILMA"))))

  for(rep.time in 1:replication){
    print(rep.time)
    set.seed(rep.time + 0527)
    # Scenario I
    if(scenario == "I"){
      message("scenario 1 ...")
      Alpha.x <- rep(0, p)
      Beta.m <- rep(0, p)
    }else if(scenario == "II"){
      message("scenario 2 ...")
      Alpha.x <- rep(runif(1, min = 1, max = 3), p)                   # X -> M (the association between exposure and mediators)
      Beta.m <- rep(0, p)	                                            # M -> Y (the association between mediators and outcome [conditional on exposure])
    }else if(scenario == "III"){
      message("scenario 3 ...")
      Alpha.x <- rep(0, p)
      Beta.m <- c(rep(runif(1, min = 1, max = 3), 50), rep(0, p-50))
    }else{
      message("scenario 4 ...")
      Alpha.x <- c(rep(0, 50), rep(runif(1, min = 1, max = 3), p-50))
      Beta.m <- c(rep(runif(1, min = 1, max = 3), 50), rep(0, p-50))
    }
    ### generating data
    epsilon <- rmvnorm(n=N, mean=rep(0, p), v)
    if(X.type == "conti"){
      exposure <- rnorm(N, mean = 0, sd = 1)
    }else{
      exposure <- c(rep(0:1, each=N/3), rep(2, N/3))
      exposure <- c(exposure, rep(0, N-length(exposure)))
    }
    conf <- c(rep(0, ceiling(N/6)), rep(1, floor(N/3)), rep(0, floor(N/3)), rep(1, ceiling(N/6)))
    conf2 <- runif(N, 0, 1)
    mediator <- 0.3 + (conf%o%Alpha.c+exposure%o%Alpha.x+epsilon)
    error <- rnorm(N)
    diz <- 0.2 + conf*Beta.c + exposure*Beta.x + mediator%*%Beta.m + error  # no interaction

    conf.mat <- matrix(cbind(conf,conf2), ncol = 2)
    exposure.mat <- matrix(exposure, ncol = 1)
    PS5.result <- PS5(M = mediator,
                      X = exposure.mat,
                      Y = diz,
                      C = conf.mat,
                      n.draw = 10000,
                      dim.reduction = T)
    # Huang and Pan 2016
    adaptive.ratio <- 0.8
    HP2016.result <- HP2016(M = mediator,
                            X = exposure,
                            Y = diz,
                            C = conf.mat,
                            n.draw = 1000,
                            adaptive = adaptive.ratio)

    # HIMA method
    HIMA.result <- tryCatch(hima(X = exposure,
                                 Y = diz,
                                 M = mediator,
                                 COV.XM = conf.mat,
                                 COV.MY = conf.mat,
                                 Y.family = "gaussian",
                                 M.family = "gaussian",
                                 penalty = "MCP"),
                            error = function(e) HIMA_result <- data.frame(Bonferroni.p = 1))
    if(is.null(HIMA.result)){
      HIMA.pval <- 1
    }else{
      HIMA.pval <- min(HIMA.result$Bonferroni.p)
    }

    # HILMA
    if(rep.time <= 100){
      hilma.result <-  hilma(Y = diz,
                             G = mediator,
                             S = exposure.mat)
    }else{
      hilma.result <- list(pvalue_beta_hat = NA)
    }


    All.PVAL[rep.time,] <- c(PS5.result$global.test,
                             HP2016.result$HP2016,
                             HIMA.pval,
                             hilma.result$pvalue_beta_hat)
    if(rep.time %% 100 == 0) print(apply(All.PVAL < 0.05, 2, function(x) mean(x, na.rm=T))*100)
  }
  return(All.PVAL)
}


cov.setting <- 0.5   # TRY 0 AND 0.5
scenario.setting <- "IV"
N <- 500
p <- 1000
type1error <- type1.sim(N = N,
                        p = p,
                        cov = cov.setting,
                        replication = 2000,
                        scenario = scenario.setting)
# type1error
apply(type1error<0.05, 2, function(x) mean(x, na.rm=T))*100

