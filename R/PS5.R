#' High-dimensional mediation analysis by single sample splitting
#'
#' \code{PS5} is used to estimate and test global mediation effect and individual mediation contribution for high-dimensional causal mediation analysis.
#'
#' @param M a matrix of high-dimensional mediators. Row is samples, and column is variables.
#' @param X a vector of exposure.
#' @param Y a vector of outcome.
#' @param C a matrix of confounding variables
#' @param n.draw number of draws for parametric bootstrap
#' @param dim.reduction a logical value. If 'TRUE' a sample splitting strategy with penalized regression regression is used for dimension reduction; if 'FALSE' dimension reduction procedure will be skipped. Default is 'TRUE'.
#' @param seed a random seed
#' @param ga value of gamma parameter. Default is 2.
#' @param Y.family either 'gaussian' (default) or 'binomial', depending on the data type of outcome (Y)
#' @param M.family either 'gaussian' (default) or 'binomial', depending on the data type of mediator (M)
#' @param penalty the penalty method. Either 'MCP' (the default), 'SCAD', or 'lasso'.
#'
#' @details XXXXXXXXXXXXXXXXXXXXXX
#'
#' @return A \code{list} of containing mediation testing result.
#' \itemize{
#'     \item{global.test: }{}
#'     \item{global.me: }{}
#'     \item{global.me.prop: }{}
#'     \item{estimation: }{}
#' }
#'
#' @examples
#'
#' @author Hung-Ching Chang
#'
#' @import ncvreg
#' @import Matrix
#' @import mvtnorm
#'
#' @export
PS5 <-function(
    M = mediators,
    X = exposure,
    Y = outcome,
    C = conf,
    n.draw = 10000,
    dim.reduction = TRUE,
    seed = 1004,
    ga = 2,
    Y.family = "gaussian",
    M.family = "gaussian",
    penalty = "MCP"
){
  set.seed(seed)
  Mediator.name <- colnames(M)
  n <- dim(M)[1]
  p <- dim(M)[2]
  ######################################################################
  ######################## dimension reduction #########################
  ##################### (sample splitting method) ######################
  ######################################################################
  if(dim.reduction){
    ss.result <- sample.splitting(M,X,Y,C,
                                  Y.family = Y.family,
                                  penalty.type = penalty)
    if(length(ss.result$mediator.selected) == 0){
      message("None of mediators are selected ...")
      result <- list(global.test = NA,
                     global.me = NA,
                     global.me.prop = NA,
                     estimation = NA)
      return(result)
    }
    if(length(ss.result$mediator.selected) > n/2){
      message("Sparsity assumption of beta might be violated ...")
    }
    M.fold2 <- as.matrix(M[-ss.result$sample.idx, ss.result$mediator.selected])
    X.fold2 <- X[-ss.result$sample.idx]
    Y.fold2 <- Y[-ss.result$sample.idx]
    C.fold2 <- as.matrix(C[-ss.result$sample.idx,])
    p <- ncol(M.fold2)
    M.selected <- as.matrix(M[,ss.result$mediator.selected])
    Mediator.name <- Mediator.name[ss.result$mediator.selected]
  }else{
    M.fold2 <- M.selected <- M
    X.fold2 <- X
    Y.fold2 <- Y
    C.fold2 <- C
  }
  message("pre-select ", p, " mediators ...")
  ######################################################################
  ############################ Estimation ##############################
  ######################################################################
  ## MY model
  fit.y <- glm(Y.fold2 ~ M.fold2 + X.fold2 + C.fold2, family = Y.family)
  cov2 <- vcov(fit.y)[(1:p)+1,(1:p)+1]
  beta <- fit.y$coef[1+(1:p)]
  NDE <- fit.y$coef[2+p]

  fit.m2 <- lm(M.selected ~ C + X)
  cov1 <- vcov(fit.m2)[(1:p)*(2+ncol(C)), (1:p)*(2+ncol(C))]  # (2+ncol(C)) to select cov of alpha
  if(p == 1){
    alpha <- fit.m2$coef["X"]
  }else{
    alpha <- fit.m2$coef["X",]
  }

  me <- alpha*beta # mediation effect
  bcov <- bdiag(cov1, cov2)

  ######################################################################
  ########################### MC sampling ##############################
  ######################################################################
  # bootstrap sampling for p-value calculation
  theta.boot <- rmvnorm(n.draw,
                        mean = c(alpha, beta),
                        sigma = as.matrix(bcov))
  alpha.boot <- theta.boot[,1:p]
  beta.boot <- theta.boot[,p+1:p]
  me.boot <- as.matrix(alpha.boot*beta.boot)   # alpha * beta
  me.boot.c <- t(t(me.boot) - apply(me.boot, 2, mean))
  me.combin <- rbind(me, me.boot.c)
  ## partial sum test
  if(p == 1){
    # if only one mediator
    PS.test <- ps.pval <- mean(abs(me)^ga < abs(me.boot.c)^ga)
  }else{
    # if multiple mediators
    # L2 norm version
    order.me.combin <- t(apply(abs(me.combin)^ga, 1,
                               function(x) sort(x, decreasing = T)))
    ps.stat <- t(apply(order.me.combin, 1, function(x){
      ps.tmp <- NULL
      for(k in 1:length(x)){
        ps.tmp <- c(ps.tmp, sum(x[1:k]))
      }
      return(ps.tmp)
    }))
    # p-value for each partial sum statistic
    ps.pval <- rep(NA, ncol(ps.stat))
    for(j in 1:ncol(ps.stat)){
      ps.pval[j] <- mean(ps.stat[1,j] < ps.stat[2:(n.draw+1),j])
    }
    ps.cauchy.stat <- 1/p*sum(tan(pi*(0.5-ps.pval)))
    PS.test <- 1 - pcauchy(ps.cauchy.stat)
  }
  message("PS5 p-value: ", PS.test)
  global.me <- sum(me)
  ######################################################################
  ############ p-value for single mediation contribution ###############
  ######################################################################
  # sobel's test
  mediator.pval <- apply(me.combin, 2, function(x) mean(abs(x[1]) < abs(x[-1])))
  mediator.pval.Bonferroni <- p.adjust(mediator.pval, "bonferroni")
  mediator.pval.BY <- p.adjust(mediator.pval, "BY")

  ######################################################################
  ############################ output result ###########################
  ######################################################################
  # total.effect <- coef(glm(Y ~ X + C, family = Y.family))[2]
  total.effect <- NDE + global.me

  total.me.proportion <- global.me/total.effect * 100
  estimation <- data.frame(`alpha` = alpha,
                           `beta` = beta,
                           `mediation contribution` = me,
                           `contribution proportion` = me/total.effect*100,
                           `p value` = mediator.pval,
                           `p value BY` = mediator.pval.BY,
                           `p value Bonf` = mediator.pval.Bonferroni)
  row.names(estimation) <- Mediator.name
  estimation <- estimation[order(estimation$p.value.BY),]
  result <- list(global.test = PS.test,
                 global.me = global.me,
                 global.me.prop = total.me.proportion,
                 estimation = estimation)
  return(result)
}