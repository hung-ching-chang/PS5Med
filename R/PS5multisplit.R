#' High-dimensional mediation analysis by multiple sample splitting
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
#' @param multi.num
#' @param cores
#'
#' @details write XXX
#'
#' @return A \code{list} of containing mediation testing result.
#' \itemize{
#'     \item{global.test: }{}
#'     \item{global.me: }{}
#'     \item{global.me.prop: }{}
#'     \item{mediation.contri: }{}
#'     \item{detail.global.test: }{}
#'     \item{detail.global.me: }{}
#'     \item{detail.global.me.prop: }{}
#'     \item{detail.estimation: }{}
#' }
#'
#' @examples
#'
#' @author Hung-Ching Chang
#'
#' @import ncvreg
#' @import Matrix
#' @import mvtnorm
#' @import parallel
#'
#' @export
PS5.multisplit <- function(
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
    penalty = "MCP",
    multi.num = 100,
    cores = detectCores() - 1
){
  # multi-split
  result <- mclapply(c(1:multi.num),
                     function(i) PS5.mediation(M = M,
                                               X = X,
                                               Y = Y,
                                               C = C,
                                               n.draw = n.draw,
                                               dim.reduction = dim.reduction,
                                               seed = (i+seed),
                                               ga = ga,
                                               Y.family = Y.family,
                                               M.family = M.family,
                                               penalty = penalty),
                     mc.cores = cores)
  # summaize result
  PS5.global.me <- PS5.global.test <- PS5.total.me.proportion <- rep(NA,multi.num)
  PS5.estimation <- NULL
  for(i in 1:multi.num){
    PS5.global.me[i] <- result[[i]]$global.me
    PS5.global.test[i] <- result[[i]]$global.test
    PS5.total.me.proportion[i] <- result[[i]]$global.me.prop
    PS5.estimation[[i]] <- result[[i]]$estimation
  }
  # p-value for global IE
  if(any(PS5.global.test == 0, rm.na = T)){PS5.global.test[PS5.global.test == 0] <- 1e-16}
  aggr.global.test <- multisplit.pval(PS5.global.test[!is.na(PS5.global.test)],
                                      min.gamma = 0.5)

  # unbiased ME estimation
  global.me <- as.numeric(quantile(PS5.global.me, .5, type = 1, na.rm = T))
  unbiased.seed <- which(PS5.global.me == global.me)[1]
  tmp.mc <- result[[unbiased.seed]]$estimation$contribution.proportion
  positive.global.me.prop <- sum(tmp.mc[tmp.mc > 0])
  negative.global.me.prop <- sum(tmp.mc[tmp.mc < 0])
  global.me.prop <- c(quantile(PS5.total.me.proportion, .5, type = 1, na.rm = T),
                      positive.global.me.prop,
                      negative.global.me.prop)
  names(global.me.prop) <- c("total", "positive","negative")

  #p-value aggregation
  mediator.list <- unlist(lapply(PS5.estimation, function(x) row.names(x)))
  all.estimation <- do.call(rbind, PS5.estimation)
  unique.mediator.name <- unique(mediator.list)
  mediation.contri <- data.frame(matrix(NA, ncol = 6, nrow = length(unique.mediator.name),
                                        dimnames=list(NULL, c("preselected.prop",
                                                              "mediation contribution",
                                                              "contribution proportion",
                                                              "pval", "pval.BY", "pval.Bonf"))))
  for(i in 1:length(unique.mediator.name)){ # top 50 genes ranked by selected proportion
    tmp.result <- all.estimation[mediator.list == unique.mediator.name[i],]  # result for single mediator
    tmp.result$p.value[tmp.result$p.value == 0] <- 1/n.draw
    tmp.result$p.value.BY[tmp.result$p.value.BY == 0] <- 1/n.draw
    tmp.result$p.value.Bonf[tmp.result$p.value.Bonf == 0] <- 1/n.draw
    # p-value
    tmp.pvalue <- c(tmp.result$p.value,
                    rep(1, multi.num - nrow(tmp.result)))
    tmp.pvalue.BY <- c(tmp.result$p.value.BY,
                       rep(1, multi.num - nrow(tmp.result)))
    tmp.pvalue.Bonf <- c(tmp.result$p.value.Bonf,
                         rep(1, multi.num - nrow(tmp.result)))
    # keep 1/2 informative result to aggregate p-value
    aggr.pvalue <- multisplit.pval(tmp.pvalue, min.gamma = nrow(tmp.result)/multi.num/2)
    aggr.pvalue.BY <- multisplit.pval(tmp.pvalue.BY, min.gamma = nrow(tmp.result)/multi.num/2)
    aggr.pvalue.Bonf <- multisplit.pval(tmp.pvalue.Bonf, min.gamma = nrow(tmp.result)/multi.num/2)

    # mediation contribution
    mediation.contri[i,] <- c(nrow(tmp.result)/multi.num,
                              median(tmp.result$mediation.contribution),
                              median(tmp.result$contribution.proportion),
                              aggr.pvalue, aggr.pvalue.BY, aggr.pvalue.Bonf)
  }
  row.names(mediation.contri) <- unique.mediator.name

  # sort by p-value (BY)
  mediation.contri <- mediation.contri[order(mediation.contri$pval.BY),]

  return(list(global.test = aggr.global.test,
              global.me = global.me,
              global.me.prop = global.me.prop,
              mediation.contri = mediation.contri,
              detail.global.test = PS5.global.test,
              detail.global.me = PS5.global.me,
              detail.global.me.prop = PS5.total.me.proportion,
              detail.estimation = PS5.estimation))
}
