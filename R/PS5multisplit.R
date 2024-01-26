#' High-dimensional mediation analysis by multiple sample splitting
#'
#' \code{PS5.multisplit} enhances the robustness through multi sample split for estimating and testing global mediation effect and individual mediation contribution.
#'
#' @param M a matrix of high-dimensional mediators. Row is samples, and column is variables.
#' @param X a vector of exposure.
#' @param Y a vector of outcome.
#' @param C a matrix of confounding variables
#' @param n.draw number of draws for parametric bootstrap. Default is 10,000.
#' @param dim.reduction a logical value. If 'TRUE' (default) a sample splitting strategy with penalized regression regression is used for dimension reduction; if 'FALSE' dimension reduction procedure will be skipped..
#' @param seed a random seed. Default is 1004.
#' @param ga value of gamma parameter. Default is 2.
#' @param Y.family either 'gaussian' (default) or 'binomial', depending on the data type of outcome (Y)
#' @param M.family either 'gaussian' (default) or 'binomial', depending on the data type of mediator (M)
#' @param penalty the penalty method. Either 'MCP' (the default), 'SCAD', or 'lasso'.
#' @param multi.num number of multi sample split. Default is 100.
#' @param cores number of cores in parallel computing. Default is 1.
#'
#'
#' @return A \code{list} containing mediation testing result.
#' \itemize{
#'     \item{global.test: }{statistical significance of global mediation effect}
#'     \item{global.me: }{estimated global mediation effect}
#'     \item{global.me.pct: }{estimated global mediation percentage}
#'     \item{mediation.contri: }{a data.frame of individual mediation contribution}
#'     \item{detail.global.test: }{a vector containing global.test from each sample split}
#'     \item{detail.global.me: }{a vector containing global.me from each sample split}
#'     \item{detail.global.me.pct: }{a vector containing global.me.pct from each sample split}
#'     \item{detail.mediation.contri: }{ A \code{list} containing mediation.contri from each sample split}
#' }
#'
#' @examples XXXXXXXXXX
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
    cores = 1
){
  # multi-split
  result <- mclapply(c(1:multi.num),
                     function(i) PS5(M = M,
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
  detailed.global.me <- detailed.global.test <- detailed.total.me.pct <- rep(NA, multi.num)
  detailed.mediation.contri <- NULL
  for(i in 1:multi.num){
    detailed.global.me[i] <- result[[i]]$global.me
    detailed.global.test[i] <- result[[i]]$global.test
    detailed.total.me.pct[i] <- result[[i]]$global.me.pct
    detailed.mediation.contri[[i]] <- result[[i]]$mediation.contri
  }
  # p-value for global IE
  if(any(detailed.global.test == 0, rm.na = T)){
    detailed.global.test[detailed.global.test == 0] <- 1e-16
  }
  aggr.global.test <- multisplit.pval(detailed.global.test[!is.na(detailed.global.test)],
                                      min.gamma = 0.5)

  # unbiased ME estimation
  global.me <- as.numeric(quantile(detailed.global.me, .5, type = 1, na.rm = T))
  unbiased.seed <- which(detailed.global.me == global.me)[1]
  tmp.mc <- result[[unbiased.seed]]$mediation.contri$contribution.pct
  positive.global.me.pct <- sum(tmp.mc[tmp.mc > 0])
  negative.global.me.pct <- sum(tmp.mc[tmp.mc < 0])
  global.me.pct <- c(quantile(detailed.total.me.pct, .5, type = 1, na.rm = T),
                      positive.global.me.pct,
                      negative.global.me.pct)
  names(global.me.pct) <- c("total", "positive", "negative")

  #p-value aggregation
  mediator.list <- unlist(lapply(detailed.mediation.contri, function(x) row.names(x)))
  all.mediation.contri <- do.call(rbind, detailed.mediation.contri)
  unique.mediator.name <- unique(mediator.list)
  mediation.contri <- data.frame(matrix(NA, ncol = 6, nrow = length(unique.mediator.name),
                                        dimnames=list(NULL, c("preselected.pct",
                                                              "mediation contribution",
                                                              "contribution pct",
                                                              "pval", "pval.BY", "pval.Bonf"))))
  for(i in 1:length(unique.mediator.name)){ # top 50 genes ranked by selected percentage
    tmp.result <- all.mediation.contri[mediator.list == unique.mediator.name[i],]  # result for single mediator
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
                              median(tmp.result$contribution.pct),
                              aggr.pvalue, aggr.pvalue.BY, aggr.pvalue.Bonf)
  }
  row.names(mediation.contri) <- unique.mediator.name

  # sort by p-value (BY)
  mediation.contri <- mediation.contri[order(mediation.contri$pval.BY),]

  return(list(global.test = aggr.global.test,
              global.me = global.me,
              global.me.pct = global.me.pct,
              mediation.contri = mediation.contri,
              detail.global.test = detailed.global.test,
              detail.global.me = detailed.global.me,
              detail.global.me.pct = detailed.total.me.pct,
              detail.mediation.contri = detailed.mediation.contri))
}
