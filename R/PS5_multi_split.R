#' @title Mediation analysis using partial sum statistic and sample splitting
#'
#' @description aggregate p-value from multi-split result
#'
#' @param M 	  # n-by-p matrix
#' @param X 	  # n      vector
#' @param Y     # n      vector
#' @param C   	# n-by-l matrix, not including intercept
#' @param n.draw 	# no. of parametric bootstrapping sample
#' @param dim.reduction, # dimension reduction
#' @param seed
#' @param ga           # gamma
#' @param Y.family
#' @param M.family
#' @param pre.selection
#' @param lambda.selection
#' @param multi.num
#' @param cores
#'
#' @details write XXX
#'
#' @return XXXXXX
#'
#' @examples
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{PS5_mediation}}
#' @export
#' @import
PS5.multi.split <- function(
    M = mediators,		  # n-by-p matrix
    X = exposure,		  	# n      vector
    Y = outcome,			  # n      vector
    C = conf,		       	# n-by-q matrix, not including intercept
    n.draw = 10000,	  	# no. of parametric bootstrapping sample
    dim.reduction = NA, # dimension reduction
    seed = 1004,
    Y.family = "gaussian",
    M.family = "gaussian",
    pre.selection = "MCP",
    multi.num = 100,
    cores = detectCores()-1,
    ga = 2
){
  # multi-split
  result <- mclapply(c(1:multi.num),
                     function(i) PS5.mediation(M = M,
                                               X = X,
                                               Y = Y,
                                               C = C,
                                               dim.reduction = dim.reduction,
                                               seed = (i+seed),
                                               n.draw = n.draw,
                                               Y.family = Y.family,
                                               M.family = M.family,
                                               pre.selection = pre.selection,
                                               ga = ga),
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
