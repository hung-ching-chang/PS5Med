#' @title Mediation analysis using partial sum statistic and sample splitting
#'
#' @description aggregate p-value from multi-split result
#'
#' @param pvalue XXX
#' @param min.gamma XXXX
#'
#' @details write XXX
#'
#' @return XXXXXX
#'
#' @examples
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{PS5_multi_split}}
#' @export
#' @import
multisplit.pval <- function(pvalue, min.gamma = 0.05){
  B <- length(pvalue)
  gamma <- seq(ceiling(min.gamma * B)/B, 1-1/B, by = 1/B)
  quant.gamma <- quantile(pvalue, gamma)/gamma
  penalty <- ifelse(length(gamma)>1, 1-log(min.gamma), 1)
  pval.final <- min(1, min(quant.gamma)*penalty)
  return(pval.final)
}
