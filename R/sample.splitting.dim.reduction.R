#' @title Mediation analysis using partial sum statistic and sample splitting
#'
#' @description aggregate p-value from multi-split result
#'
#' @param M,
#' @param X,
#' @param Y,
#' @param C,
#' @param Y.family
#' @param penalty.type
#' @param lambda.selection
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
sample.splitting.dim.reduction <- function(
    M,
    X,
    Y,
    C,
    Y.family = "gaussian",
    penalty.type = "MCP",
    lambda.selection = "BIC"
){
  n <- nrow(M)
  p <- ncol(M)
  ss <- sample(1:n, floor(n/2))
  if(Y.family == "gaussian"){
    M.res <- lm(M ~ X)$residuals
    Y.res <- lm(Y ~ X)$residuals
    D <- cbind(M.res,C)[ss,]
  }else{
    M.res <- M
    Y.res <- Y
    D <- cbind(M.res,X,C)[ss,]
  }

  # penalty regression
  fold1.fit.mcp <- ncvreg(D, Y.res[ss], family = Y.family, penalty = penalty.type,
                          returnX=FALSE, max.iter = 100000)
  selected.idx <- which.min(BIC(fold1.fit.mcp))
  lambda.mcp <- fold1.fit.mcp$lambda[selected.idx]
  fold1.fit.beta <- coef(fold1.fit.mcp, lambda = lambda.mcp)[1+1:p]
  mediator.selected <- which(fold1.fit.beta != 0)

  if(lambda.selection == "cv" &selected.idx > 95){
    # message("Use cross-validation to turn the lambda parameter ...")
    # BIC doesn't works well, then change to cross-validation
    cvfit <- cv.ncvreg(D, Y.res[ss], family = Y.family, penalty = penalty.type,
                       returnX=FALSE, max.iter = 100000)
    mediator.selected <- which(coef(cvfit)[1+1:p] != 0)
  }
  # do not keep too many mediator since it may affect the estimation process (type I error issue)
  # set a threshold: 2n/log(n)
  #num.mediator <- colSums(fold1.fit.mcp$beta!=0)
  #if(num.mediator[selected.idx] > n/log(n)){
  #  acceptable.lambda <- num.mediator < n/log(n)
  #  min.BIC <- min(BIC(fold1.fit.mcp)[acceptable.lambda])
  #  selected.idx <- which(BIC(fold1.fit.mcp) == min.BIC)
  #}

  return(list(mediator.selected = mediator.selected,
              sample.idx = ss))
}
