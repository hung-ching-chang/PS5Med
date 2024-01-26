# Internal function: sample splitting with variable selection

sample.splitting <- function(
    M,
    X,
    Y,
    C,
    Y.family = "gaussian",
    penalty.type = "MCP"
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
  fold1.fit.mcp <- ncvreg(D, Y.res[ss], family = Y.family,
                          penalty = penalty.type,
                          returnX = FALSE, max.iter = 100000)
  selected.idx <- which.min(BIC(fold1.fit.mcp))
  lambda.mcp <- fold1.fit.mcp$lambda[selected.idx]
  fold1.fit.beta <- coef(fold1.fit.mcp, lambda = lambda.mcp)[1+1:p]
  mediator.selected <- which(fold1.fit.beta != 0)

  return(list(mediator.selected = mediator.selected,
              sample.idx = ss))
}

# Internal function: calculate p-value from multi sample splitting

multisplit.pval <- function(pvalue, min.gamma = 0.05){
  B <- length(pvalue)
  gamma <- seq(ceiling(min.gamma * B)/B, 1-1/B, by = 1/B)
  quant.gamma <- quantile(pvalue, gamma)/gamma
  penalty <- ifelse(length(gamma)>1, 1-log(min.gamma), 1)
  pval.final <- min(1, min(quant.gamma)*penalty)
  return(pval.final)
}
