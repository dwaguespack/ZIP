zip.mle <- function(data){
  n <- length(data)
  xbar <- mean(data)
  posdata <- data[data>0]
  n1 <- length(posdata)
  n0 <- n - n1
  lambda0 <- xbar/(1 - n0/n)
  i <- 0
  repeat{
    lambdai <- lambda0 - (lambda0*(1-n0/n)-xbar*(1-exp(-lambda0)))/((1-n0/n)-xbar*exp(-lambda0))
    i <- i+1
    if(abs(lambdai-lambda0)<=.000001 | i>100){lambda.hat <- lambdai; break}
    lambda0 <- lambdai
  }
  xi.hat <- 1-xbar/lambda.hat
  if(n0==0){xi.hat <- 0}
  return(c(lambda.hat,xi.hat))
}