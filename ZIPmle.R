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
    if(lambda0=="NaN" | lambdai=="NaN"){print(data)}
    if(abs(lambdai-lambda0)<=.000001 | i>100){lambda.hat <- lambdai; break}
    lambda0 <- lambdai
  }
  xi.hat <- 1-xbar/lambda.hat
  if(n0==0){xi.hat <- 0}
  return(c(lambda.hat,xi.hat))
}


asymptotic.variance <- function(mles,n){
  lam.hat <- mles[1]; xi.hat <- mles[2]
  P0 <- (xi.hat + (1-xi.hat)*exp(-lam.hat))
  I11 <- ((1-exp(-lam.hat))^2)/(P0*(1-xi.hat))
  I22 <- (1-xi.hat)*(1-exp(-lam.hat))/lam.hat - (1-xi.hat)*exp(-lam.hat)*xi.hat/P0
  I12 <- -exp(-lam.hat)/P0
  I <- matrix(c(I11,I12,I12,I22),2,2)
  I.inv <- solve(I)
  gradient.mle <- c(-lam.hat,1-xi.hat)
  var <- t(gradient.mle)%*%I.inv%*%gradient.mle/n
  return(abs(var))
}

asymptotic.variance.2 <- function(mles,data){
  xbar <- mean(data)
  n <- length(data)
  n0 <- length(data[data==0])
  lam.hat <- mles[1]; xi.hat <- mles[2]
  P0 <- (xi.hat + (1-xi.hat)*exp(-lam.hat))
  I11 <- (n0*(1-exp(-lam.hat))^2)/(P0^2) + (n-n0)/((1-xi.hat)^2)
  I22 <- n*xbar/(lam.hat^2) - n0*xi.hat*(1-xi.hat)*exp(-lam.hat)/(P0^2)
  I12 <- -n0*exp(-lam.hat)/(P0^2)
  I <- matrix(c(I11,I12,I12,I22),2,2)
  I.inv <- solve(I)
  gradient.mle <- c(-lam.hat,1-xi.hat)
  var <- t(gradient.mle)%*%I.inv%*%gradient.mle
  return(abs(var))
}

asymptotic.variance.3 <- function(mles,data){
  xbar <- mean(data)
  n <- length(data)
  n0 <- length(data[data==0])
  lam.hat <- mles[1]; xi.hat <- mles[2]
  P0 <- (xi.hat + (1-xi.hat)*exp(-lam.hat))
  I11 <- n*(1-exp(-lam.hat))/P0 + n*(1-P0)/((1-xi.hat)^2)
  I22 <- n*(1-xi.hat)/lam.hat - n*xi.hat*(P0-xi.hat)/P0
  I12 <- -n*exp(-lam.hat)/P0
  I <- matrix(c(I11,I12,I12,I22),2,2)
  I.inv <- solve(I)
  gradient.mle <- c(-lam.hat,1-xi.hat)
  var <- t(gradient.mle)%*%I.inv%*%gradient.mle
  return(abs(var))
}

ci.mean.zip <- function(data,cl){
  n <- length(data)
  mles <- zip.mle(data)
  lam.hat <- mles[1]
  xi.hat <- mles[2]
  pe.mean <- (1-xi.hat)*lam.hat
  se <- sqrt(asymptotic.variance.2(mles,data))
  cv <- qt(1-(1-cl)/2,n-1)
  return(c(pe.mean-se*cv,pe.mean+se*cv))
}

ci.mean.zip2 <- function(data,cl){
  n <- length(data)
  mles <- zip.mle(data)
  lam.hat <- mles[1]
  xi.hat <- mles[2]
  pe.mean <- (1-xi.hat)*lam.hat
  se1 <- sqrt(asymptotic.variance.3(mles,data))
  se2 <- sqrt(asymptotic.variance(mles,n))
  cv <- qt(1-(1-cl)/2,n-1)
  print(c(pe.mean-se1*cv,pe.mean+se1*cv))
  return(c(pe.mean-se2*cv,pe.mean+se2*cv))
}

cov.prb.ci.mean <- function(sim,n,xi,lam,cl){
  cov <- 0
  mu <- (1-xi)*lam
  l <- c()
  for(i in 1:sim){
    ci <- ci.mean.zip(pois.rand(n,xi,lam),cl)
    if(ci[1]<mu & mu<ci[2]){cov <- cov+1}
    l[i] <- ci[2]-ci[1]
  }
  return(c(cov/sim,mean(l)))
}

