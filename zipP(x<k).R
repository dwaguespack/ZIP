pois.rand = function(nr, xi, lam){
  u = runif(nr)
  x = rep(0,nr) 
  psr = c() 
  psl = c()
  prob0 = xi+(1-xi)*exp(-lam)
  psl[1] = 0; psr[1] = prob0
  k = 2 
  repeat{
    psl[k] = psr[k-1]; psr[k] = psr[k-1]+(1-xi)*dpois(k-1,lam)
    if(psr[k] >= 1.0-1.0e-7){break}
    k = k+1
  }
  psl[k+1] <- psr[k]
  psr[k+1] <- 1
  sqx = seq(0, k+1, 1)
  for(i in 1:nr){
    ind = which(psl < u[i] & u[i] <= psr)
    x[i] = sqx[ind]
  }
  return(x)
}

cov.prob = function(nr, mr, xi, lam, n, k, cl){
  p.hat = c()
  al = (1-cl) 
  p <- xi + (1-xi)*exp(-lam)
  x = matrix(0,n,nr)
  for(h in 1:(k-1)){
    p <- p + (1-xi)*exp(-lam)*lam^(h)/(factorial(h))
  }
  covr = 0
  for(i in 1:nr){
    x[,i] = pois.rand(n, xi, lam)
  }
  for(i in 1:nr){
    for(j in 1:mr){
      xs = sample(x[,i], n, replace = T)
      p.hat[j] = length(xs[xs<k])/n
    }
    ci = quantile(p.hat, al)
    if(ci[1] <= p){covr = covr + 1}
  }
  return(covr/nr)
}