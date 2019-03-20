pois.rand = function(nr, xi, lam){
  repeat{
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
  if(sum(x)>0){break}
  }
  
  return(x)
}

cov.prob = function(nr, mr, xi, lam, n, cl){
  xb = c()
  al = .5*(1-cl) 
  x = matrix(0,n,nr)
  mu = (1-xi)*lam
  covr = 0
  for(i in 1:nr){
    x[,i] = pois.rand(n, xi, lam)
  }
  for(i in 1:nr){
    for(j in 1:mr){
      xs = sample(x[,i], n, replace = T)
      xb[j] = mean(xs)
    }
  ci = quantile(xb, c(al,1-al))
  if(ci[1] <= mu & mu <= ci[2]){covr = covr + 1}
  }
  return(covr/nr)
}

table <- function(){
  co <- c(1,3,5,7,9)
  n <- c(30,50,100)
  lam <- c(2,3,4,5,10)
  xi <- c(.1,.2,.3,.4,.5)
  t <- matrix(rep(0,3*5*5),3*5,5*2)
  for(i in 1:3){
    print(c("i ",i))
    for(j in 1:5){
      print(c("j ",j))
      for(k in 1:5){
        print(c("k ",k))
        prb <- cov.prb.ci.mean(10000,n[i],xi[j],lam[k],.95)
        t[j+5*(i-1),co[k]] <- prb[1]
        t[j+5*(i-1),co[k]+1] <- prb[2]
      }
    }
  }
  return(t)                                  
} 

