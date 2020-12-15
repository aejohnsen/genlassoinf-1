library(genlassoinf)
library(cghFLasso)
data(CGH)

# Get data and path
sigma = .46    # estimate from the end of the Hyun paper (5-fold CV)
n = length(CGH$GBM.y)
maxsteps = 25
consec = 2
D = makeDmat(n,type='tf',ord=0)
y0    = CGH$GBM.y
f0    = dualpathSvd2(y0, D=D, maxsteps, approx=T)


## Method 2: First way to /manually/ get stop-time-incorporated polyhedron.
mm = get.modelinfo(obj=f0, consec=2, sigma=sigma)
bic = mm$ic
stoptime = which.rise(bic,consec) - 1
stoptime = pmin(stoptime, n-consec - 1)
Gobj.stoprule = getGammat.with.stoprule(obj = f0, y = y0,
                                        condition.step = stoptime+consec, type ='tf',
                                        stoprule = "bic", sigma = sigma,
                                        consec = consec, maxsteps = maxsteps, D = D)
G = Gobj.stoprule$G
u = Gobj.stoprule$u


## AFTER running ONE OF Method 1-3: Form contrast and segment test p-value
locs = f0$states[[stoptime]]
pvals = rep(NA,length(locs))
vs = list()
for(ii in 1:length(locs)){
  this.sign = f0$pathobjs$s[which(f0$pathobjs$B == locs[ii])]
  my.v.lrt = make.v.tf.fp(test.knot = locs[ii],
                          adj.knot  = locs,
                          test.knot.sign = this.sign,
                          D=D)
  pval = poly.pval(y=y0, G=G, u=u, v=my.v.lrt, sigma=sigma)$pv
  cat("After fixed size model was selected, the segment test pvalue at location",
      locs[ii], "is", pval,fill=TRUE)
  pvals[ii] = pval
  vs[[ii]] = my.v.lrt
}

write.csv(cbind(locs,pvals),file="pvals.csv")

#############
# Conf. int #
#############

# Run everything below here for confidence intervals for the above contrasts

# Main p-value function

poly.pval <- function(y, G, u, v, sigma, bits=NULL) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  
  pv = tnorm.surv(z,0,sd,vlo,vup,bits)
  return(list(pv=pv,vlo=vlo,vup=vup))
}

# Main confidence interval function

poly.int <- function(y, G, u, v, sigma, alpha, gridrange=c(-100,100),
                     gridpts=100, griddepth=2, flip=FALSE, bits=NULL) {
  
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  
  xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  fun = function(x) { tnorm.surv(z,x,sd,vlo,vup,bits) }
  
  int = grid.search(xg,fun,alpha/2,1-alpha/2,gridpts,griddepth)
  tailarea = c(fun(int[1]),1-fun(int[2]))
  
  if (flip) {
    int = -int[2:1]
    tailarea = tailarea[2:1]
  }
  
  return(list(int=int,tailarea=tailarea))
}

##############################

# Assuming that grid is in sorted order from smallest to largest,
# and vals are monotonically increasing function values over the
# grid, returns the grid end points such that the corresponding
# vals are approximately equal to {val1, val2}

grid.search <- function(grid, fun, val1, val2, gridpts=100, griddepth=2) {
  n = length(grid)
  vals = fun(grid)
  
  ii = which(vals >= val1)
  jj = which(vals <= val2)
  if (length(ii)==0) return(c(grid[n],Inf))   # All vals < val1
  if (length(jj)==0) return(c(-Inf,grid[1]))  # All vals > val2
  # RJT: the above logic is correct ... but for simplicity, instead,
  # we could just return c(-Inf,Inf) 
  
  i1 = min(ii); i2 = max(jj)
  if (i1==1) lo = -Inf
  else lo = grid.bsearch(grid[i1-1],grid[i1],fun,val1,gridpts,
                         griddepth-1,below=TRUE)
  if (i2==n) hi = Inf
  else hi = grid.bsearch(grid[i2],grid[i2+1],fun,val2,gridpts,
                         griddepth-1,below=FALSE)
  return(c(lo,hi))
}

# Repeated bin search to find the point x in the interval [left, right]
# that satisfies f(x) approx equal to val. If below=TRUE, then we seek
# x such that the above holds and f(x) <= val; else we seek f(x) >= val.

grid.bsearch <- function(left, right, fun, val, gridpts=100, griddepth=1, below=TRUE) {
  n = gridpts
  depth = 1
  
  while (depth <= griddepth) {
    grid = seq(left,right,length=n)
    vals = fun(grid)
    
    if (below) {
      ii = which(vals >= val)
      if (length(ii)==0) return(grid[n])   # All vals < val (shouldn't happen)
      if ((i0=min(ii))==1) return(grid[1]) # All vals > val (shouldn't happen)
      left = grid[i0-1]
      right = grid[i0]
    }
    
    else {
      ii = which(vals <= val)
      if (length(ii)==0) return(grid[1])   # All vals > val (shouldn't happen)
      if ((i0=max(ii))==n) return(grid[n]) # All vals < val (shouldn't happen)
      left = grid[i0]
      right = grid[i0+1]
    }
    
    depth = depth+1
  }
  
  return(ifelse(below, left, right))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector

tnorm.surv <- function(z, mean, sd, a, b, bits=NULL) {
  z = max(min(z,b),a)
  
  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1
  
  # Try the multi precision floating point calculation first
  o = is.finite(mean)
  mm = mean[o]
  pp = mpfr.tnorm.surv(z,mm,sd,a,b,bits) 
  
  # If there are any NAs, then settle for an approximation
  oo = is.na(pp)
  if (any(oo)) pp[oo] = bryc.tnorm.surv(z,mm[oo],sd,a,b)
  
  p[o] = pp
  return(p)
}

##' Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, using
##' multi precision floating point calculations thanks to the Rmpfr package
mpfr.tnorm.surv <- function(z, mean=0, sd=1, a, b, bits=NULL) {
  # If bits is not NULL, then we are supposed to be using Rmpf
  # (note that this was fail if Rmpfr is not installed; but
  # by the time this function is being executed, this should
  # have been properly checked at a higher level; and if Rmpfr
  # is not installed, bits would have been previously set to NULL)
  if (!is.null(bits)) {
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)
    return(as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                        (Rmpfr::pnorm(b)-Rmpfr::pnorm(a))))
  }
  
  # Else, just use standard floating point calculations
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)
  
  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)
  
  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
           (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}

# Return Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# Riemann approximation tricks, by Max G'Sell

gsell.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  return(max.approx.frac(a/sd,b/sd,z/sd,mean/sd))
}


##############################

forwardStop <- function(pv, alpha=.10){
  if (alpha<0 || alpha>1) stop("alpha must be in [0,1]")
  if (min(pv,na.rm=T)<0 || max(pv,na.rm=T)>1) stop("pvalues must be in [0,1]")
  val=-(1/(1:length(pv)))*cumsum(log(1-pv))
  oo = which(val <= alpha)
  if (length(oo)==0) out=0
  else out = oo[length(oo)]
  return(out)
}

##############################

aicStop <- function(x, y, action, df, sigma, mult=2, ntimes=2) {
  n = length(y)
  k = length(action)
  aic = numeric(k)
  G = matrix(0,nrow=0,ncol=n)
  u = numeric(0)
  count = 0
  
  for (i in 1:k) {
    A = action[1:i]
    aic[i] = sum(lsfit(x[,A],y,intercept=F)$res^2) + mult*sigma^2*df[i]
    
    j = action[i]
    if (i==1) xtil = x[,j]
    else xtil = lsfit(x[,action[1:(i-1)]],x[,j],intercept=F)$res
    s = sign(sum(xtil*y))
    
    if (i==1 || aic[i] <= aic[i-1]) {
      G = rbind(G,s*xtil/sqrt(sum(xtil^2)))
      u = c(u,sqrt(mult)*sigma)
      count = 0
    }
    
    else {
      G = rbind(G,-s*xtil/sqrt(sum(xtil^2)))
      u = c(u,-sqrt(mult)*sigma)
      count = count+1
      if (count == ntimes) break
    }
  }
  
  if (i < k) {
    khat = i - ntimes
    aic = aic[1:i]
  }
  else khat = k
  
  return(list(khat=khat,G=G,u=u,aic=aic,stopped=(i<k)))
}

#these next two functions are used by the binomial and Cox options of fixedLassoInf

mypoly.pval.lee=
  function(y, A, b, eta, Sigma, bits=NULL) {
    # compute pvalues from poly lemma:  full version from Lee et al for full matrix Sigma
    nn=length(y)
    eta=as.vector(eta)
    temp = sum(eta*y)
    vv=as.numeric(matrix(eta,nrow=1,ncol=nn)%*%Sigma%*%eta)
    cc = Sigma%*%eta/vv
    
    z=(diag(nn)-matrix(cc,ncol=1)%*%eta)%*%y
    rho=A%*%cc
    
    vec = (b- A %*% z)/rho
    vlo = suppressWarnings(max(vec[rho<0]))
    vup = suppressWarnings(min(vec[rho>0]))
    sd=sqrt(vv)
    pv = tnorm.surv(temp,0,sd,vlo,vup,bits)
    return(list(pv=pv,vlo=vlo,vup=vup,sd=sd))
  }



mypoly.int.lee=
  function(y,eta,vlo,vup,sd, alpha, gridrange=c(-100,100),gridpts=100, griddepth=2, flip=FALSE, bits=NULL) {
    # compute sel intervals from poly lemmma, full version from Lee et al for full matrix Sigma
    
    temp = sum(eta*y)
    
    xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
    fun = function(x) { tnorm.surv(temp,x,sd,vlo,vup,bits) }
    
    int = grid.search(xg,fun,alpha/2,1-alpha/2,gridpts,griddepth)
    tailarea = c(fun(int[1]),1-fun(int[2]))
    
    if (flip) {
      int = -int[2:1]
      tailarea = tailarea[2:1]
    }
    
    return(list(int=int,tailarea=tailarea))
  }



mydiag=function(x){
  if(length(x)==1) out=x
  if(length(x)>1) out=diag(x)
  return(out)
}


for (i in 1:length(vs)){
  print(poly.int(y0,G,u,vs[[i]],sigma,.05)$int)
}


### Naive p-values and CIs

locsplus <- c(1,sort(locs),length(y0))
naive.pvals <- rep(0,length(locs))
cis <- matrix(nrow=length(locs),ncol=2)
estim <- rep(0,length(locs))
for (i in 1:length(locs)){
  st <- locsplus[i]
  I_j <- locsplus[i+1]
  en <- locsplus[i+2]
  print(c(st,I_j,en))
  x <- y0[st:I_j]
  y <- y0[(I_j+1):en]
  testobj <- t.test(x,y,var.equal=T)
  naive.pvals[i] <- testobj$p.value
  estim[i] <- testobj$estimate[2] - testobj$estimate[1]
  cis[i,] <- testobj$conf.int
}
round(estim,2)
round(naive.pvals,3)
cis

# Plot of chosen model

lassofit <- fusedlasso1d(y0)
plot(lassofit,lambda=lassofit$lambda[10],col="lightblue",main="Selected model",pch=16,xlab="Location",ylab="CGH")
abline(v=locs,col="green",lty=2)


# Sample splitting

train.ind <- seq(1,length(y0),by=2)
test.ind <- seq(2,length(y0),by=2)
train2.ind <- train.ind[seq(1,length(train.ind),by=2)]
test2.ind <- train.ind[seq(2,length(train.ind),by=2)]

train.fit <- fusedlasso1d(y0[train2.ind])
test <- y0[test2.ind]
fits <- rep(0,length(test))
MSE <- rep(0,length(train.fit$lambda))
for (i in 1:length(train.fit$lambda)){
  betas <- train.fit$beta[,i]
  fitting.f <- function(x){
    betas[ceiling(x/4)] 
  }
  for (j in 1:length(test)){
    x <- test2.ind[j]
    fits[j] <- fitting.f(x)
  }
  MSE[i] <- mean((fits-test)**2)
}
plot(MSE)
abline(h=min(MSE)+sd(MSE))

# pick 15th step
sigma = .46    # estimate from the end of the Hyun paper (5-fold CV)

n = length(train2.ind)
maxsteps = 25
D = makeDmat(n,type='tf',ord=0)
f0 = dualpathSvd2(y0[train2.ind], D=D, maxsteps, approx=T)
locs = train2.ind[f0$states[[15]]]

#
y1 = y0[test.ind]
locsplus <- c(1,sort(locs[-12]),length(y1))
naive.pvals <- rep(0,length(locs[-12]))
cis <- matrix(nrow=length(locs[-12]),ncol=2)
estim <- rep(0,length(locs[-12]))
for (i in 1:length(locs[-12])){
  print("step")
  print(i)
  st <- locsplus[i]
  I_j <- locsplus[i+1]
  en <- locsplus[i+2]
  print(c(st,I_j,en))
  x <- y1[st:I_j]
  y <- y1[(I_j+1):en]
  testobj <- t.test(x,y,var.equal=T)
  naive.pvals[i] <- testobj$p.value
  estim[i] <- testobj$estimate[2] - testobj$estimate[1]
  cis[i,] <- testobj$conf.int
}
round(estim,2)
cbind(sort(locs),round(naive.pvals,3))




