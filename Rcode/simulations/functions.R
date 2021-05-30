#### Bai and Ng (2002) PC3 ####
BNtest <- function (X,  kmax){
  T <- dim(X)[1]
  N <- dim(X)[2]
  rmax <- kmax
  v<-rep(0,rmax)
  kfactor<-1:rmax
  cNT<-min(N,T)
  Wy <- cov(X)
  bev<-eigen((t(X)%*%X)/T)$values
  for(k in 1:rmax){
    v[k]<-sum(bev[(k+1):N])
  }
  
  
  PC3<-v+v[rmax]*log(cNT)/cNT*kfactor
  r1 <-which.min(PC3)
  
  results <- list(r1=r1, ratio = v, Wy=Wy)
  return(results)
}

#### Onatski (2010)  #####
ONtest <- function (X,  kmax){
  rmax <- kmax
  nols<-4
  Wy = cov(X)
  noev<-eigen(Wy)$values
  oJ<-rep(1,(nols+1))
  ed<-noev[1:rmax]-noev[2:(rmax+1)]
  noj<-rmax+1
  for(j in 1:4){
    noy<-noev[noj:(noj+nols)]
    nox<-(noj+seq(-1,(nols-1),1))^(2/3)
    nobeta<-sum((nox-mean(nox))*(noy-mean(noy)))/sum((nox-mean(nox))^2)
    nodelte<-2*abs(nobeta)
    noer<-ifelse(max(ed-nodelte)<0,0,max(which(ed>=nodelte)))
    noj<-noer+1
  }
  r1 <- noer
  results <- list(r1=r1, Wy=Wy)
  return(results)
}


#### Lam and Yao (2011)  ####
LamYaotest <- function (X, K0, kmax){
  N <- dim(X)[1]
  P <- dim(X)[2]
  
  if (K0==1){
    dW=cov(X[2:N,], X[1:(N-1),])
    Wy=dW%*%t(dW)
    eA=eigen(Wy, symmetric=T)
    ev=sort(eA$values, decreasing = TRUE)
  } else {
    dW=cov(X[2:N,], X[1:(N-1),])
    Wy=dW%*%t(dW)
    for(k in 2:K0) { dW=cov(X[(k+1):N,], X[1:(N-k),])
    Wy=Wy + dW%*%t(dW)
    }
    eA=eigen(Wy, symmetric=T)
    ev=sort(eA$values, decreasing = TRUE)
  }
  ratio <- rep(0,kmax)
  for (j in 1:kmax){
    ratio[j] <- ev[j]/ev[(j+1)]
  }
  r1 = which.max(ratio)
  
  results <- list(r1=r1, ratio=ratio, Wy=Wy)
  return(results)
}
#### Ahn and Horenstein (2013)  ####

AhnHorensteintest <- function (X,kmax){
  N <- dim(X)[1]
  P <- dim(X)[2]
  
  Wy = var(X)
  eA = eigen(Wy, symmetric=T)
  ev = sort(eA$values, decreasing = TRUE)
  
  ratio <- rep(0,kmax)
  for (j in 1:kmax){
    ratio[j] <- ev[j]/ev[(j+1)]
  }
  r1 = which.max(ratio)
  results <- list(r1=r1, ratio=ratio, Wy=Wy)
  return(results)
}

#### Fan et al. (2020) ####
Fanetal <- function(X, kmax){
  pZ <- ncol(X)
  n <- nrow(X)
  rmax <- kmax
  sn<-cov(X);
  hatRR=cov2cor(sn); 
  lambdaZ<-eigen(hatRR)$values; 
  DD=NULL; 
  lambdaLY=lambdaZ; 
  p=pZ
  pp=rmax+2; 
  mz=rep(0,pp); 
  dmz=mz; 
  tmz=mz
  
  for (kk in 1:pp){qu=3/4
  lambdaZ1=lambdaZ[-seq(max(0, 1),kk,1)]; 
  z0=qu*lambdaZ[kk]+(1-qu)*lambdaZ[kk+1] 
  ttt=c(lambdaZ1-lambdaZ[kk], z0-lambdaZ[kk]); 
  y0=(length(lambdaZ1))/(n-1) 
  mz[kk]=-(1-y0)/lambdaZ[kk]+y0*mean(na.omit(1/ttt))
  }

  temp2018=(-1/mz)[-1]-1-sqrt(pZ/(n-1)); 
  temp1=seq(1,rmax,1); 
  temp2=cbind(temp1, temp2018[1:rmax])
  k00new=max((temp2[,1][temp2[,2]>0]), 0)+1 # estimated number of factors
  
  r1 = k00new
  results <- list(r1=r1, ratio=temp2, Wy = hatRR)
}

#### Caro and Pena (2021) ####
CaroPenatestCOR <- function (X, K0, kmax){
  N <- dim(X)[1]
  P <- dim(X)[2]
  
  wk0 = N/((K0+1)*(N-(K0/2)))
  Wy = wk0 * cor(X) %*%t (cor(X))
  dW = acf(X, lag.max = K0, type = "correlation", plot = FALSE)
  for(k in 1:K0) {
    dWk = dW$acf[(k+1),,]
    wk = (N-k)/((K0+1)*(N-(K0/2)))
    Wy = Wy + wk * dWk %*% t(dWk)
  }
  eA=eigen(Wy, symmetric=T)
  ev=sort(eA$values, decreasing = TRUE)
  
  ratio <- rep(0,kmax)
  for (j in 1:kmax){
    ratio[j] <- ev[j]/ev[(j+1)]
  }
  r1 = which.max(ratio)
  results <- list(r1=r1, ratio=ratio, Wy=Wy)
  return(results)
}



