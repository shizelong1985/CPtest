library(parallel)
library(Metrics)
library(matrixStats)
library(doParallel)
library(forecast)
library(MASS)
library(ggplot2)
library(fBasics)

# Monte Carlo exercise: Test Number of Factors
# DATA GENERATING PROCESS: sigma_E = idiosyncratic var-cov matrix
# THREE COMMON FACTORS r = 3
# DGP4_A -> sigma_E is diagonal with sigma_ei = 1
# DGP5_A -> sigma_E is diagonal with sigma_ei = 1 if i is odd and sigma_ei = 2 if i is even
# DGP6_A -> sigma_E is not diagonal (cross-section dependence) and diagonal elements sigma_ei = 1 if i is odd and sigma_ei = 2 if i is even
#
# DGP_B -> eit = theta * ei(t-1) + uit for i = 1,..,P theta N(0.5, 0.05)
#
# DGP4, DGP5, DGP6 two common factors 
# A1, B1 three common factors with autoregressive coeff = 0.9 for F1, 0.8 for F2 and 0.7 for F3
# A2, B2 three common factors with autoregressive coeff = 0.6 for F1, 0.5 for F2 and 0.4 for F3
# A3, B3 three common factors with autoregressive coeff = 0.6 for F1, 0.5 for F2 and 0.4 for F3 which errors N(0,sd=0.5)

setwd("set your path here/simulations")
start_time <- Sys.time()
############################ DGP4_A1 ############################

# register cluster and call the main.R 
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

# generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #
N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.9 # AR coefficient of factor process
a2 <- 0.8
a3 <- 0.7
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    r<-NTR[1]
    P<-NTR[2]
    N<-NTR[3]
    
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT    ####
      
      Gamma_E <- diag(1,P,P) # pos-def sym mtx specifying the cov mtx of the variables
      E <- mvrnorm(n = N, rep(0, P), Gamma_E)
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
    
    scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
    sink()
    
    bn <- sum(media[,1] == 3)
    on <- sum(media[,2] == 3)
    ly <- sum(media[,3] == 3)
    ah <- sum(media[,4] == 3)
    fan <- sum(media[,5] == 3)
    cp <- sum(media[,6] == 3)
    meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP4_A1.RDS"))
}

stopCluster(cl)

############################ DGP5_A1 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

# generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.9 # AR coefficient of factor process
a2 <- 0.8
a3 <- 0.7
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_e <- 4 + rnorm(P, 0, 0.5) 
      sigma_e[seq(from = 1, to = P, by = 2)] <- 1 
      Gamma_E <- diag(sigma_e[1:P],P,P) # pos-def sym mtx specifying the cov mtx of the variables
      E <- mvrnorm(n = N, rep(0, P), Gamma_E)
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP5_A1.RDS"))
}

stopCluster(cl)

############################ DGP6_A1 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.9 # AR coefficient of factor process
a2 <- 0.8
a3 <- 0.7
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_e <- 4 + rnorm(P, 0, 0.5) 
      sigma_e[seq(from = 1, to = P, by = 2)] <- 1
       
      
      CovMatrix <- matrix(0,P,P)
      for (a in 1:P){
        for (b in 1:P){
          CovMatrix[a,b] = 0.7^(abs(a-b))*sqrt(sigma_e[a])* sqrt(sigma_e[b])
          CovMatrix[b,a] = CovMatrix[a,b]
        }
      }
      
      E <- mvrnorm(n = N, rep(0, P), CovMatrix)
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP6_A1.RDS"))
}

stopCluster(cl)


############################ DGP4_B1 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.9 # AR coefficient of factor process
a2 <- 0.8
a3 <- 0.7

iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    r<-NTR[1]
    P<-NTR[2]
    N<-NTR[3]
    
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      CovMatrix <- diag(1,P,P)
      u <- mvrnorm(n = N, rep(0, P), CovMatrix) 
      E <- matrix(0,N,P)
      for (m in 1:(P/2)){
        for (n in 2:N){
          phi <- 0.5 + rnorm(1, 0, 0.05)
          E[n,m] <- phi*E[(n-1),m] + u[n,m]
        }
      }
      for (m in (P/2+1):P){
        E[,m] <- u[,m]
      }
      
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP4_B1.RDS"))
}

stopCluster(cl)

############################ DGP5_B1 ############################
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.9 # AR coefficient of factor process
a2 <- 0.8
a3 <- 0.7

iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_u <- 4 + rnorm(P, 0, 0.5) 
      sigma_u[seq(from = 1, to = P, by = 2)] <- 1 
      CovMatrix <- diag(sigma_u[1:P],P,P)
      
      u <- mvrnorm(n = N, rep(0, P), CovMatrix) # normal (0,1)
      E <- matrix(0,N,P)
      for (m in 1:(P/2)){
        for (n in 2:N){
          phi <- 0.5 + rnorm(1, 0, 0.05)
          E[n,m] <- phi*E[(n-1),m] + u[n,m]
        }
      }
      for (m in (P/2+1):P){
        E[,m] <- u[,m]
      }
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP5_B1.RDS"))
}

stopCluster(cl)

############################ DGP6_B1 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.9 # AR coefficient of factor process
a2 <- 0.8
a3 <- 0.7

iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_u <- 4 + rnorm(P, 0, 0.5) 
      sigma_u[seq(from = 1, to = P, by = 2)] <- 1
      
      CovMatrix <- matrix(0,P,P)
      for (a in 1:P){
        for (b in 1:P){
          CovMatrix[a,b] = 0.7^(abs(a-b))*sqrt(sigma_u[a])* sqrt(sigma_u[b]) 
          CovMatrix[b,a] = CovMatrix[a,b]
        }
      }
      
      u <- mvrnorm(n = N, rep(0, P), CovMatrix)
      
      E <- matrix(0,N,P)
      for (m in 1:(P/2)){
        for (n in 2:N){
          phi <- 0.5 + rnorm(1, 0, 0.05)
          E[n,m] <- phi*E[(n-1),m] + u[n,m]
        }
      }
      for (m in (P/2+1):P){
        E[,m] <- u[,m]
      }
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP6_B1.RDS"))
}

stopCluster(cl)











############################ DGP4_A2 ############################

# register cluster and call the main.R 
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

# generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #
N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    r<-NTR[1]
    P<-NTR[2]
    N<-NTR[3]
    
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT    ####
      
      Gamma_E <- diag(1,P,P) # pos-def sym mtx specifying the cov mtx of the variables
      E <- mvrnorm(n = N, rep(0, P), Gamma_E)
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP4_A2.RDS"))
}

stopCluster(cl)

############################ DGP5_A2 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

# generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_e <- 4 + rnorm(P, 0, 0.5) 
      sigma_e[seq(from = 1, to = P, by = 2)] <- 1
      Gamma_E <- diag(sigma_e[1:P],P,P) # pos-def sym mtx specifying the cov mtx of the variables
      E <- mvrnorm(n = N, rep(0, P), Gamma_E)
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP5_A2.RDS"))
}

stopCluster(cl)

############################ DGP6_A2 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      
      sigma_e <- 4 + rnorm(P, 0, 0.5) 
      sigma_e[seq(from = 1, to = P, by = 2)] <- 1
       
      
      CovMatrix <- matrix(0,P,P)
      for (a in 1:P){
        for (b in 1:P){
          CovMatrix[a,b] = 0.7^(abs(a-b))*sqrt(sigma_e[a])* sqrt(sigma_e[b]) 
          CovMatrix[b,a] = CovMatrix[a,b]
        }
      }
      
      E <- mvrnorm(n = N, rep(0, P), CovMatrix)
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP6_A2.RDS"))
}

stopCluster(cl)


############################ DGP4_B2 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4

iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    r<-NTR[1]
    P<-NTR[2]
    N<-NTR[3]
    
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      CovMatrix <- diag(1,P,P)
      u <- mvrnorm(n = N, rep(0, P), CovMatrix) 
      E <- matrix(0,N,P)
      for (m in 1:(P/2)){
        for (n in 2:N){
          phi <- 0.5 + rnorm(1, 0, 0.05)
          E[n,m] <- phi*E[(n-1),m] + u[n,m]
        }
      }
      for (m in (P/2+1):P){
        E[,m] <- u[,m]
      }
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP4_B2.RDS"))
}

stopCluster(cl)

############################ DGP5_B2 ############################
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4

iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_u <- 4 + rnorm(P, 0, 0.5) 
      sigma_u[seq(from = 1, to = P, by = 2)] <- 1 
      CovMatrix <- diag(sigma_u[1:P],P,P)
      
      u <- mvrnorm(n = N, rep(0, P), CovMatrix) # normal (0,1)
      E <- matrix(0,N,P)
      for (m in 1:(P/2)){
        for (n in 2:N){
          phi <- 0.5 + rnorm(1, 0, 0.05)
          E[n,m] <- phi*E[(n-1),m] + u[n,m]
        }
      }
      for (m in (P/2+1):P){
        E[,m] <- u[,m]
      }
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP5_B2.RDS"))
}

stopCluster(cl)

############################ DGP6_B2 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4

iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_u <- 4 + rnorm(P, 0, 0.5) 
      sigma_u[seq(from = 1, to = P, by = 2)] <- 1
      
      CovMatrix <- matrix(0,P,P)
      for (a in 1:P){
        for (b in 1:P){
          CovMatrix[a,b] = 0.7^(abs(a-b))*sqrt(sigma_u[a])* sqrt(sigma_u[b]) 
          CovMatrix[b,a] = CovMatrix[a,b]
        }
      }
      
      u <- mvrnorm(n = N, rep(0, P), CovMatrix)
      
      E <- matrix(0,N,P)
      for (m in 1:(P/2)){
        for (n in 2:N){
          phi <- 0.5 + rnorm(1, 0, 0.05)
          E[n,m] <- phi*E[(n-1),m] + u[n,m]
        }
      }
      for (m in (P/2+1):P){
        E[,m] <- u[,m]
      }
      
      ####            FACTORS            ####
      
      F1=arima.sim(model=list(ar=a1), n=N)
      F2=arima.sim(model=list(ar=a2), n=N)
      F3=arima.sim(model=list(ar=a3), n=N)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP6_B2.RDS"))
}

stopCluster(cl)



############################ DGP4_A3 ############################

# register cluster and call the main.R 
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

# generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #
N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    r<-NTR[1]
    P<-NTR[2]
    N<-NTR[3]
    
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT    ####
      
      Gamma_E <- diag(1,P,P) # pos-def sym mtx specifying the cov mtx of the variables
      E <- mvrnorm(n = N, rep(0, P), Gamma_E)
      
      ####            FACTORS            ####
      
      F1=arima.sim(list(order = c(1,0,0), ar = a1), n = N, rand.gen= rnorm, sd = 0.5)
      F2=arima.sim(list(order = c(1,0,0), ar = a2), n = N, rand.gen= rnorm, sd = 0.5)
      F3=arima.sim(list(order = c(1,0,0), ar = a3), n = N, rand.gen= rnorm, sd = 0.5)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP4_A3.RDS"))
}

stopCluster(cl)

############################ DGP5_A3 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

# generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_e <- 4 + rnorm(P, 0, 0.5) 
      sigma_e[seq(from = 1, to = P, by = 2)] <- 1
      Gamma_E <- diag(sigma_e[1:P],P,P) # pos-def sym mtx specifying the cov mtx of the variables
      E <- mvrnorm(n = N, rep(0, P), Gamma_E)
      
      ####            FACTORS            ####
      
      F1=arima.sim(list(order = c(1,0,0), ar = a1), n = N, rand.gen= rnorm, sd = 0.5)
      F2=arima.sim(list(order = c(1,0,0), ar = a2), n = N, rand.gen= rnorm, sd = 0.5)
      F3=arima.sim(list(order = c(1,0,0), ar = a3), n = N, rand.gen= rnorm, sd = 0.5)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP5_A3.RDS"))
}

stopCluster(cl)

############################ DGP6_A3 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      
      sigma_e <- 4 + rnorm(P, 0, 0.5) 
      sigma_e[seq(from = 1, to = P, by = 2)] <- 1
       
      
      CovMatrix <- matrix(0,P,P)
      for (a in 1:P){
        for (b in 1:P){
          CovMatrix[a,b] = 0.7^(abs(a-b))*sqrt(sigma_e[a])* sqrt(sigma_e[b]) 
          CovMatrix[b,a] = CovMatrix[a,b]
        }
      }
      
      E <- mvrnorm(n = N, rep(0, P), CovMatrix)
      
      ####            FACTORS            ####
      
      F1=arima.sim(list(order = c(1,0,0), ar = a1), n = N, rand.gen= rnorm, sd = 0.5)
      F2=arima.sim(list(order = c(1,0,0), ar = a2), n = N, rand.gen= rnorm, sd = 0.5)
      F3=arima.sim(list(order = c(1,0,0), ar = a3), n = N, rand.gen= rnorm, sd = 0.5)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP6_A3.RDS"))
}

stopCluster(cl)


############################ DGP4_B3 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    r<-NTR[1]
    P<-NTR[2]
    N<-NTR[3]
    
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      CovMatrix <- diag(1,P,P)
      u <- mvrnorm(n = N, rep(0, P), CovMatrix) 
      E <- matrix(0,N,P)
      for (m in 1:(P/2)){
        for (n in 2:N){
          phi <- 0.5 + rnorm(1, 0, 0.05)
          E[n,m] <- phi*E[(n-1),m] + u[n,m]
        }
      }
      for (m in (P/2+1):P){
        E[,m] <- u[,m]
      }
      
      ####            FACTORS            ####
      
      F1=arima.sim(list(order = c(1,0,0), ar = a1), n = N, rand.gen= rnorm, sd = 0.5)
      F2=arima.sim(list(order = c(1,0,0), ar = a2), n = N, rand.gen= rnorm, sd = 0.5)
      F3=arima.sim(list(order = c(1,0,0), ar = a3), n = N, rand.gen= rnorm, sd = 0.5)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP4_B3.RDS"))
}

stopCluster(cl)

############################ DGP5_B3 ############################
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_u <- 4 + rnorm(P, 0, 0.5) 
      sigma_u[seq(from = 1, to = P, by = 2)] <- 1 
      CovMatrix <- diag(sigma_u[1:P],P,P)
      
      u <- mvrnorm(n = N, rep(0, P), CovMatrix) # normal (0,1)
      E <- matrix(0,N,P)
      for (m in 1:(P/2)){
        for (n in 2:N){
          phi <- 0.5 + rnorm(1, 0, 0.05)
          E[n,m] <- phi*E[(n-1),m] + u[n,m]
        }
      }
      for (m in (P/2+1):P){
        E[,m] <- u[,m]
      }
      
      ####            FACTORS            ####
      
      F1=arima.sim(list(order = c(1,0,0), ar = a1), n = N, rand.gen= rnorm, sd = 0.5)
      F2=arima.sim(list(order = c(1,0,0), ar = a2), n = N, rand.gen= rnorm, sd = 0.5)
      F3=arima.sim(list(order = c(1,0,0), ar = a3), n = N, rand.gen= rnorm, sd = 0.5)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP5_B3.RDS"))
}

stopCluster(cl)

############################ DGP6_B3 ############################

cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)
clusterCall(cl, function() source("functions.R"))

#generate a log file for having track of the iteratations
writeLines(c("Starting"), "log.txt")

#           PARAMETERS         #

N <- c(125,250,500,1250) # number of time observations
R <- r <- 3 # number of static factors (q=1 dynamic factor)
K0 <- 2 # lags in L estimator
a1 <- 0.6 # AR coefficient of factor process
a2 <- 0.5
a3 <- 0.4
iteration <- 200

Pvector <- c(50,100,200,300)
for (P in Pvector){
  sample<-expand.grid(R,P,N) # posibles escenarios
  
  
  resultsParallel <- foreach(i = 1:nrow(sample), .verbose = TRUE, .inorder = TRUE) %dopar% {
    sink("log.txt", append = TRUE)
    cat(paste("Iteration",i,"\n"))
    
    NTR<-as.numeric(sample[i,])
    P<-NTR[2]
    N<-NTR[3]
    r<-NTR[1]
    
     media <- matrix(0,nrow = iteration, ncol = 6)
     
    for (j in 1:iteration){
      
      ####   IDIOSYNCRATIC COMPONENT     ####
      sigma_u <- 4 + rnorm(P, 0, 0.5) 
      sigma_u[seq(from = 1, to = P, by = 2)] <- 1
      
      CovMatrix <- matrix(0,P,P)
      for (a in 1:P){
        for (b in 1:P){
          CovMatrix[a,b] = 0.7^(abs(a-b))*sqrt(sigma_u[a])* sqrt(sigma_u[b])
          CovMatrix[b,a] = CovMatrix[a,b]
        }
      }
      
      u <- mvrnorm(n = N, rep(0, P), CovMatrix)
      
      E <- matrix(0,N,P)
      for (m in 1:(P/2)){
        for (n in 2:N){
          phi <- 0.5 + rnorm(1, 0, 0.05)
          E[n,m] <- phi*E[(n-1),m] + u[n,m]
        }
      }
      for (m in (P/2+1):P){
        E[,m] <- u[,m]
      }
      
      ####            FACTORS            ####
      
      F1=arima.sim(list(order = c(1,0,0), ar = a1), n = N, rand.gen= rnorm, sd = 0.5)
      F2=arima.sim(list(order = c(1,0,0), ar = a2), n = N, rand.gen= rnorm, sd = 0.5)
      F3=arima.sim(list(order = c(1,0,0), ar = a3), n = N, rand.gen= rnorm, sd = 0.5)
      F <- cbind(F1,F2,F3)
      
      ####            LOADINGS           ####
      
      L=matrix(runif(P*r, -0.5, 0.5), nrow=P, ncol=r)
      
      ####       DATA MATRIX (T,N)       ####
      
      X <- F%*%t(L) + E # data matrix
      
      X <- scale(X, center = TRUE, scale = FALSE)
      
      
      #          ESTIMATION           #
      
      theta<-rep(0,6)
      kmax <- 10
      #### Ahn and Horestein ####
      
      theta[4] <- AhnHorensteintest(X,kmax)$r1
      
      #### Lagged K0=2 ACOV estimation ####
      
      theta[3] <- LamYaotest(X,2,kmax)$r1
      
      #### Lagged including K0=0 ACOR estimation  ####
      
      theta[6] <-  CaroPenatestCOR(X,2,kmax)$r1 # weights
      
      #### Fan et al ###
      
      theta[5] <- Fanetal(X, kmax)$r1
      
      #### BN(2002) PC3 ####
      
      theta[1] <- BNtest(X, kmax)$r1
      
      #### Onatski(2005) ON2 ####
      
      theta[2] <- ONtest(X, kmax)$r1
      
      media[j,] <- theta
      
      print(c(i,j))
    }
     
     scenario <- paste(paste0(N,"N"), paste0(P, "P"),sep = "_")
     sink()
     
     bn <- sum(media[,1] == 3)
     on <- sum(media[,2] == 3)
     ly <- sum(media[,3] == 3)
     ah <- sum(media[,4] == 3)
     fan <- sum(media[,5] == 3)
     cp <- sum(media[,6] == 3)
     meantheta <- c(bn/iteration, on/iteration, ly/iteration, ah/iteration, fan/iteration, cp/iteration)
    return(list(meantheta=meantheta,scenario=scenario, out=media))
  }
  saveRDS(resultsParallel, file = paste0("P", P, "_resultsDGP6_B3.RDS"))
}

stopCluster(cl)

##### END TIME ####

end_time <- Sys.time()
end_time - start_time
