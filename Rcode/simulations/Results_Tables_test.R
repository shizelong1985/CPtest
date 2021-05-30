# Program that summarizes results from SimulationTest_TwoFactors.R and SimulationTest_ThreeFactors.R
# Latex Table Generation
# P : number of cross-section variables
# N : number of time series observations
#
# results by columns are the Test for the number of factors
# "BN", "ON", "LY","AH","ACT", "CP"


library(xtable)

############# Scenarios 1, 2, 3 #############
setwd("set your path here/simulations/ResultsTwoFactors")

files <- list.files()
comb_model <- paste0(rep(1:3, each = 6), rep(c("_A1", "_B1", "_A2", "_B2", "_A3", "_B3")))
comb_P <- c("50_", "100_", "200_", "300_")

scenarios <- seq(1,length(comb_model), 2)
for(comb_i in scenarios){
  files_model <- files[grepl(comb_model[comb_i], files)]
  files_modelB <- files[grepl(comb_model[(comb_i+1)], files)]
  CanCor <- matrix(NA,21,14)
  
  for(P in comb_P){
    results <- readRDS(files_model[grepl(P, files_model)])
    resultsB <- readRDS(files_modelB[grepl(P, files_modelB)])
    if(P == comb_P[1]){
      for (i in 1:4){CanCor[(i+1), 1:6] <- results[[i]]$meantheta}
      for (i in 1:4){CanCor[(i+1), 8:13] <- resultsB[[i]]$meantheta}
      
    }
    if(P == comb_P[2]){
      for (i in 1:4){CanCor[(i+6), 1:6] <- results[[i]]$meantheta}
      for (i in 1:4){CanCor[(i+6), 8:13] <- resultsB[[i]]$meantheta}
    }
    if(P == comb_P[3]){
      for (i in 1:4){CanCor[(i+11), 1:6] <- results[[i]]$meantheta}
      for (i in 1:4){CanCor[(i+11), 8:13] <- resultsB[[i]]$meantheta}
    }
    if(P == comb_P[4]){
      for (i in 1:4){CanCor[(i+16), 1:6] <- results[[i]]$meantheta}
      for (i in 1:4){CanCor[(i+16), 8:13] <- resultsB[[i]]$meantheta}
    }
  }
  
  cols <- c(1:6,8:13)
  rows <- c(2:5, 7:10, 12:15, 17:20)
  meanCol <- colMeans(CanCor[,cols], na.rm = TRUE)
  meanRow <- rowMeans(CanCor[rows,], na.rm = TRUE)
  CanCor[21,cols] <- round(meanCol, digits = 2)
  CanCor[rows,14] <- round(meanRow, digits = 2)
  CanCor[,seq(1,14,1)] <- CanCor[,c(5,6,2,1,4,3,7,12,13,9,8,11,10,14)]
  colnames(CanCor) <- c("BN", "ON", "LY","AH","ACT", "CP", " ","BN", "ON", "LY","AH","ACT", "CP", "Mean")  
  col1 <- c("N=50","T=125", "T=250", "T=500", "T=1250","N=100","T=125", "T=250", "T=500", "T=1250","N=200","T=125", "T=250", "T=500", "T=1250","N=300","T=125", "T=250", "T=500", "T=1250", "Mean")
  CanCor <- round(CanCor, digits = 2)
  angela <- cbind(col1,CanCor)
  tab <- xtable(angela, caption = paste0("tab_", comb_model[comb_i]))
  assign(paste0("tab_", comb_model[comb_i]), tab)
  
  print.xtable(tab, include.rownames = getOption("xtable.include.rownames", FALSE))
}






############# Scenarios 4, 5, 6 #############
setwd("set your path here/simulations/ResultsThreeFactors")

files <- list.files()
comb_model <- paste0(rep(4:6, each = 6), rep(c("_A1", "_B1", "_A2", "_B2", "_A3", "_B3")))
comb_P <- c("50_", "100_", "200_", "300_")

scenarios <- seq(1,length(comb_model), 2)
for(comb_i in scenarios){
  files_model <- files[grepl(comb_model[comb_i], files)]
  files_modelB <- files[grepl(comb_model[(comb_i+1)], files)]
  CanCor <- matrix(NA,21,14)
  
  for(P in comb_P){
    results <- readRDS(files_model[grepl(P, files_model)])
    resultsB <- readRDS(files_modelB[grepl(P, files_modelB)])
    if(P == comb_P[1]){
      for (i in 1:4){CanCor[(i+1), 1:6] <- results[[i]]$meantheta}
      for (i in 1:4){CanCor[(i+1), 8:13] <- resultsB[[i]]$meantheta}
      
    }
    if(P == comb_P[2]){
      for (i in 1:4){CanCor[(i+6), 1:6] <- results[[i]]$meantheta}
      for (i in 1:4){CanCor[(i+6), 8:13] <- resultsB[[i]]$meantheta}
    }
    if(P == comb_P[3]){
      for (i in 1:4){CanCor[(i+11), 1:6] <- results[[i]]$meantheta}
      for (i in 1:4){CanCor[(i+11), 8:13] <- resultsB[[i]]$meantheta}
    }
    if(P == comb_P[4]){
      for (i in 1:4){CanCor[(i+16), 1:6] <- results[[i]]$meantheta}
      for (i in 1:4){CanCor[(i+16), 8:13] <- resultsB[[i]]$meantheta}
    }
  }
  
  cols <- c(1:6,8:13)
  rows <- c(2:5, 7:10, 12:15, 17:20)
  meanCol <- colMeans(CanCor[,cols], na.rm = TRUE)
  meanRow <- rowMeans(CanCor[rows,], na.rm = TRUE)
  CanCor[21,cols] <- round(meanCol, digits = 2)
  CanCor[rows,14] <- round(meanRow, digits = 2)
  CanCor[,seq(1,14,1)] <- CanCor[,c(5,6,2,1,4,3,7,12,13,9,8,11,10,14)]
  colnames(CanCor) <- c("BN", "ON", "LY","AH","ACT", "CP", " ","BN", "ON", "LY","AH","ACT", "CP", "Mean")  
  col1 <- c("N=50","T=125", "T=250", "T=500", "T=1250","N=100","T=125", "T=250", "T=500", "T=1250","N=200","T=125", "T=250", "T=500", "T=1250","N=300","T=125", "T=250", "T=500", "T=1250", "Mean")
  CanCor <- round(CanCor, digits = 2)
  angela <- cbind(col1,CanCor)
  tab <- xtable(angela, caption = paste0("tab_", comb_model[comb_i]))
  assign(paste0("tab_", comb_model[comb_i]), tab)
  
  print.xtable(tab, include.rownames = getOption("xtable.include.rownames", FALSE))
}
