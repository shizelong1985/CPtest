#### OECD STATISTICS ####
setwd("set your path here")
library(readxl)
library(ggplot2)
library(plotly)
library(patchwork)
library(lmtest)

# Functions needed
source("functions.R")

GDP_CON_INV <- read_excel("set your path here/OECDexample/GDP_CON_INV.xlsx", 
                          col_names = FALSE)

quart <- c("Q1", "Q2", "Q3", "Q4")
years <- c("1999", "2000", "2001", "2002", "2003",
           "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020")
dates <- paste0(rep(years,each=4),quart)

countries <- c("Australia", "Austria", "Belgium", "Canada", "Chile","Czech Republic", "Denmark", "Estonia",      
"Finland","France","Germany","Hungary", "Iceland","Ireland","Israel","Italy","Japan","Korea","Latvia","Lithuania",
"Luxembourg","Mexico","Netherlands" ,"New Zealand", "Norway","Poland","Portugal","Slovak Republic",
"Slovenia","Spain","Sweden","Switzerland", "Turkey","United Kingdom","United States")
 
data <- as.matrix(GDP_CON_INV[,2:106], nrow = 88, ncol = 105)
colnames(data) <- paste0(rep(countries, each =3), c("_gdp", "_cons", "_inv"))

Ydata <- ts(data,start=1999, frequency = 4) # frequency: quarterly
logY <- log(Ydata)
difY <- diff(logY, lag = 1, differences = 1)
Y <- scale(difY, center = TRUE, scale = FALSE)

P <- ncol(Y)
N <- nrow(Y)
K0 <- 2 # max number of lags in the cummulative sum
kmax <- 10 # max number of factors rmax

dates <- dates[2:88]

names <- colnames(Y)
labs <- rep("",105)
gdps <- seq(1,105,3)
for (i in gdps){
  labs[i]<-names[i]
}

#### Ahn and Horestein ####
AH <- AhnHorensteintest(Y,kmax)
AH$r1
ratio <- AH$ratio
plot_data <- data.frame(eigvalues = seq(1,10), ratio = ratio)
ratio1_AH <- ggplot(plot_data, aes(x = eigvalues, y = ratio)) + 
  geom_line() +  geom_point() + 
  scale_x_continuous(breaks = seq(1,10, by = 3),
                     labels = as.character(seq(1,10, by = 3))) +  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5), 
        legend.title = element_blank(),axis.title.x = element_blank(), 
        axis.title.y = element_blank(),legend.position = "none")

ratio1_AH

Ev <- eigen(AH$Wy)
# loadings
Lhat <- Ev$vectors[,1]
# loadings plot
xtick<-seq(1, 105, by = 1)
#axis(side=1, at = xtick, labels = colnames(Y))

plot_data <- data.frame(Loadings = -Lhat, name = colnames(Y))
plot_data$name <- factor(plot_data$name, levels = plot_data$name)

bar_plot <- ggplot(plot_data) + geom_bar(aes(x = name, y = Loadings), stat = "identity") +  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(),axis.title.y = element_blank())+
  scale_x_discrete(labels=c(labs)) 

ggplotly(bar_plot)

# factors
fhat <- Y%*%Lhat
# factor plot 
plot_data <- data.frame(time = rep(dates, ncol(fhat)), factor = c(-fhat), color = as.factor(rep(c("factor.1"), each = nrow(fhat))))
plot <- ggplot(plot_data, aes(x = time, y = factor, group = color, color = color)) + geom_line(color='black') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_classic() 
plot <- plot + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) 
plot <- plot + scale_x_discrete(labels=c("1999", "", "","2000", "", "", "","2001", "", "", "","2002","", "", "", "","2003", "", "", "","2004", "", "", "", "2005", "","","","2006", "", "", "", "2007", "", "", "", "2008", "", "", "", "2009", "","","", "2010", "", "", "","2011", "", "", "", "2012", "", "", "", "2013", "", "", "","2014", "", "", "","2015", "", "", "","2016", "", "", "","2017", "", "", "","2018", "", "", "","2019","","","", "2020","", "", ""))
plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplotly(plot) 

eigenvaluesM<-Ev$values 

pvarM<-eigenvaluesM*100/sum(eigenvaluesM) #proportion of variability explained by the components
cumvarM<-cumsum(pvarM)#proportion of accumulated variability for the first r PCs
eig.data<-data.frame(eig=round(eigenvaluesM, digits = 3), variance=round(pvarM, digits = 3), cumvariance=round(cumvarM, digits = 3))
eig.data[1:10,]

barplot(eig.data[1:12,2], names.arg=1:12, main="variances", xlab="PC",ylab="% of variances", col="steelblue")
lines(x=1:12, eig.data[1:12,2],type="b", pch=19, col="red")

#### Caro and Peña 
## STEP 1 ####
CPW <- CaroPenatestCOR(Y, K0, kmax)
CPW$r1
ratio <- CPW$ratio
plot_data <- data.frame(eigvalues = seq(1,10),  ratio = (ratio))
ratio1_CP <- ggplot(plot_data, aes(x = eigvalues, y = ratio)) + 
  geom_line() +  geom_point() + theme_classic() +
  scale_x_continuous(breaks = seq(1,10, by = 3),
                     labels = as.character(seq(1,10, by = 3)))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5), 
        legend.title = element_blank(),axis.title.x = element_blank(), 
        axis.title.y = element_blank(),legend.position = "none")

ratio1_CP

Ev <- eigen(CPW$Wy)
# loadings
Lhat <- Ev$vectors[,1]
# loadings plot
xtick<-seq(1, 105, by = 1)
axis(side=1, at = xtick, labels = colnames(Y))

plot_data <- data.frame(Loadings = -Lhat, name = colnames(Y))
plot_data$name <- factor(plot_data$name, levels = plot_data$name)

bar_plot <- ggplot(plot_data) + geom_bar(aes(x = name, y = Loadings), stat = "identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_x_discrete(labels=c(labs))

ggplotly(bar_plot)

# factors
fhat <- Y%*%Lhat
# factor plot 
plot_data <- data.frame(time = rep(dates, ncol(fhat)), factor = c(-fhat), color = as.factor(rep(c("factor.1"), each = nrow(fhat))))
plot <- ggplot(plot_data, aes(x = time, y = factor, group = color, color = color)) + geom_line(color='black') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_classic()
plot <- plot + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) 
plot <- plot + scale_x_discrete(labels=c("1999", "", "","2000", "", "", "","2001", "", "", "","2002","", "", "", "","2003", "", "", "","2004", "", "", "", "2005", "","","","2006", "", "", "", "2007", "", "", "", "2008", "", "", "", "2009", "","","", "2010", "", "", "","2011", "", "", "", "2012", "", "", "", "2013", "", "", "","2014", "", "", "","2015", "", "", "","2016", "", "", "","2017", "", "", "","2018", "", "", "","2019","","","","2020","", "", ""))
plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplotly(plot) 

eigenvaluesM<-Ev$values 

pvarM<-eigenvaluesM*100/sum(eigenvaluesM) #proportion of variability explained by the components
cumvarM<-cumsum(pvarM)#proportion of accumulated variability for the first r PCs
eig.data<-data.frame(eig=round(eigenvaluesM, digits = 3), variance=round(pvarM, digits = 3), cumvariance=round(cumvarM, digits = 3))
eig.data[1:10,]

## STEP 2 ####
Y2 <- Y - fhat%*%Lhat 

fit <- acf(Y2, lag.max = 1, type="correlation", plot=FALSE)
cormtx <- fit$acf[1,,]
hist(cormtx, main = "Cross correlation coefficients", xlab = "", ylab="", xlim=c(-1,1))

coef <- rep(0,105)
for (i in 1:105){
  value <- arima(Y2[,i], order = c(1, 0, 0))
  coef[i] <- value$coef[1]
}
 hist(coef, main = "Autoregressive coefficients", xlab = "", ylab="", xlim=c(-0.75,0.6))  

#### Fan et al. 2020 ####
F1 <- Fanetal(Y,kmax)
F1$r1
Ev <- eigen(F1$Wy)

eigenvaluesM<-Ev$values 

pvarM<-eigenvaluesM*100/sum(eigenvaluesM) #proportion of variability explained by the components
cumvarM<-cumsum(pvarM)#proportion of accumulated variability for the first r PCs
eig.data<-data.frame(eig=round(eigenvaluesM, digits = 3), variance=round(pvarM, digits = 3), cumvariance=round(cumvarM, digits = 3))
eig.data[1:10,]

# loadings
Lhat <- Ev$vectors[,1:3]

# factors
fhat <- Y%*%Lhat

# factor plot 
plot_data <- data.frame(time = rep(dates, ncol(fhat)), factor = c(-fhat), color = as.factor(rep(c("factor.1","factor.2","factor.3"), each = nrow(fhat))))
plot <- ggplot(plot_data, aes(x = time, y = factor, group = color, color = color)) + geom_line()+ theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
plot <- plot + theme(legend.title = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank())
plot <- plot + theme(legend.position = c(0.12, 0.2))
plot <- plot + scale_x_discrete(labels=c("1999", "", "","2000", "", "", "","2001", "", "", "","2002","", "", "", "","2003", "", "", "","2004", "", "", "", "2005", "","","","2006", "", "", "", "2007", "", "", "", "2008", "", "", "", "2009", "","","", "2010", "", "", "","2011", "", "", "", "2012", "", "", "", "2013", "", "", "","2014", "", "", "","2015", "", "", "","2016", "", "", "","2017", "", "", "","2018", "", "", "","2019","","","","2020","", "", ""))

ggplotly(plot) 

# only factor 2 and 3
plot_data <- data.frame(time = rep(dates, ncol(fhat[,2:3])), factor = c(-fhat[,2:3]), color = as.factor(rep(c("factor.2","factor.3"), each = nrow(fhat))))
plot <- ggplot(plot_data, aes(x = time, y = factor, group = color)) + geom_line(aes(linetype=color))+ theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
plot <- plot + theme(legend.title = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank())
plot <- plot + theme(legend.position = c(0.12, 0.2))
plot <- plot + scale_x_discrete(labels=c("1999", "", "","2000", "", "", "","2001", "", "", "","2002","", "", "", "","2003", "", "", "","2004", "", "", "", "2005", "","","","2006", "", "", "", "2007", "", "", "", "2008", "", "", "", "2009", "","","", "2010", "", "", "","2011", "", "", "", "2012", "", "", "", "2013", "", "", "","2014", "", "", "","2015", "", "", "","2016", "", "", "","2017", "", "", "","2018", "", "", "","2019","","","","2020","", "", ""))

ggplotly(plot) 

