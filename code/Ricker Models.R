#Original code was from Steve Fleischman and Xinxian Zhang. code was modified
#and converted to JAGS by Sara Miller on 28 June 2017 and modified and
#maintained by Rich Brenner.

#load----
library(arm)
library(lmtest)
library(rjags)
library(R2OpenBUGS)
library(gdata)
library(tidyverse)

#data----
brood <- read_csv("data/Chilkoot_Sock.csv")

# cleanup data
brood %>% 
  dplyr::select(S=Spawn, R=Recruit) %>% 
  mutate(lnRS = log(R/S), n()) -> sr

n <- nrow(sr) #calculates the number of years of data.
sr.data <- list(n = n, S = sr$S, lnRS = sr$lnRS)

#analysis----
#check AR(1)
mylm <- lm(lnRS ~ S, sr)
dwtest(mylm)
pacf(residuals(mylm))


# JAGS CODE NOT R CODE
# Ricker w/o autocorrelation ----
# Ricker model for stock-recruitment analysis 
# Created by Steve Fleischman  

Ricker = function(){
  
  lnalpha ~ dnorm(0,1.0E-6)%_%T(0,3)#uninformative;*
  beta ~ dnorm(0,1.0E-6)%_%T(0,)   #uninformative; normal distrib; constrained to be >0 *
  phi <- 0                #this model does not account for autocorrelation so phi is not used (thus, phi =0)
  sigma.white ~ dnorm(0,1.0E-6)%_%T(0,)  
  resid.red.0 ~ dnorm(0,tau.red)
  
  for(y in 1:n) {lnRS[y] ~ dnorm(mean2.lnRS[y],tau.white) }
  
  mean2.lnRS[1]  <- mean1.lnRS[1] + phi * resid.red.0  
  for (y in 2:n) { mean2.lnRS[y] <- mean1.lnRS[y] + phi * resid.red[y-1] }  #NO autocorrelation model
  
  for(y in 1:n) {  mean1.lnRS[y] <- lnalpha - beta * S[y]  }
  for(y in 1:n) {  resid.red[y]  <- lnRS[y] - mean1.lnRS[y]  }
  for(y in 1:n) {  resid.white[y] <- lnRS[y] - mean2.lnRS[y]  }
  
  tau.white <- 1 / sigma.white / sigma.white        
  tau.red <- tau.white * (1-phi*phi)
  sigma.red <- 1 / sqrt(tau.red) 
  
  lnalpha.c <- lnalpha + (sigma.red * sigma.red / 2)  #adjust for calculating means of R.msy, S.msy etc.
  alpha <- exp(lnalpha)  #exponentiate to solve for alpha
  S.max <- 1 / beta
  S.eq <- S.max * lnalpha.c 
  S.msy <- S.eq * (0.5 - 0.07*lnalpha.c) #Hilborn approximation of Smsy...could use Scheuerell solution too....
  U.msy <- lnalpha.c * (0.5 - 0.07*lnalpha.c)
  R.msy <- S.msy * exp(lnalpha.c - beta * S.msy)  #Solves for recruits at Smsy
  }


#write the non-AR model to a text file to be called by WinBUGS
model_file_loc=paste("code/Chilkoot_Sockeye.txt", sep="")
write.model(Ricker, paste("code/Chilkoot_Sockeye.txt", sep=""))



# Ricker model WITH autocorrelation  ----
AR = function(){
  #lnalpha ~ dunif(0, 10) 
  lnalpha ~ dnorm(0,1.0E-6)%_%T(0,3)#uninformative;*
  beta ~ dnorm(0,1.0E-6)%_%T(0,)   #uninformative; normal distrib; constrained to be >0 *
  #beta ~ dunif(0, 10)  
  phi ~ dnorm(0,1.0E-6)%_%T(-0.98,0.98)#uninformative; btw -1 and 1
  sigma.white ~ dnorm(0,1.0E-6)%_%T(0,)   #uninformative; normal distrib; constrained to be >0 *
  resid.red.0 ~ dnorm(0,tau.red)
  
  for(y in 1:n) {lnRS[y] ~ dnorm(mean2.lnRS[y],tau.white) }
  
  mean2.lnRS[1] <- mean1.lnRS[1] + phi * resid.red.0  
  for (y in 2:n) { mean2.lnRS[y] <- mean1.lnRS[y] + phi * resid.red[y-1] }   #AR1
  
  for(y in 1:n) {  mean1.lnRS[y] <- lnalpha - beta * S[y]  } #Ricker model
  for(y in 1:n) {  resid.red[y]     <- lnRS[y] - mean1.lnRS[y]  }
  for(y in 1:n) {  resid.white[y] <- lnRS[y] - mean2.lnRS[y]  }
  
  tau.white <- 1 / sigma.white / sigma.white        
  tau.red <- tau.white * (1-phi*phi)
  sigma.red <- 1 / sqrt(tau.red)
  
  lnalpha.c <- lnalpha + (sigma.red * sigma.red / 2)  #adjust for calculating means of R.msy, S.msy etc.
  alpha <- exp(lnalpha) #exponentiate to solve for alpha 
  S.max <- 1 / beta
  S.eq <- S.max * lnalpha.c 
  S.msy <- S.eq * (0.5 - 0.07*lnalpha.c)
  U.msy <- lnalpha.c * (0.5 - 0.07*lnalpha.c)
  R.msy <- S.msy * exp(lnalpha.c - beta * S.msy)
  }


#write the AR model to a text file to be called by WinBUGS or JAGS
model_file_loc=paste("code/Chilkoot_Sockeye_AR.txt", sep="")
write.model(AR, paste("code/Chilkoot_Sockeye_AR.txt", sep=""))

# NOW BACK TO R CODE

# Results ----
# 100000 iterations, 3 chains, 10000 burn-in period, thin by 100
# Run the Ricker model that does NOT have autocorrelation

inits1 <- list(lnalpha=1.5, beta=0.0005, sigma.white=0.7, resid.red.0= 0)
inits2 <- list(lnalpha=2.0, beta=0.0010, sigma.white=0.5, resid.red.0=-1)
inits3 <- list(lnalpha=2.5, beta=0.0020, sigma.white=0.3, resid.red.0= 1)
inits <- list(inits1, inits2, inits3)

parameters <- c("lnalpha.c","lnalpha", "alpha","beta", "sigma.red", "phi", "S.max", "S.eq", "S.msy", "U.msy", "R.msy",  "resid.red.0", "sigma.white")
ptm = proc.time()
jmod <- jags.model(file='code/Chilkoot_Sockeye.txt', data=sr.data, n.chains=3, inits=inits, n.adapt=1000) 
x <- update(jmod, n.iter=100000, by=100, progress.bar='text', DIC=T, n.burnin=10000) 
post <- coda.samples(jmod, parameters, n.iter=100000, thin=100, n.burnin=10000)
post.samp <- post

#Numerical summary of each parameter (mean, median, quantiles)
summary <- summary(post)                     
stats <- summary$statistics;  colnames(stats)
quants <- summary$quantiles;  colnames(quants)
statsquants <- cbind(stats,quants)
statsquants <- statsquants[,c(1,2,4,5,7,9)] #select statquant columns of interest: Mean, SD, Time-series SE,...

data.frame(statsquants) %>%  #converts statsquants to a dataframe
  tibble::rownames_to_column() %>%  #preserves the row names (otherwise, tibble  drops these)
  mutate(CV = SD/Mean) %>% #calculates the CV
  rename(Item=rowname, "Time series SE" = Time.series.SE, "2.5%"=X2.5.) -> statsquants #converting to a tibble messed
#the column names, need
#to rename these

write.csv(statsquants, file= paste("results/Ricker.csv") ) 

# Plots ----
# Density and time series plots
post.samp <- post
nvars <- dim(post.samp[[1]])[2]
nsamps <- dim(post.samp[[1]])[1]
int <- 25

pdf("figures/Ricker_profiles.pdf",height=6, width=8)

for(j in seq(1,nvars,int)){
  par(mfrow=c(5,4),mai=c(0.3,0.3,0.2,0.2))
  for(i in 0:(int-1)){
    mindat=min(c(post.samp[[1]][,i+j],post.samp[[1]][,i+j]))
    maxdat=max(c(post.samp[[1]][,i+j],post.samp[[1]][,i+j]))
    plot(density(post.samp[[1]][,i+j]),col='blue',main=colnames(post.samp[[1]])[i+j],xlim=c(mindat,maxdat))
    lines(density(post.samp[[2]][,i+j]),col='red')
    lines(density(post.samp[[3]][,i+j]),col='green')
    
    plot(as.numeric(post.samp[[1]][,i+j]),col='blue',main=colnames(post.samp[[1]])[i+j],ylim=c(mindat,maxdat),type='l')
    lines(as.numeric(post.samp[[2]][,i+j]),col='red')
    lines(density(post.samp[[3]][,i+j]),col='green')
  }}

dev.off()


# Gelman ----
gel <- as.data.frame(gelman.diag(post, multivariate=F)[[1]])
poor.threshold=1.10#values less than 1.2 are generally considered converged
poor <- gel[gel[,1] > poor.threshold, ]

write_csv(poor, "results/Ricker_Gelman.csv")     

# DIC ----
dic.pD  <-dic.samples(jmod,n.iter=100000, thin=100,"pD",  n.burnin=10000)
dic.popt<-dic.samples(jmod,n.iter=100000, thin=100,"popt",n.burnin=10000)
dev1 <- sum(dic.pD[[1]])
pD   <- sum(dic.pD[[2]])
dic.pD <- dev1 + pD
dic.pD.summary <- data.frame(dev1, pD, dic.pD)

write_csv(dic.pD.summary, "results/Ricker_DIC.csv")  

#Create coda samples for horsetail plots and probability plots
post2 <- coda.samples(jmod, c("lnalpha", "beta", "lnalpha.c"), n.iter=100000, thin=100,n.burnin=10000) 
x <- as.array(post2)
x <- data.frame(x)
coda <- x[,1:3]
coda <- rename.vars(coda, from=c("beta.1","lnalpha.1","lnalpha.c.1"), to=c("beta","lnalpha", "lnalpha.c"))

write_csv(coda, "results/Ricker_coda.csv")  # writes csv file

#Create coda samples for lambert calc
library(gsl)
post2a <- coda.samples(jmod, c("lnalpha", "beta", "lnalpha.c"), n.iter=100000, thin=100,n.burnin=10000) 
x <- as.array(post2a)
x <- data.frame(x)
coda1 <- x[,1:3]
coda2 <- x[,4:6]
coda3 <- x[,7:9]
coda1<- rename.vars(coda1, from=c("beta.1","lnalpha.1","lnalpha.c.1"), to=c("beta","lnalpha", "lnalpha.c"))
coda2<- rename.vars(coda2, from=c("beta.2","lnalpha.2","lnalpha.c.2"), to=c("beta","lnalpha", "lnalpha.c"))
coda3<- rename.vars(coda3, from=c("beta.3","lnalpha.3","lnalpha.c.3"), to=c("beta","lnalpha", "lnalpha.c"))
coda<-rbind(coda1,coda2,coda3)
coda$Smsy_lambert <- (1-lambert_W0(exp(1-coda$lnalpha.c)))/coda$beta 
coda$Umsy_lambert <- (1-lambert_W0(exp(1-coda$lnalpha.c))) 
coda<-as.data.frame(coda)
summary<-summary(coda) 
q1<-apply(coda,2,quantile,probs=c(0,0.025,0.5,0.975,1))#percentiles
write.csv(q1, file= paste("results/quantiles_lambert.csv") )    
write.csv(summary, file= paste("results/lambert.csv") ) 
write.csv(coda, file= paste("results/coda_lambert.csv") ) 


# Run Ricker WITH AR(1) ----  
# 100000 iterations, 3 chains, 10000 burn-in period, thin by 100
inits1 <- list(lnalpha=1.5, beta=0.0005, phi= 0.3, sigma.white=0.7, resid.red.0= 0)
inits2 <- list(lnalpha=2.0, beta=0.0010, phi=-0.1, sigma.white=0.5, resid.red.0=-1)
inits3 <- list(lnalpha=2.5, beta=0.0020, phi= 0.2, sigma.white=0.3, resid.red.0= 1)
inits <- list(inits1, inits2, inits3)

parameters <- c("sigma.white", "lnalpha.c","beta", "sigma.red","alpha","lnalpha","S.msy", "phi", "S.max", "S.eq",  "U.msy", "R.msy",  "resid.red.0")
jmod <- jags.model(file='code/Chilkoot_Sockeye_AR.txt', data=sr.data, n.chains=3, inits=inits, n.adapt=1000) 
x <- update(jmod, n.iter=100000, by=100, progress.bar='text', DIC=T, n.burnin=10000) 
post <- coda.samples(jmod, parameters, n.iter=100000, thin=100, n.burnin=10000)
post.samp <- post

# Numerical summary of each parameter (mean, median, quantiles)
summary <- summary(post)                     
stats <- summary$statistics;  colnames(stats)
quants <- summary$quantiles;  colnames(quants)
statsquants <- cbind(stats,quants)
statsquants <- statsquants[,c(1,2,4,5,7,9)] #select statquant columns of interest: Mean, SD, Time-series SE,...

data.frame(statsquants) %>%  #converts statsquants to a dataframe
  tibble::rownames_to_column() %>%  #preserves the row names (otherwise, tibble  drops these)
  mutate(CV = SD/Mean) %>% #calculates the CV
  rename(Item=rowname, "Time series SE" = Time.series.SE, "2.5%"=X2.5.) -> statsquants #converting to a tibble messed
                                                                                      #the column names, need
                                                                                      #to rename these

write.csv(statsquants, file= paste("results/Ricker_AR.csv") ) 

# Density and time series plots
post.samp <- post
nvars <- dim(post.samp[[1]])[2]
nsamps <- dim(post.samp[[1]])[1]
int <- 25
pdf("figures/Ricker_AR_profiles.pdf",height=6, width=8)
for(j in seq(1,nvars,int)){
  par(mfrow=c(5,4),mai=c(0.3,0.3,0.2,0.2))
  for(i in 0:(int-1)){
    mindat=min(c(post.samp[[1]][,i+j],post.samp[[1]][,i+j]))
    maxdat=max(c(post.samp[[1]][,i+j],post.samp[[1]][,i+j]))
    plot(density(post.samp[[1]][,i+j]),col='blue',main=colnames(post.samp[[1]])[i+j],xlim=c(mindat,maxdat))
    lines(density(post.samp[[2]][,i+j]),col='red')
    lines(density(post.samp[[3]][,i+j]),col='green')
    
    plot(as.numeric(post.samp[[1]][,i+j]),col='blue',main=colnames(post.samp[[1]])[i+j],ylim=c(mindat,maxdat),type='l')
    lines(as.numeric(post.samp[[2]][,i+j]),col='red')
    lines(density(post.samp[[3]][,i+j]),col='green')
  }}
dev.off()

# Gelman ----
gel <- as.data.frame(gelman.diag(post, multivariate=F)[[1]])
poor.threshold=1.10 #values less than 1.2 are generally considered converged
poor <- gel[gel[,1] > poor.threshold,]

write_csv(poor, "results/Ricker_AR_Gelman.csv")     

# DIC ----
dic.pD  <- dic.samples(jmod,n.iter=100000, thin=100,"pD",  n.burnin=10000)
dic.popt <- dic.samples(jmod,n.iter=100000, thin=100,"popt",n.burnin=10000)
dev1 <- sum(dic.pD[[1]])
pD   <- sum(dic.pD[[2]])
dic.pD <- dev1 + pD
dic.pD.summary <- data.frame(dev1, pD, dic.pD)
write_csv(dic.pD.summary, "results/Ricker_AR_DIC.csv")  

#Create coda samples for horsetail plots and probability plots for the AR model
post2 <- coda.samples(jmod, c("lnalpha", "beta", "lnalpha.c"), n.iter=100000, thin=100,n.burnin=10000) 
x <- as.array(post2)
x <- data.frame(x)
coda <- x[,1:3] 
coda <- rename.vars(coda, from=c("beta.1","lnalpha.1","lnalpha.c.1"), to=c("beta","lnalpha", "lnalpha.c"))
write_csv(coda, "results/Ricker_AR_coda.csv") # writes csv file

#Create coda samples for lambert calc
library(gsl)
post2a <- coda.samples(jmod, c("lnalpha", "beta", "lnalpha.c"), n.iter=100000, thin=100,n.burnin=10000) 
x <- as.array(post2a)
x <- data.frame(x)
coda1 <- x[,1:3]
coda2 <- x[,4:6]
coda3 <- x[,7:9]
coda1<- rename.vars(coda1, from=c("beta.1","lnalpha.1","lnalpha.c.1"), to=c("beta","lnalpha", "lnalpha.c"))
coda2<- rename.vars(coda2, from=c("beta.2","lnalpha.2","lnalpha.c.2"), to=c("beta","lnalpha", "lnalpha.c"))
coda3<- rename.vars(coda3, from=c("beta.3","lnalpha.3","lnalpha.c.3"), to=c("beta","lnalpha", "lnalpha.c"))
coda<-rbind(coda1,coda2,coda3)
coda$Smsy_lambert <- (1-lambert_W0(exp(1-coda$lnalpha.c)))/coda$beta 
coda$Umsy_lambert <- (1-lambert_W0(exp(1-coda$lnalpha.c))) 
coda<-as.data.frame(coda)
summary<-summary(coda) 
q1<-apply(coda,2,quantile,probs=c(0,0.025,0.5,0.975,1))#percentiles
write.csv(q1, file= paste("results/AR_quantiles_lambert.csv") )    
write.csv(summary, file= paste("results/AR_lambert.csv") ) 
write.csv(coda, file= paste("results/AR_coda_lambert.csv") ) 

