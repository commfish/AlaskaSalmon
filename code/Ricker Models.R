
#Original code was from Steve Fleischman and Xinxian Zhang. code was modified
#and converted to JAGS by Sara Miller on 28 June 2017 and maintained
#Rich Brenner.

#load----
library(arm)
library(lmtest)
library(rjags)
library(R2OpenBUGS)

#data----
brood<-read.table("data/Chilkoot_Sock.csv", header=TRUE, sep=",")
sr<-brood[,2:3] #retrieve spawner recruit data for one stock
colnames(sr)<-c("S","R") #rename the column names
sr<-na.omit(sr)  #remove missing rows with missing values(na)
sr$lnRS<-log(sr$R/sr$S) #add one column,lnRS
n<-nrow(sr)
sr.data<-list(n=n,S=sr$S, lnRS=sr$lnRS)
#analysis----

#check AR(1)
dwtest(sr$lnRS ~ sr$S)
mylm <- lm(lnRS ~ S, sr)
pacf(residuals(mylm))

##Ricker model for stock-recruitment analysis 
#Created by Steve Fleischman 

#Ricker without autocorrelation
Ricker=function(){
  
  lnalpha ~ dunif(0, 10) 
  beta ~ dunif(0, 10)                 
  phi <- 0                #get rid of AR(1)
  sigma.white ~ dunif(0,10)
  resid.red.0 ~ dnorm(0,tau.red)
  
  for(y in 1:n) {lnRS[y] ~ dnorm(mean2.lnRS[y],tau.white) }
  
  mean2.lnRS[1]     <- mean1.lnRS[1] + phi * resid.red.0  
  for (y in 2:n) { mean2.lnRS[y] <- mean1.lnRS[y] + phi * resid.red[y-1] }  #no autocorrelation model
  
  for(y in 1:n) {  mean1.lnRS[y] <- lnalpha - beta * S[y]  }
  for(y in 1:n) {  resid.red[y]     <- lnRS[y] - mean1.lnRS[y]  }
  for(y in 1:n) {  resid.white[y] <- lnRS[y] - mean2.lnRS[y]  }
  
  tau.white <- 1 / sigma.white / sigma.white        
  tau.red <- tau.white * (1-phi*phi)
  sigma.red <- 1 / sqrt(tau.red)
  
  lnalpha.c <- lnalpha + (sigma.red * sigma.red / 2)  #adjust for calculating means of R.msy, S.msy etc.
  #lnalpha.c <- lnalpha
  alpha<-exp(lnalpha)
  S.max <- 1 / beta
  S.eq <- S.max * lnalpha.c 
  S.msy <- S.eq * (0.5 - 0.07*lnalpha.c)
  U.msy <- lnalpha.c * (0.5 - 0.07*lnalpha.c)
  R.msy <- S.msy * exp(lnalpha.c - beta * S.msy)
  MSY<-step(R.msy-S.msy)*(R.msy-S.msy) #if R.msy< S.msy then MSY=0.
  #step(x) = 1 if x>=0; otherwise =0 if x<0
  
  
  S.star[1]<-0
  step<-400
  for (i in 2:501) {                      #LOOP TO FIND Pr(SY>90%MSY)
    S.star[i] <- S.star[i-1]+step
    R.star[i] <- S.star[i] * exp(lnalpha.c - beta * S.star[i]) 
    SY[i] <- R.star[i] - S.star[i]
    I90[i] <- step(SY[i] - 0.9 * MSY)  
  }
}
#write the model to a text file to be called by WinBUGS
model_file_loc=paste("code/Chilkoot_Sockeye.txt", sep="")
write.model(Ricker, paste("code/Chilkoot_Sockeye.txt", sep=""))

#Ricker with autocorrelation
AR=function(){
  
  lnalpha ~ dunif(0, 10)
  beta ~ dunif(0, 10) 
  phi ~ dunif(-1,1)                #AR(1)
  sigma.white ~ dunif(0,10)
  resid.red.0 ~ dnorm(0,tau.red)
  
  for(y in 1:n) {lnRS[y] ~ dnorm(mean2.lnRS[y],tau.white) }
  
  mean2.lnRS[1]     <- mean1.lnRS[1] + phi * resid.red.0  
  for (y in 2:n) { mean2.lnRS[y] <- mean1.lnRS[y] + phi * resid.red[y-1] }   #AR1
  
  for(y in 1:n) {  mean1.lnRS[y] <- lnalpha - beta * S[y]  } #Ricker model
  for(y in 1:n) {  resid.red[y]     <- lnRS[y] - mean1.lnRS[y]  }
  for(y in 1:n) {  resid.white[y] <- lnRS[y] - mean2.lnRS[y]  }
  
  tau.white <- 1 / sigma.white / sigma.white        
  tau.red <- tau.white * (1-phi*phi)
  sigma.red <- 1 / sqrt(tau.red)
  #sigma.white<-1/sqrt(tau.white)
  #sigma<-sigma.red
  
  lnalpha.c <- lnalpha + (sigma.red * sigma.red / 2)  #adjust for calculating means of R.msy, S.msy etc.
  #lnalpha.c <- lnalpha
  alpha<-exp(lnalpha)
  S.max <- 1 / beta
  S.eq <- S.max * lnalpha.c 
  S.msy <- S.eq * (0.5 - 0.07*lnalpha.c)
  U.msy <- lnalpha.c * (0.5 - 0.07*lnalpha.c)
  R.msy <- S.msy * exp(lnalpha.c - beta * S.msy)
  
  MSY<-step(R.msy-S.msy)*(R.msy-S.msy) #if R.msy< S.msy then MSY=0.
  #step(x) = 1 if x>=0; otherwise =0 if x<0
  
  S.star[1]<-0
  step<-400
  for (i in 2:501) {                      #LOOP TO FIND Pr(SY>90%MSY)
    S.star[i] <- S.star[i-1]+step
    R.star[i] <- S.star[i] * exp(lnalpha.c - beta * S.star[i]) 
    SY[i] <- R.star[i] - S.star[i]
    I90[i] <- step(SY[i] - 0.9 * MSY)  
  }
}

#write the model to a text file to be called by WinBUGS
model_file_loc=paste("code/Chilkoot_Sockeye_AR.txt", sep="")
write.model(AR, paste("code/Chilkoot_Sockeye_AR.txt", sep=""))

#results----
#Run Ricker without AR(1)
inits1<-list(lnalpha=1.5, beta=0.0005, sigma.white=0.7, resid.red.0= 0)
inits2<-list(lnalpha=2.0, beta=0.0010, sigma.white=0.5, resid.red.0=-1)
inits3<-list(lnalpha=2.5, beta=0.0020, sigma.white=0.3, resid.red.0= 1)
inits<-list(inits1, inits2, inits3)
#parameters<-c("lnalpha","beta", "sigma.red","S.msy","MSY", "I90")
parameters<-c("lnalpha","beta", "sigma.red","S.msy","MSY")
ptm = proc.time()
jmod <- jags.model(file='code/Chilkoot_Sockeye.txt', data=sr.data, n.chains=3, inits=inits, n.adapt=1000) 
x<-update(jmod, n.iter=100000, by=5, progress.bar='text', DIC=T, n.burnin=1000) 
post <- coda.samples(jmod, parameters, n.iter=100000, thin=5, n.burnin=1000)
post.samp <- post

#Numerical summary of each parameter (mean, median, quantiles)
summary<-summary(post)                     
stats<-summary$statistics;  colnames(stats)
quants<-summary$quantiles;  colnames(quants)
statsquants <- cbind(stats,quants) 
statsquants <- statsquants[,c(1,2,4,5,7,9)] #select columns of interest
write.csv(statsquants, file= paste("results/Ricker.csv") )    

#Density and time series plots
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


#Gelman
gel <- as.data.frame(gelman.diag(post, multivariate=F)[[1]])
poor.threshold=1.10#values less than 1.2 are generally considered converged
poor <- gel[gel[,1] > poor.threshold, ]
write.csv(poor, file= paste("results/Ricker_Gelman.csv") )    

#DIC
dic.pD  <-dic.samples(jmod,n.iter=100000, thin=100,"pD",  n.burnin=1000)
dic.popt<-dic.samples(jmod,n.iter=100000, thin=100,"popt",n.burnin=1000)
dev1 <- sum(dic.pD[[1]])
pD   <- sum(dic.pD[[2]])
dic.pD <- dev1 + pD
dic.pD.summary <- data.frame(dev1, pD, dic.pD)
write.csv(dic.pD.summary, file=paste("results/Ricker_DIC.csv") ) 

#Create coda samples for horsetail plots and probability plots
post2 <- coda.samples(jmod, c("lnalpha", "beta", "lnalpha.c"), n.iter=100000, thin=100,n.burnin=1000) #n.iter=10000 & thin=10 is 1000 samples
x <- as.array(post2)
x <- data.frame(x)
coda <- x[,1:3]
coda<- rename.vars(coda, from=c("beta.1","lnalpha.1","lnalpha.c.1"), to=c("beta","lnalpha", "lnalpha.c"))
write.csv(coda, file= paste("results/Ricker_coda.csv") ,row.names=FALSE)    # writes csv file


#Run Ricker with AR(1):  
inits1<-list(lnalpha=1.5, beta=0.0005, phi= 0.3, sigma.white=0.7, resid.red.0= 0)
inits2<-list(lnalpha=2.0, beta=0.0010, phi=-0.1, sigma.white=0.5, resid.red.0=-1)
inits3<-list(lnalpha=2.5, beta=0.0020, phi= 0.2, sigma.white=0.3, resid.red.0= 1)
inits<-list(inits1, inits2, inits3)
#parameters<-c("lnalpha.c","beta", "sigma.red","S.msy", "MSY", "I90" )
parameters<-c("lnalpha.c","beta", "sigma.red","S.msy", "MSY")
jmod <- jags.model(file='code/Chilkoot_Sockeye_AR.txt', data=sr.data, n.chains=3, inits=inits, n.adapt=1000) 
x<-update(jmod, n.iter=100000, by=5, progress.bar='text', DIC=T, n.burnin=1000) 
post <- coda.samples(jmod, parameters, n.iter=100000, thin=5, n.burnin=1000)
post.samp <- post

#Numerical summary of each parameter (mean, median, quantiles)
summary<-summary(post)                     
stats<-summary$statistics;  colnames(stats)
quants<-summary$quantiles;  colnames(quants)
statsquants <- cbind(stats,quants) 
statsquants <- statsquants[,c(1,2,4,5,7,9)] #select columns of interest
write.csv(statsquants, file= paste("results/Ricker_AR.csv") ) 

#Density and time series plots
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

#Gelman
gel <- as.data.frame(gelman.diag(post, multivariate=F)[[1]])
poor.threshold=1.10#values less than 1.2 are generally considered converged
poor <- gel[gel[,1] > poor.threshold, ]
write.csv(poor, file= paste("results/Ricker_AR_Gelman.csv") )    

#DIC
dic.pD  <-dic.samples(jmod,n.iter=100000, thin=100,"pD",  n.burnin=1000)
dic.popt<-dic.samples(jmod,n.iter=100000, thin=100,"popt",n.burnin=1000)
dev1 <- sum(dic.pD[[1]])
pD   <- sum(dic.pD[[2]])
dic.pD <- dev1 + pD
dic.pD.summary <- data.frame(dev1, pD, dic.pD)
write.csv(dic.pD.summary, file=paste("results/Ricker_AR_DIC.csv") ) 

#Create coda samples for horsetail plots and probability plots
post2 <- coda.samples(jmod, c("lnalpha", "beta", "lnalpha.c"), n.iter=100000, thin=100,n.burnin=1000) #n.iter=10000 & thin=10 is 1000 samples
x <- as.array(post2)
x <- data.frame(x)
coda <- x[,1:3]
coda<- rename.vars(coda, from=c("beta.1","lnalpha.1","lnalpha.c.1"), to=c("beta","lnalpha", "lnalpha.c"))
write.csv(coda, file= paste("results/Ricker_AR_coda.csv") ,row.names=FALSE)    # writes csv file
