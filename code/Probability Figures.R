#R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
#point estimate plots, yield plots)
#created by Ben Williams (Ben.Williams@alaska.gov);Nov 3, 2016
#Changes made by Sara Miller (Sara.Miller@alaska.gov); April 2017
#results are in data/processed/..
#i and z act as ways to change range of escapement based on stock size

rm(list=ls(all=T))#Remove previous variables.
LowerB<-65000 #lower bound of recommended escapement goal range
UpperB<-140000 #upper bound of recommended escapement goal range
#load----
#Load Packages
library(plyr)
library(reshape2)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(MASS)
library(survival)
library(scatterplot3d)
library(vcd)
library(grid)
library(calibrate)
library(scales)
library(extrafont)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(cowplot)
#data----
#loadfonts(device="win") #only need to do this once; takes awhile to run!
coda <- read.csv("results/Ricker_coda.csv") #Load Data File

#data clean----
#Create profile parameters
coda %>% mutate(S.eq.c = lnalpha.c/beta, 
                S.msy.c = S.eq.c*(0.5-(((0.65*lnalpha.c)^1.27)/((8.7+lnalpha.c)^1.27))),
					      R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c), 
					      MSY.c = R.msy.c-S.msy.c, 
					      Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda

attach(coda)
#analysis----
#create function for probability profiles and figures
profile <-function(i,z,xa.start, xa.end,lnalpha.c, beta){ 
xa = seq(xa.start, xa.end, by=i) 
x =(xa+i)*z
# empty dataframes
dat <- data.frame(S0=rep(1, length(coda[,1])))
dat1 <- data.frame(S0=rep(0, length(coda[,1])))
dat2 <- data.frame(S0=rep(0, length(coda[,1])))
dat3 <- data.frame(S0=rep(1, length(coda[,1])))
dat4 <- data.frame(S0=rep(0, length(coda[,1])))
dat5 <- data.frame(S0=rep(0, length(coda[,1])))
dat6 <- data.frame(S0=rep(1, length(coda[,1])))
dat7 <- data.frame(S0=rep(0, length(coda[,1])))
dat8 <- data.frame(S0=rep(0, length(coda[,1])))
dat9 <- data.frame(S0=rep(0, length(coda[,1])))
for (i in 1:length(xa)){
  dat[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.7*coda$MSY.c), 0, ifelse(dat[,i]==0, 0,1))
  dat1[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.7*coda$MSY.c), 1,0)
  dat2[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i]))>(0.7*coda$Rmax), 1,0)
  dat3[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.8*coda$MSY.c), 0, ifelse(dat3[,i]==0, 0,1))
  dat4[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.8*coda$MSY.c), 1,0)
  dat5[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i]))>(0.8*coda$Rmax), 1,0)
  dat6[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.9*coda$MSY.c), 0, ifelse(dat6[,i]==0, 0,1))
  dat7[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.9*coda$MSY.c), 1,0)
  dat8[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i]))>(0.9*coda$Rmax), 1,0)
  dat9[,i+1] = x[i]*exp(coda$lnalpha.c-coda$beta*x[i])-x[i]
  }
# Overfishing estimate ----
dat %>% filter(complete.cases(.)) %>% summarise_each(funs(mean)) %>% gather() %>% select(value)-> of_0.7 
dat3 %>% filter(complete.cases(.)) %>% summarise_each(funs(mean)) %>% gather() %>% select(value)-> of_0.8  
dat6 %>% filter(complete.cases(.)) %>% summarise_each(funs(mean)) %>% gather() %>% select(value)-> of_0.9  
# Optimal yield estimate ----
dat1 %>% filter(complete.cases(.)) %>% summarise_each(funs(mean)) %>% gather() %>% select(value)-> oy_0.7 
dat4 %>% filter(complete.cases(.)) %>% summarise_each(funs(mean)) %>% gather() %>% select(value)-> oy_0.8 
dat7 %>% filter(complete.cases(.)) %>% summarise_each(funs(mean)) %>% gather() %>% select(value)-> oy_0.9 
# Optimal recruitment ----
dat2 %>% filter(complete.cases(.)) %>% summarise_each(funs(mean)) %>% gather() %>% select(value) -> or_0.7 
dat5 %>% filter(complete.cases(.)) %>% summarise_each(funs(mean)) %>% gather() %>% select(value) -> or_0.8 
dat8 %>% filter(complete.cases(.)) %>% summarise_each(funs(mean)) %>% gather() %>% select(value) -> or_0.9 
#Bind dataframes together
Y <- cbind(of_0.7,oy_0.7,or_0.7,of_0.8,oy_0.8,or_0.8,of_0.9,oy_0.9,or_0.9, c(0, x))
names(Y) <- c('of_0.7','oy_0.7','or_0.7','of_0.8','oy_0.8','or_0.8','of_0.9','oy_0.9',
              'or_0.9','Escapement')
#Quantiles and Medians ----
summarise_each(dat9, funs(median, q95=quantile(., 0.95, na.rm=T), q90=quantile(., 0.90, na.rm=T),
                          q10=quantile(., 0.10, na.rm=T),q5=quantile(., 0.05, na.rm=T))) -> mq

names(mq) <- c(rep(('Median'),length(x)+1), rep(('q95'),length(x)+1), rep(('q90'),length(x)+1), rep(('q10'),length(x)+1), rep(('q5'),length(x)+1))

qm <- data.frame(measure = names(mq), value = as.numeric(mq[1,]), Escapement=rep(c(0,x), length(unique(names(mq)))))
qm <- spread(qm, measure, value)
qm <- qm[c("q95", "q90", "Median","q10", "q5", "Escapement")]
Y <- Y[c("oy_0.9", "oy_0.8", "or_0.9","or_0.8", "of_0.9", "of_0.8", "oy_0.7","or_0.7","of_0.7","Escapement")]
write.csv(qm,("data/processed/QM.csv"), row.names=FALSE)
write.csv(Y,("data/processed/Y.csv"), row.names=FALSE)

#create probability profile plots (0.7)
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+ 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
Y <- read.csv("data/processed/Y.csv")
Y["Optimal_Yield0.7"] <-Y$oy_0.7 
Y["Overfishing0.7"] <-Y$of_0.7 
Y["Optimal_Recruitment0.7"] <-Y$or_0.7 
Y<-subset(Y, select=c(Escapement, Optimal_Yield0.7, Overfishing0.7, Optimal_Recruitment0.7)) 
mY <- melt(Y, id.vars='Escapement')
ggplot(mY, aes(Escapement/1000, value, lty=variable))+geom_line()+xlab('Escapement / 1,000')+ylab('Probability')+
  theme(legend.justification=c(1,0), legend.position=c(1,.5), legend.key = element_blank(),legend.title=element_blank())
ggsave("figures/0.7.png", dpi=200, dev='png', width=8, height=5, units='in')

Y <- read.csv("data/processed/Y.csv")
Y["Optimal_Yield0.8"] <-Y$oy_0.8 
Y["Overfishing0.8"] <-Y$of_0.8 
Y["Optimal_Recruitment0.8"] <-Y$or_0.8 
Y<-subset(Y, select=c(Escapement, Optimal_Yield0.8, Overfishing0.8, Optimal_Recruitment0.8)) 
mY <- melt(Y, id.vars='Escapement')
ggplot(mY, aes(Escapement/1000, value, lty=variable))+geom_line()+xlab('Escapement / 1,000')+ylab('Probability')+
  theme(legend.justification=c(1,0), legend.position=c(1,.5), legend.key = element_blank(),legend.title=element_blank())
ggsave("figures/0.8.png", dpi=200, dev='png', width=8, height=5, units='in')

Y <- read.csv("data/processed/Y.csv")
Y["Optimal_Yield0.9"] <-Y$oy_0.9 
Y["Overfishing0.9"] <-Y$of_0.9 
Y["Optimal_Recruitment0.9"] <-Y$or_0.9 
Y<-subset(Y, select=c(Escapement, Optimal_Yield0.9, Overfishing0.9, Optimal_Recruitment0.9)) 
mY <- melt(Y, id.vars='Escapement')
ggplot(mY, aes(Escapement/1000, value, lty=variable))+geom_line()+xlab('Escapement / 1,000')+ylab('Probability')+
  theme(legend.justification=c(1,0), legend.position=c(1,.5), legend.key = element_blank(),legend.title=element_blank())
ggsave("figures/0.9.png", dpi=200, dev='png', width=8, height=5, units='in')
theme_set(theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
options(scipen=99999)

Y <- read.csv("data/processed/Y.csv")
Y["OY0.7"] <-Y$oy_0.7 
Y["OY0.8"] <-Y$oy_0.8 
Y<-subset(Y, select=c(Escapement, OY0.8, OY0.7))
mY1 <- melt(Y, id.vars='Escapement')
mY1["sra"] <-"Optimal Yield Profile"
mY1["max_pct"] <- ifelse(grepl("OY0.8",mY1$variable), 
                         0.8,0.7)

Y <- read.csv("data/processed/Y.csv")
Y["OF0.7"] <-Y$of_0.7 
Y["OF0.8"] <-Y$of_0.8 
Y<-subset(Y, select=c(Escapement, OF0.7,OF0.8))
mY2 <- melt(Y, id.vars='Escapement')
mY2["sra"] <-"Overfishing Profile"
mY2["max_pct"] <- ifelse(grepl("OF0.8",mY2$variable), 
                         0.8,0.7)

Y <- read.csv("data/processed/Y.csv")
Y["OR0.7"] <-Y$or_0.7 
Y["OR0.8"] <-Y$or_0.8 
Y<-subset(Y, select=c(Escapement, OR0.7,OR0.8))
mY3 <- melt(Y, id.vars='Escapement')
mY3["sra"] <-"Optimal Recruitment Profile"
mY3["max_pct"] <- ifelse(grepl("OR0.8",mY3$variable), 
                         0.8,0.7)
mY4<-rbind(mY1,mY2, mY3)
mY4<-subset(mY4, select=c(Escapement, value, sra, max_pct))
mY4$Escapement<-as.numeric(mY4$Escapement)
mY4$value<-as.numeric(mY4$value)
mY4$max_pct<-as.factor(mY4$max_pct)
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
Fig1<-ggplot(mY4, aes(x = Escapement, y = value, linetype = max_pct)) 
Fig1<-Fig1+geom_rect(aes(xmin = LowerB, xmax = UpperB, ymin = 0, ymax = 1),
inherit.aes = FALSE, fill = "grey80", alpha = 0.3)+
geom_line()+xlab('Escapement (S)')+
scale_x_continuous(labels = comma, breaks = seq(0, 350000, 50000), limits = c(0, 350000))+
scale_linetype_discrete(name = "Percent of Max.")+
facet_grid(sra ~ .) +
theme_bw()+ theme(legend.key = element_blank())+scale_y_continuous("Probability", breaks = seq(0, 1, 0.2), limits = c(0, 1))+
theme(text=element_text(family="Times New Roman"))
ggsave("figures/0.8_0.7.png", dpi=200, dev='png', width=7, height=6, units='in')
theme_set(theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
options(scipen=99999)

QM <- read.csv("data/processed/QM.csv")
mQM <- melt(QM, id.vars='Escapement')
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
ggplot(QM, aes(Escapement, Median))+geom_line(size=1)+
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.15)+
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.15)+ xlab('Escapement (S)')+
  ylab('Expected Yield')+scale_y_continuous(labels = comma)+
 scale_x_continuous(labels = comma,breaks = seq(0, 300000, 50000), limits = c(0,300000))+
  geom_vline(xintercept = LowerB,linetype = "longdash" )+geom_vline(xintercept = UpperB ,linetype = "longdash")

ggsave("figures/Expected Sustained Yield (QM).png", dpi=200, dev='png', width=8, height=5, units='in')
}
#Run function
profile(i=10,z=500,xa.start=0, xa.end=700,lnalpha.c, beta)#can change i,z, xa.start, xa.end


