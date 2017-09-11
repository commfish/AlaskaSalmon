# R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
# point estimate plots, yield plots)
# created by Ben Williams (Ben.Williams@alaska.gov);Nov 3, 2016; 2017-7-4
# Changes made by Sara Miller (Sara.Miller@alaska.gov); April 2017

rm(list=ls(all=T))#Remove previous variables.
LowerB <- 38000 #lower bound of recommended escapement goal range
UpperB <- 86000 #upper bound of recommended escapement goal range
SMSY<-53280#Lambert W from AR_quantiles_lambert
#load----

library(tidyverse)
library(reshape2)
library(extrafont)
library(grid)
library(plyr)
library(gsl)
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+ 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))


library(scales)
# library(extrafont)
# loadfonts(device="win") #only need to do this once; takes awhile to run!

# data----

coda <- read.csv("results/Ricker_AR_coda.csv") #Load Data File
Parameters <- read.csv("data/Parameters.csv") #Load Data File (make sure this file is updated)

#data clean----
#Create profile parameters
#coda %>% 
#  mutate(S.eq.c = lnalpha.c/beta, 
#         S.msy.c = S.eq.c*(0.5-(((0.65*lnalpha.c)^1.27)/((8.7+lnalpha.c)^1.27))),
#         R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c),
#         MSY.c = R.msy.c-S.msy.c,
#         Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda
coda %>% mutate(S.eq.c = lnalpha.c/beta, 
                S.msy.c = (1-lambert_W0(exp(1-lnalpha.c)))/beta, #Lambert W
                R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c), 
                MSY.c = R.msy.c-S.msy.c, 
                Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda
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
  dat10 <- data.frame(S0=rep(0, length(coda[,1])))
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
    dat10[,i+1] = x[i]*exp(coda$lnalpha.c-coda$beta*x[i])
  }
  # Overfishing estimate ----
  f.over <- function(x){
    x %>% 
      filter(complete.cases(.)) %>% 
      summarise_all(funs(mean)) %>% 
      gather() %>% 
      select(value)
  }
  
  of_0.7 <- f.over(dat)
  of_0.8 <- f.over(dat3)
  of_0.9 <- f.over(dat6)
  
  # Optimal yield estimate ----
  oy_0.7 <- f.over(dat1)
  oy_0.8 <- f.over(dat4)
  oy_0.9 <- f.over(dat7)
  
  # Optimal recruitment ----
  or_0.7 <- f.over(dat2)
  or_0.8 <- f.over(dat5)
  or_0.9 <- f.over(dat8)
  
  #Bind dataframes together
  Y <- cbind(of_0.7,oy_0.7,or_0.7,of_0.8,oy_0.8,or_0.8,of_0.9,oy_0.9,or_0.9, c(0, x))
  names(Y) <- c('of_0.7','oy_0.7','or_0.7','of_0.8','oy_0.8','or_0.8','of_0.9','oy_0.9',
                'or_0.9','Escapement')
  
  #Quantiles and Medians ----
  summarise_all(dat9, funs(median, 
                           q95=quantile(., 0.95, na.rm=T), 
                           q90=quantile(., 0.90, na.rm=T),
                           q10=quantile(., 0.10, na.rm=T),
                           q5=quantile(., 0.05, na.rm=T))) -> mq
  names(mq) <- c(rep(('Median'),length(x)+1), 
                 rep(('q95'),length(x)+1), 
                 rep(('q90'),length(x)+1), 
                 rep(('q10'),length(x)+1), 
                 rep(('q5'),length(x)+1))
  
  qm <- data.frame(measure = names(mq), value = as.numeric(mq[1,]), Escapement=rep(c(0,x), length(unique(names(mq)))))
  qm <- spread(qm, measure, value)
  qm <- qm[c("q95", "q90", "Median","q10", "q5", "Escapement")]
  Y <- Y[c("oy_0.9", "oy_0.8", "or_0.9","or_0.8", "of_0.9", "of_0.8", "oy_0.7","or_0.7","of_0.7","Escapement")]
  write.csv(qm,("data/processed/QM.csv"), row.names=FALSE)
  write.csv(Y,("data/processed/Y.csv"), row.names=FALSE)
  
  #confidence intervals ----
  summarise_all(dat10, funs(median, 
                            q95=quantile(., 0.95, na.rm=T), 
                            q90=quantile(., 0.90, na.rm=T),
                            q10=quantile(., 0.10, na.rm=T),
                            q5=quantile(., 0.05, na.rm=T))) -> mq
  names(mq) <- c(rep(('Median'),length(x)+1), 
                 rep(('q95'),length(x)+1), 
                 rep(('q90'),length(x)+1), 
                 rep(('q10'),length(x)+1), 
                 rep(('q5'),length(x)+1))
  
  CI <- data.frame(measure = names(mq), value = as.numeric(mq[1,]), Escapement=rep(c(0,x), length(unique(names(mq)))))
  CI <- spread(CI, measure, value)
  CI <- CI[c("q95", "q90", "Median","q10", "q5", "Escapement")]
  write.csv(CI,("data/processed/CI.csv"), row.names=FALSE) #confidence intervals around S-R relationship
  #create probability profile plots (0.7, 0.8, 0.9, 0.8 & 0.9)
  Y %>% 
    dplyr::select(Escapement, oy_0.7, of_0.7,or_0.7) %>% 
    melt(., id.vars = 'Escapement') %>% 
    ggplot( aes(Escapement/1000, value, lty=variable))+geom_line()+
    xlab('Escapement (1,000)')+ylab('Probability')+
    theme(legend.justification=c(1,0), legend.position=c(1,.5), 
          legend.key = element_blank(),legend.title=element_blank())
  ggsave("figures/0.7.AR.png", dpi=200, width=8, height=5, units='in')
  
  Y %>% 
    dplyr::select(Escapement, oy_0.8, of_0.8, or_0.8) %>% 
    melt(., id.vars = 'Escapement') %>% 
    ggplot(aes(Escapement/1000, value, lty=variable))+geom_line()+
    xlab('Escapement (1,000)')+ylab('Probability')+
    theme(legend.justification=c(1,0), legend.position=c(1,.5), 
          legend.key = element_blank(),legend.title=element_blank())
  ggsave("figures/0.8.AR.png", dpi=200, width=8, height=5, units='in')
  
  Y %>% 
    dplyr::select(Escapement, oy_0.9, of_0.9, or_0.9) %>% 
    melt(., id.vars = 'Escapement')  %>% 
    ggplot(aes(Escapement/1000, value, lty=variable))+geom_line()+
    xlab('Escapement (1,000)')+ylab('Probability')+
    theme(legend.justification=c(1,0), legend.position=c(1,.5), 
          legend.key = element_blank(),legend.title=element_blank())
  ggsave("figures/0.9.AR.png", dpi=200, width=8, height=5, units='in')
  
  
  Y <- read.csv("data/processed/Y.csv")
  Y["OY0.9"] <-Y$oy_0.9 
  Y["OY0.8"] <-Y$oy_0.8 
  Y<-subset(Y, select=c(Escapement, OY0.9,OY0.8))
  mY1 <- melt(Y, id.vars='Escapement')
  mY1["sra"] <-"Optimal Yield Profile"
  mY1["max_pct"] <- ifelse(grepl("OY0.8",mY1$variable), 
                           0.8,0.9)
  
  Y <- read.csv("data/processed/Y.csv")
  Y["OF0.9"] <-Y$of_0.9 
  Y["OF0.8"] <-Y$of_0.8 
  Y<-subset(Y, select=c(Escapement, OF0.9,OF0.8))
  mY2 <- melt(Y, id.vars='Escapement')
  mY2["sra"] <-"Overfishing Profile"
  mY2["max_pct"] <- ifelse(grepl("OF0.8",mY2$variable), 
                           0.8,0.9)
  
  Y <- read.csv("data/processed/Y.csv")
  Y["OR0.9"] <-Y$or_0.9 
  Y["OR0.8"] <-Y$or_0.8 
  Y<-subset(Y, select=c(Escapement, OR0.9,OR0.8))
  mY3 <- melt(Y, id.vars='Escapement')
  mY3["sra"] <-"Optimal Recruitment Profile"
  mY3["max_pct"] <- ifelse(grepl("OR0.8",mY3$variable), 
                           0.8,0.9)
  mY4<-rbind(mY1,mY2, mY3)
  mY4<-subset(mY4, select=c(Escapement, value, sra, max_pct))
  mY4$Escapement<-as.numeric(mY4$Escapement)
  mY4$value<-as.numeric(mY4$value)
  mY4$max_pct<-as.factor(mY4$max_pct)
  colnames(mY4)[2] <- "Probability"
  windowsFonts(Times=windowsFont("TT Times New Roman"))
  theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  Fig1<-ggplot(mY4, aes(x = Escapement, y = Probability, linetype = max_pct)) 
  Fig1<-Fig1+geom_rect(aes(xmin = LowerB, xmax = UpperB, ymin = 0, ymax = 1),
                       inherit.aes = FALSE, fill = "grey80", alpha = 0.3)+geom_line()+xlab('Escapement (S)')+
    scale_x_continuous(labels = comma, breaks = seq(0, 350000, 50000), limits = c(0, 350000))+
    scale_linetype_discrete(name = "Percent of Maximum")+
    facet_grid(sra ~ .) +geom_vline(xintercept=SMSY, lwd=1.25)
  theme_bw()+ theme(legend.key = element_blank())+
    theme(text=element_text(family="Times New Roman"))
  ggsave("figures/0.8_0.9.png", dpi=200, dev='png', width=7, height=6, units='in')
  theme_set(theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  options(scipen=99999)
  
  ggplot(qm, aes(Escapement, Median))+geom_line(size=1)+
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.15)+
    geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.15)+ xlab('Escapement (S)')+
    ylab('Expected Yield')+scale_y_continuous(labels = comma)+
    scale_x_continuous(labels = comma,breaks = seq(0, 300000, 50000), limits = c(0,300000))+
    geom_vline(xintercept = LowerB,linetype = "longdash" )+geom_vline(xintercept = UpperB ,linetype = "longdash")
  ggsave("figures/expected_sustained_yield.png", dpi=200, width=8, height=5, units='in')
}
#Run function
profile(i=10,z=500,xa.start=0, xa.end=700,lnalpha.c, beta)#can change i,z, xa.start, xa.end
#profile(i=2,z=500,xa.start=0, xa.end=700,lnalpha.c, beta)#can change i,z, xa.start, xa.end
  ################################################################################################
  #Horesetail Plots
  #Figure x.- Graphical summary of knowledge of spawner-recruitment relationship for Chilkat Lake 
  #sockeye salmon as derived from age-structured state space model fitted to abundance, harvest, 
  #and age data 1976-2015. Symbols are posterior medians of R and S; error bars bracket 95% credibility 
  #intervals. Heavy dashed line is Ricker relationship constructed from ln(a) and b posterior medians,  
  #Ricker relationships are also plotted for 50 paired values of ln(a) and b sampled from the posterior 
  #probability distribution, representing plausible Ricker relationships that could have generated the 
  #observed data.  Diagonal line is replacement line (R=S).
  ################################################################################################
  coda <- read.csv("results/Ricker_AR_coda.csv") #Load Data File
  Parameters <- read.csv("data/Parameters.csv")
  QM <- read.csv("data/processed/QM.csv")
  CI<- read.csv("data/processed/CI.csv")
  coda <- subset(coda, select = c(lnalpha.c, beta))
  coda<-as.data.frame(coda[1:50,])#select first 50 rows of dataframe
  QM<-subset(QM, select = c(Escapement))
  num<-nrow(QM)
  coda<-data.frame(lnalpha.c=rep(coda$lnalpha.c,each=num), beta=rep(coda$beta, each=num))
  dataset<-cbind(coda,QM) #lnalpha.c, beta, and S
  dataset['Recruitment']<-dataset$Escapement*exp(dataset$lnalpha.c-dataset$beta*dataset$Escapement)
  dataset['Variable']<-data.frame(dataset=rep(1:50,each=num))
  
  myvars <- c("lnalpha.c", "beta")
  Parameters<-Parameters[myvars]
  Parameters<-data.frame(lnalpha.c=rep(Parameters$lnalpha.c,each=num), beta=rep(Parameters$beta, each=num))
  row.has.na <- apply(Parameters, 1, function(x){any(is.na(x))})
  sum(row.has.na)
  final.filtered <- Parameters[!row.has.na,]
  final.filtered <-cbind(final.filtered,QM)
  final.filtered['Recruitment']<-final.filtered$Escapement*exp(final.filtered$lnalpha.c-final.filtered$beta*final.filtered$Escapement)
  final.filtered['Variable']<-data.frame(final.filtered=rep(51,each=num))
  dataset<-rbind(dataset, final.filtered)
  dataset['Year']<-'NA'
  Parameters <- read.csv('./data/Parameters.csv') #Load Data File
  Parameters<-subset(Parameters, select = c(S_median, R_median,Year))
  Parameters['lnalpha.c']<-'NA'
  Parameters['beta']<-'NA'
  Parameters['Escapement']<-Parameters$spawn
  Parameters['Recruitment']<-Parameters$recruit
  Parameters['Variable']<-52
  Parameters<-subset(Parameters, select=c(lnalpha.c, beta, Escapement, Recruitment, Variable,Year))
  dataset<-rbind(dataset, Parameters)
  dataset <- dataset[order(dataset$Variable, dataset$Escapement),] 
  windowsFonts(Times=windowsFont("TT Times New Roman"))
  theme_set(theme_bw(base_size=12,base_family='Times New Roman')+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  Fig<-ggplot(data=dataset, aes(x=Escapement, y=Recruitment, group=Variable))+
    geom_line(data=subset(dataset,dataset$Variable<52),linetype="solid", size=0.5, color="grey80")+
    scale_y_continuous(labels = comma,breaks = seq(0, 1000000, 50000), limits = c(0, 1000000))+
    scale_x_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0, 350000))+ylab("Recruits (R)")+xlab("Spawners (S)")+
    geom_line(data=dataset, aes(x=Escapement, y=Escapement, group=1),linetype="solid", size=1)+#replacement line
    geom_line(data=subset(dataset,dataset$Variable==51),colour = "black", lty=2, size=2)+
    geom_text(data=subset(dataset,dataset$Variable==52), aes(x=Escapement, y=Recruitment, label=Year,family="Times"))+
    #geom_errorbar(aes(ymax = R_val97.5pc, ymin=R_val2.5pc), width=0.20,linetype = 2, colour="grey50")+
    #geom_errorbarh(aes(xmax = S_val97.5pc, xmin=S_val2.5pc), height=0.20,linetype = 2,colour="grey50")
  png(file='figures/Horsetail_Plot.png', res=200, width=6, height=4, units ="in")  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1,1)))
  vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
  print(Fig,vp=vplayout(1,1:1)) 
  dev.off()
  
  dataset <- subset(dataset,Variable == 52) 
  dataset['Escapement1']<-dataset$Escapement
  dataset<-subset(dataset, select = -c(lnalpha.c, beta, Escapement))
  dataset['Escapement']<-'NA'
  CI['Year']<-'NA'
  CI['Variable']<-51
  CI['Recruitment']<-'NA'
  CI['Escapement1']<-'NA'
  dataset['Median']<-'NA'
  dataset['q95']<-'NA'
  dataset['q90']<-'NA'
  dataset['q10']<-'NA'
  dataset['q5']<-'NA'
  dataset1<-rbind(dataset, CI)
  dataset1$Escapement<-as.numeric(dataset1$Escapement)
  dataset1$Median<-as.numeric(dataset1$Median)
  dataset1$q5<-as.numeric(dataset1$q5)
  dataset1$q95<-as.numeric(dataset1$q95)
  dataset1$q10<-as.numeric(dataset1$q10)
  dataset1$q90<-as.numeric(dataset1$q90)
  dataset1$Recruitment<-as.numeric(dataset1$Recruitment)
  dataset1$Escapement1<-as.numeric(dataset1$Escapement1)
  
  Fig1<-ggplot(data=dataset1, aes(x=Escapement, y=Median, group=Variable))+geom_line(size=2, lty=2, group=51)+
    geom_ribbon(aes(ymin = q5, ymax = q95, group=51), alpha=.15)+
    geom_ribbon(aes(ymin = q10, ymax = q90, group=51), alpha=.15)+
    xlab('Spawners (S)')+
    ylab('Recruits (R)')+scale_y_continuous(labels = comma)+
    scale_x_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0,350000))+
    scale_y_continuous(labels = comma,breaks = seq(0, 600000, 100000), limits = c(0,600000))+
    geom_line(aes(x=Escapement, y=Escapement, group=51),linetype="solid", size=1)+
    geom_text(size=3, data=dataset1, aes(x=Escapement1, y=Recruitment, group=52, label=Year,family="Times", 
                                         hjust = -0.1, vjust= -0.4))
  Fig1<-Fig1+
    geom_point(data=dataset1, aes(x=Escapement1, y=Recruitment, group=52),pch=16, size=1)
    #geom_errorbar(data=dataset1, aes(x=Escapement1, ymax=R_val97.5pc, ymin=R_val2.5pc, group=52), width=0.2,linetype = 1, colour="grey50")
  Fig1<-Fig1+
    geom_point(data=dataset1, aes(x=Escapement1, y=Recruitment, group=52),pch=16, size=1)
    #geom_errorbarh(data=dataset1, aes(x=Escapement1, y=Recruitment, xmax=S_val97.5pc, xmin=S_val2.5pc),na.rm=T,  linetype = 1, colour="grey50")
  
  png(file='figures/Horsetail_Plot_Reconfig.png', res=200, width=6, height=4, units ="in")  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1,1)))
  vplayout<-function(x,y) viewport (layout.pos.row=x, layout.pos.col=y)
  print(Fig1,vp=vplayout(1,1:1)) 
  dev.off()
  
 