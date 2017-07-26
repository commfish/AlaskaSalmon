# R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
# point estimate plots, yield plots)
# created by Ben Williams (Ben.Williams@alaska.gov);Nov 3, 2016; 2017-7-4
# Changes made by Sara Miller (Sara.Miller@alaska.gov); April 2017
# Step#1: Create three csv files from OpenBUGS or JAGS output
# a)lnalpha, beta, and lnalpha.c called 'coda' with 1000 values of each in columns with variable names
# b)Parameters.csv
# c)p_q_Nya.csv
# Step#2: Set working directory to correct location.
# Step#3: Run code.
# i and z act as ways to change range of escapement based on stock size

rm(list=ls(all=T))#Remove previous variables.
LowerB <- 65000 #lower bound of recommended escapement goal range
UpperB <- 140000 #upper bound of recommended escapement goal range

#load----

library(tidyverse)
library(reshape2)
library(extrafont)
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=12,base_family='Times New Roman')+ 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))


# library(plyr)
# library(reshape2)
# library(lattice)
# library(latticeExtra)
# library(gridExtra)
# library(ggplot2)
# library(MASS)
# library(survival)
# library(scatterplot3d)
# library(vcd)
# library(grid)
# library(calibrate)
library(scales)
# library(extrafont)
# loadfonts(device="win") #only need to do this once; takes awhile to run!

# library(RColorBrewer)
# library(tidyverse)
# library(reshape2)
# library(cowplot)

# data----

coda <- read_csv("data/Coda.csv") #Load Data File
Parameters <- read_csv("data/Parameters.csv") #Load Data File
p_q_Nya <- read_csv("data/p_q_Nya.csv") #Load Data File

#data clean----
#Create profile parameters
coda %>% 
  mutate(S.eq.c = lnalpha.c/beta, 
         S.msy.c = S.eq.c*(0.5-(((0.65*lnalpha.c)^1.27)/((8.7+lnalpha.c)^1.27))),
         R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c),
         MSY.c = R.msy.c-S.msy.c,
         Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda

# analysis----
# create function for probability profiles and figures

profile <-function(i,z,xa.start, xa.end, data){ 
  xa = seq(xa.start, xa.end, by=i) 
  x =(xa + i) * z
  
  # create empty dataframes
  dat <- data.frame(S0 = rep(0, nrow(data)))
  dat1 <- data.frame(S0 = rep(1, nrow(data)))
  
  dat2 <- dat4 <- dat5 <- dat7 <- dat8 <- dat9 <- dat
  dat3 <- dat6 <- dat1
  
  for (i in 1:length(xa)){
    dat [,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.7*coda$MSY.c), 0, ifelse(dat[,i]==0, 0,1))
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
  
  
  
  write.csv(qm,("data/processed/QM.csv"), row.names=FALSE)
  write.csv(Y,("data/processed/Y.csv"), row.names=FALSE)
  
  
  #create probability profile plots (0.7)
  
  Y %>% 
    dplyr::select(Escapement, Optimal_Yield0.7 = oy_0.7, 
                  Overfishing0.7 = of_0.7, Optimal_Recruitment0.7 = or_0.7) %>% 
    melt(., id.vars = 'Escapement') %>% 
    ggplot( aes(Escapement/1000, value, lty=variable))+geom_line()+
    xlab('Escapement (1,000)')+ylab('Probability')+
    theme(legend.justification=c(1,0), legend.position=c(1,.5), 
          legend.key = element_blank(),legend.title=element_blank())
  
  ggsave("figures/0.7.AR.png", dpi=200, width=8, height=5, units='in')
  
  
  Y %>% 
    dplyr::select(Escapement, Optimal_Yield0.8 = oy_0.8, 
                  Overfishing0.8 = of_0.8, Optimal_Recruitment0.8 = or_0.8) %>% 
    melt(., id.vars = 'Escapement') %>% 
    ggplot(aes(Escapement/1000, value, lty=variable))+geom_line()+
    xlab('Escapement (1,000)')+ylab('Probability')+
    theme(legend.justification=c(1,0), legend.position=c(1,.5), 
          legend.key = element_blank(),legend.title=element_blank())
  
  ggsave("figures/0.8.AR.png", dpi=200, width=8, height=5, units='in')
  
  
  Y %>% 
    dplyr::select(Escapement, Optimal_Yield0.9 = oy_0.9, 
                  Overfishing0.9 = of_0.9, Optimal_Recruitment0.9 = or_0.9) %>% 
    melt(., id.vars = 'Escapement')  %>% 
    ggplot(aes(Escapement/1000, value, lty=variable))+geom_line()+
    xlab('Escapement (1,000)')+ylab('Probability')+
    theme(legend.justification=c(1,0), legend.position=c(1,.5), 
          legend.key = element_blank(),legend.title=element_blank())
  
  ggsave("figures/0.9.AR.png", dpi=200, width=8, height=5, units='in')
  
  
  
  Y %>% 
    dplyr::select(Escapement, OY0.9 = oy_0.9, OY0.8 = of_0.8, OR0.9=or_0.9, 
                  OR0.8 = or_0.8, OF0.9 = of_0.9, OF0.8 = of_0.8) %>% 
    melt(., id.vars = 'Escapement')  %>% 
    mutate(sra = ifelse(grepl("OY0",variable), "Yield Profile",
                        ifelse(grepl("OR0",variable), "Recruitment Profile", "Overfishing Profile")),
           max_pct = ifelse(grepl('0.8', variable), 0.8,0.9)) %>% 
    ggplot(aes(Escapement, value, linetype = factor(max_pct)))+ 
    geom_rect(aes(xmin = LowerB, xmax = UpperB, ymin = 0, ymax = 1),
              inherit.aes = FALSE, fill = "grey80", alpha = 0.3)+geom_line()+
    xlab('Escapement (S)')+
    scale_x_continuous(labels = comma, breaks = seq(0, 350000, 100000))+
    scale_linetype_discrete(name = "Percent of Max.")+
    facet_grid(sra ~ .) + 
    theme(legend.key = element_blank(),legend.justification=c(0,0), legend.position=c(.65,.35),
          legend.background = element_rect(fill=alpha('white', 0.0)))+
    scale_y_continuous("Probability", breaks = seq(0, 1, 0.2), limits = c(0, 1))
  
  ggsave("figures/Model_3/0.8_0.9.png", dpi=200, width=7, height=6, units='in')						  
  
  
  
  
  ggplot(qm, aes(Escapement, Median))+geom_line(size=1)+
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.15)+
    geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.15)+ xlab('Escapement (S)')+
    ylab('Expected Yield')+scale_y_continuous(labels = comma)+
    scale_x_continuous(labels = comma,breaks = seq(0, 300000, 50000), limits = c(0,300000))+
    geom_vline(xintercept = LowerB,linetype = "longdash" )+geom_vline(xintercept = UpperB ,linetype = "longdash")
  
  ggsave("figures/Model_3/expected_sustained_yield.png", dpi=200, width=8, height=5, units='in')
}
#Run function
profile(i=10,z=500,xa.start=0, xa.end=700, coda)#can change i,z, xa.start, xa.end



