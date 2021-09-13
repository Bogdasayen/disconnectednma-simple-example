# Code to implement reference prediction and ALM
# Howard Thom 13th September 2021

# Load necessary libraries
require(R2OpenBUGS)

# Load necessary data
source("load.data.R")

# Global option to use random study effects or not
# Paper establishes that random study effects are preferred as base case
random.effects <- TRUE

# BUGS options
n.chains <- 2 # Code below (setting initial values) only set up for 2 but can be extended
num.sims <- 1000 #  Recommend at least 10000 and ideally use n.thin = 2 or more
burn.in <- 1000 * n.chains # Recommend at least 30000
do.debug = FALSE # Set this to TRUE to see OpenBUGS running and check errors

##############################################################################
## Load models    ############################################################
##############################################################################

# Standard independent baseline models with fixed and random study effects
# The following are the NICE TSD2 1c (RE) and 1d (FE) functions
source("independent.baselines.model.R")

# Reference prediction with fixed study effects to include single-arm studies
# The following use fixed effects with (up to) 3 covariates on baseline to include single arm studies and keeps RCT network otherwise separate from single-arm studies
source("reference.prediction.model.single.1cov.fe.R")


# Reference prediction with fixed study effects to include disconnected RCTs
# The following uses fixed effects with (up to) 1 covariate on baseline 
source("reference.prediction.model.disc.1cov.fe.R")

# Random effects models for reference prediction and ALM have to use a model file
# format in order to use the dnorm(,)I(,) syntax to truncate the prior on sd
# in disconnected evidence

# Reference prediction with random study effects to include single-arm studies
model.file.reference.prediction.single.1cov.re <- "reference.prediction.model.single.1cov.re.txt"
# Reference prediction with random study effectst to include disconnected RCTs
model.file.reference.prediction.disc.1cov.re <- "reference.prediction.model.disc.1cov.re.txt"


# ALM models fixed study effects to include single-arm studies
source("alm.model.single.fe.R")
# ALM models fixed study effects to include disconnected RCTs
source("alm.model.disc.fe.R")

# ALM models random study effects to include single-arm studies
model.file.alm.single.re <- "alm.model.single.re.txt"
# ALM models random study effects to include disconnected RCTs
model.file.alm.disc.re <- "alm.model.disc.re.txt"

##############################################################################
## Initial values ############################################################
##############################################################################

# Initial values are shared across all models
# Some values are simply unused by some models
inits1<-list(d = c(NA, rep(0.5,b.data$nt-1)), mu=rep(0.5,b.data$ns), sd = 1,
             mu.single = rep(-0.5, b.data$ns.single), m = 0.1, sd.m = 1, beta = 0.1,
             mu.base = rep(-0.5, b.data$ns.base), beta.base = 0.1,
             mu.disc = matrix(-0.5, nrow=b.data$ns.disc, ncol=max(b.data$na.disc)),
             sd = 1) # mu.disc is a vector in ALM and matrix in RP
inits2<-list(d = c(NA, rep(-0.5,b.data$nt-1)), mu=rep(-0.5,b.data$ns), sd = 0.5,
             mu.single = rep(0.5, b.data$ns.single), m = 0.5, sd.m = 0.5, beta = 0.25,
             mu.base = rep(0.5, b.data$ns.base), beta.base = 0.25,
             mu.disc = matrix(0.5, nrow=b.data$ns.disc, ncol=max(b.data$na.disc)),
             sd = 0.5) # mu.disc is a vector in ALM and matrix in RP
bugs.inits<-list(inits1,inits2)


##############################################################################
## Run independent baseline models ###########################################
##############################################################################


# Independent baseline fixed study effects
if(!random.effects) {
  bugs.object.independent.baselines<- 
    bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("mu","d"),model=model.independent.baseline.fe,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=do.debug)
    
}

# Independent baseline random study effects
if(random.effects) {
  bugs.object.independent.baselines <-
    bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("mu","d","sd"),model=model.independent.baseline.re,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=2,debug=do.debug)
    
}

##############################################################################
## Match studies for ALM #####################################################
##############################################################################


# Vectors of RCTs matched to the disconnected RCTs and single-arm studies
# Used for plug-in models
matched.rct.disc <- rep(NA,b.data$ns.disc)
matched.rct.single <- rep(NA, b.data$ns.single)

# For ALM, need to extract the plug in for the closest matching study
# Need separate data structures for single-arm and disconnected analyses for the plug-in models
b.data.disc <- b.data.single <- b.data
for(i.disc in 1:b.data$ns.disc)
{
  # Use euclidean distance to choose closest arm match
  distance <- rep(NA, b.data$ns)
  for(i.rct in 1:b.data$ns)
  {
    # Distance between unweighted average of arms of each trial
    distance[i.rct] <- dist(rbind(
      mean(b.data$x.disc[i.disc,1:b.data$na.disc[i.disc]]),
      mean(b.data$x[i.rct,1:b.data$na[i.rct]])
    ))
  }
  min.rct<-which.min(distance)
  matched.rct.disc[i.disc]<-min.rct
}
for(i.single in 1:b.data$ns.single)
{
  # Use euclidean distance to choose closest arm match
  distance <- rep(NA, b.data$ns)
  for(i.rct in 1:b.data$ns)
  {
    # Distance between unweighted average of arms of each trial
    distance[i.rct] <- dist(rbind(
      b.data$x.single[i.single],
      mean(b.data$x[i.rct,1:b.data$na[i.rct]])
    ))
  }
  min.rct<-which.min(distance)
  matched.rct.single[i.single]<-min.rct
}
b.data.single$matched.rct <- matched.rct.single
b.data.disc$matched.rct <- matched.rct.disc

# Take the mu from the matched RCT
# Adding NA at end to avoid confusion between vectors and scalars
b.data.single$mu.plugin.mean<-
  c(bugs.object.independent.baselines$summary[matched.rct.single,"mean"],NA)
b.data.single$mu.plugin.prec<-
  1/c(bugs.object.independent.baselines$summary[matched.rct.single,"sd"],NA)^2
b.data.disc$mu.plugin.mean<-
  c(bugs.object.independent.baselines$summary[matched.rct.disc,"mean"],NA)
b.data.disc$mu.plugin.prec<-
  1/c(bugs.object.independent.baselines$summary[matched.rct.disc,"sd"],NA)^2

##############################################################################
## Run reference prediction and ALM ##########################################
##############################################################################

# And informative priors for sd.disc from the connected RCTs
# This is for both plug-in and reference prediction models
if(random.effects) {
  b.data.disc$sd.connected.mean<-b.data.single$sd.connected.mean<-b.data$sd.connected.mean<-
    bugs.object.independent.baselines$summary["sd","mean"]
  b.data.disc$sd.connected.tau<-b.data.single$sd.connected.tau<-b.data$sd.connected.tau<-
    1/bugs.object.independent.baselines$summary["sd","sd"]^2
}

# Reference prediction with 1 covariate for disconnected studies
if(!random.effects) {
  bugs.object.reference.prediction.disc.fe <- 
    bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("d"),model=model.reference.prediction.disc.1cov.fe,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=do.debug)
    
} else {
  bugs.object.reference.prediction.disc.re <-
    bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("d"),model.file= model.file.reference.prediction.disc.1cov.re,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=2,debug=do.debug)
    
}

# Reference prediction with 1 covariate for single-arm studies
if(!random.effects) {
  bugs.object.reference.prediction.single.fe <- 
    bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("d"),model=model.reference.prediction.single.1cov.fe,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=2,n.thin=1,debug=do.debug)

} else {
  bugs.object.reference.prediction.single.re <- 
    bugs(data=b.data,inits=bugs.inits,parameters.to.save=c("d"),model.file=model.file.reference.prediction.single.1cov.re,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=2,n.thin=2,debug=do.debug)
    
}


# ALM disconnected RCTs
if(!random.effects) {
  bugs.object.alm.disc.fe <- 
    bugs(data=b.data.disc,inits=bugs.inits,parameters.to.save=c("d"),model=model.alm.disc.fe,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=do.debug)
    
} else {
  bugs.object.alm.disc.re <-
    bugs(data=b.data.disc,inits=bugs.inits,parameters.to.save=c("d"),model.file=model.file.alm.disc.re,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=2,debug=do.debug)
    
}

# ALM single-arm studies
if(!random.effects) {
  bugs.object.alm.single.fe <- 
    bugs(data=b.data.single,inits=bugs.inits,parameters.to.save=c("d"),model=model.alm.single.fe,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=do.debug)
    
} else {
  bugs.object.alm.single.re <- 
    bugs(data=b.data.single,inits=bugs.inits,parameters.to.save=c("d"),model.file=model.file.alm.single.re,
         clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=2,debug=do.debug)
    
}


