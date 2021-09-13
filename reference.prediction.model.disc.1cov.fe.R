# Reference prediction models with fixed study effects to include disconnected RCTs


# Data are same as TSD format: ns, nt, na, r, n, t
# The data on baseline arms (from RCTs) are
# ns.base : number of baseline arms
# r.base : number of events in baseline arms
# n.base : number of patients in baseline arms
# x.base : Matrix of covariates with one row per covariate
# The data on disconnected networks are (as in standard TSD):
# ns.disc, nt.disc, na.disc, r.disc, n.disc, t.disc
# x.disc is matrix of covariates for disconnected RCTs.
# cov.index is the index of the covariates to use in x.disc and x.base

# NOTE: All x.base and x.disc must be defined. Suggest using mean covariate value when not reported.
# NOTE: Naively pools disconnected and connected RCT evidence on treatment effects (not sure enough evidence to fit hierarchical model)

# Reference prediction with fixed study effect model to include disconnected RCTs
# Binomial likelihood, logit link 
# Simultaneous baseline and treat effects model for multi-arm trials 
model.reference.prediction.disc.1cov.fe<-function()
{
	# Model for RCTs ###################################
	for(i in 1:ns){ # LOOP THROUGH STUDIES 
		mu[i] ~ dnorm(0,0.01) # random effect on baselines
		for (k in 1:na[i]) { # LOOP THROUGH ARMS 
			r[i,k] ~ dbin(p[i,k],n[i,k]) # Binomial likelihood 
			logit(p[i,k]) <- mu[i] + delta[i,k]
			delta[i,k]<-d[t[i,k]] - d[t[i,1]] # model for linear predictor 
			rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
			dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) 
			+ (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k]))) #Deviance contribution 
		} 
	resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial 
	} 

	# Model for baseline effects ###################################
	# Adapted from Program 1 of NICE DSU TSD 5
	for (i in 1:ns.base){ # LOOP THROUGH STUDIES
		r.base[i] ~ dbin(p.base[i],n.base[i]) # Likelihood
		logit(p.base[i]) <- mu.base[i] 
		mu.base[i] ~ dnorm(mu.base.mean[i],tau.m) # Random effects model
		mu.base.mean[i] <- m + (x.base[i]-x.base.mean)*beta.base # Prediction of mean effect # Log-odds of response
	}
	beta.base~dnorm(0,0.01) #0.298 #  vague prior for covariate effects
	m ~ dnorm(0,.01) # vague prior for mean
	var.m <- 1/tau.m # between-trial variance
	tau.m <- pow(sd.m,-2) # between-trial precision = (1/between-trial variance)
	sd.m ~ dunif(0,2) # vague prior for between-trial SD


	# Model for disconnected RCTs ###################################
	# Prevent feedback to baseline model (these two may not be necessary as mu.single is cut below
	m.cut<-cut(m)
	tau.m.cut<-cut(tau.m)
	beta.base.cut<-cut(beta.base)

	for(i in 1:ns.disc){ # LOOP THROUGH STUDIES 
		for (k in 1:na.disc[i]) { # LOOP THROUGH ARMS 
			r.disc[i,k] ~ dbin(p.disc[i,k],n.disc[i,k]) # binomial likelihood 
			logit(p.disc[i,k]) <- mu.disc[i,k]
			mu.disc[i,k] ~ dnorm(mu.disc.mean[i,k], tau.m.cut)
			mu.disc.mean[i,k] <- m.cut + delta.disc[i,k] + (x.disc[i,k]-x.base.mean)*beta.base.cut #+ x.disc[cov.index[2],i,k]*beta.base.cut[2] #+ x.disc[cov.index[3],i,k]*beta.base.cut[3]  # model for linear predictor 
			delta.disc[i,k]<-d[t.disc[i,k]] - d[t.disc[i,1]] # model for linear predictor 
			rhat.disc[i,k] <- p.disc[i,k] * n.disc[i,k] # expected value of the numerators 
			dev.disc[i,k] <- 2 * (r.disc[i,k] * (log(r.disc[i,k])-log(rhat.disc[i,k])) 
			+ (n.disc[i,k]-r.disc[i,k]) * (log(n.disc[i,k]-r.disc[i,k]) - log(n.disc[i,k]-rhat.disc[i,k]))) #Deviance contribution 
		} 
	resdev.disc[i] <- sum(dev.disc[i,1:na.disc[i]]) # summed residual deviance contribution for this trial 
	} 

	# Calculate deviance ###################################
	totresdev.disc<-sum(resdev.disc[]) # Total residual deviance for disconnected RCTs
	totresdev.rct <- sum(resdev[]) # Total residual deviance for RCTs
	totresdev<-totresdev.disc+totresdev.rct #Total Residual Deviance 
	
	# Priors for remaining parameters
	d[1]<-0 # treatment effect is zero for reference treatment 
	for (k in 2:nt){ d[k] ~ dnorm(0,.01) } # vague priors for treatment effects 
}



