# Reference prediction fixed study effects to include single-arm studies

# Data are same as TSD format: ns, nt, na, r, n, t
# The data on baseline arms (from RCTs) are
# ns.base : number of baseline arms
# r.base : number of events in baseline arms
# n.base : number of patients in baseline arms
# x.base : Matrix of covariates with one row per covariate
# The data on single arm studies are:
# ns.single : number of single arm studies
# r.single : number of events in single arm studies
# n.single : number of patients in single arm studies
# t.single : treatment in single arm study
# x.single : matrix of covariates with one row per covariate
# cov.index : vector of numbers indicating which row of x.single/x.base to use for each covariate

# NOTE: All x.base and x.single must be defined. Suggest using mean covariate value when not reported.
# NOTE: Naively pools single arm and RCT evidence on treatment effects (not sure enough evidence to fit hierarchical model)



# Reference prediction fixed study effects with 1 covariate including single-arm studies
# Binomial likelihood, logit link 
# Simultaneous baseline and treatment effects model for multi-arm trials 
model.reference.prediction.single.1cov.fe<-function()
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
	beta.base~dnorm(0,0.01) #  vague prior for covariate effects
	m ~ dnorm(0,.01) # vague prior for mean
	var.m <- 1/tau.m # between-trial variance
	tau.m <- pow(sd.m,-2) # between-trial precision = (1/between-trial variance)
	sd.m ~ dunif(0,5) # vague prior for between-trial SD


	# Model for single-arm studies ###################################
	# Prevent feedback to baseline model (these two may not be necessary as mu.single is cut below
	m.cut<-cut(m)
	tau.m.cut<-cut(tau.m)
	beta.base.cut<-cut(beta.base)

	for(i in 1:ns.single)
	{
		r.single[i]~dbin(p.single[i],n.single[i])
		logit(p.single[i]) <- mu.single[i]
		mu.single[i] ~ dnorm(mu.single.mean[i], tau.m.cut)
		mu.single.mean[i] <- m.cut + delta.single[i] + (x.single[i]-x.base.mean)*beta.base.cut # model for linear predictor 
		delta.single[i]<-d[t.single[i]]
		rhat.single[i] <- p.single[i] * n.single[i] # expected value of the numerators 
		dev.single[i] <- 2 * (r.single[i] * (log(r.single[i])-log(rhat.single[i])) 
		+ (n.single[i]-r.single[i]) * (log(n.single[i]-r.single[i]) - log(n.single[i]-rhat.single[i]))) #Deviance contribution 
	}

	# Calculate deviance ###################################
	totresdev.single<-sum(dev.single[]) # Total residual deviance for single-arm studies
	totresdev.rct <- sum(resdev[]) # Total residual deviance for RCTs
	totresdev<-totresdev.single+totresdev.rct #Total Residual Deviance 
	
	# Priors for remaining parameters
	d[1]<-0 # treatment effect is zero for reference treatment 
	for (k in 2:nt){ d[k] ~ dnorm(0,.01) } # vague priors for treatment effects 
}



