#!/usr/bin/env Rscript
library(rstan)

args = commandArgs(trailingOnly=TRUE)
# 1. File containing summary statistics.
data.fn = args[1]
# 2. File to save output to.
out.fn = args[2]
# 3. Save to HDF5. If FALSE, save as RData file.
hdf5 = args[3]
# 4. If mcmc=TRUE, then do an MCMC run. Otherwise, do optimization.
mcmc = as.logical(args[4])
# 5. Number of MCMC iterations or number of draws for optimization.
niter = as.numeric(args[5])
# 6. Number of MCMC warmup iterations. Ignored for optimization.
warmup = as.numeric(args[6])
# 7. Number of MCMC chains. Ignored for optimization.
chains = as.numeric(args[7])
# 8. Number of cores (should probably be the same as the number of chains).
# Ignored for optimization.
cores = as.numeric(args[8])
# 9. mvpmm.stan file to use.
mvpmm.fn = args[9]

seed = 20171001

# Create data for Stan
data = read.table(gzfile(data.fn), header=TRUE)
numphenos = length(levels(data$code))
numvars = dim(data)[1] / numphenos
phenos = matrix(data$code, nrow=numvars, ncol=numphenos)[1,]
betas = matrix(data$BETA, nrow=numvars, ncol=numphenos)
ses = matrix(data$SE, nrow=numvars, ncol=numphenos)

stan.data <- list(
    N = numvars,
    M = numphenos,
    B = betas,
    SE = ses,
    K = 2
)

if (mcmc) {
	fit <- stan(
	  file = mvpmm.fn,
	  data = stan.data,
	  chains = chains,
	  warmup = warmup,
	  iter = niter,
	  cores = cores,
	  refresh = round(niter / 8)
	)
	# fit.summary = summary(fit, probs=c(0.025, 0.5, 0.975), digits_summary=5) 
} else {
	sm <- stan_model(file = mvpmm.fn)
	fit = optimizing(sm, data=stan.data, hessian=TRUE, as_vector=FALSE,
			 draws=niter, seed=seed)
	
	se = sqrt(diag(solve(-fit$hessian)))
	drawz = abs(fit$par$Omegacor[2,1] / se["L_Omega.1"])
	drawp = pnorm(drawz, lower.tail=F)
}

# Remove large data objects to save space and speed up loading of .RData file.
rm(data)
rm(stan.data)
rm(betas)
rm(ses)
save.image(out.fn)
