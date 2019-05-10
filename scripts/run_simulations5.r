library(tryx)
library(parallel)
library(tidyverse)

param <- bind_rows(
	expand.grid(
		nid = c(5000),
		bxy = 0.1,
		bu3y = 0.3,
		bxu3 = 0.3,
		nu1 = c(1, 5, 10, 14),
		nu2 = c(1, 5, 10, 14),
		outliers_known = c(FALSE),
		simr = c(1:2000)
	),
	expand.grid(
		nid = c(5000),
		bxy = 0,
		bu3y = 0.45,
		bxu3 = 0.45,
		nu1 = c(1, 5, 10, 14),
		nu2 = c(1, 5, 10, 14),
		outliers_known = c(FALSE),
		simr = c(1:2000)
	),
	expand.grid(
		nid = c(5000),
		bxy = c(0, 0.2),
		bu3y = 0,
		bxu3 = 0,
		nu1 = c(1, 5, 10, 14),
		nu2 = c(1, 5, 10, 14),
		outliers_known = c(FALSE),
		simr = c(1:2000)
	)
)
param$sim <- 1:nrow(param)
dim(param)


args <- commandArgs(T)
chunk <- as.numeric(args[1])
chunksize <- as.numeric(args[2])
ncores <- as.numeric(args[3])
a <- (chunk - 1) * chunksize + 1
b <- min(chunk * chunksize, nrow(param))
message("Running ", a, " to ", b)
param <- param[a:b, ]


l <- mclapply(1:nrow(param), function(i)
{
	set.seed(param$sim[i])
	out <- try({
		sim <- tryx.simulate(param$nid[i], param$nu1[i], param$nu2[i], param$bu3y[i], param$bxu3[i], param$bxy[i], outliers_known=param$outliers_known[i]) %>% tryx.sig()
		res <- sim %>% 
			tryx.analyse()
		res$mvres <- sim$mvres
		res
	})
	if(class(out) == "try-error") return(NULL)
	return(out)
}, mc.cores=ncores)

save(l, param, file=paste0("../results/sim5_", chunk, ".rdata"))
