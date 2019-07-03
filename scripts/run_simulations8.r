library(tryx)
library(parallel)
library(tidyverse)
library(simulateGP)


param <- bind_rows(
	# U1
	expand.grid(
		set = "U1",
		nid = 10000,
		ngx = c(5, 10, 15, 20, 25, 30),
		ngu1 = c(5, 10, 15, 20, 25, 30),
		ngu2 = 30,
		nu2 = 1,
		ngu3 = 30,
		vgx = 0.2,
		vgu1 = 0.6,
		vgu2 = 0.2,
		vgu3 = 0.2,
		bxy = c(0, 0.1),
		bu1x = 0.6,
		bu1y = 0.4,
		bxu3 = 0.3,
		bu3y = 0,
		vgxu2 = 0,
		vu2y = 0,
		ngxu3 = 0,
		vgxu3 = 0,
		mininum_instruments = 7,
		instrument_threshold = "bonferroni",
		outlier_threshold = 0.05,
		outliers_known = c("detected", "all"),
		directional_bias = FALSE,
		simr = c(1:1000)
	),

	# U2
	expand.grid(
		set = "U2",
		nid = 10000,
		ngx = 30,
		ngu1 = 1,
		ngu2 = 30,
		nu2 = c(5, 10, 15, 20, 25, 30),
		ngu3 = 30,
		vgx = 0.2,
		vgu1 = 0.6,
		vgu2 = 0.2,
		vgu3 = 0.2,
		bxy = c(0, 0.1),
		bu1x = 0.6,
		bu1y = 0,
		bxu3 = 0.3,
		bu3y = 0,
		vgxu2 = 0.3,
		vu2y = 0.4,
		ngxu3 = 0,
		vgxu3 = 0,
		mininum_instruments = 10,
		instrument_threshold = "bonferroni",
		outlier_threshold = 0.05,
		outliers_known = c("detected", "all"),
		directional_bias = c(FALSE, TRUE),
		simr = c(1:1000)
	),

	# U3
	expand.grid(
		set = "U3",
		nid = 10000,
		ngx = 30,
		ngu1 = 30,
		ngu2 = 30,
		nu2 = 1,
		ngu3 = 30,
		vgx = 0.2,
		vgu1 = 0.6,
		vgu2 = 0.2,
		vgu3 = 0.2,
		bxy = c(0, 0.1),
		bu1x = 0.6,
		bu1y = 0,
		bxu3 = 0.3,
		bu3y = 0.4,
		vgxu2 = 0.2,
		vu2y = 0.1,
		ngxu3 = c(1, 10, 20, 25),
		vgxu3 = 0.2,
		mininum_instruments = 10,
		instrument_threshold = "bonferroni",
		outlier_threshold = 0.05,
		outliers_known = c("detected", "all"),
		directional_bias = FALSE,
		simr = c(1:1000)
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
		sim <- tryx.simulate(
			nid = param$nid[i],
			ngx = param$ngx[i],
			ngu1 = param$ngu1[i],
			ngu2 = param$ngu2[i],
			nu2 = param$nu2[i],
			ngu3 = param$ngu3[i],
			vgx = param$vgx[i],
			vgu1 = param$vgu1[i],
			vgu2 = param$vgu2[i],
			vgu3 = param$vgu3[i],
			bxy = param$bxy[i],
			bu1x = param$bu1x[i],
			bu1y = param$bu1y[i],
			bxu3 = param$bxu3[i],
			bu3y = param$bu3y[i],
			vgxu2 = param$vgxu2[i],
			vu2y = param$vu2y[i],
			mininum_instruments = param$mininum_instruments[i],
			instrument_threshold = param$instrument_threshold[i],
			outlier_threshold = param$outlier_threshold[i],
			outliers_known = param$outliers_known[i],
			directional_bias = param$directional_bias[i]
		) %>% tryx.sig()
		res <- sim %>% 
			tryx.analyse()
		res$mvres <- sim$mvres
		res$simulation <- sim$simulation
		res
	})
	if(class(out) == "try-error") return(NULL)
	return(out)
}, mc.cores=ncores)

save(l, param, file=paste0("../results/scratch/sim8_", chunk, ".rdata"))
