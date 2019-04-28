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
		simr = c(1:1000)
	),
	expand.grid(
		nid = c(5000),
		bxy = 0,
		bu3y = 0.45,
		bxu3 = 0.45,
		nu1 = c(1, 5, 10, 14),
		nu2 = c(1, 5, 10, 14),
		outliers_known = c(FALSE),
		simr = c(1:1000)
	)
)
param$sim <- 1:nrow(param)
dim(param)



l <- mclapply(1:nrow(param), function(i)
{
	set.seed(i)
	out <- try({
		sim <- tryx.simulate(param$nid[i], param$nu1[i], param$nu2[i], param$bu3y[i], param$bxu3[i], param$bxy[i], outliers_known=param$outliers_known[i]) %>% tryx.sig()
		res <- sim %>% 
			tryx.analyse()
		res$mvres <- sim$mvres
		res
	})
	if(class(out) == "try-error") return(NULL)
	return(out)
}, mc.cores=16)

save(l, param, file="../results/sim5.rdata")

## Summarise simulations

simres <- lapply(1:length(l), function(x){
	a <- l[[x]]$estimates
	a$sim <- x
	a$no_outlier_flag <- l[[x]]$detection$no_outlier_flag
	return(a)
}) %>% bind_rows %>% inner_join(param) %>% filter(!is.na(est))
simres$no_outlier_flag[is.na(simres$no_outlier_flag)] <- FALSE
simres <- group_by(simres, sim) %>%
	do({
		x <- .
		if(x$no_outlier_flag[1])
		{
			x <- rbind(
				subset(x, est == "Raw"),
				subset(x, est == "Raw"),
				subset(x, est == "Raw")
			)
			x$est <- c("Outliers removed", "Raw", "Outliers adjusted")
		}
		x
	})

simres$nu <- simres$nu1 + simres$nu2

dim(simres)
table(simres$est)

# Rename the method
simres$outliers_known2 <- "Pleiotropy known"
simres$outliers_known2[!simres$outliers_known] <- "Pleiotropy detected"
simres$method <- paste0(simres$est, " (", simres$outliers_known2, ")")
simres$method[simres$est == "Raw"] <- "Raw"

# how different from the simulated causal effect is the estimated causal effect?
simres$diff <- simres$b - simres$bxy
simres$pdiff <- pt(abs(simres$diff)/simres$se, simres$nsnp, lower.tail=FALSE)


mvres <- lapply(1:length(l), function(x)){
	a <- l[[x]]$mvres
	if(class(a) == "list")
	{
		a$result$sim <- x
		return(a$result)
	} else {
		return(NULL)
	}
} %>% bind_rows()

save(simres, mvres file="../results/sim5_summary.rdata")