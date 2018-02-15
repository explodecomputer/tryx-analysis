library(tryx)
library(parallel)
library(tidyverse)

param <- expand.grid(
	nid = c(5000),
	bxy = c(0.2, 0),
	nu1 = c(1, 5, 10, 14),
	nu2 = c(1, 5, 10, 14),
	outliers_known = c(TRUE, FALSE),
	simr = c(1:1000)
)
param$sim <- 1:nrow(param)



l <- mclapply(1:nrow(param), function(i)
{
	set.seed(i)
	out <- simulate.tryx(param$nid[i], param$nu1[i], param$nu2[i], param$bxy[i], outliers_known=param$outliers_known[i])
	out <- outlier_sig(out)
	return(tryx.analyse(out, plot=FALSE))
}, mc.cores=16)

save(l, param, file="../results/sim2.rdata")

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


save(simres, file="../results/sim2_summary.rdata")
