library(dplyr)

load("../results/SIM5.rdata")

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
				subset(x, est == "Raw"),
				subset(x, est == "Raw")
			)
			x$est <- c("Outliers removed (all)", "Outliers removed (candidates)", "Raw", "Outliers adjusted")
		}
		x
	})

simres$nu <- simres$nu1 + simres$nu2
simres$mediator <- simres$bu3y * simres$bxu3
dim(simres)
table(simres$est)

# Rename the method
simres$outliers_known2 <- "Pleiotropy known"
simres$outliers_known2[!simres$outliers_known] <- "Pleiotropy detected"
simres$method <- paste0(simres$est, " (", simres$outliers_known2, ")")
simres$method[simres$est == "Raw"] <- "Raw"
table(simres$method)

# how different from the simulated causal effect is the estimated causal effect?
simres$diff <- simres$b - simres$bxy
simres$pdiff <- pt(abs(simres$diff)/simres$se, simres$nsnp, lower.tail=FALSE)


mvres <- lapply(1:length(l), function(x){
	a <- l[[x]]$mvres
	if(class(a) == "list")
	{
		a$result$sim <- x
		return(a$result)
	} else {
		return(NULL)
	}
}) %>% bind_rows()

mvres <- inner_join(mvres, param, by="sim")

mvresx <- subset(mvres, exposure == "X")
mvresx$diff <- mvresx$b - mvresx$bxy
mvresx$pdiff <- pt(abs(mvresx$diff)/mvresx$se, mvresx$nsnp, lower.tail=FALSE)

mvresx$est <- "Multivariable MR"
mvresx$method <- "Multivariable MR (Pleiotropy detected)"
mvresx$outliers_known <- FALSE
mvresx$outliers_known2 <- "Pleiotropy detected"

simr <- bind_rows(simres, mvresx)

table(simr$est)
table(simr$method)
save(simres, mvres, mvresx, simr, file="../results/sim5_summary.rdata")


library(dplyr)
library(ggplot2)
load("../results/sim5_summary.rdata")


