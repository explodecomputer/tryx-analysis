
## Summarise simulations

load("../results/sim5.rdata")

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
