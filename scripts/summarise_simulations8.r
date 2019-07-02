library(dplyr)

load("../results/sim8.rdata")

simres <- lapply(1:length(l), function(x){
	a <- l[[x]]$estimates
	a$sim <- x
	a$no_outlier_flag <- l[[x]]$simulation$no_outlier_flag
	a$n_instruments <- l[[x]]$simulation$n_instruments
	a$n_valid_instruments <- l[[x]]$simulation$n_valid_instruments
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
				subset(x, est == "Raw"),
				subset(x, est == "Raw"),
				subset(x, est == "Raw")
			)
			x$est <- c("Outliers removed (all)", "Outliers removed (candidates)", "Raw", "Outliers adjusted", "Mixed", "Multivariable MR")
		}
		x
	})

simres$mediator <- simres$bu3y * simres$bxu3
simres$prop_invalid <- (simres$n_instruments - simres$n_valid_instruments) / simres$n_instruments

dim(simres)
table(simres$est)



# Rename the method
simres$method <- simres$est
simres$method[simres$method == "Outliers removed (all)"] <- "Removed (detected outliers)"
simres$method[simres$method == "Outliers removed (candidates)"] <- "Removed (detected candidates)"
simres$method[simres$method == "Outliers adjusted" & simres$outliers_known == "all"] <- "Adjusted (all variants)"
simres$method[simres$method == "Outliers adjusted" & simres$outliers_known == "detected"] <- "Adjusted (outliers)"
simres$method[simres$method == "Multivariable MR" & simres$outliers_known == "all"] <- "MVMR (all variants)"
simres$method[simres$method == "Multivariable MR" & simres$outliers_known == "detected"] <- "MVMR (outliers)"
simres$method[simres$method == "Mixed" & simres$outliers_known == "all"] <- "Hybrid (all variants)"
simres$method[simres$method == "Mixed" & simres$outliers_known == "detected"] <- "Hybrid (outliers)"

table(simres$method)

# how different from the simulated causal effect is the estimated causal effect?
simres$diff <- simres$b - simres$bxy
simres$pdiff <- pt(abs(simres$diff)/simres$se, simres$nsnp, lower.tail=FALSE)
save(simres, file="../results/sim8_summary.rdata")

