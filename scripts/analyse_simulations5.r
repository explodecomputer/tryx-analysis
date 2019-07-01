library(dplyr)
library(ggplot2)
load("../results/sim5_summary.rdata")
simr5 <- simr
simr5$sim <- paste0(5, simr5$sim)
load("../results/sim6_summary.rdata")
simr6 <- simr
simr6$sim <- paste0(6, simr6$sim)
load("../results/sim7_summary.rdata")
simr7 <- simr
simr7$sim <- paste0(7, simr7$sim)

simr <- simr7
# simr <- bind_rows(simr6, simr7)
# simr$directional_bias[is.na(simr$directional_bias)] <- FALSE

simr$nu <- simr$nu1 + simr$nu2
simr$mediator <- simr$bu3y * simr$bxu3
table(simr$mediator, simr$bxy)


temp2 <- filter(simr, bxy >= 0) %>%
	group_by(method, nu, bxy, mediator, nid, outliers_known) %>%
	summarise(
		n=n(),
		psig=sum(pval < 0.05, na.rm=TRUE) / n(),
		bias=sum(pdiff < 0.05, na.rm=TRUE) / n(),
		eff=mean(b, na.rm=T),
		isq=mean(Isq, na.rm=TRUE),
		est=first(est)
	)

ggplot(temp2 %>% filter(mediator == 0), aes(y=psig, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates with p < 0.05", colour="Method", shape="Outlier detection", linetype="Outlier detection")

ggplot(temp2 %>% filter(mediator == 0), aes(y=bias, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")



ggplot(temp2 %>% filter(mediator != 0), aes(y=psig, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates with p < 0.05", colour="Method", shape="Outlier detection", linetype="Outlier detection")

ggplot(temp2, aes(y=bias, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known)) +
facet_grid(mediator == 0 ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")

ggplot(temp2, aes(y=eff, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known)) +
facet_grid(mediator == 0 ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")


ggplot(simr %>% filter(est %in% c("Multivariable MR", "Outliers adjusted")), aes(x=nu, y=b)) +
geom_boxplot(alpha=0.4, aes(fill=est)) +
facet_grid(mediator == 0 ~ bxy) 