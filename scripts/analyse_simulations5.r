library(dplyr)
library(ggplot2)
load("../results/sim5_summary.rdata")

simr$nu <- simr$nu1 + simr$nu2
simr$mediator <- simr$bu3y * simr$bxu3
table(simr$mediator, simr$bxy)


temp2 <- filter(simr, bxy >= 0) %>%
	group_by(method, nu, bxy, mediator, nid) %>%
	summarise(
		n=n(),
		psig=sum(pval < 0.05, na.rm=TRUE) / n(),
		bias=sum(pdiff < 0.05, na.rm=TRUE) / n(),
		eff=mean(b, na.rm=T),
		isq=mean(Isq, na.rm=TRUE),
		outliers_known2=first(outliers_known2),
		est=first(est)
	)

ggplot(temp2 %>% filter(mediator == 0), aes(y=psig, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates with p < 0.05", colour="Method", shape="Outlier detection", linetype="Outlier detection")

ggplot(temp2 %>% filter(mediator == 0), aes(y=bias, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")



ggplot(temp2 %>% filter(mediator != 0), aes(y=psig, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates with p < 0.05", colour="Method", shape="Outlier detection", linetype="Outlier detection")

ggplot(temp2, aes(y=bias, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
facet_grid(mediator == 0 ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")

ggplot(temp2, aes(y=eff, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
facet_grid(mediator == 0 ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")


ggplot(simr %>% filter(est %in% c("Multivariable MR", "Outliers adjusted")), aes(x=nu, y=b)) +
geom_boxplot(alpha=0.4, aes(fill=est)) +
facet_grid(mediator == 0 ~ bxy) 