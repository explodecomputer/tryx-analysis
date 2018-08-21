library(tidyverse)
load("../results/sim3_summary.rdata")
simres$est[simres$est == "Outliers removed"] <- "Outliers removed (candidates)"

## Main plots

temp2 <- filter(simres, bxy >= 0) %>%
	group_by(method, nu, bxy, nid) %>%
	summarise(
		n=n(),
		psig=sum(pval < 0.05) / n(),
		bias=sum(pdiff < 0.05) / n(),
		isq=mean(Isq),
		outliers_known2=first(outliers_known2),
		est=first(est)
	)

ggplot(temp2, aes(y=psig, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates with p < 0.05", colour="Method", shape="Outlier detection", linetype="Outlier detection")
ggsave("../images/sig1.pdf")

ggplot(temp2, aes(y=bias, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")
ggsave("../images/bias1.pdf")


temp3 <- filter(simres, bxy >= 0) %>%
	group_by(method, nu, bxy) %>%
	summarise(pow=sum(pval < 0.05)/n(), est=first(est), outliers_known2=first(outliers_known2))
ggplot(temp3, aes(y=pow, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
facet_grid(. ~ bxy) +
scale_colour_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates with p < 0.05", colour="Method", shape="Outlier detection", linetype="Outlier detection")
ggsave("../images/sig2.pdf")


temp3 <- filter(simres, bxy >= 0) %>%
	group_by(method, nu) %>%
	summarise(bias=sum(pdiff < 0.05)/n(), est=first(est), outliers_known2=first(outliers_known2))
ggplot(temp3, aes(y=bias, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
scale_colour_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")
ggsave("../images/bias2.pdf")


temp3 <- filter(simres, bxy >= 0) %>%
	group_by(method, bxy) %>%
	summarise(pow=sum(pval < 0.05)/n())
ggplot(temp3, aes(y=pow, x=method)) +
geom_bar(aes(fill=as.factor(bxy)), stat="identity", position="dodge") +
scale_fill_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
labs(x="", y="Proportion of estimates with p < 0.05", fill="Simulated\ncausal effect") +
geom_hline(yintercept=0.05, linetype="dotted")
ggsave("../images/sig3.pdf")


temp3 <- filter(simres, bxy >= 0) %>%
	group_by(method) %>%
	summarise(bias=sum(pdiff < 0.05)/n())
ggplot(temp3, aes(y=bias, x=method)) +
geom_bar(aes(fill=method), stat="identity") +
scale_fill_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position="none") +
geom_hline(yintercept=0.05, linetype="dotted") +
labs(y="Proportion of estimates that are biased", x="")
ggsave("../images/bias3.pdf")





## Other plots

ggplot(simres, aes(y=diff, x=method)) +
geom_boxplot() +
facet_grid(nid ~ bxy) +
theme(axis.text.x=element_text(angle=90))

ggplot(simres, aes(y=pdiff, x=method)) +
geom_boxplot() +
facet_grid(nid ~ bxy) +
theme(axis.text.x=element_text(angle=90))

ggplot(simres, aes(y=-log10(pdiff), x=method)) +
geom_boxplot() +
facet_grid(nid ~ bxy) +
theme(axis.text.x=element_text(angle=90))



temp2 <- group_by(simres, method, bxy, nid, nu) %>%
summarise(pow=sum(pval < 1e-5)/n(), isq=mean(Isq))

ggplot(temp2, aes(y=pow, x=nu)) +
geom_line(aes(colour=method)) +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual")

temp3 <- group_by(simres, method, bxy, nid) %>%
summarise(pow=sum(pval < 0.05)/n(), isq=mean(Isq))

ggplot(temp3, aes(y=pow, x=method)) +
geom_point() +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
geom_hline(yintercept = 0.05, linetype="dotted") +
theme(axis.text.x=element_text(angle=90))

ggplot(temp3, aes(y=isq, x=method)) +
geom_point() +
facet_grid(nid ~ bxy) +
scale_colour_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90))


