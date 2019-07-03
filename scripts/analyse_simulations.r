---
title: Analysing simulations
author: Gibran Hemani
---

Simulations run using this model [https://photos.google.com/photo/AF1QipNbxPRZyr7iODKRwuwrEmYALTpxnyudq7xovC7v](https://photos.google.com/photo/AF1QipNbxPRZyr7iODKRwuwrEmYALTpxnyudq7xovC7v)



library(tidyverse)
load("../results/sim8_summary.rdata")
# simres <- subset(simres, est != "Outliers removed")
# simres$est[simres$est == "Outliers removed"] <- "Outliers removed (candidates)"

## Main plots


simres$prop_invalid_bin <- cut(simres$prop_invalid, breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), ordered=TRUE)
hist(simres$prop_invalid)
table(simres$prop_invalid_bin)


simresu1 <- subset(simres, set == "U1")
dim(simresu1)

sumu1 <- simresu1 %>% 
	group_by(method, prop_invalid_bin, set) %>%
	summarise(
		n=n(),
		psig=sum(pval < 1e-3, na.rm=TRUE) / n(),
		bias=sum(pdiff < 0.05, na.rm=TRUE) / n(),
		isq=mean(Isq)
	)

sumu1b <- simresu1 %>% 
	group_by(method, prop_invalid_bin, set, bxy) %>%
	summarise(
		n=n(),
		psig=sum(pval < 1e-3, na.rm=TRUE) / n(),
		bias=sum(pdiff < 0.05, na.rm=TRUE) / n(),
		isq=mean(Isq)
	)

sumu1 %>% 
filter(!grepl("Hybrid", method)) %>%
ggplot(aes(y=bias, x=prop_invalid_bin, group=as.factor(method)), colour=as.factor(method)) +
geom_point(aes(colour=as.factor(method))) +
geom_line(aes(colour=as.factor(method))) +
scale_colour_brewer(type="qual") +
labs(x="Proportion of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection") +
ylim(c(0,1))


sumu1b %>% 
filter(!grepl("Hybrid", method)) %>%
ggplot(aes(y=psig, x=prop_invalid_bin, group=as.factor(method)), colour=as.factor(method)) +
geom_point(aes(colour=as.factor(method))) +
geom_line(aes(colour=as.factor(method))) +
facet_grid(. ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Proportion of SNPs with pleiotropic effects", y="False discovery rate / power", colour="Method", shape="Outlier detection", linetype="Outlier detection")




simresu2 <- subset(simres, set == "U2")
dim(simresu2)

sumu2 <- simresu2 %>% 
	group_by(method, prop_invalid_bin, set, directional_bias) %>%
	summarise(
		n=n(),
		psig=sum(pval < 1e-6, na.rm=TRUE) / n(),
		bias=sum(pdiff < 0.05, na.rm=TRUE) / n(),
		isq=mean(Isq)
	)

sumu2b1 <- simresu2 %>% 
	subset(bxy == 0) %>%
	group_by(method, prop_invalid_bin, set, bxy, directional_bias) %>%
	summarise(
		n=n(),
		psig=sum(pval < 0.05, na.rm=TRUE) / n(),
		bias=sum(pdiff < 0.05, na.rm=TRUE) / n(),
		isq=mean(Isq)
	)
sumu2b2 <- simresu2 %>% 
	subset(bxy != 0) %>%
	group_by(method, prop_invalid_bin, set, bxy, directional_bias) %>%
	summarise(
		n=n(),
		psig=sum(pval < 1e-6, na.rm=TRUE) / n(),
		bias=sum(pdiff < 0.05, na.rm=TRUE) / n(),
		isq=mean(Isq)
	)
sumu2b <- bind_rows(sumu2b1, sumu2b2)

sumu2 %>% 
filter(!grepl("Hybrid", method), !is.na(prop_invalid_bin)) %>%
ggplot(aes(y=bias, x=prop_invalid_bin, group=as.factor(method)), colour=as.factor(method)) +
geom_point(aes(colour=as.factor(method))) +
geom_line(aes(colour=as.factor(method))) +
facet_grid(. ~ directional_bias) +
scale_colour_brewer(type="qual") +
labs(x="Proportion of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")

sumu2b %>% 
filter(!grepl("Hybrid", method), !is.na(prop_invalid_bin)) %>%
ggplot(aes(y=psig, x=prop_invalid_bin, group=as.factor(method)), colour=as.factor(method)) +
geom_point(aes(colour=as.factor(method))) +
geom_line(aes(colour=as.factor(method))) +
facet_grid(directional_bias ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")


simresu3 <- subset(simres, set == "U3")
dim(simresu3)

sumu3 <- simresu3 %>% 
	group_by(method, prop_invalid_bin, set, bxy) %>%
	summarise(
		n=n(),
		psig=sum(pval < 0.05, na.rm=TRUE) / n(),
		bias=sum(pdiff < 0.05, na.rm=TRUE) / n(),
		isq=mean(Isq)
	)

sumu3 %>% 
filter(!grepl("Hybrid", method), !is.na(prop_invalid_bin)) %>%
ggplot(aes(y=bias, x=prop_invalid_bin, group=as.factor(method)), colour=as.factor(method)) +
geom_point(aes(colour=as.factor(method))) +
geom_line(aes(colour=as.factor(method))) +
facet_grid(. ~ bxy) +
scale_colour_brewer(type="qual") +
labs(x="Proportion of SNPs with pleiotropic effects", y="Proportion of estimates that are biased", colour="Method", shape="Outlier detection", linetype="Outlier detection")








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


temp2$lab1 <- "Null model"
temp2$lab1[temp2$bxy == 0.2] <- "Causal model"

templ <- gather(subset(temp2, select=c(lab1, est, nu, psig, bias, outliers_known2)), "key", "value", psig, bias)
templ$key[templ$key == "bias"] <- "Proportion with bias"
templ$key[templ$key == "psig"] <- "Proportion detected at p < 0.05"
ggplot(templ, aes(y=value, x=nu)) +
geom_point(aes(colour=as.factor(est), shape=outliers_known2)) +
geom_line(aes(colour=as.factor(est), linetype=outliers_known2)) +
# facet_grid(key ~ lab1, scale="free_y") +
facet_wrap(~ key + lab1, scale="free_y") +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of simulations", colour="Method", shape="Outlier detection", linetype="Outlier detection")
ggsave("../images/allsims.pdf")

templs <- subset(templ, (outliers_known2 == "Pleiotropy detected" & est %in% c("Outliers adjusted", "Outliers removed (all)")) | est == "Raw")
templs
templs$est[templs$est == "Outliers removed (all)"] <- "Outliers removed"
ggplot(templs, aes(y=value, x=nu)) +
geom_point(aes(colour=as.factor(est))) +
geom_line(aes(colour=as.factor(est))) +
# facet_grid(key ~ lab1, scale="free_y") +
facet_wrap(~ key + lab1, scale="free_y") +
scale_colour_brewer(type="qual") +
labs(x="Number of SNPs with pleiotropic effects", y="Proportion of simulations", colour="Method")
ggsave("../images/allsims_simple.pdf")

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


