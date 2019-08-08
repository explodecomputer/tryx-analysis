library(tidyverse)

load("../results/mvmr_simulations.rdata")
m <- subset(m, !is.na(nprox))

s <- m %>% group_by(method, sim, plei, nprox) %>%
summarise(n=sum(b!=0), b=sum(b), model=ifelse(first(plei), "Pleiotropic", "Simple"))

tab <- group_by(s, model, method, nprox) %>% summarise(est=mean(b), se=sd(b)/sqrt(n()), bias_perc=mean((b-0.4) / 0.4)*100, n_variables=mean(n), n_sims=n()) 

tab %>%
knitr::kable()

write.csv(tab, file="../results/mvmr_simulations_summary.csv")


