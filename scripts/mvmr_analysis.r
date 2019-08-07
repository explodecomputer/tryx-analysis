library(tidyverse)

load("../results/mvmr_simulations.rdata")

s <- m %>% group_by(method, sim, plei, nprox) %>%
summarise(n=sum(b!=0), b=sum(b), model=ifelse(first(plei), "Pleiotropic", "Simple"))

group_by(s, model, method, nprox) %>% summarise(est=mean(b), se=sd(b)/sqrt(n()), bias_perc=mean((b-0.4) / 0.4)*100, n_variables=mean(n), n_sims=n()) %>%
knitr::kable()



