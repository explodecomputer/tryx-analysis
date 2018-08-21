library(TwoSampleMR)
library(dplyr)
ao <- available_outcomes()
a <- extract_instruments(subset(ao, grepl("Usual walking", trait))$id)
b <- extract_outcome_data(a$SNP, 2)
dat <- harmonise_data(a,b)
dat <- steiger_filtering(dat)
mr(dat)
mr_heterogeneity(dat)
dat %>% filter(steiger_dir) %>% mr
dat %>% filter(steiger_dir) %>% mr_heterogeneity

mr_forest_plot(mr_singlesnp(dat))
mr_leaveoneout_plot(mr_leaveoneout(dat))
filter(dat, SNP != "rs9972653") %>% mr

subset(dat, SNP == "rs9972653")
filter(dat, st)



a <- extract_instruments(subset(ao, grepl("Alcohol intake fr", trait))$id)
b <- extract_outcome_data(a$SNP, 2)
dat <- harmonise_data(a,b)
dat <- steiger_filtering(dat)
mr(dat)
mr_heterogeneity(dat)
dat %>% filter(steiger_dir) %>% mr
dat %>% filter(steiger_dir) %>% mr_heterogeneity

