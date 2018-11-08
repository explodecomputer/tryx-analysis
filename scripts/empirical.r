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
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################


rm(list=ls())


devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("WSpiller/RadialMR")
devtools::install_github("explodecomputer/tryx")


library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(TwoSampleMR)
library(MRInstruments)
library(RadialMR)
library(tryx)


#########################################################################################################################
#Example 1: systolic blood pressure and coronary artery disease
#########################################################################################################################

#Exposure
ao <- available_outcomes()
sbp <- extract_instruments(outcomes="UKB-a:360")

#Outcome
chd <-  extract_outcome_data(snps=sbp$SNP, outcomes = 7)

dat_sbp <- harmonise_data(sbp, chd)

############################################Radial MR#######################################################

dat_sbp2 <- dat_sbp[dat_sbp$mr_keep, ]
radial_mr <-format_radial(dat_sbp2$beta.exposure, dat_sbp2$beta.outcome, dat_sbp2$se.exposure, dat_sbp2$se.outcome, dat_sbp2$SNP)
radial_res_ivw <-ivw_radial(radial_mr, 0.05/nrow(dat_sbp),3,0.00001, summary = FALSE)
radial_plot <-plot_radial(radial_res_ivw,TRUE,FALSE,TRUE)
radial_funnel <- funnel_radial(radial_res_ivw,TRUE)



#############################################RESULT##########################################################
#TRYX
tryxscan_sbp <- tryx.scan(dat_sbp, use_proxies = TRUE)
tryxscan_sbp <- tryx.sig(tryxscan_sbp)
#tryx.network(tryxscan_sbp)
adj_sbp <- tryx.analyse(tryxscan_sbp)


#Multivariable MR

id_remove <- c("UKB-a:22", "UKB-a:24", "UKB-a:61", "UKB-a:62", "UKB-a:63", "UKB-a:131", "UKB-a:201", "UKB-a:218", 
               "UKB-a:222", "UKB-a:310", "UKB-a:359", "UKB-a:360", "UKB-a:434", "UKB-a:435", "UKB-a:436", "UKB-a:437", 
               "UKB-a:449", "UKB-a:450", "UKB-a:489", "UKB-a:490", "UKB-a:533", "UKB-a:534", "798")

ta_sbp <- tryx.analyse.mv(tryxscan_sbp, id_remove=id_remove)    


#MR
result<- mr(dat_sbp, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger <- mr_pleiotropy_test(dat_sbp)

#MR-remove outliers (all)
outlier <- tryxscan_ur$outliers
dat_outout <- dat_ur[!(dat_ur$SNP %in% outlier),]
res_outout <- mr(dat_outout, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger_ro <- mr_pleiotropy_test(dat_outout)

#MR-remove outlierS (candidate)
outlier <- tryxscan_ur$outliers
dat_ca_out <- dat_ur[!(dat_ur$SNP %in% adj_ur$adj$SNP),]
res_ca_out <- mr(dat_ca_out, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger_ca_out <- mr_pleiotropy_test(dat_ca_out)



#########################################################################################################################
#Example 2: urate and coronary heart disease
#########################################################################################################################

#Exposure
ur <- extract_instruments(outcomes=1055)

#Outcome
chd_ur <-  extract_outcome_data(snps=ur$SNP, outcomes = 7)

dat_ur <- harmonise_data(ur, chd_ur)

############################################Radial MR#######################################################

dat_ur2 <- dat_ur[dat_ur$mr_keep, ]
radial_mr <-format_radial(dat_ur2$beta.exposure, dat_ur2$beta.outcome, dat_ur2$se.exposure, dat_ur2$se.outcome, dat_ur2$SNP)
radial_res_ivw <-ivw_radial(radial_mr, 0.05/nrow(dat_ur2), summary = FALSE)
radial_plot <-plot_radial(radial_res_ivw,TRUE,FALSE,TRUE)
radial_funnel <- funnel_radial(radial_res_ivw,TRUE)


#############################################RESULT##########################################################
#TRYX
tryxscan_ur <- tryx.scan(dat_ur, use_proxies = TRUE)
tryxscan_ur <- tryx.sig(tryxscan_ur)
adj_ur <- tryx.analyse(tryxscan_ur)


#Multivariable MR
id_remove <- c("UKB-a:24", "UKB-a:218", "UKB-a:490", "UKB-a:222", "UKB-a:450", "UKB-a:437", "UKB-a:449", "UKB-a:310", 
               "UKB-a:435", "UKB-a:489", "UKB-a:61", "UKB-a:131", "UKB-a:359", "UKB-a:533", "UKB-a:22", "UKB-a:107", 
               "UKB-a:360", "UKB-a:434", "798", "1105", "UKB-a:201", "UKB-a:63", "id:1055")


ta_ur <- tryx.analyse.mv(tryxscan_ur, id_remove=id_remove)               


#MR
result<- mr(dat_ur, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger <- mr_pleiotropy_test(dat_ur)

#MR-remove outliers (all)
outlier <- tryxscan_ur$outliers
dat_outout <- dat_ur[!(dat_ur$SNP %in% outlier),]
res_outout <- mr(dat_outout, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger_ro <- mr_pleiotropy_test(dat_outout)

#MR-remove outlierS (candidate)
outlier <- tryxscan_ur$outliers
dat_ca_out <- dat_ur[!(dat_ur$SNP %in% adj_ur$adj$SNP),]
res_ca_out <- mr(dat_ca_out, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger_ca_out <- mr_pleiotropy_test(dat_ca_out)


#########################################################################################################################
#Example 3: BMI and education
#########################################################################################################################

#Exposure
#bmi <- extract_instruments(outcomes=2)
edu <- extract_instruments(outcomes = 1001)

#Outcome
#edu <-  extract_outcome_data(snps=bmi$SNP, outcomes = 1001)
bmi <- extract_outcome_data(snps=edu$SNP, outcomes = 2)

#dat_bmi <- harmonise_data(bmi, edu)
dat_edu <- harmonise_data(edu, bmi)

############################################Radial MR#######################################################

dat_edu2 <- dat_edu[dat_edu$mr_keep, ]
radial_mr <-format_radial(dat_edu2$beta.exposure, dat_edu2$beta.outcome, dat_edu2$se.exposure, dat_edu2$se.outcome, dat_edu2$SNP)
radial_res_ivw <-ivw_radial(radial_mr, 0.05/nrow(dat_edu2), summary = FALSE)
radial_plot <-plot_radial(radial_res_ivw,TRUE,FALSE,TRUE)
radial_funnel <- funnel_radial(radial_res_ivw,TRUE)


#############################################RESULT##########################################################
#TRYX
tryxscan_edu <- tryx.scan(dat_edu, use_proxies = TRUE)
#tryxscan_bmi <- tryx.scan(dat_bmi, use_proxies = TRUE)
tryxscan_edu <- tryx.sig(tryxscan_edu)
#tryxscan_bmi <- tryx.sig(tryxscan_bmi)
adj_edu <- tryx.analyse(tryxscan_edu)



#Multivariable MR
#id_remove <- c("UKB-a:34", "UKB-a:35", "UKB-a:249", "UKB-a:248", "UKB-a:278", "UKB-a:266", "UKB-a:265", "UKB-a:282", 
#               "UKB-a:267", "UKB-a:275", "UKB-a:264", "UKB-a:270", "UKB-a:274", "UKB-a:273", "UKB-a:276", "UKB-a:277", "UKB-a:280", 
#               "UKB-a:287", "UKB-a:288", "UKB-a:291", "UKB-a:286", "UKB-a:284", "UKB-a:283", "UKB-a:289", "UKB-a:388", "UKB-a:35", 
#               "UKB-a:28", "UKB-a:271", "UKB-a:269", "UKB-a:248", "UKB-a:279", "UKB-a:281", "UKB-a:285", "UKB-a:292", "UKB-a:382", 
#               "UKB-a:293", "UKB-a:290", "UKB-a:34", "UKB-a:272", "UKB-a:268", "UKB-a:249", "UKB-a:299", "UKB-a:397", "UKB-a:398", "UKB-a:196",
#               90, 1096, 85, 91, 999, 92, 93, 72, 2)


id_remove <- c("UKB-a:35", "UKB-a:279", "UKB-a:276", "UKB-a:268", "UKB-a:278", "UKB-a:280", "UKB-a:275",
               "UKB-a:282", "UKB-a:284", "UKB-a:281", "UKB-a:286", "UKB-a:285", "UKB-a:289", "UKB-a:292", 
               "UKB-a:288", "UKB-a:397", "UKB-a:248", "UKB-a:277", "UKB-a:267", "UKB-a:265", "UKB-a:274", 
               "UKB-a:388", "UKB-a:287", "UKB-a:293", "UKB-a:283", "UKB-a:290", "UKB-a:264", "UKB-a:196", 
               "UKB-a:291", "UKB-a:249", "UKB-a:382", "UKB-a:389", "UKB-a:399", "UKB-a:34", "UKB-a:266", 
               "UKB-a:270", "UKB-a:271", "UKB-a:398")

ta_edu <- tryx.analyse.mv(tryxscan_edu, id_remove=id_remove)   
#ta_bmi <- tryx.analyse.mv(tryxscan_bmi, id_remove=id_remove)  


#MR
result<- mr(dat_edu, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger <- mr_pleiotropy_test(dat_edu)

#MR-remove outliers (all)
outlier <- tryxscan_edu$outliers
dat_outout <- dat_edu[!(dat_edu$SNP %in% outlier),]
res_outout <- mr(dat_outout, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger_ro <- mr_pleiotropy_test(dat_outout)

#MR-remove outlierS (candidate)
outlier <- tryxscan_edu$outliers
dat_ca_out <- dat_edu[!(dat_edu$SNP %in% adj_edu$adj$SNP),]
res_ca_out <- mr(dat_ca_out, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger_ca_out <- mr_pleiotropy_test(dat_ca_out)



#########################################################################################################################
#Example 4: sleep and schizophrenia
#########################################################################################################################

#Exposure
slp <- extract_instruments(outcomes="UKB-a:9")

#Outcome
scz <-  extract_outcome_data(snps=slp$SNP, outcomes = 22)
dat_slp <- harmonise_data(slp, scz)


############################################Radial MR#######################################################

dat_slp2 <- dat_slp[dat_slp$mr_keep, ]
radial_mr <-format_radial(dat_slp2$beta.exposure, dat_slp2$beta.outcome, dat_slp2$se.exposure, dat_slp2$se.outcome, dat_edu2$SNP)
radial_res_ivw <-ivw_radial(radial_mr, 0.05/nrow(dat_edu2), summary = FALSE)
radial_plot <-plot_radial(radial_res_ivw,TRUE,FALSE,TRUE)
radial_funnel <- funnel_radial(radial_res_ivw,TRUE)


#############################################RESULT##########################################################
#TRYX
tryxscan_slp <- tryx.scan(dat_slp, use_proxies = TRUE)
tryxscan_slp <- tryx.sig(tryxscan_slp)
adj_slp <- tryx.analyse(tryxscan_slp)


#Multivariable MR
id_remove <- c("UKB-a:12", "UKB-a:15")
ta_slp <- tryx.analyse.mv(tryxscan_slp, id_remove=id_remove)   


#MR
result<- mr(dat_slp, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger <- mr_pleiotropy_test(dat_slp)

#MR-remove outliers (all)
outlier <- tryxscan_slp$outliers
dat_outout <- dat_slp[!(dat_slp$SNP %in% outlier),]
res_outout <- mr(dat_outout, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger_ro <- mr_pleiotropy_test(dat_outout)

#MR-remove outlierS (candidate)
outlier <- tryxscan_slp$outliers
dat_ca_out <- dat_slp[!(dat_slp$SNP %in% adj_slp$adj$SNP),]
res_ca_out <- mr(dat_ca_out, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
egger_ca_out <- mr_pleiotropy_test(dat_ca_out)
