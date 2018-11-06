library(data.table)

ao <- available_outcomes()

bw.child <- extract_instruments(outcomes="UKB-a:318")
chd <-  extract_outcome_data(snps=bw.child$SNP, outcomes = 7)
dat_bw <- harmonise_data(bw.child, chd)
bwch <- mr(dat_bw, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


st.ht <- extract_instruments(outcomes="UKB-a:389")
chd <- extract_outcome_data(snps=st.ht$SNP, outcomes = 7)
dat_st.ht <- harmonise_data(st.ht, chd)
st.htch <-mr(dat_st.ht, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


ldlc <- extract_instruments(outcomes="300")
chd <- extract_outcome_data(snps=ldlc$SNP, outcomes = 7)
dat_ldl <- harmonise_data(ldlc, chd)
ldlch <-mr(dat_ldl, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


hdlc <- extract_instruments(outcomes="299")
chd <- extract_outcome_data(snps=hdlc$SNP, outcomes = 7)
dat_hdl <- harmonise_data(hdlc, chd)
hdlch <-mr(dat_hdl, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


tc <- extract_instruments(outcomes="301")
chd <- extract_outcome_data(snps=tc$SNP, outcomes = 7)
dat_tc <- harmonise_data(tc, chd)
tcch <-mr(dat_tc, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


ibu <- extract_instruments(outcomes="UKB-a:136")
chd <- extract_outcome_data(snps=ibu$SNP, outcomes = 7)
dat_ibu <- harmonise_data(ibu, chd)
ibuch <-mr(dat_ibu, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

#########################################################################################################################

alco <- extract_instruments(outcomes="UKB-a:25")
bmi <- extract_outcome_data(snps=alco$SNP, outcomes = 2)
dat_alco <- harmonise_data(alco, bmi)
alcbmi <-mr(dat_alco, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


walk <- extract_instruments(outcomes="UKB-a:513")
bmi <- extract_outcome_data(snps=walk$SNP, outcomes = 2)
dat_walk <- harmonise_data(walk, bmi)
walkch <-mr(dat_walk, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


#########################################################################################################################

bw <- extract_instruments(outcomes="UKB-a:198")
chd <-  extract_outcome_data(snps=bw$SNP, outcomes = 7)
dat_bww <- harmonise_data(bw, chd)
bwwch <- mr(dat_bww, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


htat10 <- extract_instruments(outcomes="UKB-a:35")
chd <-  extract_outcome_data(snps=htat10$SNP, outcomes = 7)
dat_htat10 <- harmonise_data(htat10, chd)
htat10ch <- mr(dat_htat10, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


hc <- extract_instruments(outcomes="UKB-a:388")
chd <-  extract_outcome_data(snps=hc$SNP, outcomes = 7)
dat_hc <- harmonise_data(hc, chd)
hcch <- mr(dat_hc, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


imparm <- extract_instruments(outcomes="UKB-a:273")
chd <-  extract_outcome_data(snps=imparm$SNP, outcomes = 7)
dat_imparm <- harmonise_data(imparm, chd)
imparmch <- mr(dat_imparm, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


hypo <- extract_instruments(outcomes="UKB-a:77")
chd <-  extract_outcome_data(snps=hypo$SNP, outcomes = 7)
dat_hypo <- harmonise_data(hypo, chd)
hypoch <- mr(dat_hypo, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


smk <- extract_instruments(outcomes="UKB-a:17")
chd <-  extract_outcome_data(snps=smk$SNP, outcomes = 7)
dat_smk <- harmonise_data(smk, chd)
smkch <- mr(dat_smk, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


levo <- extract_instruments(outcomes="UKB-a:190")
chd <-  extract_outcome_data(snps=levo$SNP, outcomes = 7)
dat_levo <- harmonise_data(levo, chd)
levoch <- mr(dat_levo)


wc <- extract_instruments(outcomes="UKB-a:382")
chd <-  extract_outcome_data(snps=wc$SNP, outcomes = 7)
dat_wc <- harmonise_data(wc, chd)
wcch <- mr(dat_wc, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


######################################################################################################

malab <- extract_instruments(outcomes="UKB-a:101")
sch <-  extract_outcome_data(snps=malab$SNP, outcomes = 22)
dat_malab <- harmonise_data(malab, sch)
malabsch <- mr(dat_malab, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


impleg <- extract_instruments(outcomes="UKB-a:271")
sch <-  extract_outcome_data(snps=impleg$SNP, outcomes = 22)
dat_impleg <- harmonise_data(impleg, sch)
implegsch <- mr(dat_impleg, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


memo <- extract_instruments(outcomes="UKB-a:197")
sch <-  extract_outcome_data(snps=memo$SNP, outcomes = 22)
dat_memo <- harmonise_data(memo, sch)
memosch <- mr(dat_memo, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


sink(file = "tryx_example_re_310718.csv")
cat("Empirical analysis 1.  Systolic blood pressure (mmHg) and coronary heart disease (log odds)\n")
write.csv(bwch, row.names = F, quote= F)
cat("--------\n")
write.csv(st.htch, row.names = F, quote= F)        
cat("--------\n")
write.csv(ldlch, row.names = F, quote= F)
cat("--------\n")
write.csv(hdlch, row.names = F, quote= F)      
cat("--------\n")
write.csv(tcch, row.names = F, quote= F)      
cat("--------\n")
write.csv(ibuch, row.names = F, quote= F) 
cat("--------\n")
cat("Empirical analysis 2.  Years of schooling (years) and body mass index (kg/m2)\n")
write.csv(alcbmi, row.names = F, quote= F)     
cat("--------\n")
write.csv(walkch, row.names = F, quote= F)            
cat("--------\n")
cat("Empirical analysis 3.  Urate (mg/dl) and coronary heart disease (log odds)\n")
write.csv(bwwch, row.names = F, quote= F)           
cat("--------\n")
write.csv(htat10ch, row.names = F, quote= F)           
cat("--------\n")       
write.csv(hcch, row.names = F, quote= F)           
cat("--------\n")               
write.csv(imparmch, row.names = F, quote= F)           
cat("--------\n")          
write.csv(hypoch, row.names = F, quote= F)           
cat("--------\n")          
write.csv(smkch, row.names = F, quote= F)           
cat("--------\n")           
write.csv(levoch, row.names = F, quote= F)           
cat("--------\n")           
write.csv(wcch, row.names = F, quote= F)           
cat("--------\n")
cat("Empirical analysis 4.  Sleep duration (hour/night) and schizophrenia (log odds)\n")          
write.csv(malabsch, row.names = F, quote= F)            
cat("--------\n")           
write.csv(implegsch, row.names = F, quote= F)           
cat("--------\n")           
write.csv(memosch, row.names = F, quote= F) 
sink()





library(data.table)

ao <- available_outcomes()

egger.bwch <- mr_pleiotropy_test(dat_bw)
egger.st.htch <- mr_pleiotropy_test(dat_st.ht)
egger.ldlch <- mr_pleiotropy_test(dat_ldl)
egger.hdlch <- mr_pleiotropy_test(dat_hdl)
egger.tcch <- mr_pleiotropy_test(dat_tc)
egger.ibuch <- mr_pleiotropy_test(dat_ibu)
#########################################################################################################################
egger.alcbmi <- mr_pleiotropy_test(dat_alco)
egger.walkch <- mr_pleiotropy_test(dat_walk)
#########################################################################################################################
egger.bwwch <- mr_pleiotropy_test(dat_bww)
egger.htat10ch <- mr_pleiotropy_test(dat_htat10)
egger.hcch <- mr_pleiotropy_test(dat_hc)
egger.imparmch <- mr_pleiotropy_test(dat_imparm)
egger.hypoch  <- mr_pleiotropy_test(dat_hypo)
egger.smkch <- mr_pleiotropy_test(dat_smk)
egger.levoch <- mr_pleiotropy_test(dat_levo)
egger.wcch <- mr_pleiotropy_test(dat_wc)
######################################################################################################
egger.malabsch <- mr_pleiotropy_test(dat_malab)
egger.implegsch <- mr_pleiotropy_test(dat_impleg)
egger.memosch <- mr_pleiotropy_test(dat_memo)


######################################################################################################


sink(file = "tryx_example_egger_310718.csv")
cat("Empirical analysis 1.  Systolic blood pressure (mmHg) and coronary heart disease (log odds)\n")
write.csv(egger.bwch, row.names = F, quote= F)
cat("--------\n")
write.csv(egger.st.htch, row.names = F, quote= F)        
cat("--------\n")
write.csv(egger.ldlch, row.names = F, quote= F)
cat("--------\n")
write.csv(egger.hdlch, row.names = F, quote= F)      
cat("--------\n")
write.csv(egger.tcch, row.names = F, quote= F)      
cat("--------\n")
write.csv(egger.ibuch, row.names = F, quote= F) 
cat("--------\n")
cat("Empirical analysis 2.  Years of schooling (years) and body mass index (kg/m2)\n")
write.csv(egger.alcbmi, row.names = F, quote= F)     
cat("--------\n")
write.csv(egger.walkch, row.names = F, quote= F)            
cat("--------\n")
cat("Empirical analysis 3.  Urate (mg/dl) and coronary heart disease (log odds)\n")
write.csv(egger.bwwch, row.names = F, quote= F)           
cat("--------\n")
write.csv(egger.htat10ch, row.names = F, quote= F)           
cat("--------\n")       
write.csv(egger.hcch, row.names = F, quote= F)           
cat("--------\n")               
write.csv(egger.imparmch, row.names = F, quote= F)           
cat("--------\n")          
write.csv(egger.hypoch, row.names = F, quote= F)           
cat("--------\n")          
write.csv(egger.smkch, row.names = F, quote= F)           
cat("--------\n")           
write.csv(egger.levoch, row.names = F, quote= F)           
cat("--------\n")           
write.csv(egger.wcch, row.names = F, quote= F)           
cat("--------\n")
cat("Empirical analysis 4.  Sleep duration (hour/night) and schizophrenia (log odds)\n")          
write.csv(egger.malabsch, row.names = F, quote= F)            
cat("--------\n")           
write.csv(egger.implegsch, row.names = F, quote= F)           
cat("--------\n")           
write.csv(egger.memosch, row.names = F, quote= F) 
sink()

