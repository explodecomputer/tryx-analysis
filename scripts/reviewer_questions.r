# 5) It would be interesting to see how the multivariable MR with the hypothesized exposure + the candidate exposure(s) (meaning without adjusting for outliers) compares to the tested methods. Is it perfectly equivalent to adjusting for outliers in MR-TRYX?

library(TwoSampleMR)
library(simulateGP)
g <- make_geno(10000, 19, 0.5)
e1 <- choose_effects(10, 0.05)
e2 <- choose_effects(10, 0.05)
x1 <- make_phen(e1, g[,1:10])
x2 <- make_phen(e2, g[,10:19])

y <- make_phen(c(0.2, 0.2), cbind(x1, x2))

mvdat <- make_mvdat(list(x1,x2), y, g)
mv_multiple(mvdat)

est1 <- get_effs(x1,y,g)
est2 <- get_effs(x2,y,g)

# path from g10 -> x2 -> y

x2inst <- subset(est2, SNP !=10 & pval.exposure < 5e-8)
x2_y <- mr(x2inst, )$b

path <- est2$beta.exposure[est2$SNP == 10] * x2_y

est1adj <- est1
est1adj$beta.outcome[est1adj$SNP == 10] <-  est1adj$beta.outcome[est1adj$SNP == 10] - path


mr(subset(est1adj, pval.exposure < 5e-8))

