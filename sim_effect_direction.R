library(tryx)
library(simulateGP)
library(parallel)
library(tidyverse)


mv_lasso <- function(mvdat)
{
  b <- glmnet::cv.glmnet(x=mvdat$exposure_beta, y=mvdat$outcome_beta, weight=1/mvdat$outcome_se^2, intercept=0)
  c <- coef(b, s = "lambda.min")
  x <- unique(c(rownames(c)[!c[,1] == 0]))
  return(x)
}


l <- mclapply(1:1000, function(i)
  {
  set.seed(i)
  # Genotypes for 80 SNPs and 100000 individuals, where the SNPs have MAF = 0.5:
  nid <- 100000
  g <- make_geno(nid, 80, 0.5)
  
  # Choose effect sizes for instruments for each trait
  # These SNPs instrument some exposure, and together explain 5% of the variance
  effs1 <- choose_effects(50, 0.05)
  effs2 <- choose_effects(50, 0.05)
  
  # Create X1 and X2, where they overlap some variants
  x1 <- make_phen(effs1, g[,1:50])
  x2 <- make_phen(effs2, g[,31:80])

  # Create Y - x1 has a -0.3 influence on it and x2 has a +0.3 influence on it
  y <- make_phen(c(-0.3, 0.3), cbind(x1, x2))

  # Create proxy phenotypes for x1 and x2
  nprox <- 5
  x1p <- matrix(0, nid, nprox)
  x2p <- matrix(0, nid, nprox)
  for(j in 1:nprox)
  {
    x1p[,j] <- make_phen(runif(1), x1)
    x2p[,j] <- make_phen(runif(1), x2)
  }
    
  # Perform separate MR on each
  dat1 <- get_effs(x1, y, g)
  dat2 <- get_effs(x2, y, g)
  
  
  library(TwoSampleMR)
  TwoSampleMR::mr(subset(dat1, pval.exposure < 5e-8))
  TwoSampleMR::mr(subset(dat2, pval.exposure < 5e-8))
  
  # Do multivariable MR
  # First get the effects for x1, x2 and y, and put them in mv format
  mvdat <- make_mvdat(list(x1, x2), y, g)
  

  # Perform LASSO MV MR

  exposure_names <- c("x1", "x2")
  x <- c(list(x1, x2))
  names(x) <- exposure_names
  colnames(mvdat$exposure_beta) <- colnames(mvdat$exposure_se) <- colnames(mvdat$exposure_se) <- exposure_names

  out <- list()
  out$mvmr_full <- mv_multiple(mvdat)$result
  out$lasso <- mv_lasso(mvdat)
    
  mvdat <- make_mvdat(x[out$lasso], y, g)
  colnames(mvdat$exposure_beta) <- colnames(mvdat$exposure_se) <- colnames(mvdat$exposure_se) <- out$lasso
  out$mvmr_selected <- mv_multiple(mvdat)$result
  return(out)
  
})


save(l, file="/newhome/yc16575/simulation/effect_direction_130319.rdata")


library(dplyr)
out <- lapply(1:length(l), function(i) {
  x <- l[[i]]
  x[[1]]$what <- "mvmr_full"
  x[[3]]$what <- "mvmr_selected"
  x <- rbind(x[[1]],x[[3]])
  x$sim <- i
  return(x)
  }) %>% bind_rows

out %>%
  group_by(what, exposure) %>% #I found that this command doesn’t work for list format.
  summarise(
    n=n(),
    nsnp=mean(nsnp),
    b=mean(b),
    se=mean(se),
    p=mean(pval)
  )

