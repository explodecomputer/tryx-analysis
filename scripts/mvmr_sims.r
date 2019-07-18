library(tryx)
library(simulateGP)
library(parallel)
library(tidyverse)

mv_lasso <- function(mvdat)
{
  b <- glmnet::cv.glmnet(x=mvdat$exposure_beta, y=mvdat$outcome_beta, weight=1/mvdat$outcome_se^2, intercept=0)
  c <- coef(b, s = "lambda.min")
  i <- !c[,1] == 0
  d <- data_frame(exposure=rownames(c)[i], b=c[i,])
  return(d)
}

runsim <- function(snps1 = 1:50, snps2 = 31:80, nprox=5, effx1=0.2, effx2=0.2, nid=50000)
{
  g <- make_geno(nid, length(unique(c(snps1, snps2))), 0.5)    
  effs1 <- choose_effects(length(snps1), 0.05)
  effs2 <- choose_effects(length(snps2), 0.05)
  x1 <- make_phen(effs1, g[,snps1])
  x2 <- make_phen(effs2, g[,snps2])
  y <- make_phen(c(effx1, effx2), cbind(x1, x2))
  x1p <- list(x1)
  x2p <- list(x2)
  for(j in 1:nprox)
  {
    x1p[[j+1]] <- make_phen(runif(1), x1)
    x2p[[j+1]] <- make_phen(runif(1), x2)
  }
  names(x1p) <- paste0("a", 1:(nprox+1))
  names(x2p) <- paste0("b", 1:(nprox+1))
  mvdat <- make_mvdat(c(x1p, x2p), y, g)
  out <- list()
  out$mvmr_full <- mv_multiple(mvdat)$result
  out$lasso <- mv_lasso(mvdat)

  mvdat <- make_mvdat(c(x1p[names(x1p) %in% out$lasso$exposure], x2p[names(x2p) %in% out$lasso$exposure]), y, g)
  colnames(mvdat$exposure_beta) <- colnames(mvdat$exposure_se) <- colnames(mvdat$exposure_se) <- out$lasso$exposure
  out$mvmr_selected <- mv_multiple(mvdat)$result
  return(out)
}

l1 <- mclapply(1:500, function(S)
  {
    message(S)
    set.seed(S)
    runsim()
  }, mc.cores=16)

l2 <- mclapply(1:500, function(S)
  {
    message(S)
    set.seed(S)
    runsim(snps1=1:40, snps2=41:80)
  }, mc.cores=16)

m1 <- lapply(1:length(l1), function(x)
{
  bind_rows(
    mutate(l1[[x]][[1]], method="full"),
    mutate(l1[[x]][[1]] %>% filter(pval < 0.05), method="pval"),
    mutate(l1[[x]][[2]], method="lasso1"),
    mutate(l1[[x]][[3]], method="lasso2"),
  ) %>% mutate(sim=x, plei=TRUE)
}) %>% bind_rows %>% as_tibble

m2 <- lapply(1:length(l2), function(x)
{
  bind_rows(
    mutate(l2[[x]][[1]], method="full"),
    mutate(l2[[x]][[1]] %>% filter(pval < 0.05), method="pval"),
    mutate(l2[[x]][[2]], method="lasso1"),
    mutate(l2[[x]][[3]], method="lasso2")
  ) %>% mutate(sim=x, plei=FALSE)
}) %>% bind_rows %>% as_tibble

m <- bind_rows(m1, m2) %>% 
group_by(sim, plei) %>%
do({
  x <- .
  for(exp in c("a", "b"))
  {
    for(meth in c("full", "pval", "lasso1", "lasso2"))
    {
      if(sum(grepl(exp, x$exposure) & x$method == meth) == 0)
      {
        x <- bind_rows(x, tibble(exposure=exp, b=0, method=meth, sim=x$sim[1], plei=x$plei[1]))
      }
    }
  }
  x
})

save(m, file="../results/mvmr_simulations.rdata")

s <- m %>% group_by(expa=grepl("a",exposure), method, sim, plei) %>%
summarise(n=sum(b!=0), b=sum(b))

group_by(s, expa, method, plei) %>% summarise(bias=mean(b)-0.2, sd=sd(b), n=mean(n), count=n())
