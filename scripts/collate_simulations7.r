library(dplyr)

a <- list()
b <- list()

for(i in 1:640)
{
	message(i)
	load(paste0("../results/scratch/sim7_", i, ".rdata"))
	a[[i]] <- lapply(l, function(x) { x$plot <- NULL; return(x)})
	b[[i]] <- param
}

l <- unlist(a, recursive=FALSE)
param <- bind_rows(b)
save(l, param, file="../results/sim7.rdata")
