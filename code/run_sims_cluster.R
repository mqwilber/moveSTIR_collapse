
# Register a parallel backend
library(doParallel)
# Detect the number of available cores and create cluster
cl <- makeCluster(length(future::availableWorkers()))
# Activate cluster for foreach library
registerDoParallel(cl)

# Load functions
source("code/functions.R")
# Load functions in each node
clusterEvalQ(cl, source("code/functions.R"))


social <- c(0,0.5, 0.7, 0.9,0.93, 0.96, 1)
nus <- 1/(24*c(1/12, 1/3, 1, 3, 7))
gridres <- c(5, 10, 20)
iterations <- 20

trajs <- (replicate(iterations, lapply(social, \(x) simulate_tracks(tau = 5, dp = 2, social = x))))
### New version, edited to have only 3 grid sizes for the same sim.
UDout <- foreach(i = seq_along(trajs), .packages = c("move","ctmm")) %:% foreach(gr = gridres) %dopar% {
  list(SIM = i,
       TRAJS = trajs[[i]],
       UDS = getUDs(trajs[[i]],res = gr))
}

UDout<-unlist(UDout,recursive=FALSE)

out <- foreach(i = seq_along(UDout), .packages = c("move", "ctmm")) %:% 
  foreach(nu = nus, .combine = 'rbind', .inorder = FALSE) %dopar% {
    A <- UDout[[i]]$TRAJS
    UDS <- UDout[[i]]$UDS
    SIM <- UDout[[i]]$SIM
    prods <- getUDprod(UDS)
    # cors <- getCorrs(A, uds[[i]][[2]])
    fois <- getFOI(A, UDS, nu = nu)
    cellarea <- prod(res(fois[[1]][[1]]))
    totareai <- summary(UDS[[2]][[1]])$CI[2]
    totareaj <- summary(UDS[[2]][[2]])$CI[2]
    
    c(SIM,
      nu, 
      cellarea,
      totareai,
      totareaj,
      cellStats(fois[[1]][[3]], sum), 
      cellStats(fois[[1]][[1]],sum), 
      cellStats(fois[[1]][[2]],sum), 
      overlap(UDS[[2]])$CI[,,2][1,2])
    
  }
stopCluster(cl)
outdf <- as.data.frame(do.call(rbind,out))
names(outdf) <- c("sim", "nu", "Ax","Atoti","Atotj", "foi_ud", "foi_full1", "foi_full2", "overlap")
outdf$social <- rep(social, each = length(nus)*length(gridres))
outname <- paste0("outputs/sim_res_", format(Sys.time(), "%y%m%d"), ".csv")
write.csv(outdf,outname, quote = F, row.names = F)