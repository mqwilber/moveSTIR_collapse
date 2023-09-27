
# Register a parallel backend
library(parallel)
library(doParallel)
# Detect the number of available cores and create cluster
cl <- makeCluster(future::availableWorkers())
# Activate cluster for foreach library
registerDoParallel(cl)

# Load functions
source("functions.R")
# Load functions in each node
clusterCall(cl, source,"./functions.R")


social <- 0.5
nus <- 1/2
gridres <- 10
iterations <- 5

trajs <- (replicate(iterations, lapply(social, \(x) simulate_tracks(tau = 5, dp = 2, social = x))))
### New version, edited to have only 3 grid sizes for the same sim.
UDout <- foreach(i = seq_along(trajs), .packages = c("move","ctmm")) %:% foreach(gr = gridres) %dopar% {
  list(SIM = i,
       TRAJS = trajs[[i]],
       UDS = getUDs(trajs[[i]],res = gr))
}

out <- foreach(i = seq_along(uds), .packages = c("move", "ctmm"), .combine = 'rbind', .inorder = FALSE) %:% 
  foreach(nu = nus, .combine = 'rbind', .inorder = FALSE) %dopar% {
    A <- UDout[[i]]$TRAJS
    UDS <- UDout[[i]]$UDS
    SIM <- UDout[[i]]$SIM
    prods <- getUDprod(UDS)
    # cors <- getCorrs(A, uds[[i]][[2]])
    fois <- getFOI(A, prods, nu = nu)
    cellarea <- prod(res(fois[[1]]))
    totareai <- summary(UDS[[2]][[1]])$CI[2]
    totareaj <- summary(UDS[[2]][[2]])$CI[2]
    
    c(SIM,
      nu, 
      cellarea,
      totareai,
      totareaj,
      cellStats(1*1/24/cellarea*1/nu*prods[[1]][[1]], sum), 
      cellStats(fois[[1]],sum), 
      cellStats(fois[[2]],sum), 
      overlap(uds[[2]])$CI[,,2][1,2])
  }
### Delete if the above works as expected
# out <- list()
# for (i in seq_along(trajs)) {
#   A <- trajs[[i]]
#   out[[i]] <- foreach(gr = gridres, 
#                       .combine = 'rbind', .packages = c("move", "ctmm"), .inorder = FALSE) %:%
#     foreach(nu = nus, .combine = 'rbind', .inorder = FALSE) %dopar% {
#       uds <- getUDs(A, gr)
#       prods <- getUDprod(uds)
#       cors <- getCorrs(A, uds[[2]])
#       fois <- getFOI(A,prods, nu = nu)
#       cellarea <- prod(res(fois[[1]]))
#       totareai <- summary(uds[[2]][[1]])$CI[2]
#       totareaj <- summary(uds[[2]][[2]])$CI[2]
#       
#       c(i,
#         nu, 
#         cellarea,
#         totareai,
#         totareaj,
#         cellStats(1*1/24/cellarea*1/nu*prods[[1]][[1]], sum), 
#         cellStats(fois[[1]],sum), 
#         cellStats(fois[[2]],sum), 
#         overlap(uds[[2]])$CI[,,2][1,2])
#     }
# } 

stopCluster(cl)
outdf <- as.data.frame(do.call(rbind, out))
names(outdf) <- c("sim", "nu", "Ax","Atoti","Atotj", "foi_ud", "foi_full1", "foi_full2", "overlap")
outdf$social <- rep(social, each = length(nus)*length(gridres))
outname <- paste0("sim_res_", format(Sys.time(), "%y%m%d-%H%M"), "_test.csv")
write.csv(outdf,outname, quote = F, row.names = F)