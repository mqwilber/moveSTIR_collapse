
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


social <- c(0,0.5, 0.7, 0.9,0.93, 0.96, 1)
nus <- 1/(24*c(1/12, 1/3, 1, 3, 7))
gridres <- c(5, 10, 20)
iterations <- 10

trajs <- (replicate(iterations, lapply(social, \(x) simulate_tracks(tau = 5, dp = 2, social = x))))
out <- list()
for (i in seq_along(trajs)) {
  A <- trajs[[i]]
  out[[i]] <- foreach(gr = gridres, 
                      .combine = 'rbind', .packages = c("move", "ctmm"), .inorder = FALSE) %:%
    foreach(nu = nus, .combine = 'rbind', .inorder = FALSE) %dopar% {
      uds <- getUDs(A, gr)
      prods <- getUDprod(uds)
      cors <- getCorrs(A, uds[[2]])
      fois <- getFOI(A,prods, nu = nu)
      cellarea <- prod(res(fois[[1]]))
      totarea <- length(fois[[1]])*cellarea
      
      c(i,
        nu, 
        cellarea,
        totarea,
        cellStats(1*1/24/cellarea*1/nu*prods[[1]][[1]], sum), 
        cellStats(fois[[1]],sum), 
        cellStats(fois[[2]],sum), 
        overlap(uds[[2]])$CI[,,2][1,2])
    }
} 

stopCluster(cl)
outdf <- as.data.frame(do.call(rbind, out))
names(outdf) <- c("sim", "nu", "Ax","Atot" "foi_ud", "foi_full1", "foi_full2", "overlap")
outdf$social <- rep(social, each = length(nus)*length(gridres))
outname <- paste0("sim_res_", format(Sys.time(), "%y%m%d-%H%M"), ".csv")
write.csv(outdf,outname, quote = F, row.names = F)