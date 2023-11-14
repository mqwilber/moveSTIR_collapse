
# Register a parallel backend
library(parallel)
library(doParallel)
# Detect the number of available cores and create cluster
cl <- makeCluster(4)
# Activate cluster for foreach library
registerDoParallel(cl)

# Load functions
source("functions.R")
# Load functions in each node
clusterCall(cl, source,"./functions.R")

# pars <- expand.grid(#nu = 1/(24*c(1/12, 1/6, 1/3, 1,3,7)), 
#                     social = c(0.7,0.9, 0.925,0.95, 0.975,1), 
#                     replicate = 1:10)
social <- c(0,0.5, 0.7, 0.9,0.93, 0.96, 1)
nus <- 1/(24*3)#1/(24*c(1/12, 1/6, 1/3, 1, 3, 7))
gridres <- 10# c(5, 10, 20)
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
      test <- getFOI(A,prods, nu = nu)
      cellarea <- prod(res(test[[1]]))
      
      c(nu, 
        gr,
        cellStats(1/nu*prods[[1]][[1]], sum), 
        sapply(test, cellStats,sum)[1], 
        sapply(test, cellStats,sum)[2], 
        overlap(uds[[2]])$CI[,,2][1,2] )
    }
} 

stopCluster(cl)
outdf <- as.data.frame(do.call(rbind, out))
names(outdf) <- c("nu", "res", "foi_ud", "foi_full1", "foi_full2", "overlap")
outdf$social <- rep(social, each = length(nus)*length(gridres))
# outname <- paste0("sim_res_nu", Sys.Date(), ".csv")
write.csv(outdf,"sim_res_soc_3d_10m.csv", quote = F, row.names = F)