#### FUNCTIONS ####
getProds <- function(x,y) {
  # Check whether the grids overlap/extend to a common grid for both, or all
  r1 <- extend(x, extent(merge(x,y)), value=0)
  r2 <- extend(y, extent(merge(x,y)), value=0)
  # UD and SD products
  udprod <- r1*r2
  sdprod <- sqrt(r1*(1-r1))*sqrt(r2*(1-r2))
  list(UD = udprod, SD = sdprod)
}
# function to get vector of positions (cell index) for two individuals at a
# time. The grid is common to both, obtained from the getProds function
getPositions <- function(X,R) {
  grd <- R$UD
  pos1 <- cellFromXY(grd, xy = as.matrix(X[[1]][,c("x","y")], ncol = 2))
  pos2 <- cellFromXY(grd, xy = as.matrix(X[[2]][,c("x","y")], ncol = 2))
  list(pos1,pos2)
}
# function to get the correlation matrices for a single pair. Outputs are the
# correlation from a->b, from b->a, and the proportion of "significant"
# correlations within the cell, assessed by confidence intervals. There is an
# option to prewhiten, or to estimate CI using bootstrapping
pairCorrs <- function(X, prewt = TRUE, ci = c("reg", "bs"), export = F, niter = 20) {
  pos1 <- X[[1]]
  pos2 <- X[[2]]
  ovlpcells <- unique(pos1)[unique(pos1) %in% unique(pos2)]
  sigcells <- numeric(length = length(ovlpcells))
  if(length(ovlpcells)==0) {
    cat("\nThere are no overlap cells")
    NA
  } else {
    nsteps <- length(pos1)
    maxlag <- ceiling(nsteps/2)
    cormat_ab <- cormat_ba <- matrix(0, nrow = maxlag+1, ncol = length(ovlpcells))
    for (j in seq_along(ovlpcells)) {
      cell <- ovlpcells[j]
      a <- b <- numeric(nsteps)
      a[cell==pos1] <- 1
      b[cell==pos2] <- 1
      xcorr <- if(prewt) TSA::prewhiten(a,b,lag.max = maxlag, plot = F)$ccf else ccf(a, b, lag.max = maxlag, plot = F)
      xcorr_vals <- as.numeric(xcorr$acf)
      sigcells[j] <- switch(ci[1], 
                            reg = mean(abs(xcorr_vals)>(1.96/sqrt(xcorr$n.used))),
                            bs = corboot(AB = cbind(a,b), cors = xcorr,n = niter))
      cormat_ab[,j] <- xcorr_vals[(maxlag+1):1]
      cormat_ba[,j] <- xcorr_vals[(maxlag+1):length(xcorr_vals)]
    }
    dimnames(cormat_ab) <- dimnames(cormat_ba) <- list(lag = 0:maxlag, cell = ovlpcells)
    # Export
    if(export) {
      write.csv(cormat_ab, paste0("outputs/correlations_10min_",ids[ind1],"-",ids[ind2],"_",format(Sys.Date(), "%m%d"), ".csv"))
      write.csv(cormat_ba, paste0("outputs/correlations_10min_",ids[ind2],"-",ids[ind1],"_",format(Sys.Date(), "%m%d"), ".csv"))
    }
    return(list(CAB = cormat_ab, CBA = cormat_ba, psig = sigcells))
  }
}


interpTrajs <- function(x,y, lag = "10 min") {
  tr1 <- range(x$timestamp)
  tr2 <- range(y$timestamp)
  # Check if there is temporal overlap
  if(max(tr1[1],tr2[1])>min(min(tr1[2],tr2[2]))) {
    cat("There is no temporal overlap \n")
    # foi_ab <- foi_ba <- beta/Area*lam*(1/nu*udprod)
  } else {
    tseq <- seq(max(tr1[1],tr2[1]),min(tr1[2],tr2[2]), lag)
    lags <- as.numeric(tseq-min(tseq))
    nsteps <- length(tseq)
    # Interpolate trajectories
    interp_traj_1 <- predict(telemetries[[ind1]], FITS[[ind1]], t = tseq)
    interp_traj_2 <- predict(telemetries[[ind2]], FITS[[ind2]], t = tseq)
    list(interp_traj_1,interp_traj_2)
  }
}

corboot <- function(AB, cors, n = 1000) {
  maxlag <- (length(cors$acf)-1)/2
  M <- replicate(n, {
    a2 <- sample(AB[,1])
    b2 <- sample(AB[,2])
    ccf(a2,b2,lag.max = maxlag, plot = F)$acf
  })
  Q <- apply(M, 1, quantile,probs = c(0.025, 0.975))
  mean(cors$acf<Q[1,] | cors$acf>Q[2,])
}
calcFOI <- function() {
  
}