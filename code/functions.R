##### --------------- FUNCTION DEFINITIONS --------------#####
simulate_tracks <- function(tau = 10, g = 1e3, steps = 500, dp = 0.2, social = NULL, smooth = TRUE) {
  nind = 2
  p.lambdas = matrix(rpois(nind*2, dp*g)-dp*g, ncol = 2)# home range centers
  p.x0 = p.lambdas # initial positions
  p.phi = 1 # trajectory smoothing parameter
  nsteps = steps
  
  xs = replicate(nind, numeric(nsteps))
  ys = replicate(nind, numeric(nsteps))
  
  for (s in seq_len(nind)) {
    x <- y <- numeric(nsteps)
    x[1] <- mux <- p.x0[s,1]
    y[1] <- muy <- p.x0[s,2]
    xnoise = rnorm(nsteps)
    ynoise = rnorm(nsteps)
    # take successive steps
    for(i in seq_len(nsteps-1)) {
      x[i+1] = x[i]-1/tau*(x[i]-p.lambdas[s,1])+sqrt(g)*xnoise[i]
      y[i+1] = y[i]-1/tau*(y[i]-p.lambdas[s,2])+sqrt(g)*ynoise[i]
    }
    xs[,s] = x
    ys[,s] = y
  }
  
  if(is.numeric(social)) {
    W <- replicate(nsteps, diag(nrow = nind, ncol = nind))
    W[1,2,] <- social
    # make symmetric
    W <- apply(W, 3, \(x) x*upper.tri(x, diag = T)+t(x*upper.tri(x)), simplify = FALSE)
    # convolution kernel hij(tau) is the weights matrix, divided by the row sums
    h_soc <- lapply(W, \(x) x/rowSums(x))
    # create empty matrices with the same dimensions as the position matrices
    xsoc <- matrix(nrow = nrow(xs), ncol = ncol(xs))
    ysoc <- matrix(nrow = nrow(ys), ncol = ncol(ys))
    # fill the array using dot product of the step*ind x and y matrices, multiplied by the rows of h_soc
    for (t in 1:nsteps) {
      xs[t,] <- h_soc[[t]]%*%xs[t,]
      ys[t,] <- h_soc[[t]]%*%ys[t,]
    }
  }
  
  
  if(smooth) {
    # create smoother matrix for convolution
    H_inl <- outer(seq_len(nsteps),seq_len(nsteps), \(x,y) abs(x-y)/p.phi*besselK(abs(x-y)/p.phi, nu = 1))
    diag(H_inl) <- 1
    H_inl <- H_inl/rowSums(H_inl)
    xs <- apply(xs, 2, "%*%", H_inl)
    ys <- apply(ys, 2, "%*%", H_inl)
  }
  list(x=xs,y=ys)
}

# Function to get the utilization distributions. Output is a list, first element
# are the ctmm fits, second are the individual UDs
getUDs <- function(X, res = 10, ctmm = TRUE) {
  xs <- X[[1]]
  ys <- X[[2]]
  nind <- ncol(xs)
  nsteps <- nrow(xs)
  if (ctmm) {
    # create move objects
    moveobjs <- lapply(seq_len(nind), \(c) move(x = xs[,c],y = ys[,c],time = as.POSIXct(3600*seq_len(nsteps),origin = Sys.time()), animal = paste0("ind",c), proj = "+proj=tmerc"))
    # convert to telemetry
    telemetries <- lapply(moveobjs, as.telemetry)
    # ctmm fit
    GUESS <- lapply(telemetries, \(x) ctmm.guess(x, interactive = F))
    FITS <- lapply(seq_along(telemetries), \(x) ctmm.select(telemetries[[x]], GUESS[[x]]))
    # AKDE
    UDS <- akde(telemetries, FITS, res = res)
    return(list(FITS, UDS, method = "AKDE"))
  } else {
    # Create SpatialPoints object
    dat_sp <- sp::SpatialPointsDataFrame(coords = data.frame(x=xs,y=ys), 
                                                     data = data.frame(id = rep(seq_len(nind),each = nsteps)))
    # Create grid on which UD will be estimated
    kde_grid <- sp::SpatialGrid(grid = GridTopology(cellcentre.offset =bbox(dat_sp)[,1]-bbox(dat_sp)[,1]%%10, cellsize = c(contact_dist,contact_dist), cells.dim = ceiling(apply(bbox(dat_sp), 1, diff)/contact_dist)))
    kde_grid_px <- sp::SpatialPixels(SpatialPoints(coordinates(kde_grid)))
    # Estimate UD
    UDS <- adehabitatHR::kernelUD(dat_sp, grid = kde_grid_px)
    return(list(NA,UDS,method = "KDE"))
  }
}

getUDprod <- function(X) {
  uds <- X[[2]]
  nind <- length(uds)
  ### UD and SD products ###
  rs <- lapply(uds, raster, DF = "PMF")
  cellarea <- prod(res(rs[[1]]))
  
  # extend all rasters to have the same extent. 
  rs <- lapply(rs, extend, y = extent(do.call(merge,rs)), value = 0)
  # Create lists to store objects
  # UDprods <- UDsds <- list()
  OUT <- list()
  # Possible combinations
  combs <- combn(length(uds), 2)
  # Calculate UD and SD products from UD pair values
  for (i in seq_len(ncol(combs))) {
    ind1 <- combs[1,i]
    ind2 <- combs[2,i]
    r1 <- rs[[ind1]]
    r2 <- rs[[ind2]]
    # product of probabilities, divided by area
    udprob <- r1*r2
    # UDprods[[paste(ind1,ind2,sep = "-")]] <- UDprods[[paste(ind2,ind1,sep = "-")]] <- udprob
    sdprob <- sqrt(r1*(1-r1))*sqrt(r2*(1-r2))
    # UDsds[[paste(ind1,ind2,sep = "-")]] <- UDsds[[paste(ind2,ind1,sep = "-")]] <- sdprob
    OUT[[paste(ind1,ind2,sep = "-")]] <- list(UD = udprob, SD = sdprob)
  }
  # return(list(UD = UDprods, SD = UDsds))
  OUT
}
# function to calculate the correlations. Output is a list, where every element
# is a lags by cells matrix of correlation between two individuals
getCorrs <- function(xy, prods, prewt = TRUE, fltr = "none") {
  gridcors2 <- list()
  xs <- xy[[1]]
  ys <- xy[[2]]
  nind <- ncol(xs)
  nsteps <- nrow(xs)
  combs <- combn(nind,2)
  for (i in 1:ncol(combs)) {
    r1 <- raster(prods[[i]]$UD)
    ind1 <- combs[1,i]
    ind2 <- combs[2,i]
    # get position histories (i.e. which cell was each individual in at every time
    # point). 
    pos1 <- cellFromXY(r1, xy = cbind(xs[,ind1],ys[,ind1]))
    pos2 <- cellFromXY(r1, xy = cbind(xs[,ind2],ys[,ind2]))
    
    # keep only cells that both visited at some point
    ovlpcells <- unique(pos1)[unique(pos1) %in% unique(pos2)]
    sigcells <- numeric(length = length(ovlpcells))
    if(length(ovlpcells)==0) {
      # Still fill in an item in the list, but write just NA. 
      gridcors2[[paste(ind1,ind2,sep = "-")]] <- list(CAB = NA, CBA = NA, psig = NA)
    } else {
      maxlag <- ceiling(nsteps/2)
      cormat_ab <- cormat_ba <- matrix(0,nrow = maxlag+1, ncol = length(ovlpcells))
      for (j in seq_along(ovlpcells)) {
        cell <- ovlpcells[j]
        a <- b <- numeric(nsteps)
        a[cell==pos1] <- b[cell==pos2]<- 1
        xcorr <- if(prewt) TSA::prewhiten(a,b,lag.max = maxlag, plot = F)$ccf else ccf(a, b, lag.max = maxlag, plot = F)
        xcorr_vals <- as.numeric(xcorr$acf)
        FILT <- fltr
        if(FILT == "bs") {
          BSC <- funtimes::ccf_boot(a,b,lag.max = maxlag, plot = "none")
          xcorr_vals <- BSC$r_P
          xcorr_vals <- xcorr_vals*(xcorr_vals<BSC$lower_P | xcorr_vals>BSC$upper_P)
        } else if(FILT == "reg") {
          xcorr_vals <- xcorr_vals*(abs(xcorr_vals)>1.96/sqrt(xcorr$n.used))
        }

        # sigcells[j] <- switch(ci.method[1], 
        #                       reg = mean(abs(xcorr_vals)>(1.96/sqrt(xcorr$n.used))),
        #                       bs = do.call(function(a,b, cors, n = 1000) {
        #                         maxlag = ceiling(nsteps/2)
        #                         M <- replicate(n, {
        #                           a2 <- sample(a, length(a), T)
        #                           b2 <- sample(b, length(b), T)
        #                           ccf(a2,b2,lag.max = maxlag, plot = F)$acf
        #                         })
        #                         Q <- apply(M, 1, quantile,probs = c(0.025, 0.975))
        #                         mean(cors$acf<Q[1,] | cors$acf>Q[2,])
        #                       }, list(a = a, b = b, cors = xcorr)))

        cormat_ab[,j] <- xcorr_vals[(maxlag+1):1]
        cormat_ba[,j] <- xcorr_vals[(maxlag+1):length(xcorr_vals)]
      }
      dimnames(cormat_ab) <- dimnames(cormat_ba) <- list(NULL, cell = ovlpcells)
      gridcors2[[paste(ind1,ind2,sep = "-")]] <- list(CAB = cormat_ab, CBA = cormat_ba, N = nsteps)
    }
  }
  return(gridcors2)
}
# function to calculate the per-cell FOI. Calculates correlation within. Output
# is list with each element 2 FOI rasters (one for each direction)
getFOI <- function(xy, uds = NULL, spr.rm = TRUE, beta = 1, lambda = 1, nu = 1/(24*7)) {
  foirasts2 <- list()
  UDS <- if(is.null(uds)) getUDs(xy) else uds
  PRODS <- getUDprod(UDS)
  gridcors2 <- getCorrs(xy,PRODS)
  for (gi in seq_along(gridcors2)) {
    udp <- PRODS[[gi]]$UD
    sdp <- PRODS[[gi]]$SD
    CORRS <- gridcors2[[gi]]
    cellarea <- prod(res(udp))
    
    if (all(sapply(CORRS[c(1,2)],is.array))) {
      corcells <- as.numeric(colnames(CORRS[[1]]))
      corrastAB <- corrastBA <- udp
      values(corrastAB) <- values(corrastBA) <- 0
      # get lags
      lags <- 0:(nrow(CORRS[[1]])-1)
      dtau <- unique(diff(lags))
      # scale and integrate correlation at every cell. In simulations we will
      # use a step function of 10 time steps
      if(spr.rm) {
        CORRS$CAB <- CORRS$CAB*(abs(CORRS$CAB)>(1.96/sqrt(CORRS$N)))
        CORRS$CBA <- CORRS$CBA*(abs(CORRS$CBA)>(1.96/sqrt(CORRS$N)))
        }
      corrastAB[corcells] <- if(dim(CORRS$CAB)[2]>1) colSums(CORRS$CAB[1:10,]*dtau) else sum(CORRS$CAB[1:10,]*dtau)
      corrastBA[corcells] <- if(dim(CORRS$CBA)[2]>1) colSums(CORRS$CBA[1:10,]*dtau) else sum(CORRS$CBA[1:10,]*dtau)
      FOIAB <- beta/cellarea*lambda*(1/nu*udp+sdp*corrastAB)
      FOIAB <- FOIAB*(FOIAB>=0)
      FOIBA <- beta/cellarea*lambda*(1/nu*udp+sdp*corrastBA)
      FOIBA <- FOIBA*(FOIBA>=0)
      foirasts2[[gi]] <- list(FAB = FOIAB, FBA = FOIBA, FUD = beta/cellarea*lambda*(1/nu*udp))
      
      # for (j in 1:2) {
      #   values(corrast) <- 0
      #   corrast[corcells] <- if (sum(sig.indx)>1) colSums(gridcors2[[gi]][[j]][,sig.indx]*exp(-nu*lags)*dtau) else sum(gridcors2[[gi]][[j]][,sig.indx]*exp(-nu*lags)*dtau)
      #   FOI <- beta/cellarea*lambda*(1/nu*udp+sdp*corrast)
      #   foirasts2[[gi]][[j]] <- FOI*(FOI>=0)
      # }
      # foi <- beta/cellarea*lambda*(1/nu*udp+sdp*corrast)
    } else {
      foirasts2[[gi]] <- replicate(3, beta/cellarea*lambda*(1/nu*udp), simplify = F)
      names(foirasts2[[gi]]) <- c("FAB", "FBA", "FUD")
    }
    # add UD-only output to all FOIs
    foirasts2[[gi]][[3]] <- beta/cellarea*lambda*(1/nu*udp)
    # to delete
    # foi <- foi*(foi>=0)
  }
  return(foirasts2)
}
