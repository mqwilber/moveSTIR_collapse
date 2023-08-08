# Analysis of deer data using the pMoveSTIR framework
# Author: Juan S. Vargas Soto
# Date: July 2023

# Intro 
# We have developed a framework to estimate the expected FOI from movement data.
# GPS tracking data is described using ctmm, and from that we build utilization
# distributions. These are the basis for calculating the overlap (as product of
# UDs) and the product of the SDs at every cell. We then calculate the location
# histories and the cross-correlation for every cell. 

# Load libraries
library(move)
library(ctmm)

# Set epidemiological parameters: threshold contact distance and parasite decay rate
contact_dist <- 10
nu <- 1/(3600*24*7)# 7 days in seconds

# Import data
dat_raw <- read.csv("../data/cleaned_movement_data_07132023.csv")
# the data has 13 columns: animal id, x, y, time, dop (error), capture method,
# lat and long, age, sex, duplicate animal_id x_ y_ t_ dop method_of_capture
# lat_capture long_capture age sex duplicate_ fast_step_ fast_roundtrip_ Some of
# the first steps are fast steps (animals going back to their home range from
# capture site?). I'll remove them since they seem to be different from the
# actual home range in some cases, and they expand the grid excessively
dat <- dat_raw[!(dat_raw$fast_step_ | dat_raw$fast_roundtrip_),]

# filter to keep only a subset of individuals
ids <- c("151600", "151571","151589","151599","151592")
dat <- dat[dat$animal_id %in% ids,]


#### Create move and telemetry objects. ####

# For the deer data I have the data as a csv,
# so I need to import the individual columns. If it were a Move object you can
# use that directly.
dat_move <- move(x = dat$x_,y = dat$y_, 
                 time = as.POSIXct(dat$t_, format = "%Y-%m-%dT%H:%M:%SZ", tz="Etc/GMT-6"), 
                 data = data.frame(HDOP=dat$dop), 
                 proj = "+proj=utm +zone=16N", 
                 animal = dat$animal_id)
plot(dat_move, col = hcl.colors(length(ids), "Dark 3"),pch=16)
legend("bottomright", col = hcl.colors(length(ids), "Dark 3"), legend = ids, pch=16)

telemetries <- as.telemetry(dat_move)

#### CTMM and AKDE ####
# Fit a CTMM model to the trajectories
GUESS <- lapply(telemetries, \(i) ctmm.guess(i, interactive = F))
FITS <- lapply(seq_along(telemetries), \(i) ctmm.select(telemetries[[i]], GUESS[[i]])) 
names(FITS) <- names(telemetries)

# USE CTMM to build UDs with AKDE. 
# At this step you use the threshold distance set at the beginning as the resolution of the grid.
# epigrid <- raster(extent(dat_move), resolution = contact_dist, crs = crs(dat_move))
for (i in seq_along(FITS)) {
  # file output name
  outname <- paste0("../outputs/", names(FITS)[i], ".tif")
  # create UD
  UD <- akde(data = telemetries[[i]], CTMM = FITS[[i]], 
             grid = list(dr = c(contact_dist, contact_dist), 
                         align.to.origin = TRUE))
  # save to disk. This saves the probability mass, not the density 
  writeRaster(UD, outname, DF = "PMF")
  rm(UD)
}



#### Combine UDs and trajectories to calculate pairwise FOI ####
fnames <- paste0("../outputs/X", ids, ".tif")

# Possible combinations of individuals
combs <- combn(length(ids), 2)

# Calculate UD and SD products from UD pair values, the correlations, and
# corresponding FOIs
for (i in seq_len(ncol(combs))) {
  # read in UDs, two at a time
  ind1 <- combs[1,i]
  ind2 <- combs[2,i]
  
  r1 <- raster(fnames[ind1]) 
  r2 <- raster(fnames[ind2])
  
  # Check whether the grids overlap/extend to a common grid for both, or all
  r1 <- extend(r1, extent(merge(r1,r2)))
  r2 <- extend(r2, extent(merge(r1,r2)))
  
  # cell area
  Area <- prod(res(r1))
  # get the UD and sd products
  udprod <- r1*r2/Area
  sdprod <- sqrt(r1*(1-r1))*sqrt(r2*(1-r2))/Area
  # export outputs
  writeRaster(udprod, paste0("../outputs/UDprod_",ids[ind1],"-",ids[ind2],".tif"))
  writeRaster(sdprod, paste0("../outputs/SDprod_",ids[ind1],"-",ids[ind2],".tif"))
  
  ### CORRELATIONS
  # to estimate the correlation, I have to put the tracks in a common time
  # frame. For this, I interpolate the positions for a regular set of times
  tr1 <- range(telemetries[[ind1]]$timestamp)
  tr2 <- range(telemetries[[ind2]]$timestamp)
  
  if(max(tr1[1],tr2[1])>min(min(tr1[2],tr2[2]))) {# Check if there is temporal overlap
    cat("There is no temporal overlap between", ids[ind1], "and", ids[ind2])
    foi_ab <- foi_ba <- beta/Area*lam*(1/nu*udprod*Area)
  } else {
    tseq <- seq(max(tr1[1],tr2[1]),min(tr1[2],tr2[2]), "30 mins")
    lags <- as.numeric(tseq-min(tseq))
    nsteps <- length(tseq)
    # Interpolate trajectories
    interp_traj_1 <- predict(telemetries[[ind1]], FITS[[ind1]], t = tseq) 
    interp_traj_2 <- predict(telemetries[[ind2]], FITS[[ind2]], t = tseq)
    
    # get position histories 
    pos1 <- cellFromXY(r1, xy = as.matrix(interp_traj_1[,c("x","y")], ncol = 2))
    pos2 <- cellFromXY(r2, xy = as.matrix(interp_traj_2[,c("x","y")], ncol = 2))
    
    # keep only cells that both visited at some point
    ovlpcells <- pos1[pos1 %in% pos2]
    if(length(ovlpcells)==0) {
      cat("There are no overlap cells between", ids[ind1], "and", ids[ind2])
      foi_ab <- foi_ba <- beta/Area*lam*(1/nu*udprod*Area)
    } else {
      maxlag <- nsteps-1
      cormat_ab <- cormat_ba <- matrix(0, nrow = nsteps, ncol = length(ovlpcells))
      for (j in seq_along(ovlpcells)) {
        cell <- ovlpcells[j]
        cellx <- xFromCell(r1,cell)
        celly <- yFromCell(r1,cell)
        a <- b <- numeric(nsteps)
        a[match(cell, pos1)] <- b[match(cell, pos2)]<- 1
        xcorr <- ccf(a,b,lag.max = maxlag, plot = F)
        xcorr_vals <- as.numeric(xcorr$acf)
        cormat_ab[,j] <- rev(xcorr_vals[1:nsteps])
        cormat_ba[,j] <- xcorr_vals[nsteps:length(xcorr_vals)]
      }
      dimnames(cormat_ab) <- dimnames(cormat_ba) <- list(lag = tseq-min(tseq), cell = ovlpcells)
      # Export 
      write.csv(cormat_ab, paste0("../outputs/correlations_",ids[ind1],"-",ids[ind2],".csv"))
      write.csv(cormat_ba, paste0("../outputs/correlations_",ids[ind2],"-",ids[ind1],".csv"))
      
      ## FOI
      # scale and integrate correlation at every cell
      corcells <- ovlpcells
      corvals_ab <- corvals_ba <- numeric(length(udprod))
      corvals_ab[corcells] <- colSums(cormat_ab*exp(-nu*lags))
      corvals_ba[corcells] <- colSums(cormat_ba*exp(-nu*lags))
      corrast_ab <- corrast_ba <- udprod
      values(corrast_ab) <- corvals_ab
      values(corrast_ba) <- corvals_ba
      # Calculate FOI
      foi_ab <- beta/Area*lam*(1/nu*udprod*Area+sdprod*Area*corrast_ab)
      foi_ba <- beta/Area*lam*(1/nu*udprod*Area+sdprod*Area*corrast_ba)
      
      # Keep only positive values
      foi_ab <- foi_ab*(foi_ab>=0)
      foi_ba <- foi_ba*(foi_ba>=0)
    }
  }
  # export
  writeRaster(foi_ab, paste0("../outputs/FOI_",ids[ind1],"-",ids[ind2],".tif"))
  writeRaster(foi_ba, paste0("../outputs/FOI_",ids[ind2],"-",ids[ind1],".tif"))
}

#### Correlation analysis ####
cors <- lapply(list.files("../outputs/","corr", full.names = T), read.csv)
# Check the length of all correlations
sapply(cors, ncol)
# There is a wide variation, between just 151 rows to more than 5000 rows.
# Similarly, there is variation in spatial overlap, between only 3 cells with
# overlap, and up to 2500 cells

# Find the maximum correlation and corresponding lag
sapply(cors, \(x) apply(x, 2, which.max)) |> sapply(min)


#### Visualize total pairwise FOI ####
totfois <- lapply(list.files("../outputs/", "FOI(.*)tif$", full.names = T), raster)|>sapply(cellStats,sum)
totfoipairs <- list.files("../outputs/", "FOI(.*)tif$")|>strsplit("[[:punct:]]+")|>sapply("[",c(2,3))
totfoisdf <- data.frame(t(rbind(totfoipairs,totfois)))
names(totfoisdf) <- c("ind1","ind2","FOI")

# compare to overlap
totfoisdf$overlap <- overlap(FITS)$CI[,,2][lower.tri(overlap(FITS)$CI[,,2]) | upper.tri(overlap(FITS)$CI[,,2])]
totfoisdf

plot(totfoisdf$overlap, totfoisdf$FOI, pch = 16, las = 1, xlab = "Home range overlap", ylab = "Total FOI")
sapply(cors,nrow)
plot(sapply(cors,nrow),totfoisdf$FOI)

# Pairwise FOI plot
library(ggplot2)
ggplot(totfoisdf)+geom_raster(aes(ind1,ind2,fill=as.numeric(FOI)))+
  coord_equal()+
  theme_minimal(base_size = 16)+
  labs(x = "", y="", fill = "FOI")+
  scale_fill_gradientn(colors = hcl.colors(10, "Plasma"))

# Plot correlated pair 
plot(dat[dat$animal_id %in% c("151571","151589"),c("x_","y_")], type = "n", asp = 1, las = 1, xlab = "", ylab = "")
lines(dat[dat$animal_id =="151571",c("x_","y_")], col = hcl.colors(5, "Dark 3")[2])
lines(dat[dat$animal_id =="151589",c("x_","y_")], col = hcl.colors(5, "Dark 3")[3])

#### FOI with and without correlation ####
corpairs <- data.frame(ind1 = list.files("../outputs/", "corr") |> substr(14,19),
           ind2 = list.files("../outputs/", "corr") |> substr(21,26),
           foi_cor = 0,
           foi_ud = 0,
           dif = 0)
uds <- list.files("../outputs/", "UDprod(.*).tif$", full.names = T) 

for (i in seq_len(nrow(corpairs))) {
  ind1 <- corpairs[i,1]
  ind2 <- corpairs[i,2]
  f <- raster(paste0("../outputs/FOI_",ind1,"-",ind2,".tif"))
  ud <- beta*lam/nu*raster(uds[grepl(ind1, uds) & grepl(ind2, uds)])
  corpairs[i,"foi_cor"] <- cellStats(f,sum)
  corpairs[i,"foi_ud"] <- cellStats(ud,sum)
  corpairs[i,"dif"] <- corpairs[i,"foi_cor"]-corpairs[i,"foi_ud"]
}
corpairs <- dplyr::left_join(totfoisdf,corpairs)
corpairs
corpairs$foi_cor/corpairs$foi_ud
