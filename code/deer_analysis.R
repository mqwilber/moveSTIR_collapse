# Analysis of deer data using the pMoveSTIR framework
# Author: Juan S. Vargas Soto
# Date: July 2023

#### Introduction ####

# We have developed a framework to estimate the expected FOI from movement data.
# GPS tracking data is described using ctmm, and from that we build utilization
# distributions. These are the basis for calculating the overlap (as product of
# UDs) and the product of the SDs at every cell. We then calculate the location
# histories and the cross-correlation for every cell.

# Load libraries
library(tidyverse)
library(ggpubr)
library(cowplot)
library(move)
library(ctmm)

#### Parameters ####
# Set epidemiological parameters: threshold contact distance and parasite decay rate
contact_dist <- 10 # meters
nu <- 1/(3600*24*7)# 7 days, in seconds
beta <- 1 # search efficiency, m^2/s
lam <- 1/(3600*24) # deposition rate, s^-1

#### Load and prep data ####
# Import data
dat_raw <- read.csv("../data/cleaned_movement_data_07132023.csv")
# the data has 13 columns: animal id, x, y, time, dop (error), capture method,
# lat and long, age, sex, duplicate animal_id x_ y_ t_ dop method_of_capture
# lat_capture long_capture age sex duplicate_ fast_step_ fast_roundtrip_
# I'll remove locations flagged as unrealistically fast steps
dat1 <- dat_raw[!(dat_raw$fast_step_ | dat_raw$fast_roundtrip_),]

# filter to keep only a subset of individuals
ids <- c("151571","151575","151583","151589","151599")
dat1 <- dat1[dat1$animal_id %in% ids,]

# Remove also some initial locations that were collected before the collars were
# on the animals (during storage and transport), these show up as points outside
# of the core areas on the road and in a town, east of ~ 306000
dat1 <- dat1[dat1$x_<306000,]

# Transform the time to POSIXct, and filter out the first positions, so that we
# are only using the core home range area.
dat1$datetime <- as.POSIXct(dat1$t_, format = "%Y-%m-%dT%H:%M:%SZ", tz = "Etc/GMT-6")

# Check the range of dates for each animal

# Plot wrt to positions, to see when animals get to their home ranges and which
# dates should be deleted
ggplot(dat1, aes(x_,y_,color=months(datetime)))+geom_point()+facet_wrap(~animal_id)

# I will keep only May and June
dat1 <- subset(dat1, dat1$datetime>="2023-05-01"& dat1$datetime<"2023-07-01")

#### Create move and telemetry objects. ####

# For the deer data I have the data as a csv,
# so I need to import the individual columns. If it were a Move object you can
# use that directly.
dat_move <- move(x = dat1$x_,y = dat1$y_, 
                 time = dat1$datetime, 
                 data = data.frame(HDOP=dat1$dop), 
                 proj = "+proj=utm +zone=16N", 
                 animal = dat1$animal_id)

# Create telemetry object and plot positions
telemetries <- as.telemetry(dat_move)
plot(telemetries, col = hcl.colors(5, "Dark 3"), ann=F, units=F,xaxt='n',yaxt='n')
grid()
legend("bottomright", col = hcl.colors(length(ids), "Dark 3"), legend = paste("Ov",1:5), pch=1, cex=0.8)
axis(1, cex.axis=0.8)
axis(2, cex.axis=0.8)

#### CTMM and AKDE ####
# Fit a CTMM model to the trajectories
GUESS <- lapply(telemetries, ctmm.guess, interactive = F)
FITS <- mapply(ctmm.select, telemetries, GUESS, SIMPLIFY = FALSE)
# export
saveRDS(FITS, "outputs/deer_ctmm_fits.rds")

# USE CTMM to build UDs with AKDE. 
# At this step you use the threshold distance set at the beginning as the resolution of the grid.
# FITS <- readRDS("outputs/deer_ctmm_fits.rds")

deer_UDs <- akde(telemetries, FITS,grid = list(dr = c(contact_dist,contact_dist), align.to.origin = T))
# clusterMap(cl, akde, telemetries,FITS, grid = list(dr = c(contact_dist,contact_dist), align.to.origin = T))
outnames <- paste0("outputs/UD_", names(FITS), ".tif")
# Export UD PMFs
mapply(writeRaster,UD,outnames, DF = "PMF", overwrite = TRUE)
# Export 95% home range polygons
outnames <- paste0("outputs/HRpoly_", names(FITS), ".shp")
mapply(writeShapefile,UD,outnames)
hr95 <- lapply(UD, SpatialPolygonsDataFrame.UD)

#### Correlations ####
# get names of files with UD grids
fnames <- paste0("outputs/UD_X", ids, ".tif")

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
  r1 <- extend(r1, extent(merge(r1,r2)), value=0)
  r2 <- extend(r2, extent(merge(r1,r2)), value=0)
  
  # cell area
  Area <- prod(res(r1))
  # get the UD and sd products
  udprod <- r1*r2
  sdprod <- sqrt(r1*(1-r1))*sqrt(r2*(1-r2))
  # if(export) {
  #   # export outputs
  #   writeRaster(udprod, paste0("outputs/UDprod_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)
  #   writeRaster(sdprod, paste0("outputs/SDprod_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)
  # }

  ### CORRELATIONS
  
  # to estimate the correlation, I have to put the tracks in a common time
  # frame. For this, I interpolate the positions for a regular set of times
# new code... works
  interp_trajs <- interpTrajs(telemetries[[ind1]], telemetries[[ind2]])

  # get position histories ... works
  positions <- getPositions(X = interp_trajs, R = list(r1,r2))
  
  # new code ... works
  cormats <- getCorrs(X = positions,) 
# old code to delete ###
        # ovlpcells <- unique(pos1)[unique(pos1) %in% unique(pos2)]
    # if(length(ovlpcells)==0) {
    #   cat("There are no overlap cells between", ids[ind1], "and", ids[ind2], "\n")
    #   # foi_ab <- foi_ba <- beta/Area*lam*(1/nu*udprod)
    # } else {
    #   maxlag <- nsteps
    #   cormat_ab <- cormat_ba <- matrix(0, nrow = nsteps, ncol = length(ovlpcells))
    #   for (j in seq_along(ovlpcells)) {
    #     cell <- ovlpcells[j]
    #     a <- b <- numeric(nsteps)
    #     a[cell==pos1] <- b[cell==pos2]<- 1
    #     xcorr <- ccf(a,b,lag.max = maxlag, plot = F)
    #     xcorr_vals <- as.numeric(xcorr$acf)
    #     cormat_ab[,j] <- xcorr_vals[nsteps:1]
    #     cormat_ba[,j] <- xcorr_vals[nsteps:length(xcorr_vals)]
    #   }
    #   dimnames(cormat_ab) <- dimnames(cormat_ba) <- list(lag = lags, cell = ovlpcells)
    #   # Export 
    #   write.csv(cormat_ab, paste0("outputs/correlations_10min_",ids[ind1],"-",ids[ind2],".csv"))
    #   write.csv(cormat_ba, paste0("outputs/correlations_10min_",ids[ind2],"-",ids[ind1],".csv"))
      
  #     # ## FOI
  #     # # scale and integrate correlation at every cell
  #     # corcells <- ovlpcells
  #     # corvals_ab <- corvals_ba <- numeric(length(udprod))
  #     # corvals_ab[corcells] <- colSums(cormat_ab*exp(-nu*lags)*unique(diff(lags)))
  #     # corvals_ba[corcells] <- colSums(cormat_ba*exp(-nu*lags)*unique(diff(lags)))
  #     # corrast_ab <- corrast_ba <- udprod
  #     # values(corrast_ab) <- corvals_ab
  #     # values(corrast_ba) <- corvals_ba
  #     # # Export scaled cross-corr rasters
  #     # writeRaster(corrast_ab, paste0("outputs/Corrast_10min_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)
  #     # writeRaster(corrast_ba, paste0("outputs/Corrast_10min_",ids[ind2],"-",ids[ind1],".tif"), overwrite = T)
  #     # 
  #     # Calculate FOI
  #     foi_ab <- beta/Area*lam*(1/nu*udprod+sdprod*corrast_ab)
  #     foi_ba <- beta/Area*lam*(1/nu*udprod+sdprod*corrast_ba)
  #     
  #     # Keep only positive values
  #     foi_ab <- foi_ab*(foi_ab>=0)
  #     foi_ba <- foi_ba*(foi_ba>=0)
    }
  # }
  # # export
  # writeRaster(foi_ab, paste0("outputs/FOI_10min_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)
  # writeRaster(foi_ba, paste0("outputs/FOI_10min_",ids[ind2],"-",ids[ind1],".tif"), overwrite = T)
}

# There is no spatial overlap (at the cell scale) between 151583 and 151571, 83
# and 75, 83 and 89, and 71 with 75

#### FOI #### 

#The effect of correlation on FOI is highly dependent on the decay parameter.
#Here we calculate FOI for a set of different decay rates. Deer in
#our study area have two parasites of interest, SARS-CoV-2, which is transmitted
#in close proximity over short periods, and Chronic Wasting Disease, which is a
#transmitted through a prion that can survive long times in the environment. We
#will recalculate the FOIs for the deer in our study for these two cases. For
#SARS-CoV-2 we assume transmission occurs within one hour, through saliva
#droplets. For CWD, the probability of persistence from one year to the next has
#been estimated as 0.91 (Cook et al. 2022), which translates to a decay rate of
#nu = 0.096/yr, or 3.04e-9/s. Analyzing this requires to use the exact expression
#of the decay function rather than assuming integrating over infinite lags. In
#practice it is virtually the same as assuming no decay over the study period

udprodfiles <- list.files("outputs/", "UDprod(.*)[[:digit:]].tif$", full.names = T) 
sdprodfiles <- list.files("outputs/", "SDprod(.*)[[:digit:]].tif$", full.names = T) 
deer_FOIs_nus_prewt <- data.frame()
for (i in 1:ncol(combs)) {
  ind1 = combs[1,i]
  ind2 = combs[2,i]
  id1 = ids[ind1]
  id2 = ids[ind2]
  
  udprod <- raster(udprodfiles[grepl(id1, udprodfiles) & grepl(id2, udprodfiles)])
  sdprod <- raster(sdprodfiles[grepl(id1, sdprodfiles) & grepl(id2, sdprodfiles)])
  corfiles <- c(paste0("outputs/correlations_10min_",id1,"-",id2,"_1109.csv"),
                paste0("outputs/correlations_10min_",id2,"-",id1,"_1109.csv"))
  nus <- c(SARS = 1/3600, CWD = 3.04e-9) # decay rates in seconds^-1
  Area <- prod(res(udprod)) # cell area in m^2
  lagt <- 60*24*3600 # period to consider for exact integral scaling UD product
  foiudsars <- beta*lam/Area*1/nus[1]*udprod
  foiudcwd <- beta*lam/Area*(1-exp(-nus[2]*lagt))/nus[2]*udprod
  foidirect <- beta*lam/Area*udprod
  for (j in 1:2) {
    if (file.exists(corfiles[j])) {
      cat("using correlations for", id1, id2,"\n")
      corrs <- read.csv(corfiles[j], row.names = 1)
      # keep only correlations at cells with non-spurious correlations. Which
      # ones are spurious is determined based on a CI threshold equal to
      # 0+-1.96/sqrt(n), where n is the length of the series. If the number of
      # correlations that exceed the threshold is greater than 5% of the series
      # then the correlations are kept for the whole cell. Alternatively, only
      # those specific values are kept, and the rest discarded.
      nlags <- nrow(corrs*2) # times 2 because I only used up to half the maximum lag
      CI_THRESH <- 1.96/sqrt(nlags) 
      # index of which cells have more "significant" correlations than expected
      # by chance if series were independent
      sigCells <- colSums(abs(corrs)>CI_THRESH)>=(0.05*nlags) # 0.05 is the signif. threshold corresponding to 1.96
      cat(sum(sigCells),"/",length(sigCells), "(",100*sum(sigCells)/length(sigCells), "%) cells have non spurious correlations\n")
      # replace values in cells that did not have more than 5% of correlations beyond threshold.
      # corrs[,!sigCells] <- 0
      # lags <- (as.numeric(row.names(corrs))-1)*600
      lags <- as.numeric(row.names(corrs))
      dtau <- unique(diff(lags))
      corrintSARS <- colSums(exp(-nus[1]*lags)*corrs*dtau)
      corrintCWD <- colSums(exp(-nus[2]*lags)*corrs*dtau)
      corrastCWD <- corrastSARS <- udprod
      values(corrastCWD) <- values(corrastSARS) <- 0
      corcells <- as.numeric(substring(names(corrs),2))
      corrastCWD[corcells] <- corrintCWD
      corrastSARS[corcells] <- corrintSARS
      foiCWDrast <- beta*lam/Area*(udprod*(1-exp(-nus[2]*lagt))/nus[2]+sdprod*corrastCWD)
      foiSARSrast <- beta*lam/Area*(udprod*1/nus[1]+sdprod*corrastSARS)
      raster::plot(foiSARSrast)
    }else {
      cat("using UD only for",id1, "and", id2,"\n")
      foiCWDrast <- foiudcwd
      foiSARSrast <- foiudsars
    }
    foiCWDrast <- max(foiCWDrast,0)
    foiSARSrast <- max(foiSARSrast,0)
    idi <- ifelse(j == 1, id1,id2)
    idj <- ifelse(j == 1, id2,id1)
    deer_FOIs_nus_prewt <- rbind(deer_FOIs_nus_prewt, c(idi,idj, file.exists(corfiles[j]),
                                            cellStats(foiCWDrast, sum),
                                            cellStats(foiSARSrast, sum),
                                            cellStats(foiudcwd, sum),
                                            cellStats(foiudsars, sum),
                                            cellStats(foidirect, sum)
                                            )
                           )
  }
}

names(deer_FOIs_nus_prewt) <- c("ind1", "ind2", "correlation","FOI_CWD","FOI_SARS", 
                          "FOI_UD_CWD", "FOI_UD_SARS", "FOI_direct")
deer_FOIs_nus_prewt

# The FOI experienced by deer differs orders of magnitude depending on the
# parasite of interest and its respective decay rate. For parasites that persist
# for longer in the environment like CWD, the FOI is three orders of magnitude
# higher than for a parasite with fast decay rates. 

####  FOI for different distances ####
# Get new UDs using a different contact distance, for example 20 m
FITS <- readRDS("../outputs/deer_ctmm_fits.rds")
UDSdeer20 <- list()
for (i in seq_along(FITS)) {
  # file output name
  outname <- paste0("outputs/UD_", names(FITS)[i], "_20m.tif")
  # create UD
  UD <- UDSdeer20[[i]] <- akde(data = telemetries[[i]], CTMM = FITS[[i]], 
             grid = list(dr = c(2*contact_dist, 2*contact_dist), 
                         align.to.origin = TRUE))
  # save to disk. This saves the probability mass, not the density 
  writeRaster(UD, outname, DF = "PMF")
}; rm(UD, outname)
fnames <- paste0("outputs/UD_X", ids, "_20m.tif")
# Get UD prods, SD prods, and correlations at 20 m
for (i in seq_len(ncol(combs))) {
  # read in UDs, two at a time
  ind1 <- combs[1,i]
  ind2 <- combs[2,i]
  
  r1 <- raster(fnames[ind1]) 
  r2 <- raster(fnames[ind2])
  
  # Check whether the grids overlap/extend to a common grid for both, or all
  r1 <- extend(r1, extent(merge(r1,r2)), value = 0)
  r2 <- extend(r2, extent(merge(r1,r2)), value = 0)
  
  # cell area
  Area <- prod(res(r1))
  # get the UD and sd products
  udprod <- r1*r2
  sdprod <- sqrt(r1*(1-r1))*sqrt(r2*(1-r2))
  # export outputs
  writeRaster(udprod, paste0("outputs/UDprod_",ids[ind1],"-",ids[ind2],"_20m.tif"))
  writeRaster(sdprod, paste0("outputs/SDprod_",ids[ind1],"-",ids[ind2],"_20m.tif"))
  
  ### CORRELATIONS
  # to estimate the correlation, I have to put the tracks in a common time
  # frame. For this, I interpolate the positions for a regular set of times
  tr1 <- range(telemetries[[ind1]]$timestamp)
  tr2 <- range(telemetries[[ind2]]$timestamp)
  # Check if there is temporal overlap
  if(max(tr1[1],tr2[1])>min(min(tr1[2],tr2[2]))) {
    cat("There is no temporal overlap between", ids[ind1], "and", ids[ind2])
    foi_ab <- foi_ba <- beta/Area*lam*(1/nu*udprod)
  } else {
    tseq <- seq(max(tr1[1],tr2[1]),min(tr1[2],tr2[2]), "10 mins")
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
    } else {
      maxlag <- nsteps
      cormat_ab <- cormat_ba <- matrix(0, nrow = nsteps, ncol = length(ovlpcells))
      for (j in seq_along(ovlpcells)) {
        cell <- ovlpcells[j]
        a <- b <- numeric(nsteps)
        xcorr <- ccf(a,b,lag.max = maxlag, plot = F)
        xcorr_vals <- as.numeric(xcorr$acf)
        cormat_ab[,j] <- (xcorr_vals[nsteps:1])
        cormat_ba[,j] <- xcorr_vals[nsteps:length(xcorr_vals)]
      }
      dimnames(cormat_ab) <- dimnames(cormat_ba) <- list(lag = lags, cell = ovlpcells)
      # Export 
      write.csv(cormat_ab, paste0("outputs/correlations_",ids[ind1],"-",ids[ind2],"_20m.csv"))
      write.csv(cormat_ba, paste0("outputs/correlations_",ids[ind2],"-",ids[ind1],"_20m.csv"))
    }
  }
}; rm(fnames,cormat_ab,cormat_ba,nsteps,maxlag,a,b,Area,interp_traj_1,interp_traj_2)
# Calculate the FOIs for the nu of SARS
udprodfiles <- list.files("outputs", "UDprod(.*)20m", full.names = T)
sdprodfiles <- list.files("outputs", "SDprod(.*)20m", full.names = T)
deer_FOIs_20m <- data.frame()
for (i in 1:ncol(combs)) {
  ind1 = combs[1,i]
  ind2 = combs[2,i]
  id1 = ids[ind1]
  id2 = ids[ind2]
  
  udprod <- raster(udprodfiles[grepl(id1, udprodfiles) & grepl(id2, udprodfiles)])
  sdprod <- raster(sdprodfiles[grepl(id1, sdprodfiles) & grepl(id2, sdprodfiles)])
  corfiles <- c(paste0("outputs/correlations_",id1,"-",id2,"_20m.csv"),
                paste0("outputs/correlations_",id2,"-",id1,"_20m.csv"))
  nus <- 1/3600 # decay rates in seconds^-1
  Area <- prod(res(udprod))
  lagt <- 60*24*3600
  foiud <- beta*lam/Area*1/nus*udprod
  for (j in 1:2) {
    if (file.exists(corfiles[j])) {
      corrs <- read.csv(corfiles[j], row.names = 1)
      lags <- as.numeric(row.names(corrs))
      dtau <- unique(diff(lags))
      corrint <- colSums(exp(-nus*lags)*corrs*dtau)
      corrast <- udprod
      values(corrast) <- 0
      corcells <- as.numeric(substring(names(corrs),2))
      corrast[corcells] <- corrint
      foirast <- beta*lam/Area*(udprod*1/nus+sdprod*corrast)
    }else {
      foirast <- foiud
    }
    foirast <- max(foirast,0)
    idi <- ifelse(j == 1, id1,id2)
    idj <- ifelse(j == 1, id2,id1)
    deer_FOIs_20m <- rbind(deer_FOIs_20m, c(idi,idj, file.exists(corfiles[j]),
                                            cellStats(foirast, sum),
                                            cellStats(foiud, sum)
    )
    )
  }
}; rm(idi,idj,foirast,foiud, corrs, lags)
names(deer_FOIs_20m) <- c("ind1", "ind2", "correlation", "FOI_SARS_20", "FOI_UD_SARS_20")


# Doubling the contact distance from 10 to 20 m would mean that Ax increases
# 4-fold, from 100 to 400, and FOI should decrease by 75% at the local scale.
# This would nonetheless be offset by the corresponding reduction in the number
# of cells. For the deer, we observe some increases (at least assuming SARS
# biology), ranging from 0.2% to 108% wrt FOI calculated at 10 m distance.
# Decreases ranged between less than 1% to 99.4%. Greater base FOIs are
# correlated with a greater decrease in FOI. We observed the greatest decrease
# for the pair that had the greatest FOI.
left_join(deer_FOIs,deer_FOIs_20m, by = c("ind1","ind2")) %>% 
  mutate(across(starts_with("FOI"),as.numeric)) %>% 
  mutate(dratio = FOI_SARS_20/FOI_SARS,
         dratioUD = FOI_UD_SARS_20/FOI_UD_SARS) %>% 
  ggplot(aes(FOI_SARS,dratio))+geom_point()+scale_x_log10()

####---- Home Range overlap FOI ----####

# We calculate the FOI for Home Range overlap using equation 8 in the paper. For
# this, we need to calculate the HR areas of each individual, as well as the
# area of overlap
# hr95 <- lapply(readRDS("outputs/deer_HR_95.rds"), polygons)
# get areas
hr95_areas <- sapply(hr95, \(x) x@polygons[[2]]@area)
# create overlap areas matrix
ovlpAreas <- matrix(numeric(length(hr95)^2),ncol = length(hr95))
# dimnames(ovlpAreas) <- list(names(hr95_areas),names(hr95_areas))
# fill in the matrix
for(i in 1:5) {
  for(j in 1:5){
    hr1 <- polygons(hr95[[i]])[2]
    hr2 <- polygons(hr95[[j]])[2]
    if(rgeos::gIntersects(hr1,hr2)) {
      intersection <- raster::intersect(hr1,hr2)
      ovlpAreas[i,j] <- intersection@polygons[[1]]@area
    }
  }
}

# this can be used, multiplied by the epi parameters, to estimate an FOI in the
# same units (equation 8)
nus <- c(SARS = 1/3600, CWD = 3.04e-9) # decay rates in seconds^-1
HR_FOI_SARS <- beta*lam*1/nus[1]*ovlpAreas/outer(hr95_areas,hr95_areas)
HR_FOI_CWD <- beta*lam*(1-exp(-nus[2]*3600*24*60))/nus[2]*ovlpAreas/outer(hr95_areas,hr95_areas)
diag(HR_FOI_SARS) <- NA
diag(HR_FOI_CWD) <- NA

# The total FOI with the covariance term is obtained by the sum of the cell FOIs
# obtained previously. We put it all together in a single dataframe
deer_FOIs <- deer_FOIs_nus %>% arrange(ind1,ind2) %>% 
  add_column(FOI_hr_cwd = HR_FOI_CWD[complete.cases(as.numeric(HR_FOI_CWD))],
             FOI_hr_sars = HR_FOI_SARS[complete.cases(as.numeric(HR_FOI_SARS))],
             overlap = as.numeric(overlap(deer_UDs)$CI[,,2])[-c(1,7,13,19,25)]) %>% 
  mutate(across(starts_with("FOI"),as.numeric))
deer_FOIs

#For a parasite with short persistence like SARS, the FOI can be three orders of
#magnitude higher than for CWD which persists longer in the environment. This is
#mainly due to the effect of covariance. For a few pairs, the covariance
#component makes up more than 99% of the whole FOI. The covariance also
#decreases FOI in some cases. For CWD, the covariance always decreased the FOI
#estimation, with respect to UD only. It seems also covariance inverts the
#expectation with respect to which is more likely to transmit


          
####---- Network analysis ----####

# We can use the estimated FOIs to assign edge weight to a contact network among
# the deer. We can compare this network with a network built based on Home Range
# overlap, and compare the networks with and without a covariance term, as well
# as for direct and indirect transmission

## Visualize networks
library(igraph)
grf <- graph_from_data_frame(deer_FOIs[,c(1,2)])

# Specify layout customly
tkplot(grf)
custom_grf_coords <- tk_coords(1)
custom_grf_coords

co <- layout_with_fr(grf,coords = custom_grf_coords,niter = 1)
# pdf("../docs/figures/networks.pdf", width = 8,height = 4, family = "sans")
# x11(width = 11, height = 8.5,family = "HersheySans")

####---- VISUALIZATION ----####
pdf(file.path(tempfile(fileext = ".pdf")))

# plot FOI surface
udprod <- raster("outputs/UDprod_151571-151589.tif")
sdprod <- raster("outputs/SDprod_151571-151589.tif")
corrast <- udprod
values(corrast) <- 0
cors <- read.csv("outputs/correlations_10min_151571-151589.csv", row.names = 1)
lags <- as.numeric(rownames(cors))
intcors <- colSums(cors*exp(-1/3600*lags)*600)
corcells <- as.numeric(substring(names(cors), 2))
corrast[corcells] <- intcors
foirast <- beta*lam/100*(udprod/nus[1]+sdprod*corrast)
plot(dat1[dat1$animal_id %in% c(151571,151589),c("x_","y_")], cex.axis=0.8,ann = F, asp = 1, type = 'n', xaxt='n',yaxt='n')
points(dat1[dat1$animal_id==151571,c('x_','y_')], col = mycols[1], cex = 0.5)
points(dat1[dat1$animal_id==151589,c('x_','y_')], col = mycols[4], cex = 0.5)
raster::plot(foirast, alpha = 0.8, add=T, col=hcl.colors(25,'Rocket',rev=T))
grid()
raster::plot(udprod,add=T,alpha=0.8, col)
# link between UD product and correlation



# Pairwise FOI plot
ggplot(deer_FOIs)+
  geom_raster(aes(ind1,ind2,fill=FOI_SARS))+
  coord_equal()+
  scale_x_discrete(labels = paste("Ov",1:5))+
  scale_y_discrete(labels = paste("Ov",1:5))+
  theme_classic(base_size = 14)+
  theme(axis.line = element_blank())+
  labs(x = "", y="", fill = "FOI")+
  scale_fill_gradientn(colors = hcl.colors(20, "Plasma"), trans = "log10")

# FOI w/cov/FOI/without for 2 pathogens
cols2 <- hcl.colors(2,"Dynamic")
ggplot(deer_FOIs)+
  geom_boxplot(aes(x="CWD",y = (FOI_CWD-FOI_UD_CWD)/FOI_CWD), color = cols2[1], outlier.shape = NA)+
  geom_boxplot(aes(x="SARS",y = (FOI_SARS-FOI_UD_SARS)/FOI_SARS), color = cols2[2], outlier.shape = NA)+
  geom_jitter(aes(x="CWD",y = (FOI_CWD-FOI_UD_CWD)/FOI_CWD), col = cols2[1])+
  geom_jitter(aes(x="SARS",y = (FOI_SARS-FOI_UD_SARS)/FOI_SARS), col = cols2[2])+
  labs(x = "Pathogen", y = "FOI ratio")+
  # scale_y_log10()+
  theme_pubr()+
  scale_x_discrete(labels = c("CWD","SARS-CoV-2"))

ggplot(deer_FOIs)+
  geom_boxplot(aes(x="CWD",y = FOI_CWD/FOI_UD_CWD), color = cols2[1])+
  geom_boxplot(aes(x="SARS",y = FOI_SARS/FOI_UD_SARS), color = cols2[2])+
  # geom_jitter(aes(x="CWD",y = (FOI_CWD-FOI_UD_CWD)/FOI_CWD), col = hcl.colors(2,"BluGrn")[1])+
  # geom_jitter(aes(x="SARS",y = (FOI_SARS-FOI_UD_SARS)/FOI_SARS), col = hcl.colors(2,"BluGrn")[2])+
  labs(x = "Pathogen", y = "FOI ratio")+
  scale_y_log10()+
  theme_pubr()+
  scale_x_discrete(labels = c("CWD","SARS-CoV-2"))

# FOI with and without correlation vs overlap
deer_FOIs %>% ggplot()+
  geom_point(aes(overlap, FOI_SARS), col="darkred")+
  geom_point(aes(overlap, FOI_UD_SARS), col="steelblue")+
  scale_y_log10()+
  theme_minimal(base_size = 14)+
  labs(x = "Home Range Overlap", y="FOI")

# FOI ratio wrt overlap
pathcols <- hcl.colors(2, "Dynamic")
deer_FOIs %>% ggplot()+
  geom_point(aes(overlap, FOI_SARS/FOI_UD_SARS), color = pathcols[1])+
  geom_point(aes(overlap, FOI_CWD/FOI_UD_CWD), color = pathcols[2])+
  scale_y_log10()+
  theme_minimal(base_size = 14)+
  labs(x = "Home Range Overlap", y="FOI ratio")


# Combined figure
# tracks
p1 <- ggplot(dat1)+geom_point(aes(x_,y_, color = factor(animal_id)), size=0.5)+
  coord_equal()+
  theme_minimal(base_size = 10)+
  labs(x = "Easting (m)", y = "Northing (m)", color = "Animal ID")+
  theme(legend.position = c(0.95,0.05), legend.justification = c(0.9,0.1))
# matrix plot
p2 <- ggplot(deer_FOIs)+geom_raster(aes(ind1,ind2,fill=foi))+
  coord_equal()+
  theme_minimal(base_size = 10)+
  labs(x = "", y="", fill = "FOI")+
  scale_fill_gradientn(colors = hcl.colors(10, "Plasma"))+
  theme(legend.key.size = unit(1/8,"inch"), legend.text = element_text(size = 6))

pdf("../docs/figures/deer_results.pdf", width = 8,height = 6)
plot_grid(p1,plot_grid(p3,p2, labels = c("b", "c"), ncol=1, axis = "lr", align = "v", rel_heights = c(0.8,1)), labels = "a", rel_widths = c(1,0.7))
dev.off()

# Plot combined figure
x11(width = 9, height = 6)
pdf("../docs/figures/deer_results.pdf", width = 9,height = 6)
mycols <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77")
par(omi = c(3,0.5,0.5,5.5),mai = c(0.5,0.5,0,0))
#raster::plot(dat_move,type='n', xlab = "Easting (m)", ylab = "Northing (m)", cex.axis = 0.8, cex.lab = 0.8)
plot(deer_UDs, col.level=NA, col.DF=NA,units=F, ann = F,xaxt = 'n',yaxt = 'n')
grid()
axis(1,cex.axis = 0.6, col = 'gray')
axis(2,cex.axis = 0.6, col = 'gray')
mtext("Easting (m)", 1, line = 2, cex = 0.8)
mtext("Northing (m)", 2, line = 2, cex=0.8)
plot(telemetries, UD=deer_UDs, 
     error=F, col=mycols,col.DF = mycols,
     col.grid=NA,
     col.level = mycols,
     level=NA,lwd.level = 1.5, labels=NA,
     add=T)
mtext("a", 2,line=3,cex=2,padj=0,las=1,at = 3890000)

par(mar=c(2,2,2,2), omi = c(0.5,0.5,3,6),new=T)
plot(delete.edges(grf,which(deer_FOIs$FOI_SARS==0)), layout = co, 
     vertex.label.dist = 2, vertex.label.degree = pi/2, vertex.label.color = "black",vertex.label.family = "sans",
     edge.width = scale(deer_FOIs[,5]/median(deer_FOIs[,5]), center = F), edge.arrow.size=0, 
     vertex.color=mycols, vertex.label = NA)
mtext("FOI with correlation", 3, 0.5)
# text(-1,1,"c", cex=2)

par(mar=c(2,2,2,2), omi = c(0.5,3,3,3),new=T)
plot(delete.edges(grf,which(deer_FOIs$FOI_UD_SARS==0)), layout = co, 
     vertex.label.dist = 2, vertex.label.degree = pi/2, vertex.label.color = "black",vertex.label.family = "sans",
     edge.width = scale(deer_FOIs[,7]/min(deer_FOIs[,5]), center = F), edge.arrow.size=0,
     vertex.color=mycols, vertex.label = NA)
mtext("FOI without correlation", 3, 0.5)
# text(-1,1,"d", cex=2)

par(mar=c(2,2,2,2), omi = c(0.5,6,3,0.5),new=T)
plot(delete.edges(grf,which(deer_FOIs$foi_hr_sars==0)), layout = co, 
     vertex.label.dist = 2, vertex.label.degree = pi/2, vertex.label.color = "black",vertex.label.family = "sans",
     edge.width = scale(deer_FOIs[,10]/min(deer_FOIs[,5]), center = F), edge.arrow.size=0,
     vertex.color=mycols, vertex.label = NA)
mtext("FOI HR overlap", 3, 0.5)
# text(-1,1,"e", cex=2)

par(mar=c(4,4,1,0), omi = c(3,4,0.5,3),new=T)
plot(deer_FOIs$overlap, deer_FOIs$FOI_SARS/deer_FOIs$FOI_UD_SARS, log = "y", cex.axis = 0.8, xlab = "Home Range Overlap", ylab = "FOI ratio", col = pathcols[1], pch=16, cex.axis=0.8, xlim = c(0,1), las=1, cex.lab=0.8)
points(deer_FOIs$overlap, deer_FOIs$FOI_CWD/deer_FOIs$FOI_UD_CWD, col = pathcols[2], pch=16)
grid()

dev.off()

#### FUNCTIONS ####
r1 <- raster(fnames[ind1]) 
r2 <- raster(fnames[ind2])
getProds <- function(r1,r2) {
  # Check whether the grids overlap/extend to a common grid for both, or all
  r1 <- extend(r1, extent(merge(r1,r2)), value=0)
  r2 <- extend(r2, extent(merge(r1,r2)), value=0)
  # UD and SD products
  udprod <- r1*r2
  sdprod <- sqrt(r1*(1-r1))*sqrt(r2*(1-r2))
  list(UD = udprod, SD = sdprod)
}
# function to get vector of positions (cell index) for two individuals at a
# time. The grid is common to both
getPositions <- function(X,R) {
  pos1 <- cellFromXY(R[[1]], xy = as.matrix(X[[1]][,c("x","y")], ncol = 2))
  pos2 <- cellFromXY(R[[2]], xy = as.matrix(X[[2]][,c("x","y")], ncol = 2))
  list(pos1,pos2)
}
# function to get the correlation matrices for a single pair. Outputs are the
# correlation from a->b, from b->a, and the proportion of "significant"
# correlations within the cell, assessed by confidence intervals. There is an
# option to prewhiten, or to estimate CI using bootstrapping
getCorrs <- function(X, export = T, ci = c("none", "prewt", "bs")) {
  OUT <- list()
  pos1 <- X[[1]]
  pos2 <- X[[2]]
  ovlpcells <- unique(pos1)[unique(pos1) %in% unique(pos2)]
  if(length(ovlpcells)==0) {
    cat("There are no overlap cells\n")
  } else {
    nsteps <- length(pos1)
    lagmax <- ceiling(nsteps/2)
    cormat_ab <- cormat_ba <- matrix(0, nrow = lagmax+1, ncol = length(ovlpcells))
    for (j in seq_along(ovlpcells)) {
      cell <- ovlpcells[j]
      a <- b <- numeric(nsteps)
      a[cell==pos1] <- 1
      b[cell==pos2] <- 1
      xcorr <- switch (ci[1],
                       none = ccf(a, b, lag.max = maxlag, plot = F),
                       pw = TSA::prewhiten(a,b,lag.max = maxlag, plot = F)$ccf,
                       bs = ccf(a, b, lag.max = maxlag, plot = F)
      )
      xcorr_vals <- as.numeric(xcorr$acf)
      sigcells[j] <- switch(ci[1], 
                            none = 1,
                            pw = mean(abs(xcorr_vals)>(1.96/sqrt(xcorr$n.used))),
                            bs = corboot(AB = cbind(a,b), cors = xcorr))
      
      cormat_ab[,j] <- xcorr_vals[(lagmax+1):1]
      cormat_ba[,j] <- xcorr_vals[(lagmax+1):length(xcorr_vals)]
    }
    # dimnames(cormat_ab) <- dimnames(cormat_ba) <- list(lag = 1:(lagmax+1), cell = ovlpcells)
    # Export
    if(export) {
      write.csv(cormat_ab, paste0("outputs/correlations_10min_",ids[ind1],"-",ids[ind2],"_",format(Sys.Date(), "%m%d"), ".csv"))
      write.csv(cormat_ba, paste0("outputs/correlations_10min_",ids[ind2],"-",ids[ind1],"_",format(Sys.Date(), "%m%d"), ".csv"))
    }
    return(list(CAB = cormat_ab, CBA = cormat_ba, psig = sigcells))
  }
}

# Old function
# getCorrs <- function(XY, grd, write = FALSE) {
#   # read in UDs, two at a time
#   ind1 <- combs[1,i]
#   ind2 <- combs[2,i]
#   
#   # export outputs
#   writeRaster(udprod, paste0("outputs/UDprod_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)
#   writeRaster(sdprod, paste0("outputs/SDprod_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)
#   ### CORRELATIONS
#   
#   # to estimate the correlation, I have to put the tracks in a common time
#   # frame. For this, I interpolate the positions for a regular set of times
#   interptrajs <- interpTrajs(telemetries, ind1,ind2)
#   
#   # get position histories 
#   pos1 <- cellFromXY(r1, xy = as.matrix(interptrajs[[1]][,c("x","y")], ncol = 2))
#   pos2 <- cellFromXY(r2, xy = as.matrix(interptrajs[[2]][,c("x","y")], ncol = 2))
#   
#   # keep only cells that both visited at some point
#   ovlpcells <- unique(pos1)[unique(pos1) %in% unique(pos2)]
#   if(length(ovlpcells)==0) {
#     cat("There are no overlap cells between", ids[ind1], "and", ids[ind2], "\n")
#     foi_ab <- foi_ba <- beta/Area*lam*(1/nu*udprod)
#   } else {
#     maxlag <- nsteps
#     cormat_ab <- cormat_ba <- matrix(0, nrow = nsteps, ncol = length(ovlpcells))
#     for (j in seq_along(ovlpcells)) {
#       cell <- ovlpcells[j]
#       a <- b <- numeric(nsteps)
#       a[cell==pos1] <- b[cell==pos2]<- 1
#       xcorr <- ccf(a,b,lag.max = maxlag, plot = F)
#       xcorr_vals <- as.numeric(xcorr$acf)
#       cormat_ab[,j] <- xcorr_vals[nsteps:1]
#       cormat_ba[,j] <- xcorr_vals[nsteps:length(xcorr_vals)]
#     }
#     dimnames(cormat_ab) <- dimnames(cormat_ba) <- list(lag = tseq-min(tseq), cell = ovlpcells)
#     # Export 
#     if(write) {
#       write.csv(cormat_ab, paste0("outputs/correlations_10min_",ids[ind1],"-",ids[ind2],".csv"))
#       write.csv(cormat_ba, paste0("outputs/correlations_10min_",ids[ind2],"-",ids[ind1],".csv"))
#     }
#   }
# }

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