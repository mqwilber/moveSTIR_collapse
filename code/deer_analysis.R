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
deer_corrs_pw_flt <- list()
for (i in seq_len(ncol(combs))) {
  # read in UDs, two at a time
  ind1 <- combs[1,i]
  ind2 <- combs[2,i]
  
  r1 <- raster(fnames[ind1]) 
  r2 <- raster(fnames[ind2])
  

  # get UD and SD products
  prods <- getProds(r1,r2)
  # cell area
  Area <- prod(res(r1))
  
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
  positions <- getPositions(X = interp_trajs, R = prods)
  
  # new code ... works
  cormats <- pairCorrs(X = positions, prewt = T, fltr = "reg")
  deer_corrs_pw_flt[[paste0(ind1,"-",ind2)]] <- cormats
}

saveRDS(deer_corrs_pw_flt,"outputs/deer_corrs_prewht_flt.rds")


# There is no spatial overlap (at the cell scale) between 151583 and 151571, 83
# and 75, 83 and 89, and 71 with 75

#### Expected correlation ####
# Correlations can arise spuriously due to chance. In order to estimate the
# expected correlation at different lags one should fit another model.
# Correlation is a random variable, with random effects of the pair, and the lag
# We can fit a mixed effects model to estimate the expected correlation.

# Import all correlations
comb_corrs <- sapply(deer_corrs_pw_flt, function(x) {
  if(!any(is.na(x))) {
    data.frame(lag = rownames(x[[1]])[1:1008],
               cor = c(as.numeric(x[[1]][1:1008,]), as.numeric(x[[2]][1:1008,])))
  }
})
comb_corrs <- comb_corrs[!sapply(comb_corrs,is.null)]
comb_corrs_df <- do.call(rbind,comb_corrs)
comb_corrs_df$pair <- factor(substr(rownames(comb_corrs_df),1,3))
rm(comb_corrs)

# create mixed effect model
library(lme4)
LMM <- lmer(cor~(1|lag)+(1|pair),comb_corrs_df)

# We can see if there is a relationship between the number of visits and the
# estimated correlation. I am assuming that few visits can result in
# artificially high correlations that do not represent the expected correlation.
corvisits <- sapply(seq_along(deer_corrs_pw_flt), function(i) {
  X <- deer_corrs_pw_flt[[i]]
  UDP <- raster(udprodfiles[i])
  if(!any(is.na(X))) {
    cab <- X$CAB
    cellUD <- values(UDP)[as.numeric(colnames(cab))]
    cumcor <- colSums(cab)
    vmin <- apply(X$nvisits, 1, min)
    vtot <- colSums(X$nvisits)
    vmax <- apply(X$nvisits,1,max)
    vprod <- apply(X$nvisits,1,prod)
    data.frame(cumcor, vmax, vmin, vprod,cellUD)
  } 
})

corvisitsdf <- do.call(rbind,corvisits)
corvisitsdf$pair <- rep(names(deer_corrs_pw_flt)[!sapply(corvisits,is.null)],unlist(sapply(corvisits, nrow)))
corvisitsdf$cell <- rownames(corvisitsdf)
ggplot(corvisitsdf, aes(vmin, cumcor/cellUD, color = pair))+geom_point()+
  labs(x = "Number of visits", y = "Cumulative correlation", color = "Pair")+
  scale_y_log10()+
  scale_x_log10()+
  theme_light(base_size = 18)
cor(corvisitsdf)

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
  cors <- deer_corrs_pw_flt[[i]]

  # corfiles <- c(paste0("outputs/correlations_10min_",id1,"-",id2,"_1109.csv"),
  #               paste0("outputs/correlations_10min_",id2,"-",id1,"_1109.csv"))
  nus <- c(SARS = 1/3600, CWD = 3.04e-9) # decay rates in seconds^-1
  Area <- prod(res(udprod)) # cell area in m^2
  lagt <- 60*24*3600 # period to consider for exact integral scaling UD product
  foiudsars <- beta*lam/Area*1/nus[1]*udprod
  foiudcwd <- beta*lam/Area*(1-exp(-nus[2]*lagt))/nus[2]*udprod
  foidirect <- beta*lam/Area*udprod
  for (j in 1:2) {
    # if (file.exists(corfiles[j])) {
    CorrExists <- !all(is.na(cors))
    if(CorrExists) {
      cat("using correlations for", id1, id2,"\n")
      # corrs <- read.csv(corfiles[j], row.names = 1)
      # keep only non-spurious correlations. Which
      # ones are spurious is determined based on a CI threshold equal to
      # 0+-1.96/sqrt(n), where n is the length of the series. This has already been done in the correlation calculation
      nlags <- nrow(corrs*2) # times 2 because I only used up to half the maximum lag
      CI_THRESH <- 1.96/sqrt(nlags) 
      # index of which cells have more "significant" correlations than expected
      # by chance if series were independent
      lags <- (as.numeric(row.names(cors[[j]]))-1)*600
      # lags <- as.numeric(row.names(corrs))
      dtau <- unique(diff(lags))
      corrintSARS <- colSums(exp(-nus[1]*lags)*cors[[j]]*dtau)
      corrintCWD <- colSums(exp(-nus[2]*lags)*cors[[j]]*dtau)
      corrastCWD <- corrastSARS <- udprod
      values(corrastCWD) <- values(corrastSARS) <- 0
      # corcells <- as.numeric(substring(names(corrs),2))
      corcells <- as.numeric(dimnames(cors[[j]])$cell)
      corrastCWD[corcells] <- corrintCWD
      corrastSARS[corcells] <- corrintSARS
      foiCWDrast <- beta*lam/Area*(udprod*(1-exp(-nus[2]*lagt))/nus[2]+sdprod*corrastCWD)
      foiSARSrast <- beta*lam/Area*(udprod*1/nus[1]+sdprod*corrastSARS)
      raster::plot(foiSARSrast)
    } else {
      cat("using UD only for",id1, "and", id2,"\n")
      foiCWDrast <- foiudcwd
      foiSARSrast <- foiudsars
    }
    foiCWDrast <- max(foiCWDrast,0)
    foiSARSrast <- max(foiSARSrast,0)
    idi <- ifelse(j == 1, id1,id2)
    idj <- ifelse(j == 1, id2,id1)
    deer_FOIs_nus_prewt <- rbind(deer_FOIs_nus_prewt, c(idi,idj, CorrExists,
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
write.csv(deer_FOIs_nus_prewt, "outputs/deer_totFOI_CWDSARS.csv", row.names = F, quote = F)

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
deer_FOIs_nus_prewt2 <- deer_FOIs_nus_prewt %>% arrange(ind1,ind2) %>% 
  add_column(FOI_hr_cwd = HR_FOI_CWD[complete.cases(as.numeric(HR_FOI_CWD))],
             FOI_hr_sars = HR_FOI_SARS[complete.cases(as.numeric(HR_FOI_SARS))],
             overlap = as.numeric(overlap(deer_UDs)$CI[,,2])[-c(1,7,13,19,25)]) %>% 
  mutate(across(starts_with("FOI"),as.numeric))
deer_FOIs_nus_prewt2

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
grf <- graph_from_data_frame(deer_FOIs_nus_prewt2[,c(1,2)])

# Specify layout customly
tkplot(grf)
custom_grf_coords <- tk_coords(1)
custom_grf_coords

co <- layout_with_fr(grf,coords = custom_grf_coords,niter = 1)
# pdf("../docs/figures/networks.pdf", width = 8,height = 4, family = "sans")
# x11(width = 11, height = 8.5,family = "HersheySans")



####---- R_0 ----####
# Given the total pairwise expected FOI, we can calculate the R0 as
# the spectral radius (the greatest absolute eigenvalue) of the matrix R = F*U
# F are the pairwise FOIs, and U is the inverse of a diagonal matrix with values -lambda
Fmat_SARS <- deer_FOIs_nus_prewt2 %>% select(ind1,ind2,FOI_SARS) %>% 
  pivot_wider(names_from = ind2, values_from = FOI_SARS, names_sort = T, values_fill = 0) %>% 
  column_to_rownames("ind1") %>% as.matrix()
Fmat_UD <- deer_FOIs_nus_prewt2 %>% select(ind1,ind2,FOI_UD_SARS) %>% 
  pivot_wider(names_from = ind2, values_from = FOI_UD_SARS, names_sort = T, values_fill = 0) %>% 
  column_to_rownames("ind1") %>% as.matrix()
Fmat_HR <- deer_FOIs_nus_prewt2 %>% select(ind1,ind2,FOI_hr_sars) %>% 
  pivot_wider(names_from = ind2, values_from = FOI_hr_sars, names_sort = T, values_fill = 0) %>% 
  column_to_rownames("ind1") %>% as.matrix()
gSARS <- 1/(5*24*3600) # recovery rate is 5 days in seconds

R0_SARS <- calcR0(X = Fmat_SARS, g = gSARS) 
R0_UD <- calcR0(X = Fmat_UD, g = gSARS) 
R0_HR <- calcR0(X = Fmat_HR, g = gSARS) 

####---- VISUALIZATION ----####
# Plot combined figure
x11(width = 9, height = 6)
# pdf("docs/figures/deer_results.pdf", width = 9,height = 6)

mycols <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77")
png("docs/figures/deer_results.png", 9,6,"in",res = 300)
par(omi = c(3,0.2,0.1,5.5),mar = c(1,2,0,0),new=T)
#raster::plot(dat_move,type='n', xlab = "Easting (m)", ylab = "Northing (m)", cex.axis = 0.8, cex.lab = 0.8)
plot.new()
plot.window(c(299500,302600), c(3884500,3888900), asp = 1, xaxs = "i",yaxs = "i")
plot(deer_UDs, col.level=NA, col.DF=NA,units=F, ann = F,xaxt = 'n',yaxt = 'n', bty='n')
grid()
axis(1,cex.axis = 0.6, col = 'gray', tcl = -0.3, padj = -2)
axis(2,cex.axis = 0.6, col = 'gray', tcl = -0.3, padj = 2)
box()
mtext("Easting (m)", 1, 1, cex = 0.8)
mtext("Northing (m)", 2, 1, cex = 0.8)
plot(telemetries, UD=deer_UDs, 
     error=F, col=mycols,col.DF = mycols,
     col.grid=NA,
     col.level = mycols,
     level=NA,lwd.level = 1.5, labels=NA,
     add=T)
# Network plots
refwidth <- min(deer_FOIs_nus_prewt2[deer_FOIs_nus_prewt2$FOI_hr_sars>0,"FOI_hr_sars"])
refsize <- min(tapply(deer_FOIs_nus_prewt2$FOI_SARS,deer_FOIs_nus_prewt2$ind2,sum))
par(mar=c(2,2,2,2), omi = c(0.2,0.5,3.2,6),new=T)
with(list(db = deer_FOIs_nus_prewt2, 
          Vsize = tapply(deer_FOIs_nus_prewt2$FOI_SARS, deer_FOIs_nus_prewt2$ind2,sum)),{
            plot(grf, layout = co, 
                 # vertex.label.dist = 2, vertex.label.degree = pi/2, vertex.label.color = "black",vertex.label.family = "sans",
                 edge.width = log(deer_FOIs_nus_prewt2$FOI_SARS/refwidth), edge.arrow.size=0, 
                 vertex.color=mycols, vertex.label = NA, vertex.size = 3*log(Vsize/refsize))
            mtext("FOI with correlation", 3, 0, cex = 1.3)
            mtext(bquote(R[0] == .(format(R0_SARS$R0,digits = 2))), 1,0.5)
            # mtext("c",2,at = 2, cex=2,las=1, line=2)
            })


par(mar=c(2,2,2,2), omi = c(0.2,3.25,3.2,3.25),new=T)
with(list(db = deer_FOIs_nus_prewt2, 
          Vsize = tapply(deer_FOIs_nus_prewt2$FOI_UD_SARS, deer_FOIs_nus_prewt2$ind2,sum)),{
            plot(grf, layout = co, 
                 # vertex.label.dist = 2, vertex.label.degree = pi/2, vertex.label.color = "black",vertex.label.family = "sans",
                 edge.width = log(deer_FOIs_nus_prewt2$FOI_UD_SARS/refwidth), edge.arrow.size=0, 
                 vertex.color=mycols, vertex.label = NA, vertex.size = 3*log(Vsize/refsize))
            mtext("FOI without correlation", 3, 0, cex = 1.3)
            mtext(bquote(R[0] == .(format(R0_UD$R0,digits = 2))), 1,0.5)
          })

par(mar=c(2,2,2,2), omi = c(0.2,6,3.2,0.5),new=T)
with(list(db = deer_FOIs_nus_prewt2, 
          Vsize = tapply(deer_FOIs_nus_prewt2$FOI_hr_sars, deer_FOIs_nus_prewt2$ind2,sum)),{
            plot(grf, layout = co, 
                 # vertex.label.dist = 2, vertex.label.degree = pi/2, vertex.label.color = "black",vertex.label.family = "sans",
                 edge.width = log(deer_FOIs_nus_prewt2$FOI_hr_sars/refwidth), edge.arrow.size=0, 
                 vertex.color=mycols, vertex.label = NA, vertex.size = 3*log(Vsize/refsize))
            mtext("FOI HR overlap", 3, 0,cex=1.3)
            mtext(bquote(R[0] == .(format(R0_HR$R0,digits = 2))), 1,0.5)
          })

# plot FOI surface for 71-89 and 71-99
par(omi = c(3,3.2,0.1,3),new=T, cex = 0.8)
with(list(udprod = raster("outputs/UDprod_151571-151589.tif"),
          sdprod = raster("outputs/SDprod_151571-151589.tif")), {
            corrast <- udprod
            values(corrast) <- 0
            # cors <- read.csv("outputs/correlations_10min_151571-151589.csv", row.names = 1)
            cors <- deer_corrs_pw_flt[[3]]$CAB
            lags <- as.numeric(rownames(cors))*600
            intcors <- colSums(cors*exp(-1/3600*lags)*600)
            # corcells <- as.numeric(substring(names(cors), 2))
            corcells <- as.numeric(dimnames(cors)$cell)
            corrast[corcells] <- intcors
            foirast <- max(beta*lam/100*(udprod/nus[1]+sdprod*corrast),0)
            plot(dat1[dat1$animal_id %in% c(151571,151589),c("x_","y_")], 
                 ann = F, asp = 1, type = 'n', xaxt='n', yaxt='n', bty='n')
            raster::plot(foirast*3600*24, col=hcl.colors(25,'Oranges', rev=T), legend.width = 1.5,
                         legend.shrink = 0.5, lwd = 0,
                         legend.args = list(text = expression(FOI (days^-1)), side = 3, adj=0), 
                         add=T)
            sp::plot(polygons(hr95[[1]])[2], add = T, border = mycols[1], col = NA, lwd = 3)
            sp::plot(polygons(hr95[[4]])[2], add = T, border = mycols[4], col = NA, lwd = 3)
            # points(dat1[dat1$animal_id==151571,c('x_','y_')], col = mycols[1], cex = 0.5)
            # points(dat1[dat1$animal_id==151589,c('x_','y_')], col = mycols[4], cex = 0.5)
          })


par(mar=c(2,2,0,0), omi = c(3,6.6,0.5,0.3),new=T)
with(list(db = deer_FOIs_nus_prewt2), {
  plot(db$overlap, db$FOI_SARS/db$FOI_UD_SARS, xlim = c(0,1), type = 'n',log = "y",ann = F, xaxt = 'n',yaxt = 'n', cex.axis=0.8)
  grid()
  points(db$overlap, db$FOI_CWD/db$FOI_UD_CWD, bg = pathcols[2], pch=21)
  points(db$overlap, db$FOI_SARS/db$FOI_UD_SARS, bg = pathcols[1], pch=21)
  axis(1,cex.axis = 0.8, tcl = -0.3, padj = -1.5)
  axis(2,cex.axis = 0.8,  tcl = -0.3, hadj = 0.5, las=1)
  mtext("Home Range Overlap", 1,1.1)
  mtext("FOI ratio",2, 1.2)
})

# add labels
par(oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot.new()
plot.window(c(0,1),c(0,1))
text(x = c(0,0.4,0.77,0), 
     y = c(0.97,0.97,0.97,0.43), 
     labels = c("a","b","c","d"), cex = 1.4, font = 2)
dev.off()


with(list(db = deer_FOIs_nus_prewt2), {
  refval <- median(c(db$FOI_SARS,db$FOI_CWD))
  par(mar = c(5,5,1,2)+0.1)
  plot(db$overlap, db$FOI_CWD/refval, xlab = "Home Range Overlap", ylab = "Relative FOI", col = pathcols[2], pch=16, xlim = c(0,1), las=1, cex=1.2, log='y',ann=F)
  points(db$overlap, db$FOI_SARS/refval, col = pathcols[1], pch=16, cex=1.2)
  mtext("Relative FOI",side=2,line = 3.5)
})



# link between UD product and correlation



# Pairwise FOI plot
ggplot(deer_FOIs_nus_prewt)+
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
# ggplot(deer_FOIs_nus_prewt2)+
#   geom_boxplot(aes(x="CWD",y = (FOI_CWD-FOI_UD_CWD)/FOI_CWD), color = cols2[1], outlier.shape = NA)+
#   geom_boxplot(aes(x="SARS",y = (FOI_SARS-FOI_UD_SARS)/FOI_SARS), color = cols2[2], outlier.shape = NA)+
#   geom_jitter(aes(x="CWD",y = (FOI_CWD-FOI_UD_CWD)/FOI_CWD), col = cols2[1])+
#   geom_jitter(aes(x="SARS",y = (FOI_SARS-FOI_UD_SARS)/FOI_SARS), col = cols2[2])+
#   labs(x = "Pathogen", y = "FOI ratio")+
#   # scale_y_log10()+
#   theme_pubr()+
#   scale_x_discrete(labels = c("CWD","SARS-CoV-2"))

ggplot(deer_FOIs_nus_prewt2)+
  geom_hline(yintercept = 1, linetype =2)+
  geom_boxplot(aes(x="CWD",y = FOI_CWD/FOI_UD_CWD), color = cols2[1])+
  geom_boxplot(aes(x="SARS",y = FOI_SARS/FOI_UD_SARS), color = cols2[2])+
  # geom_jitter(aes(x="CWD",y = (FOI_CWD-FOI_UD_CWD)/FOI_CWD), col = hcl.colors(2,"BluGrn")[1])+
  # geom_jitter(aes(x="SARS",y = (FOI_SARS-FOI_UD_SARS)/FOI_SARS), col = hcl.colors(2,"BluGrn")[2])+
  labs(x = "Pathogen persistence", y = "FOI ratio")+
  # scale_y_log10()+
  # scale_y_continuous(limits = c(0,12))+
  theme_pubclean(base_size = 14)+
  scale_x_discrete(labels = c("Long (CWD)","Short (SARS-CoV-2)"))

# FOI with and without correlation vs overlap
deer_FOIs_nus_prewt2 %>% ggplot()+
  geom_point(aes(overlap, FOI_SARS), col="darkred")+
  geom_point(aes(overlap, FOI_UD_SARS), col="steelblue")+
  scale_y_log10()+
  theme_pubclean(base_size = 14)+
  labs(x = "Home Range Overlap", y="FOI")

# FOI ratio wrt overlap
pathcols <- hcl.colors(2, "Dynamic")
deer_FOIs_nus_prewt2 %>% ggplot()+
  geom_point(aes(overlap, FOI_SARS/FOI_UD_SARS), color = pathcols[1])+
  geom_point(aes(overlap, FOI_CWD/FOI_UD_CWD), color = pathcols[2])+
  scale_y_log10()+
  theme_pubclean(base_size = 14)+
  labs(x = "Home Range Overlap", y="FOI ratio")


# Combined figure
# tracks
p1 <- ggplot(dat1)+geom_point(aes(x_,y_, color = factor(animal_id)), size=0.5)+
  coord_equal()+
  theme_minimal(base_size = 10)+
  labs(x = "Easting (m)", y = "Northing (m)", color = "Animal ID")+
  theme(legend.position = c(0.95,0.05), legend.justification = c(0.9,0.1))+
  scale_color_manual(values = mycols)
# matrix plot
p2 <- ggplot(deer_FOIs_nus_prewt2)+geom_raster(aes(ind1,ind2,fill=FOI_SARS/FOI_UD_SARS))+
  coord_equal()+
  theme_minimal(base_size = 10)+
  labs(x = "", y="", fill = "FOI ratio")+
  scale_fill_gradientn(colors = hcl.colors(25, "Plasma"),trans = "log10")#+
  theme(legend.key.size = unit(1/8,"inch"), legend.text = element_text(size = 6))

plot_grid(p1,plot_grid(p3,p2, labels = c("b", "c"), ncol=1, axis = "lr", align = "v", rel_heights = c(0.8,1)), labels = "a", rel_widths = c(1,0.7))
