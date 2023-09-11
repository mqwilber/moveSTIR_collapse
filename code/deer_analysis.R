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
library(tidyverse)
library(cowplot)
library(move)
library(ctmm)

#### Parameters ####
# Set epidemiological parameters: threshold contact distance and parasite decay rate
contact_dist <- 10 # meters
nu <- 1/(3600*24*7)# 7 days, in seconds
beta <- 1 # search efficiency, m^2/s
lam <- 1/(3600*2) # deposition rate, s^-1

#### Load and prep data ####
# Import data
dat_raw <- read.csv("../data/cleaned_movement_data_07132023.csv")
# the data has 13 columns: animal id, x, y, time, dop (error), capture method,
# lat and long, age, sex, duplicate animal_id x_ y_ t_ dop method_of_capture
# lat_capture long_capture age sex duplicate_ fast_step_ fast_roundtrip_
# I'll remove locations flagged as unrealistically fast steps
dat1 <- dat_raw[!(dat_raw$fast_step_ | dat_raw$fast_roundtrip_),]

# filter to keep only a subset of individuals
ids <- c("151564","151571","151589","151592","151599")
dat1 <- dat1[dat1$animal_id %in% ids,]

# Remove also some initial locations that were collected before the collars were
# on the animals (during storage and transport), these show up as points outside
# of the core areas on the road and in a town, east of ~ 306000
dat1 <- dat1[dat1$x_<306000,]


#### Create move and telemetry objects. ####

# For the deer data I have the data as a csv,
# so I need to import the individual columns. If it were a Move object you can
# use that directly.
dat_move <- move(x = dat1$x_,y = dat1$y_, 
                 time = as.POSIXct(dat1$t_, format = "%Y-%m-%dT%H:%M:%SZ", tz = "Etc/GMT-6"), 
                 data = data.frame(HDOP=dat1$dop), 
                 proj = "+proj=utm +zone=16N", 
                 animal = dat1$animal_id)
# Plot positions
plot(dat_move, type='n',xlab = "", ylab = "")
grid()
points(dat_move, col = hcl.colors(length(ids), "Dark 3"))
legend("bottomright", col = hcl.colors(length(ids), "Dark 3"), legend = ids, pch=1)

telemetries <- as.telemetry(dat_move)

#### CTMM and AKDE ####
# Fit a CTMM model to the trajectories
GUESS <- lapply(telemetries, \(i) ctmm.guess(i, interactive = F))
FITS <- lapply(seq_along(telemetries), \(i) ctmm.select(telemetries[[i]], GUESS[[i]])) 
names(FITS) <- names(telemetries)
# export
saveRDS(FITS, "../outputs/deer_ctmm_fits.rds")

# USE CTMM to build UDs with AKDE. 
# At this step you use the threshold distance set at the beginning as the resolution of the grid.
FITS <- readRDS("../outputs/deer_ctmm_fits.rds")
for (i in seq_along(FITS)) {
  # file output name
  outname <- paste0("../outputs/", names(FITS)[i], ".tif")
  # create UD
  UD <- akde(data = telemetries[[i]], CTMM = FITS[[i]], 
             grid = list(dr = c(contact_dist, contact_dist), 
                         align.to.origin = TRUE))
  # save to disk. This saves the probability mass, not the density 
  writeRaster(UD, outname, DF = "PMF", overwrite = TRUE)
  rm(UD)
}

#### FOI ####
# get names of files with UD grids
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
  udprod <- r1*r2
  sdprod <- sqrt(r1*(1-r1))*sqrt(r2*(1-r2))
  # export outputs
  writeRaster(udprod, paste0("../outputs/UDprod_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)
  writeRaster(sdprod, paste0("../outputs/SDprod_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)

  ### CORRELATIONS
  # to estimate the correlation, I have to put the tracks in a common time
  # frame. For this, I interpolate the positions for a regular set of times
  tr1 <- range(telemetries[[ind1]]$timestamp)
  tr2 <- range(telemetries[[ind2]]$timestamp)
  
  if(max(tr1[1],tr2[1])>min(min(tr1[2],tr2[2]))) {# Check if there is temporal overlap
    cat("There is no temporal overlap between", ids[ind1], "and", ids[ind2], "\n")
    foi_ab <- foi_ba <- beta/Area*lam*(1/nu*udprod)
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
    ovlpcells <- unique(pos1)[unique(pos1) %in% unique(pos2)]
    if(length(ovlpcells)==0) {
      cat("There are no overlap cells between", ids[ind1], "and", ids[ind2], "\n")
      foi_ab <- foi_ba <- beta/Area*lam*(1/nu*udprod)
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
      # Export scaled cross-corr rasters
      writeRaster(corrast_ab, paste0("../outputs/Corrast_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)
      writeRaster(corrast_ba, paste0("../outputs/Corrast_",ids[ind2],"-",ids[ind1],".tif"), overwrite = T)
      
      # Calculate FOI
      foi_ab <- beta/Area*lam*(1/nu*udprod+sdprod*corrast_ab)
      foi_ba <- beta/Area*lam*(1/nu*udprod+sdprod*corrast_ba)
      
      # Keep only positive values
      foi_ab <- foi_ab*(foi_ab>=0)
      foi_ba <- foi_ba*(foi_ba>=0)
    }
  }
  # export
  writeRaster(foi_ab, paste0("../outputs/FOI_",ids[ind1],"-",ids[ind2],".tif"), overwrite = T)
  writeRaster(foi_ba, paste0("../outputs/FOI_",ids[ind2],"-",ids[ind1],".tif"), overwrite = T)
}

#### Step randomization (Spiegel 2016) #### 

#One way to determine the contribution of the correlation in movement is to
#randomize the steps within each individual, and recalculate the correlation. We
#would expect if correlation is playing a sizeable role, randomizing the steps
#of tracks that are originally correlated tracks would reduce the overall FOI I
#can do this randomization using the spatsoc package, and then redo the process
#of calculating correlations and FOI. The steps are the same, so the UD and SD
#should stay constant 

#### Visualize total pairwise FOI ####
totfois <- lapply(list.files("../outputs/", "FOI(.*)tif$", full.names = T), raster)|>sapply(cellStats,sum)
totfoipairs <- list.files("../outputs/", "FOI(.*)tif$")|>strsplit("[[:punct:]]+")|>sapply("[",c(2,3))
totfoisdf <- data.frame(foi = totfois, ind1 = totfoipairs[1,], ind2 = totfoipairs[2,])

# Plot vs overlap
totfoisdf$overlap <- overlap(FITS)$CI[,,2][lower.tri(overlap(FITS)$CI[,,2]) | upper.tri(overlap(FITS)$CI[,,2])]
head(totfoisdf)

ggplot(totfoisdf)+geom_point(aes(overlap, foi))+
  labs(x = "Home range overlap",
       y = "Total pairwise FOI")+
  lims(x = c(0,1))


# Pairwise FOI plot
ggplot(totfoisdf)+geom_raster(aes(ind1,ind2,fill=as.numeric(foi)))+
  coord_equal()+
  theme_minimal(base_size = 16)+
  labs(x = "", y="", fill = "FOI")+
  scale_fill_gradientn(colors = hcl.colors(10, "Plasma"))

# Pairwise FOI plot (relative to minimum FOI)
ggplot(totfoisdf)+geom_raster(aes(ind1,ind2,fill=as.numeric(foi)/min(as.numeric(foi))))+
  coord_equal()+
  theme_minimal(base_size = 16)+
  labs(x = "", y="", fill = "Relative FOI")+
  scale_fill_gradientn(colors = hcl.colors(10, "Plasma"))
# The FOI for the pair with highly correlated movement is more than 9000 times
# greater than the lowest observed FOI


#### Correlation analysis ####
cors <- lapply(list.files("../outputs/","corr", full.names = T), read.csv, row.names = 1)
# Check the length of all correlations
sapply(cors, ncol)
# There is variation in spatial overlap, between only 1 cells with
# overlap, and up to 1064 cells

# Find the maximum correlation and corresponding lag
data.frame(name = list.files("../outputs/","corr"),
  maxcor = sapply(cors, max),
           maxclag = sapply(cors, \(x) apply(x, 2, which.max)) |> sapply(min)-1,
  ncells = sapply(cors, ncol)) %>% separate(1,into = c(NA,"ind1","ind2",NA)) %>% arrange(ncells)

# For 2 pairs, the maximum correlation is at lag 0, meaning they have direct
# encounters. For these, the maximum correlation is 1. This is likely two
# individuals encountering each other only once in a given cell. For others the
# maximum correlation is close to 1, but at higher lags. These cases correspond
# to spatial but not temporal overlap, with ~1 visits per cell for each
# individual. Finally, in two cases the maximum correlation is actually low, in
# the order of 1e-4. These correlations happened at high lags, and are likely
# produced by multiple visits to the same cell. 

sapply(cors, summary)

# It seems the correlation is small (e-2) in most cases, let's check.
sapply(cors, function(x) {
  lags <- (seq_len(nrow(x))-1)*1800
  summary(colSums(x*exp(-nu*lags)))
}) 

# The integral term with nu = 1/7days is close to 1 at the highest for a single
# cell. In most cases the effect is negative, reducing the FOI wrt overlap-only.
# The effect is additionally small, in the order of 1e-2. When you multiply by a
# sd product that is already two orders of magnitude lower than the UD product,
# then the contribution of the covariance term is very small
sapply(cors, function(x) {
  lags <- (seq_len(nrow(x))-1)*1800
  c(negcells = sum(colSums(x*exp(-nu*lags))<0),
    totcells = ncol(x))
})|> rbind(pair = substr(list.files("../outputs/","corr"), 14, 26)) |> t()

# In most cases, the cumulative scaled correlation is negative, and would reduce
# the FOI, but by how much?

#### FOI with and without correlation ####
corrastfiles <- list.files("../outputs/", "Corrast(.*).tif$", full.names = T)
udprodfiles <- list.files("../outputs/", "UDprod(.*).tif$", full.names = T) 
sdprodfiles <- list.files("../outputs/", "SDprod(.*).tif$", full.names = T) 

pdf("../docs/figures/Cov_contrib.pdf")
for (i in seq_along(corrastfiles)) {
  ind1 = substr(basename(corrastfiles[i]),9,14)
  ind2 = substr(basename(corrastfiles[i]),16,21)
  corrast <- raster(corrastfiles[i])
  Area <- prod(res(corrast))
  ud <- 1/nu*raster(udprodfiles[grepl(ind1, udprodfiles) & grepl(ind2, udprodfiles)])
  sd <- raster(sdprodfiles[grepl(ind1, sdprodfiles) & grepl(ind2, sdprodfiles)])
  covrast <- sd*corrast
  plot(covrast/ud, main = paste("Covariance/UD",ind1,ind2))
}
dev.off()

# The approach below calculates only total FOI, not by cell
corpairs <- data.frame(ind1 = list.files("../outputs/", "corr") |> substr(14,19),
           ind2 = list.files("../outputs/", "corr") |> substr(21,26),
           foi_cor = 0,
           foi_ud = 0,
           dif = 0)


for (i in seq_len(nrow(corpairs))) {
  ind1 <- corpairs[i,1]
  ind2 <- corpairs[i,2]
  f <- raster(paste0("../outputs/FOI_",ind1,"-",ind2,".tif"))
  ud <- beta*lam/nu*raster(uds[grepl(ind1, uds) & grepl(ind2, uds)])
  corpairs[i,"foi_cor"] <- cellStats(f,sum)
  corpairs[i,"foi_ud"] <- cellStats(ud,sum)
  corpairs[i,"dif"] <- corpairs[i,"foi_cor"]-corpairs[i,"foi_ud"]
}
corpairs <- left_join(totfoisdf,corpairs)
corpairs
corpairs$foi_cor/corpairs$foi_ud
# The covariance contribution is in the order of 2e-3 at the highest.



          
#### Recalculate FOI for different distances ####
# Get new UDs using a different contact distance, for example 20 m
FITS <- readRDS("../outputs/deer_ctmm_fits.rds")
UDSdeer20 <- list()
for (i in seq_along(FITS)) {
  # file output name
  outname <- paste0("../outputs/UD_", names(FITS)[i], "_20m.tif")
  # create UD
  UD <- UDSdeer20[[i]] <- akde(data = telemetries[[i]], CTMM = FITS[[i]], 
             grid = list(dr = c(2*contact_dist, 2*contact_dist), 
                         align.to.origin = TRUE))
  # save to disk. This saves the probability mass, not the density 
  writeRaster(UD, outname, DF = "PMF")
  rm(UD)
}
fnames <- paste0("../outputs/UD_X", ids, "_20m.tif")
# Get FOIS
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
  writeRaster(udprod, paste0("../outputs/UDprod_",ids[ind1],"-",ids[ind2],"_20m.tif"))
  writeRaster(sdprod, paste0("../outputs/SDprod_",ids[ind1],"-",ids[ind2],"_20m.tif"))
  
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
      write.csv(cormat_ab, paste0("../outputs/correlations_",ids[ind1],"-",ids[ind2],"_20m.csv"))
      write.csv(cormat_ba, paste0("../outputs/correlations_",ids[ind2],"-",ids[ind1],"_20m.csv"))
      
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
  writeRaster(foi_ab, paste0("../outputs/FOI_",ids[ind1],"-",ids[ind2],"_20m.tif"))
  writeRaster(foi_ba, paste0("../outputs/FOI_",ids[ind2],"-",ids[ind1],"_20m.tif"))
}
# compare with 10 m
foidistcomp <- data.frame(foi = sapply(list.files("../outputs/","FOI(.*).tif$", full.names = T), raster)|>sapply(cellStats,sum),
                          d = ifelse(grepl("20m",list.files("../outputs/","FOI(.*).tif$")),20,10),
                          pair = substr(list.files("../outputs/","FOI(.*).tif$"), 5,17)
)
foidistcomp %>% pivot_wider(names_from = d, values_from = foi,names_prefix = "d") %>% 
  mutate(absdif = (d20-d10),reldif = (d20-d10)/d10) %>% 
  summary()



#### Recalculate FOI with different nu #### 
# The effect of correlation on FOI
#is highly dependent on the decay parameter. Here I will recalculate FOI for a
#set of different parameter values, one shorter and one longer '
corfiles <- list.files("../outputs/","corr", full.names = T)
corfilepairs <- basename(corfiles)|>substr(14,26)
uds <- list.files("../outputs/", "UDprod(.*).tif$", full.names = T) 
sds <- list.files("../outputs/", "SDprod(.*).tif$", full.names = T) 
nus <- 1/(3600*24*c(1/24,1/12,1/6, 1/3,1/2,1,3,7))
fois <- expand.grid(pair = corfilepairs,nu = nus, foicor = 0,foiud = 0)
cnt=1
for (j in seq_along(nus)) {
  nup <- nus[j]
  for (i in seq_along(corfiles)) {
    ind1 <- basename(corfiles[i])|>substr(14,19)
    ind2 <- basename(corfiles[i])|>substr(21,26)
    # get correlation file, corresponding ud and sd products
    corrs <- read.csv(corfiles[i], row.names = 1)
    lags <- as.numeric(row.names(corrs))
    corrint <- colSums(exp(-nup*lags)*corrs)
    ud <- Area/nup*raster(uds[grepl(ind1, uds) & grepl(ind2, uds)])
    sd <- Area*raster(sds[grepl(ind1, sds) & grepl(ind2, sds)])
    corrast <- ud
    values(corrast) <- 0
    corcells <- as.numeric(substring(names(corrs),2))
    values(corrast)[corcells] <- corrint
    foi <- beta*lam/Area*(ud+sd*corrast)
    foi <- foi*(foi>=0)
    fois[cnt,3] <- cellStats(foi,sum)
    fois[cnt,4] <- cellStats(beta*lam/Area*ud,sum)
    cnt=cnt+1
  }
}
#### FIGURES ####
# Plot deer FOI
pdf("../outputs/deer_FOI.pdf", width = 11)
par(mfrow = c(2,2))
for (i in 1:ncol(combs)) {
  # get IDs
  ind1 <- ids[combs[1,i]]
  ind2 <- ids[combs[2,i]]
  ud <- if (file.exists(paste0("../outputs/UDprod_",ind1,"-",ind2,".tif"))) {
    raster(paste0("../outputs/UDprod_",ind1,"-",ind2,".tif"))
  }else {raster(paste0("../outputs/UDprod_",ind2,"-",ind1,".tif"))}
  
  sd <- if (file.exists(paste0("../outputs/SDprod_",ind1,"-",ind2,".tif"))) {
    raster(paste0("../outputs/SDprod_",ind1,"-",ind2,".tif"))
  }else{raster(paste0("../outputs/SDprod_",ind2,"-",ind1,".tif"))}
  ud <- ud*prod(res(ud))/nu
  sd <- sd*prod(res(sd))
  foi <- raster(paste0("../outputs/FOI_",ind1,"-",ind2,".tif"))
  r[is.na(r)] <- 0
  
  plot(telemetries[combs[,i]], col = hcl.colors(5,"Dark 3")[combs[,i]], main = paste("Tracks", ind1, ind2))
  plot(ud,main = paste("scaled UD product", ind1, ind2))
  plot(sd,main = paste("SD product", ind1, ind2))
  plot(foi, main = paste("FOI",ind1,ind2))
}
dev.off()

# PLot of FOI as a function of nu
fois %>% 
  ggplot()+geom_line(aes(1/(nu*3600*24),foicor, color = pair), show.legend = F)+
  theme_classic(base_size = 24)+
  labs(x = "Time to decay (days) ",
       y =  "Force of Infection")
# PLot of difference in FOI with and without correlation, for different values
# of nu. 
fois %>% mutate(dfoi = (foicor-foiud)/foiud) %>% 
  ggplot()+geom_line(aes(1/(nu*24*3600),dfoi, color = pair), show.legend = F)+
  theme_classic(base_size = 24)+
  labs(x = "Time to decay (days) ",
       y = "Relative contribution",
       color = "Individuals")
# How much does FOI vary?
fois %>% #filter(nu == 1/(3600*24)|nu == 1/(3600*24*14)) %>% # select two decay rates
  group_by(pair) %>% mutate(doi = foicor/min(foicor)) %>% view()

# Combined figure
p1 <- ggplot(dat1)+geom_point(aes(x_,y_, color = factor(animal_id)), shape=1)+
  coord_equal()+
  theme_minimal(base_size = 10)+
  labs(x = "Easting (m)", y = "Northing (m)", color = "Animal ID")+
  theme(legend.position = c(0.95,0.05), legend.justification = c(0.9,0.1))
p2 <- ggplot(totfoisdf)+geom_raster(aes(ind1,ind2,fill=foi))+
  coord_equal()+
  theme_minimal(base_size = 10)+
  labs(x = "", y="", fill = "FOI")+
  scale_fill_gradientn(colors = hcl.colors(10, "Plasma"))+
  theme(legend.key.size = unit(1/8,"inch"), legend.text = element_text(size = 6))
p3 <- ggplot(totfoisdf)+geom_point(aes(overlap, foi))+
  theme_classic(base_size = 10)+
  labs(x = "Home range overlap",
       y = "Total pairwise FOI")+
  lims(x = c(0,1))
# p4 <- fois %>% 
#   ggplot()+geom_line(aes(1/(nu*3600*24),foicor, color = pair), show.legend = F)+
#   theme_classic(base_size = 8)+
#   labs(x = "Time to decay (days) ",
#        y =  "Force of Infection")
# p5 <- fois %>% mutate(dfoi = (foicor-foiud)/foiud) %>% 
#   ggplot()+geom_line(aes(1/(nu*24*3600),dfoi, color = pair), show.legend = F)+
#   theme_classic(base_size = 8)+
#   labs(x = "Time to decay (days) ",
#        y = "Relative contribution",
#        color = "Individuals")
pdf("../docs/figures/deer_results.pdf", width = 8,height = 6)
plot_grid(p1,plot_grid(p3,p2, labels = c("b", "c"), ncol=1, axis = "lr", align = "v", rel_heights = c(0.8,1)), labels = "a", rel_widths = c(1,0.7))
dev.off()
