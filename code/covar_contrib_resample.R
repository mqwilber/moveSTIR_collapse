# Covariance contribution, resampling method

# One way to determine the covariance contribution, and whether it is higher
# than expected by chance, is to resample the detection histories. If there is
# temporal correlation, the scaled cross-corr should be different from random

# get the correlations
corfiles <- list.files("../outputs/", "corr(.*)30min", full.names = T)
corrastfiles <- list.files("../outputs/", "Corrast(.*)30min(.*).tif$", full.names = T)
udprodfiles <- list.files("../outputs/", "UDprod(.*).tif$", full.names = T) 
sdprodfiles <- list.files("../outputs/", "SDprod(.*).tif$", full.names = T) 
pdf("../outputs/covcontribresampled.pdf")
covcontrib_df <- data.frame(ind1 = numeric(), ind2 = numeric(),
                            iteration = numeric(),
                            totfoi1 = numeric(), totfoi2 = numeric(),
                            CovContr_mean1 = numeric(),
                            CovContr_mean2 = numeric(),
                            CovContr_sd1 = numeric(),
                            CovContr_sd2 = numeric())
                            
for (i in seq_along(corfiles)) {
  # read file
  d <- read.csv(corfiles[i])
  # scale and add the correlations for every column
  # corrast <- raster(corrastfiles[i])
  
  id1 = substr(basename(corfiles[i]),20,25)
  id2 = substr(basename(corfiles[i]),27,31)
  ind1 <- match(id1,ids)
  ind2 <- match(id2,ids)
  
  ud <- 1/nu*raster(udprodfiles[grepl(id1, udprodfiles) & grepl(id2, udprodfiles)])
  sd <- raster(sdprodfiles[grepl(id1, sdprodfiles) & grepl(id2, sdprodfiles)])
  # covrast <- sd*corrast
  
  # hist(max((ud+covrast)/ud,0)[(max((ud+covrast)/ud,0))!=1],main = paste("(Covariance+UD)/UD",ind1,ind2))
  
  ### CORRELATIONS
  # to estimate the correlation, I have to put the tracks in a common time
  # frame. For this, I interpolate the positions for a regular set of times
  tr1 <- range(telemetries[[ind1]]$timestamp)
  tr2 <- range(telemetries[[ind2]]$timestamp)
  
  tseq <- seq(max(tr1[1],tr2[1]),min(tr1[2],tr2[2]), "30 mins")
  lags <- as.numeric(tseq-min(tseq))
  nsteps <- length(tseq)
  # Interpolate trajectories
  interp_traj_1 <- predict(telemetries[[ind1]], FITS[[ind1]], t = tseq) 
  interp_traj_2 <- predict(telemetries[[ind2]], FITS[[ind2]], t = tseq)
  
  # get position histories 
  pos1 <- cellFromXY(ud, xy = as.matrix(interp_traj_1[,c("x","y")], ncol = 2))
  pos2 <- cellFromXY(ud, xy = as.matrix(interp_traj_2[,c("x","y")], ncol = 2))
  # Resample
  for (it in 1:40) {
    if(it>1) {
      pos1 <- sample(pos1)
      pos2 <- sample(pos2)
    }
    
    # keep only cells that both visited at some point
    ovlpcells <- unique(pos1)[unique(pos1) %in% unique(pos2)]
    
    maxlag <- nsteps-1
    cormat_ab <- cormat_ba <- matrix(0, nrow = nsteps, ncol = length(ovlpcells))
    for (j in seq_along(ovlpcells)) {
      cell <- ovlpcells[j]
      cellx <- xFromCell(ud,cell)
      celly <- yFromCell(ud,cell)
      a <- b <- numeric(nsteps)
      a[match(cell, pos1)] <- b[match(cell, pos2)]<- 1
      xcorr <- ccf(a,b,lag.max = maxlag, plot = F)
      xcorr_vals <- as.numeric(xcorr$acf)
      cormat_ab[,j] <- rev(xcorr_vals[1:nsteps])
      cormat_ba[,j] <- xcorr_vals[nsteps:length(xcorr_vals)]
      
      dimnames(cormat_ab) <- dimnames(cormat_ba) <- list(lag = tseq-min(tseq), cell = ovlpcells)
    }
    
    ## FOI
    # scale and integrate correlation at every cell
    corcells <- ovlpcells
    corvals_ab <- corvals_ba <- numeric(length(ud))
    corvals_ab[corcells] <- colSums(cormat_ab*exp(-nu*lags))
    corvals_ba[corcells] <- colSums(cormat_ba*exp(-nu*lags))
    corrast_ab <- corrast_ba <- ud
    values(corrast_ab) <- corvals_ab
    values(corrast_ba) <- corvals_ba
    
    covrast1 <- sd*corrast_ab
    covrast2 <- sd*corrast_ba
    
    # quantify the contribution: if it's positive then it's the ratio, if
    # it's negative and of greater magnitude than the ud prod it's the ud prod negative
    contrib1 <- (covrast1>=0)*covrast1-(covrast1<0&abs(covrast1)>ud)*ud+(covrast1<0&abs(covrast1)<=ud)*covrast1
    contrib2 <- (covrast2>=0)*covrast2-(covrast2<0&abs(covrast2)>ud)*ud+(covrast2<0&abs(covrast2)<=ud)*covrast2
    subs(contrib1,data.frame(NA,0))
    foi1 = max(ud+covrast1,0)
    foi2 = max(ud+covrast2,0)
    
    covcontrib_df <- rbind(covcontrib_df,c(
      ind1=id1,
      ind2=id2,
      iteration = it,
      totfoi1 = cellStats(foi1,sum),
      totfoi2 = cellStats(foi2,sum),
      CovContr_mean1 = cellStats(contrib1/ud, 'mean'),
      CovContr_mean2 = cellStats(contrib2/ud, 'mean'),
      CovContr_sd1 = cellStats(contrib1/ud, 'sd'),
      CovContr_sd2 = cellStats(contrib2/ud, 'sd')))
  }
}
names(covcontrib_df) <- c("ind1", "ind2", "iteration", "foi1","foi2","covcontrmean1", "covcontrmean2", "covcontrsd1", "covcontrsd2")
write.csv(covcontrib_df, "resampled_covariance_contribution.csv", quote = F, row.names = F)

#### Visualization ####
covcontrib_df <- read.csv("resampled_covariance_contribution.csv")
covcontrib_df %>% ggplot(aes(100*as.numeric(covcontrmean1)))+facet_grid(ind1~ind2, scales = "free")+geom_histogram()+geom_vline(aes(xintercept= 100*as.numeric(covcontrmean1)), filter(covcontrib_df, iteration==1), col = "red")+labs(x="Covariance relative contribution (%)")+geom_vline(xintercept = 0, linetype = 2)

#### Networks ####

# It seems that for some pairs the correlation influences the FOI, although the
# overall contribution may still be small. We can compare the FOI estimated
# using networks to represent the link between different individuals. We can
# assign three different values to the edge weights: the FOI ignoring the
# covariance term, i.e. assuming independent movement. This is analogous to
# using the CDE to infer contact.The second option is to use the FOI estimated
# with covariance and the original data. Finally, we can recalculate the
# contribution of covariance randomizing the time signature. The differences
# would potentially indicate the relative importance of including or not the
# covariance term for different pairs.

# First let's get the matrices
combs <- combn(1:5, 2)
ids <- c("151564","151571","151589","151592","151599")
edgeweights <- sapply(1:ncol(combs), \(i) {
  ind1 <- ids[combs[1,i]]
  ind2 <- ids[combs[2,i]]
  udprod <- raster(udprodfiles[grepl(ind1, udprodfiles) & grepl(ind2, udprodfiles)])
  sdprod <- raster(sdprodfiles[grepl(ind1, sdprodfiles) & grepl(ind2, sdprodfiles)])
  cellarea <- prod(res(udprod))
  corfile <- paste0("../outputs/Corrast_30min_", ind1,"-",ind2, ".tif")
  foiud <- beta*lam/cellarea*udprod/nu
  if(file.exists(corfile)) {
    foicovorign <- beta*lam/cellarea*(udprod/nu+sdprod*raster(corfile))
    covrast <- udprod/nu*mean(covcontrib_df[covcontrib_df$ind1==ind1 & covcontrib_df$ind2==ind2, "covcontrmean1"])
    foicovobs <- beta*lam/cellarea*(udprod/nu+covrast)
  } else {foicovorign <- foicovobs <- foiud}
  c(ind1, ind2, cellStats(foiud, sum),
    cellStats(foicovorign, sum),
    cellStats(foicovobs, sum))
  }
  )

edgeweights_df <- as.data.frame(t(edgeweights))
names(edgeweights_df) <- c("ind1", "ind2", "foi_UD", "foi_cov","foi_cov_rand")
grf <- graph_from_data_frame(edgeweights_df)
plot(grf)
#### Step randomization (Spiegel 2016) #### 

#Another way to determine the contribution of the correlation in movement is to
#randomize the steps within each individual, and recalculate the correlation. We
#would expect if correlation is playing a sizeable role, randomizing the steps
#of tracks that are originally correlated tracks would reduce the overall FOI I
#can do this randomization using the spatsoc package, and then redo the process
#of calculating correlations and FOI. The steps are the same, so the UD and SD
#should stay constant 
library(spatsoc)
library(data.table)

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
dat_dt <- data.table(x = dat1$x_,y = dat1$y_, 
                 datetime = as.POSIXct(dat1$t_, format = "%Y-%m-%dT%H:%M:%SZ", tz = "Etc/GMT-6"), 
                 id = dat1$animal_id)
randomizedtracks <- randomizations(dat_dt, 
                                   type = "trajectory", 
                                   id = "id", 
                                   coords = c("x", "y"), 
                                   datetime = "datetime", iterations = 20)
for (i in unique(randomizedtracks$iteration)) {
  
  
}