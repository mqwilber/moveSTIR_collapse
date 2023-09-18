# PMoveSTIR simulation example
# Author: Juan S. Vargas
# Date: September 18, 2023

# Simulation and analysis of movement data. We simulate movement of two
# individuals with possible attraction and temporal correlation, and analyze
# their trajectories using the PMoveSTIR framework.

set.seed(23)
# the functions are defined in the functions.R script
source("functions.R")


# simulate trajectories of two individuals
A <- simulate_tracks(tau = 5, dp = 2, social = 0)
# visualize
plot.new()
par(mar = c(0.5,0.5,0.5,0.5))
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
grid()
sapply(1:2, \(i) lines(A$x[,i],A$y[,i], col = hcl.colors(2, palette = "Dark 3")[i]))

# get the utilization distributions of both individuals
uds <- getUDs(A)

# Plot UDs together
plot(uds[[2]], col.grid=NA, DF = "PDF", level=0,
     col.level = hcl.colors(2, palette = "Dark 3"),
     col.DF = hcl.colors(2, palette = "Dark 3"), 
     xaxt="n", yaxt="n",labels=c("","0.95",""))

# Plot UDs separately
image(raster(uds[[2]][[1]], DF="PMF"), asp=1, col = colorRampPalette(c("white",hcl.colors(2,"Dark 3")[1]))(18), bty = "o", xlim = range(A[[1]]),ylim = range(A[[2]]))
image(raster(uds[[2]][[2]], DF="PMF"), asp=1, col = colorRampPalette(c("white",hcl.colors(2,"Dark 3")[2]))(18), bty = "o", xlim = range(A[[1]]),ylim = range(A[[2]]))

# Get the product of the UDs, p_i(x) and p_j(x), as well as the product of their
# standard deviations sqrt(p_i(x)(1-p_i(x)))*sqrt(p_j(x)(1-p_j(x)))
prods <- getUDprod(uds)

zrng <- range(unlist(lapply(prods,lapply,cellStats,range)))
with(list(nu = 1/(24*7)), plot(prods$UD[[2]]/nu, col = hcl.colors(255, "Plasma"), bty = "n",xaxt="n", yaxt="n", asp=1, xlab="", ylab = ""))
plot(prods$SD[[2]], col = hcl.colors(255, "Plasma"), bty = "n",xaxt="n", yaxt="n", asp=1, xlab="", ylab = "")

# get correlation vector at cells with overlap
cors <- getCorrs(A, uds[[2]])
# calculate cell FOI
test <- getFOI(A,prods)
# Plot cor raster
corrast <- prods$UD[[1]]
values(corrast) <- NA
corvals <- with(list(nu = 1/24/7), colSums(cors[[1]]*exp(-nu*(0:(nrow(cors[[1]])-1)))))
corrast[as.numeric(names(corvals))] <- corvals
zlim = c(-max(abs(cellStats(corrast,range))),max(abs(cellStats(corrast,range))))
plot(corrast, col = hcl.colors(255, "Blue-Red 3") ,  yaxt="n",xaxt="n", zlim = zlim)
grid()
# Plot FOI rasters
par(mfrow = c(2,1), mar = c(0,1,0,4), pty = "s")
plot(test[[1]], col = hcl.colors(255,"Plasma"), xaxt = "n", yaxt = "n")
lines(extent(test[[1]]))
grid()
plot(log(test[[1]]), col = hcl.colors(255,"Plasma"), xaxt = "n", yaxt = "n")


# get total FOI
sapply(test, cellStats,sum)

# Analyze how the FOI would be different under different decay rates, and how
# the contribution of covariance changes

# Compare with regular MoveSTIR results



# I can run simulations for multiple scenarios, changing the decay rate nu, the
# grid resolution d, as well as the interaction strength. I run this through the
# run_sims_cluster code in the uni HPC, and import the results here
outdf <- rbind(read.csv("C:/Users/juans/Downloads/sim_res_230914_1924.csv"),
               read.csv("C:/Users/juans/Downloads/sim_res_230914_1931.csv"))
#Correct the calculation of the UD only FOI. In the simulation I did not scale
#by epi pars beta, lambda, or divided by area
outdf$foi_ud2 <- outdf$foi_ud*1*1/24/outdf$res^2
head(outdf)                            

### FIGURES ####
outdf %>% select(-foi_ud) %>% 
  pivot_longer(cols = starts_with("foi"), names_to = "calc", values_to = "foi") %>% 
  ggplot()+geom_point(aes(1/(nu*24),foi,color=calc))+
  facet_wrap(facets = vars(social), scales = "free_y")+
  labs(x = expression(paste("Decay time ",1/nu," (days)")),
       y = "Force of infection",
       color = "Method")+
  theme_classic(base_size = 16)


# FOI vs interaction, color nu
outdf %>% select(-foi_ud) %>% 
  ggplot(aes(factor(social), foi_full1, color=factor(1/nu)))+
  geom_boxplot(outlier.shape = NA)+theme_classic(base_size = 14)+
  geom_point(position = position_jitterdodge())+
  labs(x = "Interaction strength", y = "Force of infection", color = "Mean decay time")+
  theme(legend.position = c(0.1,0.95), legend.justification = c(0.1,0.95), legend.background = element_blank())+
  scale_color_discrete(labels = c("2 h", "8 h", "1 day", "3 days", "7 days"))

# foi vs decay time
p1 <- outdf %>% 
  filter(social == 0.9) %>% 
  ggplot(aes(1/(nu*24), foi_full1))+
  # stat_smooth(method = "lm", color="black", linetype= 2)+
  geom_boxplot(aes(group = factor(1/(nu*24))))+
  # geom_point(position = position_jitter())+
  theme_classic(base_size = 14)+
  labs(x = "Mean decay time (days)", y = "Force of infection")

# covaraince contribution
p2 <- outdf %>% 
  filter(social == 0.9) %>% 
  ggplot(aes(1/(24*nu), (foi_full1-foi_ud2)/foi_full1))+
  geom_boxplot(aes(group = 1/(24*nu)))+
  theme_classic(base_size = 14)+
  labs(x = "Mean decay time (days)", y = "Relative contribution")


# FOI vs interax, with and without covar
p3 <- outdf %>% select(-foi_ud) %>% 
  pivot_longer(cols = starts_with("foi"), names_to = "calc", values_to = "foi") %>% 
  filter(near(nu, 1/24),calc!="foi_full2") %>% 
  ggplot(aes(social, foi, color = calc, group = paste(social,calc)))+
  geom_boxplot()+
  labs(x = "Interaction strength", y = "Force of infection", color = "Covariance")+
  theme_classic(base_size = 14)+
  scale_y_continuous(limits = c(0,NA))+
  scale_color_discrete(labels = c("With","Without"))+
  theme(legend.position = c(0.1, 0.95), legend.justification = c(0.1,0.95), legend.background = element_blank())+scale_x_log10()

# Ratio of with vs without covar
p4 <- outdf %>% mutate(ratio = foi_full1/foi_ud2) %>% 
  filter(near(nu, 1/(24*1))) %>% 
  ggplot(aes(social, ratio, group = social))+
  geom_boxplot()+
  labs(x = "Interaction strength", y = "Relative difference in FOI")+
  theme_classic(base_size = 14)+
  scale_y_log10()
# FOIvs overlap
outdf %>% pivot_longer(cols = starts_with("foi"), names_to = "calc", values_to = "foi") %>% 
  filter(social<1, calc == "foi_full1") %>% 
  ggplot(aes(overlap, foi, group = nu,color=factor(1/(24*nu))))+
  geom_smooth(method = "lm", se = F)+
  geom_point()+
  labs(x = "Home range overlap", y = "Force of infection", color = "Mean decay time")+
  theme_classic(base_size = 14)+
  scale_color_discrete(labels = c("2 h", "4 h", "8 h", "12 h", "1 day", "3 days", "7 days"))+
  theme(legend.position = c(0.1,0.95), legend.justification = c(0.1,0.95), legend.background = element_blank())+
  scale_y_log10()

pdf("../docs/figures/sim_results.pdf", width = 8, height = 6)
cowplot::plot_grid(p1,p2,p3,p4, labels = "auto", align = "hv")
dev.off()
# Linear model of FOI vs overlap, for moderate interaction strength and just the
# full calculation. 
outdf %>% pivot_longer(cols = starts_with("foi"), names_to = "calc", values_to = "foi") %>% 
  filter(social<1, calc == "foi_full1") %>% 
  lm(foi~overlap+nu, data = .) %>% summary()
# linear model shows relationship between overlap and foi. nu has a clear
# negative effect; faster decay rates result in lower FOI


##### --------------- FUNCTION DEFINITIONS --------------#####
# simulate_tracks <- function(tau = 10, g = 1e3, steps = 500, dp = 0.2, social = NULL, smooth = TRUE) {
#   nind = 2
#   p.lambdas = matrix(rpois(nind*2, dp*g)-dp*g, ncol = 2)# home range centers
#   p.x0 = p.lambdas # initial positions
#   p.phi = 1 # trajectory smoothing parameter
#   nsteps = steps
#   
#   xs = replicate(nind, numeric(nsteps))
#   ys = replicate(nind, numeric(nsteps))
#   
#   for (s in seq_len(nind)) {
#     x <- y <- numeric(nsteps)
#     x[1] <- mux <- p.x0[s,1]
#     y[1] <- muy <- p.x0[s,2]
#     xnoise = rnorm(nsteps)
#     ynoise = rnorm(nsteps)
#     # take successive steps
#     for(i in seq_len(nsteps-1)) {
#       x[i+1] = x[i]-1/tau*(x[i]-p.lambdas[s,1])+sqrt(g)*xnoise[i]
#       y[i+1] = y[i]-1/tau*(y[i]-p.lambdas[s,2])+sqrt(g)*ynoise[i]
#     }
#     xs[,s] = x
#     ys[,s] = y
#   }
#   
#   if(is.numeric(social)) {
#     W <- replicate(nsteps, diag(nrow = nind, ncol = nind))
#     W[1,2,] <- social
#     # make symmetric
#     W <- apply(W, 3, \(x) x*upper.tri(x, diag = T)+t(x*upper.tri(x)), simplify = FALSE)
#     # convolution kernel hij(tau) is the weights matrix, divided by the row sums
#     h_soc <- lapply(W, \(x) x/rowSums(x))
#     # create empty matrices with the same dimensions as the position matrices
#     xsoc <- matrix(nrow = nrow(xs), ncol = ncol(xs))
#     ysoc <- matrix(nrow = nrow(ys), ncol = ncol(ys))
#     # fill the array using dot product of the step*ind x and y matrices, multiplied by the rows of h_soc
#     for (t in 1:nsteps) {
#       xs[t,] <- h_soc[[t]]%*%xs[t,]
#       ys[t,] <- h_soc[[t]]%*%ys[t,]
#     }
#   }
#   
#   
#   if(smooth) {
#     # create smoother matrix for convolution
#     H_inl <- outer(seq_len(nsteps),seq_len(nsteps), \(x,y) abs(x-y)/p.phi*besselK(abs(x-y)/p.phi, nu = 1))
#     diag(H_inl) <- 1
#     H_inl <- H_inl/rowSums(H_inl)
#     xs <- apply(xs, 2, "%*%", H_inl)
#     ys <- apply(ys, 2, "%*%", H_inl)
#   }
#   list(x=xs,y=ys)
# }
# 
# # Function to get the utilization distributions. Output is a list, first element
# # are the ctmm fits, second are the individual UDs
# getUDs <- function(X, dr = NULL) {
#   xs <- X[[1]]
#   ys <- X[[2]]
#   nind <- ncol(xs)
#   nsteps <- nrow(xs)
#   # create move objects
#   moveobjs <- lapply(seq_len(nind), \(c) move(x = xs[,c],y = ys[,c],time = as.POSIXct(3600*seq_len(nsteps),origin = Sys.time()), animal = paste0("ind",c), proj = "+proj=tmerc"))
#   # convert to telemetry
#   telemetries <- lapply(moveobjs, as.telemetry)
#   # ctmm fit
#   GUESS <- lapply(telemetries, \(i) ctmm.guess(i, interactive = F))
#   FITS <- lapply(seq_along(telemetries), \(i) ctmm.select(telemetries[[i]], GUESS[[i]]))
#   # AKDE
#   uds <- if(is.null(dr)) akde(telemetries, FITS) else akde(telemetries, FITS, grid = list(dr = c(dr,dr)))
#   return(list(FITS,uds))
# }
# 
# getUDprod <- function(X) {
#   uds <- X[[2]]
#   nind <- length(uds)
#   ### UD and SD products ###
#   rs <- lapply(uds, raster, DF = "PMF")
#   cellarea <- prod(res(rs[[1]]))
#   
#   # extend all rasters to have the same extent. 
#   rs <- lapply(rs, extend, y = extent(do.call(merge,rs)), value = 0)
#   # Create lists to store objects
#   UDprods <- UDsds <- list()
#   # Possible combinations
#   combs <- combn(length(uds), 2)
#   # Calculate UD and SD products from UD pair values
#   for (i in seq_len(ncol(combs))) {
#     ind1 <- combs[1,i]
#     ind2 <- combs[2,i]
#     r1 <- rs[[ind1]]
#     r2 <- rs[[ind2]]
#     # product of probabilities, divided by area
#     udprob <- r1*r2
#     UDprods[[paste(ind1,ind2,sep = "-")]] <- UDprods[[paste(ind2,ind1,sep = "-")]] <- udprob/cellarea
#     sdprob <- sqrt(r1*(1-r1))*sqrt(r2*(1-r2))
#     UDsds[[paste(ind1,ind2,sep = "-")]] <- UDsds[[paste(ind2,ind1,sep = "-")]] <- sdprob/cellarea
#   }
#   return(list(UD = UDprods, SD = UDsds))
# }
# # function to calculate the correlations. Output is a list, where every element
# # is a lags by cells matrix of correlation between two individuals
# getCorrs <- function(xy, r) {
#   gridcors2 <- list()
#   xs <- xy[[1]]
#   ys <- xy[[2]]
#   nind <- ncol(xs)
#   nsteps <- nrow(xs)
#   combs <- combn(nind,2)
#   for (i in 1:ncol(combs)) {
#     r1 <- raster(r[[i]])
#     ind1 <- combs[1,i]
#     ind2 <- combs[2,i]
#     # get position histories (i.e. which cell was each individual in at every time
#     # point). 
#     pos1 <- cellFromXY(r1, xy = cbind(xs[,ind1],ys[,ind1]))
#     pos2 <- cellFromXY(r1, xy = cbind(xs[,ind2],ys[,ind2]))
#     
#     # keep only cells that both visited at some point
#     ovlpcells <- unique(pos1)[unique(pos1) %in% unique(pos2)]
#     if(length(ovlpcells)==0) {
#       # Still fill in an item in the list, but write just NA. 
#       gridcors2[[paste(ind1,ind2,sep = "-")]] <- NA
#       gridcors2[[paste(ind2,ind1,sep = "-")]] <- NA
#     } else {
#       maxlag <- nsteps-1
#       cormat_ab <- cormat_ba <- matrix(0,nrow = nsteps, ncol = length(ovlpcells))
#       for (j in seq_along(ovlpcells)) {
#         cell <- ovlpcells[j]
#         a <- b <- numeric(nsteps)
#         a[match(cell, pos1)] <- b[match(cell, pos2)]<- 1
#         xcorr <- ccf(a,b,lag.max = maxlag, plot = F)
#         xcorr_vals <- as.numeric(xcorr$acf)
#         cormat_ab[,j] <- rev(xcorr_vals[1:nsteps])
#         cormat_ba[,j] <- xcorr_vals[nsteps:length(xcorr_vals)]
#       }
#       dimnames(cormat_ab) <- dimnames(cormat_ba) <- list(NULL, cell = ovlpcells)
#       gridcors2[[paste(ind1,ind2,sep = "-")]] <- cormat_ab
#       gridcors2[[paste(ind2,ind1,sep = "-")]] <- cormat_ba
#     }
#   }
#   return(gridcors2)
# }
# # function to calculate the per-cell FOI. Calculates correlation within. Output
# # is list with each element a FOI raster
# getFOI <- function(xy, uds, beta = 1, lambda = 1/24, nu = 1/(24*7)) {
#   foirasts2 <- list()
#   gridcors2 <- getCorrs(xy,uds[[1]])
#   for (i in seq_along(gridcors2)) {
#     if (is.array(gridcors2[[i]])) {
#       corcells <- as.numeric(colnames(gridcors2[[i]]))
#       corvals <- numeric(length(uds[[1]][[i]]))
#       # get lags
#       lags <- 0:(nrow(gridcors2[[i]])-1)
#       # scale and integrate correlation at every cell
#       corvals[corcells] <- colSums(gridcors2[[i]]*exp(-nu*lags))
#       udp <- uds[[1]][[i]]
#       sdp <- uds[[2]][[i]]
#       corrast <- udp
#       values(corrast) <- corvals
#       cellarea <- prod(res(corrast))
#       
#       foi <- beta/cellarea*lambda*(1/nu*udp*cellarea+sdp*cellarea*corrast)
#       foi <- foi*(foi>=0)
#     } else {
#       foi <- beta/cellarea*lambda*(1/nu*udp*cellarea)
#     }
#     
#     foirasts2[[i]] <- foi
#   }
#   return(foirasts2)
# }


