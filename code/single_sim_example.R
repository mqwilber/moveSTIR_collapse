# PMoveSTIR simulation example
# Author: Juan S. Vargas
# Date: September 18, 2023

library(tidyverse)
library(raster)
library(ctmm)
library(move)
library(ggpubr)
library(ggdist)
#### Intro ####

# Simulation and analysis of movement data. We simulate movement of two
# individuals with possible attraction and temporal correlation, and analyze
# their trajectories using the PMoveSTIR framework.

set.seed(23)
# the functions are defined in the functions.R script
source("code/functions.R")


# simulate trajectories of two individuals
A <- simulate_tracks(tau = 5, dp = 2, social = 0.9)
# visualize

plot.new()
# pdf("../docs/figures/example_tracks.pdf",width = 3,height = 3)
par(mar = c(0.5,0.5,0.5,0.5))
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
grid()
sapply(1:2, \(i) lines(A$x[,i],A$y[,i], col = hcl.colors(2, palette = "Dark 3")[i]))
# dev.off()

# get the utilization distributions of both individuals
uds <- getUDs(A)

# Plot UDs together
pdf("../docs/figures/example_UDs_overlay.pdf", width = 3, height = 3)
par(mar = c(0.5,0.5,0.5,0.5))
plot(extent(uds[[2]]), type='n', xaxt="n",yaxt='n', asp=1)
grid()
plot(uds[[2]], col.grid=NA, DF = "PDF", level=0,
     col.level = hcl.colors(2, palette = "Dark 3"),
     col.DF = hcl.colors(2, palette = "Dark 3"), 
     xaxt="n", yaxt="n",labels=c("","0.95",""),add=T)
dev.off()
# # Plot UDs separately
# image(raster(uds[[2]][[1]], DF="PMF"), asp=1, col = colorRampPalette(c("white",hcl.colors(2,"Dark 3")[1]))(18), bty = "o", xlim = range(A[[1]]),ylim = range(A[[2]]))
# image(raster(uds[[2]][[2]], DF="PMF"), asp=1, col = colorRampPalette(c("white",hcl.colors(2,"Dark 3")[2]))(18), bty = "o", xlim = range(A[[1]]),ylim = range(A[[2]]))

# Get the product of the UDs, p_i(x) and p_j(x), as well as the product of their
# standard deviations sqrt(p_i(x)(1-p_i(x)))*sqrt(p_j(x)(1-p_j(x)))
prods <- getUDprod(uds)
# Plot
zrng <- range(unlist(lapply(prods,lapply,cellStats,range)))
# pdf("../docs/figures/example_UDprod.pdf", width = 4, height = 3,family = "sans")
par(mar = c(0.5,0.5,0.5,4))
with(list(nu = 1/(24*7)), raster::plot(prods[[1]]$UD/nu, col = hcl.colors(25, "Plasma"), bty = "n",ann=F))
# dev.off()

pdf("../docs/figures/example_SDprod.pdf", width = 4, height = 3,family = "sans")
par(mar = c(0.5,0.5,0.5,4))
plot(extent(uds[[2]]), type='n', xaxt="n",yaxt='n', asp=1)
raster::plot(prods[[1]]$SD, col = hcl.colors(25, "Plasma"), bty = "n",xaxt="n", yaxt="n", asp=1, xlab="", ylab = "")
dev.off()

# get correlation vector at cells with overlap
cors <- getCorrs(A, prods, prewt = T, ci.method = 'reg')

#' Most correlation values are small and would come up even if the series were independent, 
#' we need to filter out those values and keep only the ones that are actually meaningful
#' (i.e. significant). The criterion of significance at 95% for independent data is 1.96/sqrt(n), 
#' where n is the length of the series. In addition to this filtering, we need to account for
#' potential correlations that appear due to the autocorrelation in the data. 
#' One way to do this is to prewhiten the series before calculating the cross-correlation. 
#' This is, to fit an autoregressive model to one of the series, and filter
#' both of them using the autocorrelation coefficients.
length(raster(uds[[2]][[1]]))
ncol(cors$`1-2`$CBA)

#' There are 39270 cells overall, but only 36 that both individuals visited. For all others the
#' correlation component of FOI is null. 
#' Assuming a step survival function of parasites in the environment, where the parasites last for 10
#' time units, the cumulative correlation component for unfiltered data would be
DF <- data.frame(udp = prods[[1]]$UD[as.numeric(colnames(cors[[1]]$CAB))],
           sdp = prods[[1]]$SD[as.numeric(colnames(cors[[1]]$CAB))], 
           cor1 = colSums(cors[[1]]$CAB[1:10,]),
           cor2 = colSums(cors[[1]]$CBA[1:10,]))
hist((DF$udp+DF$sdp*DF$cor1)/DF$udp, main = "ratio FOI with/without correlation")

#' Most of the cumulative correlations are negative, so the effect on FOI is to decrease it. 
#' However, there are some values where there are high correlations that bring the estimate way
#' up. 
#' Now let's do the same but keeping only the correlation values that are significant
corsf <- cors[[1]]$CAB*(abs(cors[[1]]$CAB)>(1.96/sqrt(500)))
DF <- data.frame(udp = prods[[1]]$UD[as.numeric(colnames(cors[[1]]$CAB))],
                 sdp = prods[[1]]$SD[as.numeric(colnames(cors[[1]]$CAB))], 
                 cor1 = colSums(corsf[1:10,]))
hist((DF$udp+DF$sdp*DF$cor1)/DF$udp, main = "ratio FOI with/without correlation")

#' The effect here is that the small, negative correlations go away, but large 
#' correlations stay. The local effect on FOI is thus largely the same.


# Plot correlation raster
corrast <- prods[[1]]$UD
values(corrast) <- 0
# corvals <- with(list(nu = 1/24/7), colSums(cors[[1]]$CAB*exp(-nu*(0:(nrow(cors[[1]])-1)))))
corvals <- DF$cor1
corrast[as.numeric(colnames(cors[[1]]$CAB))] <- corvals
zlim = c(-max(abs(cellStats(corrast,range))),max(abs(cellStats(corrast,range))))
# pdf("../docs/figures/example_corrRast.pdf",width = 4, height = 3, family = "sans")
raster::plot(corrast, col = hcl.colors(25, "Blue-Red 3") ,  yaxt="n",xaxt="n")
# dev.off()

# There are three random cells where there is some correlation


# A better way to visualize this is to plot the locations of the cells with correlations
# get the coordinates of the cells
corrcoords <- xyFromCell(corrast[[1]],Which(corrast[[1]]!=0, T))
x11(width = 4.5,height = 3.5)
par(fin = c(4,3),pin=c(3,3),ann=F, mai=c(0,0,0,1))
corvals <- corrast[!is.na(corrast)]
pcols <- rgb(colorRamp(hcl.colors(25, "Blue-Red 3"))(scale(corvals,center=-max(abs(corvals)),scale = 2*max(abs(corvals)))), maxColorValue = 255)
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
raster::plot(corrast, col = hcl.colors(25, "Blue-Red 3") ,  yaxt="n",xaxt="n", zlim = zlim,add=T)
grid()
sapply(1:2, \(i) lines(A$x[,i],A$y[,i], col = hcl.colors(2, palette = "Pastel 1")[i]))
# points with correlation
points(corrcoords, pch = 21, bg = pcols)

# Do the same without prewhitening
cors2 <- getCorrs(A, prods, prewt = F, ci.method = 'reg')
corsf <- cors2[[1]]$CAB*(abs(cors2[[1]]$CAB)>(1.96/sqrt(500)))
corrast <- prods[[1]]$UD
values(corrast) <- 0
corvals <- DF$cor1
corrast[as.numeric(colnames(cors[[1]]$CAB))] <- colSums(corsf[1:10,])
corrcoords <- xyFromCell(corrast[[1]],Which(corrast[[1]]!=0, T))
pcols <- rgb(colorRamp(hcl.colors(25, "Blue-Red 3"))(scale(corvals,center=-max(abs(corvals)),scale = 2*max(abs(corvals)))), maxColorValue = 255)
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
raster::plot(corrast, col = hcl.colors(25, "Blue-Red 3") ,  yaxt="n",xaxt="n", zlim = zlim,add=T)
grid()
sapply(1:2, \(i) lines(A$x[,i],A$y[,i], col = hcl.colors(2, palette = "Pastel 1")[i]))
points(corrcoords, pch = 21, bg = pcols)

# The result is practically the same without prewhitening

# calculate cell FOI
testFOI <- getFOI(A,uds)

# Plot FOI rasters
x11(width = 4.5,height = 3.5)
par(fin = c(4,3),pin=c(3,3),ann=F, mai=c(0,0,0,1))
raster::plot(extent(testFOI[[1]][[1]]),type='n',xaxt='n',yaxt='n')
raster::plot(testFOI[[1]][[1]], col = hcl.colors(25,"Reds 3", rev = T), add=T)
raster::plot(1e-20+log10(testFOI[[1]]), col = hcl.colors(25,"Reds 3", rev = T), add=T)
# plot(log(test[[1]]), col = hcl.colors(255,"Plasma"), xaxt = "n", yaxt = "n")
# add points with higher values
corrcoords <- xyFromCell(testFOI[[1]][[1]],Which(testFOI[[1]][[1]]>1e-6,T))
corvals <- testFOI[[1]][[1]][testFOI[[1]][[1]]>1e-6]/cellStats(testFOI[[1]][[1]],max)
pcols <- rgb(colorRamp(hcl.colors(25, "Reds 3",rev=T))(scale(corvals,center=-max(abs(corvals)),scale = 2*max(abs(corvals)))), maxColorValue = 255)
points(corrcoords, pch = 21, bg = pcols)

# zoomed in version
par(fin = c(4,3),pin=c(3,3),ann=F, mai=c(0,0,0,1))
raster::plot(extent(corrcoords),type='n',xaxt='n',yaxt='n')
raster::plot(test[[1]], col = hcl.colors(25,"Plasma"), add=T)



# get total FOI
sapply(test, cellStats,sum)


# Compare with regular MoveSTIR results


#### Testing multiple scenarios ####
# I can run simulations for multiple scenarios, changing the decay rate nu, the
# grid resolution d, as well as the interaction strength. I run this through the
# run_sims_cluster code in the UTK HPC, and import the results here
outdf <- read.csv("outputs/231124_sims_out.txt",header = F)[,-10]
social <- c(0,0.5, 0.7, 0.9,0.93, 0.96, 1)
nus <- 1/(24*c(1/12, 1/3, 1, 3, 7))
gridres <- c(5, 10, 20)
names(outdf) <- c("sim", "nu", "Ax","Atoti","Atotj", "foi_ud", "foi_full1", "foi_full2", "overlap")
outdf$social <- rep(social, each = 15)                           
head(outdf)
### Visualization ####

# FOI relative to minimum value, with and without covariance term. This plot
# shows that greater interaction strength leads to greater FOI in general. This
# effect is much more pronounced when you consider the covariance term. The
# difference is even apparent for no interaction. Greater interaction strengths
# lead to a greater effect of covariance, with respect to UD-only grid

p1 <- outdf %>% group_by(sim) %>% 
  # slice_max(order_by = Ax,n=5) %>% 
  # filter(near(nu, 1/(24*1))) %>% 
  slice_head(n=1) %>%
  ggplot()+
  # stat_halfeye(aes(social, foi_full1/mean(foi_ud)), .width=0, point_colour=NA, fill = "darkred")+
  # stat_halfeye(aes(social, foi_ud/mean(foi_ud)), .width=0, point_colour=NA, fill = "steelblue")+
  geom_boxplot(aes(social, foi_full1/mean(foi_ud), group = social), width = 0.015, color = "darkred", position = position_nudge(0.01))+
  geom_boxplot(aes(social, foi_ud/mean(foi_ud), group = social), width=0.015, color = "steelblue", position = position_nudge(-0.01))+
  labs(x = "Interaction strength", y = "Relative FOI")+
  theme_minimal(base_size = 12)+
  # theme(panel.grid.major.y = element_line())+
  scale_y_log10()
p1

p2 <- filter(outdf, social%in%c(0,0.5)) %>% slice_head(by=sim, n=1) %>% 
  ggplot(aes(overlap, (foi_full1/foi_ud)))+
  geom_hline(yintercept = 1,linetype=2)+
  geom_point()+
  theme_minimal(base_size = 12)+
  labs(x = "Home Range Overlap", y = "FOI ratio")
  # scale_y_log10()
p2
# Effects of epi parameters: nu and contact distance

# To show how the FOI changes with respect to decay rate, I plot the change in
# FOI for the same simulation wrt decay time. The best fit lines converge at 0,
# since I'm plotting the relative change with respect to the value for the
# highest rate. This shows, in a way, the effect of correlation. If FOI depended
# solely on UD and overlap, the increase would be linear, and there would be no
# difference between different levels of interaction
outdf %>% group_by(sim,nu) %>% slice_head(n=1) %>% 
  filter(social %in% c(0,0.7,0.96,1)) %>% 
  group_by(sim) %>%   mutate(dfoi = (foi_full1-foi_full1[nu==min(nu)])/(foi_full1[nu==min(nu)]),
                             dfoiud = (foi_ud-foi_ud[nu==min(nu)])/(foi_ud[nu==min(nu)])) %>% 
  ggplot(aes(1/nu/24, dfoi,color=factor(social)))+
  geom_smooth(method='lm', se = F)+
  geom_point(aes(1/nu/24,dfoiud), color='red')+
  geom_point()+
  scale_color_discrete(type=hcl.colors(4, "BluGrn", rev=T))+
  theme_classic(base_size = 12)+
  labs(x = "Decay time (days)", y = expression(paste(Delta," FOI")), color = "Interaction")+
  theme(legend.position = c(0.9,0.2))
 

# Another way to show this is to do something similar to b with the HR overlap,
# where I plot the ratio of FOI with covariance over FOI without, wrt decay
# time.
outdf %>% group_by(sim, nu) %>% slice_max(order_by = Ax, n=1) %>% 
  ggplot()+
  geom_hline(yintercept = 1,linetype=2)+
  geom_boxplot(aes(1/nu/24,foi_full1/foi_ud, color = factor(social),group = paste(nu,social)))+
  # geom_point(aes(1/nu/24, foi_full1/foi_ud,color=factor(social)))+
  scale_y_log10()+
  theme_classic(base_size = 12)+
  labs(x= "Decay time (days)", y = "FOI ratio")

# Yet another way, boxplot with diff color for with and without (as in a.). This
# show the increase in FOI, and how the relative difference decreases with
# longer decay times. 
outdf %>% group_by(sim) %>% slice_max(order_by = Ax,n=5) %>% 
  filter(social==0.9) %>%
  ggplot()+
  geom_boxplot(aes(1/(nu*24), foi_full1/min(foi_ud), group = factor(nu)), color = "darkred", position = position_nudge(0.05))+
  geom_boxplot(aes(1/(nu*24), foi_ud/min(foi_ud), group = factor(nu)), color = "steelblue", position=position_nudge(-0.05))+
  # stat_smooth(method = "lm", linetype= 2, color="black")+
  # geom_boxplot(aes(group = factor(nu)))+
  # facet_grid(~social)+
  # geom_point()+
  theme_classic(base_size = 14)+
  labs(x = "Mean decay time (days)", y = "Relative FOI")+
  scale_y_log10()+
  theme(panel.grid.major.y = element_line())
  
# Another way change in FOI wrt change in nu for different interaction values
# outdf %>% group_by(sim) %>% slice_max(order_by = Ax,n=5) %>% arrange(sim,nu) %>% 
#   mutate(dfoicov = (lead(foi_full1)-foi_full1),
#          dfoiud = (lead(foi_ud)-foi_ud),
#          dnu = (lead(1/nu)-1/nu)) %>% 
#   ggplot()+
#   geom_hline(yintercept = 0, linetype=2)+
#   geom_smooth(aes(social,dfoiud/dnu),se=F,color="steelblue")+
#   geom_point(aes(social,dfoicov/dnu), color="darkred")+
#   geom_point(aes(social,dfoiud/dnu),color='steelblue')+
#   theme_classic(base_size = 12)+
#   labs(x = "Interaction strength", y = expression(paste("Relative effect of",nu,"\n",Delta,"FOI"/Delta, nu)))

# Finally, a basic one in two parts. One showing that FOI increases, and another
# that the contribution decreases for longer decay times. Version showing
# different interaction values
p3 <- outdf %>%  group_by(sim) %>% 
  filter(Ax == max(Ax)) %>% 
  # slice_max(order_by = Ax,n=5) %>% 
  filter(social %in% c(0,0.7,0.96,1)) %>% 
  ggplot(aes(1/nu/24,foi_full1/mean(foi_full1),color=factor(social)))+
  geom_smooth(method = "lm", se=F)+
  geom_point()+
  theme_minimal(base_size = 12)+
  labs(x = "Decay time (days)", y = "Relative FOI", color = "Interaction")+
  scale_color_discrete(type=hcl.colors(4, "BluGrn", rev=T))
p3
# covariance contribution, 4interaction values
p4 <- outdf %>%  group_by(sim) %>% 
  filter(Ax == max(Ax)) %>% 
  filter(social %in% c(0,0.7,0.96,1)) %>% 
  ggplot(aes(1/nu/24,(foi_full1)/foi_ud,color=factor(social)))+
  # geom_hline(yintercept = 1, linetype=2)+
  # geom_smooth(se=F)+
  geom_point()+
  scale_y_log10()+
  theme_minimal(base_size = 12)+
  labs(x = "Decay time (days)", y = "FOI ratio", color = "Interaction")+  
  scale_color_discrete(type=hcl.colors(4, "BluGrn", rev=T))
p4

  
# Longer decay times, or higher decay rates, lead to higher overall FOI. The FOI
# is proportional to the mean decay time. There is also a strong effect of
# increasing FOI as the interaction strength increases. This increase would not
# be captured if we considered only direct encounter, ignoring temporal
# correlation in movement
# There is a linear increase in FOI proportional to the mean decay time. The
# relative contribution of the covariance, however, decreases. At longer decay
# times the FOI estimated assuming independnet movement is similar to the FOI
# that considers correlation

# Effect of changing grid size, matching trajectories

p5 <- outdf %>% filter(near(nu, 1/2), social %in% c(0,0.7,0.96,1)) %>% 
  group_by(sim) %>% 
  ggplot(aes(Ax/Atoti/1e4,foi_full1/mean(foi_full1), color = factor(social)))+
  stat_smooth(method = "lm", se = F)+
  geom_point()+
  theme_minimal(base_size = 12)+
  labs(x = expression(paste("Relative cell size (",A[x]/A[tot],")")), y = "Relative FOI", color = "Interaction")+
  scale_color_discrete(type=hcl.colors(4, "BluGrn", rev=T))+
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=11))
p5
p6 <- outdf %>% filter(near(nu, 1/2), social %in% c(0,0.7,0.96,1)) %>% 
  group_by(sim) %>% 
  ggplot(aes(Ax/Atoti/1e4,foi_full1/foi_ud, color = factor(social)))+
  # geom_hline(yintercept = 1, linetype=2)+
  # stat_smooth(method = "lm", se = F)+
  geom_point()+
  theme_minimal(base_size = 12)+
  labs(x = expression(paste("Relative cell size (",A[x]/A[tot],")")), y = "FOI ratio", color = "Interaction")+
  scale_color_discrete(type=hcl.colors(4, "BluGrn", rev=T))+
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=11))+
  scale_y_log10()
p6
# The effect of cell size varies. For perfectly overlapping individuals larger
# cells result in lower FOI. The opposite is true for slightly lower interaction
# strength. Furthermore, low interaction strengths show no real effect of cell
# size
pdf("docs/figures/sim_results.pdf", width = 9,height = 6)
ggarrange(p1,p3,p5,p2,p4,p6, legend.grob = get_legend(p3), legend = "right", align = "hv", labels="auto")
dev.off()



# FOIvs overlap
filter(outdf, near(nu,1/24)) %>% slice_min()%>% pivot_longer(cols = starts_with("foi"), names_to = "calc", values_to = "foi") %>% 
  filter(social<1, calc == "foi_full1") %>% 
  ggplot(aes(overlap/Ax, foi, group = nu,color=factor(1/(24*nu))))+
  geom_smooth(method = "lm", se = F)+
  geom_point()+
  labs(x = "Home range overlap", y = "Force of infection", color = "Mean decay time")+
  theme_classic(base_size = 14)+
  scale_color_discrete(labels = c("2 h", "4 h", "8 h", "12 h", "1 day", "3 days", "7 days"))+
  theme(legend.position = c(0.1,0.95), legend.justification = c(0.1,0.95), legend.background = element_blank())+
  scale_y_log10()

outdf %>% filter(near(nu,1/24)) %>% slice_max(order_by = Ax, n = 1, by = sim) %>% 
  ggplot()+
  geom_point(aes(1/Atot/min(1/Atot),foi_ud/min(foi_ud)), color = "steelblue", size=2)+
  geom_point(aes(1/Atot/min(1/Atot), foi_full1/min(foi_ud)), color="darkred", size=2)+
  theme_classic(base_size = 14)+
  labs(x = "Relative FOI - HR overlap", y = "Relative FOI - PMoveSTIR")


outdf %>% filter(near(nu,1/24)) %>% 
  mutate(foi_HR = 1*1/24*1/nu*overlap^2/(Atot/7)) %>% 
  ggplot()+
  geom_abline(slope = 1, linetype=2)+
  geom_point(aes(foi_HR,foi_ud), color = "steelblue", size=2)+
  geom_point(aes(foi_HR,foi_full1), color = "darkred", size=2)+
  theme_classic(14)+
  labs(x = "FOI - HR overlap", y = "FOI - PMoveSTIR")
  
# FOI estimated using only UD is proportional to the one using the home range
# overlap, but when you include the covariance term you can get significantly
# higher FOI values, particularly for high interaction strengths
  


# Linear model of FOI vs overlap, for moderate interaction strength and just the
# full calculation. 
outdf %>% pivot_longer(cols = starts_with("foi"), names_to = "calc", values_to = "foi") %>% 
  filter(social<1, calc == "foi_full1") %>% 
  lm(foi~overlap+nu, data = .) %>% summary()
# linear model shows relationship between overlap and foi. nu has a clear
# negative effect; faster decay rates result in lower FOI

