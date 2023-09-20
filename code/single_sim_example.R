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
A <- simulate_tracks(tau = 5, dp = 2, social = 0.9)
# visualize

plot.new()
pdf("../docs/figures/example_tracks.pdf",width = 3,height = 3)
par(mar = c(0.5,0.5,0.5,0.5))
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
grid()
sapply(1:2, \(i) lines(A$x[,i],A$y[,i], col = hcl.colors(2, palette = "Dark 3")[i]))
dev.off()

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
# Plot UDs separately
image(raster(uds[[2]][[1]], DF="PMF"), asp=1, col = colorRampPalette(c("white",hcl.colors(2,"Dark 3")[1]))(18), bty = "o", xlim = range(A[[1]]),ylim = range(A[[2]]))
image(raster(uds[[2]][[2]], DF="PMF"), asp=1, col = colorRampPalette(c("white",hcl.colors(2,"Dark 3")[2]))(18), bty = "o", xlim = range(A[[1]]),ylim = range(A[[2]]))

# Get the product of the UDs, p_i(x) and p_j(x), as well as the product of their
# standard deviations sqrt(p_i(x)(1-p_i(x)))*sqrt(p_j(x)(1-p_j(x)))
prods <- getUDprod(uds)

zrng <- range(unlist(lapply(prods,lapply,cellStats,range)))
pdf("../docs/figures/example_UDprod.pdf", width = 4, height = 3,family = "sans")
par(mar = c(0.5,0.5,0.5,4))
with(list(nu = 1/(24*7)), raster::plot(prods$UD[[2]]/nu, col = hcl.colors(255, "Plasma"), bty = "n",xaxt="n", yaxt="n", asp=1, xlab="", ylab = ""))
dev.off()

pdf("../docs/figures/example_SDprod.pdf", width = 4, height = 3,family = "sans")
par(mar = c(0.5,0.5,0.5,4))
raster::plot(prods$SD[[2]], col = hcl.colors(255, "Plasma"), bty = "n",xaxt="n", yaxt="n", asp=1, xlab="", ylab = "")
dev.off()

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
pdf("../docs/figures/example_corrRast.pdf",width = 4, height = 3, family = "sans")
raster::plot(corrast, col = hcl.colors(25, "Blue-Red 3") ,  yaxt="n",xaxt="n", zlim = zlim)
grid()
dev.off()
# Plot FOI rasters
par(mfrow = c(2,1), mar = c(0,1,0,4), pty = "s")
plot(test[[1]], col = hcl.colors(255,"Plasma"), xaxt = "n", yaxt = "n")
lines(extent(test[[1]]))
grid()
plot(log(test[[1]]), col = hcl.colors(255,"Plasma"), xaxt = "n", yaxt = "n")

# A better way to visualize this is to plot the locations of the cells with correlations
# get the coordinates of the cells
corrcoords <- xyFromCell(corrast[[1]],Which(corrast[[1]]!=0, T))



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

