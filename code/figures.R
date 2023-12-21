#### Steps diagram base ####
x11(width = 9, height = 6, pointsize = 10, bg = "white")
par(fig = c(0/18,3.5/18,6.25/12, 9.75/12), mar = c(0,0,0,0))
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
sapply(1:2, \(i) lines(A$x[601:800,i],A$y[601:800,i], col = hcl.colors(2, palette = "Dark 3")[i]))

par(fig = c(4/18, 7.5/18,6.25/12,9.75/12), new=T)
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
plot(uds[[2]], col.grid=NA, DF = "PDF", level=0,
     col.level = hcl.colors(2, palette = "Dark 3"),
     col.DF = hcl.colors(2, palette = "Dark 3"), 
     xaxt="n", yaxt="n",labels=c("","0.95",""),add=T)

par(fig = c(8/18, 11/18,8.5/12,11.5/12), new=T)
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
with(list(nu = 1/(24*7)), raster::plot(prods[[1]]$UD/nu, 
                                       col = hcl.colors(25, "Plasma"), bty = "n",ann=F,xaxt = 'n', yaxt = 'n', add=T, legend.width=1.8))


par(fig = c(8/18, 11/18,4.5/12,7.5/12), new=T)
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
with(list(), raster::plot(prods[[1]]$SD, 
                          col = hcl.colors(25, "Plasma"), bty = "n",ann=F,xaxt = 'n', yaxt = 'n', add=T, legend.width=1.8))


corrast <- prods$`1-2`$UD
corrcoords <- xyFromCell(corrast,as.numeric(colnames(cors$`1-2`$CAB)))
corvals <- DF$cor2
values(corrast) <- 0
corrast[as.numeric(colnames(cors[[1]]$CAB))] <- corvals

par(fig = c(8/18, 11/18,0.5/12,3.5/12), new=T)
pcols <- rgb(colorRamp(hcl.colors(25, "Plasma"))(scale(corvals, center=F, scale = max(corvals))), maxColorValue = 255)
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
with(list(), raster::plot(corrast, col = hcl.colors(25, "Plasma") ,  
                          yaxt="n",xaxt="n",add=T, legend.width=1.8))
# sapply(1:2, \(i) lines(A$x[601:800,i],A$y[601:800,i], col = hcl.colors(2, palette = "Pastel 1")[i]))
# points with correlation
par(fig = c(8/18, 11/18,0.5/12,3.5/12), new=T)
points(corrcoords[order(corvals),], pch = 15, col = pcols[order(corvals)], cex=0.8)

DF$FOI <- DF$udp/(1/24/7)+DF$sdp*DF$cor2

foirast <- corrast
values(foirast) <- 0
foivals <- DF$FOI
foirast[as.numeric(colnames(cors[[1]]$CAB))] <- foivals
par(fig = c(12.5/18,16/18,4.25/12,7.75/12), new=T)
plot(A$x,A$y,type = 'n',asp=1,xlab = "", ylab = "",xaxt="n", yaxt="n")
raster::plot(foirast, col = hcl.colors(15, "Rocket", rev=T), 
             add=T, legend.width = 1.8)
pcols <- rgb(colorRamp(hcl.colors(25, "Rocket", rev=T))(scale(foivals, center=F, scale = max(foivals))), maxColorValue = 255)
par(fig = c(12.5/18,16/18,4.25/12,7.75/12), new=T)
# plot.new(); plot.window(c(-100,100),c(-100,100), asp = 1)
points(corrcoords[order(foivals),], pch = 15, col = pcols[order(foivals)], cex=0.9)

#### Analytical figure ####
## Code for figure for analytical results
options(scipen = 999)
library(ggplot2)
library(data.table)

prob = seq(0.01, 0.99, len=100)
relative_contrib = ((1 - prob) / (prob))
corr_vals = c(0.1, 0.25, 0.5, 0.75, 1)
res = data.table(sapply(corr_vals, function(i) relative_contrib*i))
colnames(res) = as.character(corr_vals)
res$prob = prob
res_melt = melt(res, variable.name="corr", value.name="foi", id.var="prob")

ggplot(data=res_melt) + geom_line(aes(x=prob, y=foi, color=corr)) +   
  geom_hline(aes(yintercept=1), linetype="dashed") +
  theme_classic(base_size = 12) +
  labs(colour="Correlation") +
  annotate("text", x=0.5, y=10, label = "Correlation dominates FOI") +
  annotate("text", x=0.5, y=0.01, label = "Spatial overlap dominates FOI") +
  scale_color_manual(values=c('#fee5d9','#fcae91','#fb6a4a','#de2d26','#a50f15')) +
  ylab("Relative contribution of correlation") + 
  xlab(expression(paste("Area of contact / Total area ", (A[x] / A[tot])))) +
  theme(legend.position = c(0,0.05), legend.justification = c(0,0), legend.background = element_rect(fill = NA))+
  scale_y_log10(breaks=c(0.001, 0.01, 0.1, 1, 10, 100))

ggsave("../docs/figures/correlation_analytical_figure.pdf", width=6.5, height=5)

#### example cross-correlation functions ####
x11()
par(mfrow=c(2,2))
par(omi = c(0.4,0.4,0.2,0))
# Two individuals together, going back and forth between two patches
with(list(a = rep(rep(c(1,0),each=10),20)),plot(ccf(a,a, plot = F, lag.max = 40)[0:30,], type='l', ci.col = "gray70"), ci.type = "ma")
mtext("a)", 2, line = 1.5, at=1.45,las=1, cex=1.1)
with(list(a = rep(rep(c(1,0),each=10),20),
          b = rep(rep(c(0,1),each=10),20)),plot(ccf(a,b, plot = F, lag.max = 40)[0:30,], type='l', ci.col = "gray70",xlim = c(0,30)))
mtext("b)", 2, line = 1.5, at=1.45,las=1, cex=1.1)
with(list(a = rep(rep(c(1,0),each=10),20),
          b = rep(rep(c(1,0),each=5),20)),plot(ccf(a,b, plot = F, lag.max = 40)[0:30,], type='l', ci.col = "gray70",xlim = c(0,30)))
mtext("c)", 2, line = 1.5, at=0.2,las=1, cex=1.1)
with(list(a = rbinom(400,1,0.01)),plot(ccf(a,a, plot = F, lag.max = 40)[0:30,], type='l', ci.col = "gray70",xlim = c(0,30)))
mtext("d)", 2, line = 1.5, at=1.2,las=1, cex=1.1)
mtext("Lag",1, outer = T, line = 1, cex=1.3)
mtext("Correlation",2, outer = T, line = 1, cex=1.3)

