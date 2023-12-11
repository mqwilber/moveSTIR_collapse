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
			 theme_classic() +
			 labs(colour="Correlation") +
			 annotate("text", x=0.5, y=10, label = "Correlation dominates FOI") +
			 annotate("text", x=0.5, y=0.01, label = "Spatial overlap dominates FOI") +
			 scale_color_manual(values=c('#fee5d9','#fcae91','#fb6a4a','#de2d26','#a50f15')) +
			 ylab("FOI due to correlation / FOI due to spatial overlap") + xlab("Area of contact / Total area (A_x / A_tot)") +
			 scale_y_log10(breaks=c(0.001, 0.01, 0.1, 1, 10, 100))

ggsave("../docs/figures/correlation_analytical_figure.pdf", width=6.5, height=5)