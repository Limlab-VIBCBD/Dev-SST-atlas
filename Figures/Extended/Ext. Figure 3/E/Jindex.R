install_github("crazyhottommy/scclusteval")
library(dplyr)
library(ggplot2)
library(readr)
library(scclusteval)

#load data
jc_to_plot_all<-readRDS("Jaccard_allClusters.RDS")
colors_plot<-readRDS("palette_clusters.RDS")

#plot jaccard values for all clusters
pdf('JC_ClusterStability_plot_all_clusters.pdf', width=15, height = 10)
jc_to_plot_all %>% ggplot2::ggplot(ggplot2::aes(x = cluster, y = jaccard, fill = cluster)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .7, width = 3) +
  ggplot2::geom_point(ggplot2::aes(y = jaccard, color = cluster), position = position_jitter(width = .15), size = .5, alpha = 0.7) +
  ggplot2::geom_boxplot(width = .2, outlier.shape = NA, alpha = 0.5) +
  ggplot2::theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,1) +
  ggplot2::theme(legend.position="none") +  scale_colour_manual(values=colors_plot) +  scale_fill_manual(values=colors_plot) + geom_hline(yintercept = c(0.75), linetype = 2,color="red2")
dev.off()




