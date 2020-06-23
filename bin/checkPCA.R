
library(ggplot2)
library(cowplot)

path <- "./GT"
l <- list()
for (sim in c("simEmpirical", "simPopStructure", "simUnrelated", "simRelated")){
  ve <- read.table(sprintf("%s/%s.eigenvec", path, sim))

  df <- ve[,2:3]
  colnames(df) <- c("PC1", "PC2")
  
  l[[sim]] <- ggplot(data = df) + 
    geom_point(aes(x = PC1, y = PC2)) +
    theme_classic(base_size = 22) +
    ggtitle(sim) + theme(plot.title = element_text(hjust = 0.5)) 
}

plot_grid(plotlist = l)

