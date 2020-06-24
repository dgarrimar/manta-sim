
library(ggplot2)
library(cowplot)

l <- list()
for (sim in c("simEmpirical", "simPopStructure", "simUnrelated", "simRelated")){
  va <- read.table(sprintf("simGT/%s.eigenval", sim))$V1
  ve <- read.table(sprintf("simGT/%s.eigenvec", sim))
  
  pve <- (va/sum(va))[1:2]
  df <- ve[,2:3]
  colnames(df) <- c("PC1", "PC2")
  
  l[[sim]] <- ggplot(data = df) + 
    geom_point(aes(x = PC1, y = PC2)) +
    theme_classic(base_size = 22) +
    xlab(sprintf("PC1(%.2f)", pve[1])) + 
    ylab(sprintf("PC2(%.2f)", pve[2])) +
    ggtitle(sim) + theme(plot.title = element_text(hjust = 0.5)) 
}

plot_grid(plotlist = l)

