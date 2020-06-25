
library(ggplot2)
library(cowplot)
library(data.table)

path = "GT/1k_wg"
l <- list()
for (sim in c("simEmpirical", "simPopStructure", "simUnrelated", "simRelated")){
  va <- fread(sprintf("%s/%s.eigenval", path, sim), data.table = F)$V1
  ve <- fread(sprintf("%s/%s.eigenvec", path, sim), data.table = F)
  
  pve <- (va/sum(va))[1:2] * 100
  df <- ve[,3:4]
  colnames(df) <- c("PC1", "PC2")
  
  l[[sim]] <- 
    ggplot(data = df) + 
    geom_point(aes(x = PC1, y = PC2)) +
    theme_classic(base_size = 22) +
    xlab(sprintf("PC1(%.2f)", pve[1])) + 
    ylab(sprintf("PC2(%.2f)", pve[2])) +
    ggtitle(sim) + theme(plot.title = element_text(hjust = 0.5)) 
}

plot_grid(plotlist = l)

