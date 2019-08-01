library(data.table)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(fpc)

args = commandArgs(trailingOnly = TRUE)

table <- fread(file = args[1], sep = ",", data.table = FALSE, header = TRUE)
table <- table[,colSums(table != 0) > 0]
pca <- prcomp(table[,c(-1)], scale. = TRUE)
pca_vals <- as.data.frame(pca$x)

#get 90% variance pca number
summary <- data.frame(summary(pca)$importance)
cumsum <- as.vector(summary[3,])
numpca <- length(cumsum[cumsum<.9])
pamk <- pamk(pca_vals[,c(1:numpca)], krange = 1:10, usepam = FALSE)
pca_vals$new_method <- as.factor(pamk$pamobject$clustering)

clust_melt <- read.csv(args[2], header = TRUE, sep = ",")
clust_melt <- clust_melt[order(clust_melt$frame),]
pca_vals <- pca_vals[clust_melt$frame,]
pca_vals$old_method <- factor(clust_melt$cluster)
ggsave(paste0(dirname(args[2]), "/", "old_method.png"), ggplot(pca_vals, aes(x=PC1, y=PC2, color=old_method)) + geom_point())
ggsave(paste0(dirname(args[2]), "/", "new_method.png"), ggplot(pca_vals, aes(x=PC1, y=PC2, color=new_method)) + geom_point())
p <- plot_ly(pca_vals, x=~PC1, y=~PC2, z=~PC3, color = ~old_method) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))
htmlwidgets::saveWidget(p, paste0(dirname(args[2]), "/", "3d_old.html"), selfcontained = FALSE)
p <- plot_ly(pca_vals, x=~PC1, y=~PC2, z=~PC3, color = ~new_method) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))
htmlwidgets::saveWidget(p, paste0(dirname(args[2]), "/", "3d_new.html"), selfcontained = FALSE)

