options(warn=-1)

library(data.table)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(fpc)
library(markovchain)

#Arguments: 
#1) input feature table 
#2) output directory for figures and data
#3) output file for cluster list

args = commandArgs(trailingOnly = TRUE)
if(length(args) == 0){
  table <- fread(file = "clusterres_dihedaral_pca_coords.csv", sep = ",", data.table = FALSE, header = TRUE)
} else {
  table <- fread(file = args[1], sep = ",", data.table = FALSE, header = TRUE)
}
table <- table[,colSums(table != 0) > 0]
#table <- read.csv("all_pca_coords.csv",header = TRUE, sep = ",")

if (length(args) > 0) {
  outputDir = paste(args[2],"/",sep = "")
} else {
  outputDir = ""
}


pca <- prcomp(table[,c(-1)], scale. = TRUE)

#pull pca and add cluster
pca_vals <- as.data.frame(pca$x)

#get 90% variance pca number
if(ncol(pca_vals) > 4) {
  summary <- data.frame(summary(pca)$importance)
  cumsum <- as.vector(summary[3,])
  numpca <- length(cumsum[cumsum<.9])
} else {
  numpca <- ncol(pca_vals)
}


#run k means on the 90% pcas and plot
#k <- kmeans(pca_vals[,c(1:numpca)], 3, nstart=25, iter.max=1000)
if (nrow(pca_vals) < 3) {
  pca_vals$pamk <- as.factor(rep(1, nrow(pca_vals)))
} else if (nrow(pca_vals) < 11) {
  pamk <- pamk(pca_vals[,c(1:numpca)], krange = 1:(nrow(pca_vals) - 1), usepam = TRUE)
  pca_vals$pamk <- as.factor(pamk$pamobject$clustering)
} else {
  pamk <- pamk(pca_vals[,c(1:numpca)], krange = 1:10, usepam = FALSE)
  pca_vals$pamk <- as.factor(pamk$pamobject$clustering)
}

#pca_vals$k <- as.factor(k$cluster)


#see top contributors to variance
#sort(abs(pca$rotation[,1]), decreasing = TRUE)[1:10]
#sort(abs(pca$rotation[,2]), decreasing = TRUE)[1:10]



x <- pca_vals$pamk
trans <- markovchainFit(x)
#num <- max(x)
# trans <- matrix(nrow = num, ncol = num, 0)
# for (t in 1:(length(x) - 1)) trans[x[t], x[t + 1]] <- trans[x[t], x[t + 1]] + 1
# for (i in 1:num) trans[i, ] <- trans[i, ] / sum(trans[i, ])

mc <- new("markovchain", states = trans$estimate@states, transitionMatrix = trans$estimate@transitionMatrix, name = "traj")
mcss <- as.data.frame(cbind(trans$estimate@states,t(steadyStates(mc))))
mcss[,1] <- as.factor(mcss[,1])
colnames(mcss) <- c("Cluster", "Steady State Probabilities")
#ggplot(mcss, aes(x=V1, y=V2)) + geom_bar(stat = "identity")


##Output all figures
pdf(paste(outputDir,"markov_model.pdf", sep = ""))
plot(mc)
dev.off()

ggsave(filename = paste(outputDir,"cluster_pca.jpg",sep = ""), ggplot(pca_vals, aes(x=PC1, y=PC2, color=pamk)) + geom_point())

ggsave(filename = paste(outputDir,"cluster_pca_centroids.jpg", sep = ""), ggplot(pca_vals, aes(x=PC1, y=PC2, color=pamk)) + geom_point() + 
         geom_point(data = as.data.frame(pamk$pamobject$medoids), aes(x=PC1, y=PC2), colour="black", size = 5, pch = 1))

p <- plot_ly(pca_vals, x=~PC1, y=~PC2, z=~PC3, color = ~pamk) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))
htmlwidgets::saveWidget(p, paste(outputDir,"clustering_3d.html", sep = ""), selfcontained = FALSE)



#Output data
clusters <- data.frame()
write(paste(pamk$pamobject$i.med[1], paste(which(pca_vals$pamk == 1 & !c(1:length(pca_vals$pamk)) %in% pamk$pamobject$i.med), collapse = ","), sep = ","), args[3])
if (pamk$nc >= 2) {
  for(c in c(2:pamk$nc)) {
    write(paste(pamk$pamobject$i.med[c], paste(which(pca_vals$pamk == c & !c(1:length(pca_vals$pamk)) %in% pamk$pamobject$i.med), collapse = ","), sep = ","), args[3], append = TRUE)
  }
}

write.csv(mcss, paste(outputDir,"state_probabilities.txt", sep = ""), row.names = FALSE, quote = FALSE)

#clusters with SE > .1 outputted for further sampling

weakcluster <- which(unname(apply(trans$standardError, 1, function(x){sum(x>.1)})) > 0)
sink(paste(outputDir,"weak_clusters.txt", sep = ""))
cat(paste(weakcluster,collapse=","))
sink()

clusters <- data.frame(frame = c(1:length(pamk$pamobject$clustering)), cluster = pamk$pamobject$clustering)
write.csv(clusters, paste(outputDir,"../cluster_by_frame.csv", sep = ""), row.names = FALSE, quote = FALSE, eol = "\n")

###extra
# tsne <- Rtsne(X = table[,-1], dims = 2, perplexity=50, verbose=TRUE, max_iter = 850, initial_dims = numpca)
# tsne_vals <- as.data.frame(tsne$Y)
# tsne_vals$cluster <- factor(clust_melt$cluster)
# tsne_vals <- tsne_vals[tsne_vals$cluster != 5,]
# 
# plot_ly(tsne_vals, x=~V1, y=~V2, z=~V3, color = as.factor(pamk$pamobject$clustering)) %>% add_markers() %>%
#   layout(scene = list(xaxis = list(title = 'PC1'),
#                       yaxis = list(title = 'PC2'),
#                       zaxis = list(title = 'PC3')))
# htmlwidgets::saveWidget(p, "test.html")
# ggplot(tsne_vals, aes(x=V1, y=V2, color=cluster)) + geom_point()
# 
# vartsne <- function(complexity) {
#   tsne <- Rtsne(X = table[,-1], dims = 2, perplexity=complexity, verbose=TRUE, max_iter = 1000, initial_dims = 114)
#   tsne_vals <- as.data.frame(tsne$Y)
#   x <- ggplot(tsne_vals, aes(x=V1, y=V2)) + geom_point()
#   ggsave(paste(complexity, "backbone_tsne.jpg"),x)
# }
# 
# 
# #comparing to old clustering method
# 
# #read python script parsing
# clust_melt <- read.csv("~/Documents/georgetown/summer_2019/pca/clusterlist_5.csv", header = FALSE, sep = ",")
# colnames(clust_melt) <- c("frame","cluster")
# clust_melt <- clust_melt[order(clust_melt$frame),]
# pca_vals$cluster <- factor(clust_melt$cluster)
# pca_vals <- pca_vals[pca_vals$cluster != 5,]
# ggplot(pca_vals, aes(x=PC1, y=PC2, color=cluster)) + geom_point()
# plot_ly(pca_vals, x=~PC1, y=~PC2, z=~PC3, color = ~cluster) %>% add_markers() %>%
#   layout(scene = list(xaxis = list(title = 'PC1'),
#                       yaxis = list(title = 'PC2'),
#                       zaxis = list(title = 'PC3')))
