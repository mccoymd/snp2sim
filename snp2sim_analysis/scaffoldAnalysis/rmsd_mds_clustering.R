library(data.table)
library(ggplot2)
library(fpc)
library(plotly)
library(htmlwidgets)
3library(markovchain)
#Arguments: 
#1) input feature table 
#2) output directory for figures and data
#3) output file for cluster list

args = commandArgs(trailingOnly = TRUE)
if(length(args) == 0){
  table <- fread(file = "cluster_sample_rmsd.csv", sep = ",", data.table = FALSE, header = FALSE)
} else {
  table <- fread(file = args[1], sep = ",", data.table = FALSE, header = FALSE)
}

table <- as.dist(t(table))
#table[lower.tri(table)] <- t(table)[lower.tri(t(table))]

if (length(args) > 0) {
  outputDir = paste(args[2],"/",sep = "")
} else {
  outputDir = ""
}

fit <- cmdscale(as.dist(table), 25, eig = TRUE)

pamk <- pamk(fit$points, krange = 1:10, usepam = FALSE)

x <- pamk$pamobject$clustering
# trans <- markovchainFit(x)
# mc <- new("markovchain", states = trans$estimate@states, transitionMatrix = trans$estimate@transitionMatrix, name = "traj")
# mcss <- as.data.frame(cbind(trans$estimate@states,t(steadyStates(mc))))
# mcss[,1] <- as.factor(mcss[,1])
# colnames(mcss) <- c("Cluster", "Steady State Probabilities")
# pdf(paste(outputDir,"markov_model.pdf", sep = ""))
# plot(mc)
# dev.off()

plot1 <- ggplot(as.data.frame(fit$points), aes(x = V1, y = V2, color = as.factor(pamk$pamobject$clustering))) + geom_point()
plot2 <- ggplot(as.data.frame(fit$points), aes(x=V1, y=V2, color= as.factor(pamk$pamobject$clustering))) + geom_point() + 
  geom_point(data = as.data.frame(pamk$pamobject$medoids), aes(x=V1, y=V2), colour="black", size = 5, pch = 1)
threed <- plot_ly(as.data.frame(fit$points), x=~V1, y=~V2, z=~V3, color = as.factor(pamk$pamobject$clustering)) %>% add_markers() %>%
  layout(scene = list(xaxis = list(title = 'MDS1'),
                      yaxis = list(title = 'MDS2'),
                      zaxis = list(title = 'MDS3')))

ggsave(filename = paste(outputDir,"cluster_pca.jpg",sep = ""), plot1)
ggsave(filename = paste(outputDir,"cluster_pca_centroids.jpg",sep = ""), plot2)
htmlwidgets::saveWidget(threed, paste(outputDir,"clustering_3d.html", sep = ""), selfcontained = FALSE)

clusters <- data.frame()
write(paste(pamk$pamobject$i.med[1], paste(which(x == 1 & !c(1:length(x)) %in% pamk$pamobject$i.med), collapse = ","), sep = ","), args[3])
if (pamk$nc > 2) {
  for(c in c(2:pamk$nc)) {
    write(paste(pamk$pamobject$i.med[c], paste(which(x == c & !c(1:length(x)) %in% pamk$pamobject$i.med), collapse = ","), sep = ","), args[3], append = TRUE)
  }
}

# write.csv(mcss, paste(outputDir,"state_probabilities.txt", sep = ""), row.names = FALSE, quote = FALSE)
# 
# weakcluster <- which(unname(apply(trans$standardError, 1, function(x){sum(x>.1)})) > 0)
# sink(paste(outputDir,"weak_clusters.txt", sep = ""))
# cat(paste(weakcluster,collapse=","))
# sink()

clusters <- data.frame(frame = c(1:length(pamk$pamobject$clustering)), cluster = pamk$pamobject$clustering)
write.csv(clusters, paste(outputDir,"../cluster_by_frame.csv", sep = ""), row.names = FALSE, quote = FALSE)
