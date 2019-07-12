#! /usr/local/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
path = paste0(dirname(args[1]), "/figures/")
table <- read.table(args[1], header = TRUE, sep = "")
pdbSet <- unique(table$scaffold)
ligSet <- unique(table$ligand)
varSet <- unique(table$variant)
data <- table
# data <- lapply(varSet, function(variant){
#   pdb <- sample(pdbSet, 2)
#   data1 <- read.table(paste0(paste0("resultsSummary/vinaSummary", pdb[1]), ".txt"), header=TRUE, sep="\t")
#   data2 <- read.table(paste0(paste0("resultsSummary/vinaSummary", pdb[2]), ".txt"), header=TRUE, sep="\t")
#   data1 <- data1[data1$variant == variant,]
#   data2 <- data2[data2$variant == variant,]
#   data1$pdb <- rep("scaffold1", nrow(data1))
#   data2$pdb <- rep("scaffold2", nrow(data2))
#   rbind(data1[data1$variant == variant,], data2[data2$variant == variant,])
# })

#fulldata <- do.call(rbind, data)
#fulldata <- fulldata[fulldata$rank == 1,]
fulldata <- data[data$rank == 1, ]
wtdata <- subset(fulldata, fulldata$variant == "wt")
fulldata <- fulldata[fulldata$variant != "wt",]

#wtdata <- aggregate(affinity ~ ligand, wtdata, min)
wtdata <- aggregate(weighted_affinity ~ ligand, wtdata, sum)
colnames(wtdata)[colnames(wtdata) == "weighted_affinity"] <-
  "wtAffinity"
weighted_avgs <- aggregate(weighted_affinity ~ ligand + variant, fulldata, sum)

fulldata <- subset(fulldata, select = -c(rank, rmsd_ub, rmsd_lb, weighted_affinity, scaffold, affinity))

fulldata <- merge(fulldata, wtdata, by = "ligand")
fulldata <- merge(fulldata, weighted_avgs, by = c("ligand","variant"))

fulldata <- unique(fulldata)

fulldata$relEnergy <- fulldata$weighted_affinity - fulldata$wtAffinity

fulldata$perChange <-
  (fulldata$relEnergy / abs(fulldata$wtAffinity)) * 100

old_path <- path
for(lib in unique(fulldata$library)) {
  plotdata <- fulldata[fulldata$library == lib,]
  path <- paste0(old_path, lib, "/")
  if(!dir.exists(path)) {
    dir.create(path)
  }
#ligand vs relative energy, split by variant
ggsave(
  paste0(path, "lig_relEnergy.jpg"),
  ggplot(plotdata, aes(ligand, relEnergy, fill = variant, width = .5)) +
    theme_light() +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Affinity Relative to WT (kcal/mol)") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(. ~ variant)
)

#variant vs relative energy, split by ligand
ggsave(
  paste0(path, "variant_relEnergy.jpg"),
  ggplot(plotdata, aes(variant, relEnergy, fill = variant)) +
    theme_light() +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    ) +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Affinity Relative to WT (kcal/mol)") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(. ~ ligand)
)

#ligand vs percent change, split by variant
ggsave(
  paste0(path, "lig_perChange.jpg"),
  ggplot(plotdata, aes(ligand, perChange, fill = ligand)) +
    theme_light() +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Affinity Percent Change (%)") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(. ~ variant)
)

#variant vs percent change, split by ligand, colored by variant
ggsave(
  paste0(path, "variant_relEnergy_2.jpg"),
  ggplot(plotdata, aes(variant, perChange, fill = variant)) +
    theme_light() +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    ) +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Affinity Percent Change (%)") +
    theme(axis.text.x = element_blank()) +
    facet_grid(. ~ ligand)
)

# empty_bar=20
# to_add = data.frame( matrix(NA, empty_bar*length(unique(plotdata$variant)), ncol(plotdata)) )
# colnames(to_add) = colnames(plotdata)
# to_add$variant=rep(unique(plotdata$variant), each=empty_bar)
# plotdata=rbind(plotdata, to_add)
# plotdata <- plotdata[order(plotdata$variant),]
# ggplot(fulldata, aes_string("variant", "relEnergy", fill = "ligand")) +
#   theme_minimal() +
#   theme(text=element_text(size=15)) +
#   geom_bar(stat="identity",position=position_dodge(), na.rm = FALSE) + 
#   labs(y="Binding Affinity Relative to WT (kcal/mol)") +
#   theme(axis.text.x = element_text(angle = 90), panel.grid = element_blank()) + coord_polar() 
# 
 }