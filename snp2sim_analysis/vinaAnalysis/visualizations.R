#! /usr/local/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
path = paste0(dirname(args[1]), "/")
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

wtdata <- aggregate(affinity ~ ligand, wtdata, min)

fulldata <- merge(fulldata, wtdata, by = "ligand")
fulldata$relEnergy <- fulldata$affinity.x - fulldata$affinity.y
colnames(fulldata)[colnames(fulldata) == "affinity.x"] <-
  "absAffinity"
colnames(fulldata)[colnames(fulldata) == "affinity.y"] <-
  "wtAffinity"

#ligand vs relative energy, split by variant
ggsave(
  paste0(path, "lig_relEnergy.jpg"),
  ggplot(fulldata, aes(ligand, relEnergy, fill = scaffold)) +
    theme_light() +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Affinity Relative to WT (kcal/mol)") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(. ~ variant)
)

#variant vs relative energy, split by variant and ligand
ggsave(
  paste0(path, "variant_relEnergy.jpg"),
  ggplot(fulldata, aes(variant, relEnergy, fill = scaffold)) +
    theme_light() +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    ) +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Affinity Relative to WT (kcal/mol)") +
    theme(axis.text.x = element_blank()) +
    facet_grid(variant ~ ligand)
)

#variant vs relative energy, split by ligand
ggsave(
  paste0(path, "variant_relEnergy_ligandsplit.jpg"),
  ggplot(fulldata, aes(variant, relEnergy, fill = scaffold)) +
    theme_light() +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    ) +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Affinity Relative to WT (kcal/mol)") +
    theme(axis.text.x = element_blank()) +
    facet_grid(. ~ ligand)
)

fulldata$perChange <-
  (fulldata$relEnergy / fulldata$wtAffinity) * 100

#ligand vs percent change, split by variant
ggsave(
  paste0(path, "lig_perChange.jpg"),
  ggplot(fulldata, aes(ligand, perChange, fill = scaffold)) +
    theme_light() +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Affinity Percent Change (%)") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(. ~ variant)
)

#variant vs percent change, split by variant and ligand, colored by scaffold
ggsave(
  paste0(path, "variant_relEnergy.jpg"),
  ggplot(fulldata, aes(variant, perChange, fill = scaffold)) +
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
    facet_grid(variant ~ ligand)
)

#variant vs percent change, split by ligand, colored by variant
ggsave(
  paste0(path, "variant_relEnergy_2.jpg"),
  ggplot(fulldata, aes(variant, perChange, fill = variant)) +
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

#To do: instead of average, use scaffold representation as weighting
a <- aggregate(relEnergy ~ variant + ligand, fulldata, mean)
fulldata <- merge(fulldata, a, by = c("ligand", "variant"))
colnames(fulldata)[length(fulldata)] <- "averageEnergy"
colnames(fulldata)[colnames(fulldata) == "relEnergy.x"] <-
  "relEnergy"
#variant vs average energy (over scaffolds), split by ligand
ggsave(
  paste0(path, "variant_avgEnergy.jpg"), 
  ggplot(fulldata, aes(variant, averageEnergy, fill = variant)) +
    theme_light() +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    ) +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Average Binding Affinity Relative to WT (kcal/mol)") +
    theme(axis.text.x = element_blank()) +
    facet_grid(. ~ ligand)
)

#variant vs average energy (over scaffold), no split
ggsave(
  paste0(path, "variant_avgEnergy_2.jpg"),
  ggplot(fulldata, aes(variant, averageEnergy, fill = ligand)) +
    theme_light() +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    ) +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Average Binding Affinity Relative to WT (kcal/mol)")
)

