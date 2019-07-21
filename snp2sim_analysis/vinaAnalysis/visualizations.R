#! /usr/local/bin/Rscript

library(ggplot2)
library(htmlwidgets)
library(highcharter)
library(viridis)

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

wtdata <- aggregate(affinity ~ ligand, wtdata, min)

fulldata <- merge(fulldata, wtdata, by = "ligand")
fulldata$relEnergy <- fulldata$affinity.x - fulldata$affinity.y
colnames(fulldata)[colnames(fulldata) == "affinity.x"] <-
  "absAffinity"
colnames(fulldata)[colnames(fulldata) == "affinity.y"] <-
  "wtAffinity"

fulldata$perChange <-
  (fulldata$relEnergy / abs(fulldata$wtAffinity)) * 100

if(args[2] == "True"){
  error <- aggregate(absAffinity ~ ligand + variant, fulldata, sd)
  colnames(error)[colnames(error) == "absAffinity"] <- "std_dev"
  fulldata <- merge(fulldata, error, by = c("ligand", "variant"))
  meanval <- aggregate(absAffinity ~ ligand + variant, fulldata, mean)
  colnames(meanval)[colnames(meanval) == "absAffinity"] <- "meanAffinity"
  fulldata <- merge(fulldata, meanval,  by = c("ligand", "variant"))
  
  fulldata <- fulldata[!duplicated(fulldata[,c("ligand", "variant", "meanAffinity")]),]
  fulldata$relEnergy <- fulldata$meanAffinity-fulldata$wtAffinity
  fulldata$perChange <- (fulldata$relEnergy / abs(fulldata$wtAffinity)) * 100
  fulldata$perSD <- sqrt((fulldata$std_dev^2 / abs(fulldata$wtAffinity^2)) * 10000)
  fulldata <- subset(fulldata, select = -c(rank, rmsd_ub, rmsd_lb, absAffinity, scaffold, trial))
}

old_path <- path
for(lib in unique(fulldata$library)) {
  plotdata <- fulldata[fulldata$library == lib,]
  path <- paste0(old_path, lib, "/")
  if(!dir.exists(path)) {
    dir.create(path)
  }
  
#ligand vs relative energy, split by variant
  cur <- ggplot(plotdata, aes(ligand, relEnergy, fill = ligand, width = .5)) +
    theme_light() +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Energy Relative to WT (kcal/mol)") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(. ~ variant) +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    )
  if(args[2] == "True"){
    cur <- cur + geom_errorbar(aes(ymin=relEnergy-std_dev, ymax=relEnergy+std_dev))
  }
  ggsave(
    paste0(path, "lig_relEnergy.jpg"),
    cur, width = 20
  )
  
  #variant vs relative energy, split by ligand
  cur <- ggplot(plotdata, aes(variant, relEnergy, fill = variant)) +
    theme_light() +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Binding Energy Relative to WT (kcal/mol)") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(. ~ ligand) +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    ) 
  
  if(args[2] == "True"){
    cur <- cur + geom_errorbar(aes(ymin=relEnergy-std_dev, ymax=relEnergy+std_dev))
  }
  
  
  ggsave(
    paste0(path, "variant_relEnergy.jpg"),
    cur, width = 20
  )
  
  #ligand vs percent change, split by variant
  cur <- ggplot(plotdata, aes(ligand, perChange, fill = ligand)) +
    theme_light() +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Percent Change in Binding Energy (%)") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(. ~ variant) +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    )
  if(args[2] == "True"){
    cur <- cur + geom_errorbar(aes(ymin=perChange-perSD, ymax=perChange+perSD))
  }
  ggsave(
    paste0(path, "lig_perChange.jpg"),
    cur, width = 20
  )
  
  #variant vs percent change, split by ligand, colored by variant
  cur <- ggplot(plotdata, aes(variant, perChange, fill = variant)) +
    theme_light() +
    theme(text = element_text(size = 15)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(y = "Percent Change in Binding Energy (%)") +
    theme(axis.text.x = element_blank()) +
    facet_wrap(. ~ ligand) +
    geom_hline(
      yintercept = 0,
      color = "black",
      size = 1.3
    )
  if(args[2] == "True"){
    cur <- cur + geom_errorbar(aes(ymin=perChange-perSD, ymax=perChange+perSD))
  }
  ggsave(
    paste0(path, "variant_perChange.jpg"),
    cur, width = 20
  )
  
  if(args[2] == "True"){
    
  ggsave(
    paste0(path, "ribbon_by_variant.jpg"),
    ggplot(plotdata, aes(ligand, group = variant, fill = variant)) +
      theme_light() +
      theme(text = element_text(size = 15)) +
      geom_ribbon(aes(ymin=relEnergy-std_dev, ymax=relEnergy+std_dev)) +
      geom_line(aes(y = relEnergy)) +
      labs(y = "Binding Energy Relative to WT (kcal/mol)") +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(. ~ variant) +
      geom_hline(
        yintercept = 0,
        color = "black",
        size = 1.3
      ), width = 20
  )
  ggsave(
    paste0(path, "ribbon_by_ligand.jpg"),
    ggplot(plotdata, aes(variant, group = ligand, fill = ligand)) +
      theme_light() +
      theme(text = element_text(size = 15)) +
      geom_ribbon(aes(ymin=relEnergy-std_dev, ymax=relEnergy+std_dev)) +
      geom_line(aes(y = relEnergy)) +
      labs(y = "Binding Energy Relative to WT (kcal/mol)") +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(. ~ ligand) +
      geom_hline(
        yintercept = 0,
        color = "black",
        size = 1.3
      ), width = 20
  )
  }
  h <- hchart(fulldata, "heatmap", hcaes(x = variant, y = ligand, value = perChange)) %>%
    hc_colorAxis(stops = color_stops(40, inferno(40))) %>% 
    hc_title(text = paste0("Binding energy of small molecules in ",unique(fulldata$protein)[1]," variants"))
  
  htmlwidgets::saveWidget(h, paste0(path, "heatmap.html"), selfcontained = FALSE)
}
