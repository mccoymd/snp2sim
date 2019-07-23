library(shiny)
library(ggplot2)
library(highcharter)
library(viridis)
library(plotly)

args = commandArgs(trailingOnly = TRUE)
table <- read.table(args[1], header = TRUE, sep = "")
pdbSet <- unique(table$scaffold)
ligSet <- unique(table$ligand)
varSet <- unique(table$variant)
data <- table

fulldata <- data[data$rank == 1, ]


if(args[2] == "True"){
  
  meanval <- aggregate(absAffinity ~ ligand + variant + scaffold, fulldata, mean)
  colnames(meanval)[colnames(meanval) == "absAffinity"] <- "meanAffinity"
  fulldata <- merge(fulldata, meanval,  by = c("ligand", "variant"))
  fulldata <- fulldata[!duplicated(fulldata[,c("ligand", "variant", "scaffold", "meanAffinity")]),]
  fulldata$weighted_affinity <- fulldata$meanAffinity * fulldata$weight
  
  
  
  error <- aggregate(absAffinity ~ ligand + variant + scaffold, fulldata, var)
  colnames(error)[colnames(error) == "absAffinity"] <- "var_scaff"
  fulldata <- merge(fulldata, error, by = c("ligand", "variant", "scaffold"))
  fulldata <- fulldata[!duplicated(fulldata[,c("ligand", "variant", "scaffold", "var_scaff")]),]
  fulldata$var_scaff <- fulldata$var_scaff * fulldata$weight^2
  
  std <- aggregate(var_scaff ~ ligand + variant, fulldata, sum)
  colnames(std)[colnames(std) == "var_scaff"] <- "std_dev"
  std$std_dev <- sqrt(std$std_dev)
  fulldata <- merge(fulldata, std, by = c("ligand", "variant"))
  fulldata <- fulldata[!duplicated(fulldata[,c("ligand", "variant", "std_dev")]),]
  
} else {
  fulldata$weighted_affinity <- fulldata$affinity * fulldata$weight
}

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

fulltable <- fulldata
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("snp2sim Drug Binding Analysis"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput('x', label = "X Data",
                  choices = c("ligand", "variant"),
                  selected = "variant"
      ),
      selectInput('y', label = "Y Data",
                  choices = c("relEnergy", "perChange"),
                  selected = "relEnergy"
      ),
      selectInput('fill', label = "Fill",
                  choices = c("None", "ligand", "variant", "library"),
                  selected = "ligand"
      ),
      selectInput('facety', label = "Group",
                  choices = c("None", "ligand", "variant", "library"),
                  selected = ""
      ),
      selectInput('bargene', "Variant:",
                  choices = unique(fulldata$variant), 
                  multiple=TRUE, 
                  selectize=TRUE,
                  selected =unique(fulldata$variant)),
      selectInput('barlig', "Ligand:",
                  choices = unique(fulldata$ligand), 
                  multiple=TRUE, 
                  selectize=TRUE,
                  selected =unique(fulldata$ligand)),
      uiOutput('numgraphsrow'),
      conditionalPanel(
        "args[2] == True",
        checkboxInput('errorBars', label = "Show error bars?",
                      value = FALSE)
      ),
      downloadButton('download', 'Download figure')
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        id = "drugbinding_plots",
        tabPanel("Barchart",plotlyOutput("plot")),
        tabPanel("Heatmap",highchartOutput("heat"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  plots <- reactiveValues()
  
  output$numgraphsrow <- renderUI({
    if(input$facety != "None") {
      if (input$facety == "ligand"){
        sliderInput('numgraphsrow', label = "Number of charts per row",
                    min = 1, max = length(input$barlig),
                    value = length(unique(fulldata[,input$facety])),
                    step = 1)
      } else if(input$facety == "variant"){
        sliderInput('numgraphsrow', label = "Number of charts per row",
                    min = 1, max = length(input$bargene),
                    value = length(unique(fulldata[,input$facety])),
                    step = 1)
      } else {
        sliderInput('numgraphsrow', label = "Number of charts per row",
                    min = 1, max = length(unique(fulldata[,input$facety])),
                    value = length(unique(fulldata[,input$facety])),
                    step = 1)
      }
    }
  })
  
  
  output$heat <- renderHighchart({
    h <- hchart(fulldata, "heatmap", hcaes(x = variant, y = ligand, value = perChange)) %>%
      hc_colorAxis(stops = color_stops(40, inferno(40))) %>% 
      hc_title(text = paste0("Binding energy of small molecules in ",unique(fulldata$protein)[1]," variants"))
  })
  output$plot <- renderPlotly({
    part <- fulldata[fulldata$variant %in% input$bargene & fulldata$ligand %in% input$barlig,]
    if (input$fill != "None"){
      plot <- ggplot(part, aes_string(input$x, input$y, fill = input$fill)) +
        theme_light() +
        theme(text=element_text(size=15)) +
        geom_bar(stat="identity",position=position_dodge()) + 
        theme(axis.text.x = element_text(angle = 90))
    }
    else{
      plot <- ggplot(part, aes_string(input$x, input$y)) +
        theme_light() +
        theme(text=element_text(size=15)) +
        geom_bar(stat="identity",position=position_dodge()) + 
        labs(y="Binding Affinity Relative to WT (kcal/mol)") +
        theme(axis.text.x = element_text(angle = 90))
    }
    if (input$facety != "None"){
      plot <- plot + facet_wrap(. ~ get(input$facety), ncol = input$numgraphsrow, scales="free_x")
    }
    if (input$y == "relEnergy") {
      plot = plot + labs(y="Binding Energy Relative to WT (kcal/mol)")
    }
    else if (input$y == "perChange") {
      plot = plot + labs(y="Percent Change in Binding Energy (%)")
    }
    
    if (args[2] == "True") {
      if(input$errorBars){
        plot <- plot + geom_errorbar(aes(ymin=relEnergy-std_dev, ymax=relEnergy+std_dev), position = position_dodge())
      }
    }
    plots$plot <- plot
    plots$plot
  })
  
  output$download <- downloadHandler(
    filename = function() {"figure.png"},
    content = function(file) {
      ggsave(file, plot = plots$plot, device = "png")
    }
  )
}

# Run the application 
runApp(list(ui = ui, server = server), launch.browser = TRUE)

