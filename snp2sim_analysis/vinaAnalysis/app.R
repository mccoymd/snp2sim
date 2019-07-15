library(shiny)
library(ggplot2)
library(highcharter)
library(viridis)

args = commandArgs(trailingOnly = TRUE)
table <- read.table(args[1], header = TRUE, sep = "")
pdbSet <- unique(table$scaffold)
ligSet <- unique(table$ligand)
varSet <- unique(table$variant)
data <- table

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
      downloadButton('download', 'Download figure')
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        id = "drugbinding_plots",
        tabPanel("Barchart",plotOutput("plot")),
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
      sliderInput('numgraphsrow', label = "Number of charts per row",
                  min = 1, max = length(unique(fulldata[,input$facety])),
                  value = length(unique(fulldata[,input$facety])),
                  step = 1)
    }
  })
  output$heat <- renderHighchart({
    h <- hchart(fulldata, "heatmap", hcaes(x = variant, y = ligand, value = perChange)) %>%
      hc_colorAxis(stops = color_stops(40, inferno(40))) %>% 
      hc_title(text = paste0("Binding energy of small molecules in ",unique(fulldata$protein)[1]," variants"))
  })
  output$plot <- renderPlot({
    part <- fulldata[fulldata$variant %in% input$bargene && fulldata$ligand %in% input$barlig,]
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
      plots$plot = plot + labs(y="Binding Energy Relative to WT (kcal/mol)")
    }
    else if (input$y == "perChange") {
      plots$plot = plot + labs(y="Percent Change in Binding Energy (%)")
    }
    
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

