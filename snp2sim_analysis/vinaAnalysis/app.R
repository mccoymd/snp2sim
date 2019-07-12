library(shiny)
library(ggplot2)

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
      uiOutput('numgraphsrow'),
      checkboxInput('circle', label = "Circularize?", value = FALSE),
      downloadButton('download', 'Download figure')
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("plot")
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
  
  output$plot <- renderPlot({
    if (input$fill != "None"){
      plot <- ggplot(fulldata, aes_string(input$x, input$y, fill = input$fill)) +
        theme_light() +
        theme(text=element_text(size=15)) +
        geom_bar(stat="identity",position=position_dodge()) + 
        labs(y="Binding Affinity Relative to WT (kcal/mol)") +
        theme(axis.text.x = element_text(angle = 90))
    }
    else{
      plot <- ggplot(fulldata, aes_string(input$x, input$y)) +
        theme_light() +
        theme(text=element_text(size=15)) +
        geom_bar(stat="identity",position=position_dodge()) + 
        labs(y="Binding Affinity Relative to WT (kcal/mol)") +
        theme(axis.text.x = element_text(angle = 90))
    }
    if (input$facety != "None"){
      plots$plot <- plot + facet_wrap(. ~ get(input$facety), ncol = input$numgraphsrow, scales="free_x")
    }
    else{
      plots$plot <- plot 
    }
    if (input$circle) {
      plots$plot <- plots$plot + coord_polar()
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

