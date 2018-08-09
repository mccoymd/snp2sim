#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

table <- read.table("vinaSummary.PDL1cgc.txt", header=TRUE, sep="\t")
pdbSet <- unique(table$scaffold)
ligSet <- unique(table$ligand)
varSet <- unique(table$variant)
data <- table

fulldata <- data[data$rank == 1, ]
wtdata <- subset(fulldata, fulldata$variant == "wt")
fulldata <- fulldata[fulldata$variant != "wt",]

wtdata <- aggregate(affinity ~ ligand, wtdata, min)

fulldata <- merge(fulldata, wtdata, by = "ligand")
fulldata$relEnergy <- fulldata$affinity.x - fulldata$affinity.y
colnames(fulldata)[colnames(fulldata) == "affinity.x"] <- "absAffinity"
colnames(fulldata)[colnames(fulldata) == "affinity.y"] <- "wtAffinity"

fulldata$perChange <- (fulldata$relEnergy / fulldata$wtAffinity) * 100

a <- aggregate(relEnergy ~ variant + ligand, fulldata, mean)
fulldata <- merge(fulldata, a, by = c("ligand", "variant"))
colnames(fulldata)[length(fulldata)] <- "averageEnergy"
colnames(fulldata)[colnames(fulldata) == "relEnergy.x"] <- "relEnergy"

fulltable <- fulldata
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("snp2sim Drug Binding Analysis"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
          selectInput('x', label = "X Data",
                      choices = colnames(fulltable),
                      selected = "variant"
          ),
          selectInput('y', label = "Y Data",
                      choices = colnames(fulltable),
                      selected = "relEnergy"
          ),
          selectInput('fill', label = "Fill",
                      choices = c("None", colnames(fulltable)),
                      selected = "scaffold"
          ),
          selectInput('facety', label = "Group vertically",
                      choices = c("None", colnames(fulltable)),
                      selected = ""
          ),
          selectInput('facetx', label = "Group horizontally",
                      choices = c("None", colnames(fulltable)),
                      selected = ""
          ),
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
  
  fulltable <- reactiveVal({
    table <- read.table("vinaSummary.PDL1cgc.txt", header=TRUE, sep="\t")
    pdbSet <- unique(table$scaffold)
    ligSet <- unique(table$ligand)
    varSet <- unique(table$variant)
    data <- table
    
    fulldata <- data[data$rank == 1, ]
    wtdata <- subset(fulldata, fulldata$variant == "wt")
    fulldata <- fulldata[fulldata$variant != "wt",]
    
    wtdata <- aggregate(affinity ~ ligand, wtdata, min)
    
    fulldata <- merge(fulldata, wtdata, by = "ligand")
    fulldata$relEnergy <- fulldata$affinity.x - fulldata$affinity.y
    colnames(fulldata)[colnames(fulldata) == "affinity.x"] <- "absAffinity"
    colnames(fulldata)[colnames(fulldata) == "affinity.y"] <- "wtAffinity"
    
    fulldata$perChange <- (fulldata$relEnergy / fulldata$wtAffinity) * 100
    
    a <- aggregate(relEnergy ~ variant + ligand, fulldata, mean)
    fulldata <- merge(fulldata, a, by = c("ligand", "variant"))
    colnames(fulldata)[length(fulldata)] <- "averageEnergy"
    colnames(fulldata)[colnames(fulldata) == "relEnergy.x"] <- "relEnergy"
    
    fulldata
  })
  
   plots <- reactiveValues()
    
   
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
     if (input$facetx != "None" && input$facety != "None"){
       plots$plot <- plot + facet_grid(get(input$facetx) ~ get(input$facety))
     }
     else if (input$facetx == "None" && input$facety != "None"){
       plots$plot <- plot + facet_grid(. ~ get(input$facety))
     }
     else if (input$facetx != "None" && input$facety == "None"){
       plots$plot <- plot + facet_grid(get(input$facetx) ~ .)
     }
     else if (input$facetx == "None" && input$facety == "None"){
       plots$plot <- plot 
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
shinyApp(ui = ui, server = server)

