library(shiny)

shinyUI(fluidPage(
  titlePanel("Integrative Analysis to Select Genes Regulated by Methylation"),
  sidebarLayout(
    sidebarPanel(
      h3("Genes with L-pattern", align = "left"),
      #h6("For more information visit the article:",
         #br(),
         #p(a(href = "doc.pdf", target="_blank",
             #"Alex Sanchez-Pla et al. (2016) Integrative Analysis to Select Genes Regulated by Methylation in a Cancer Colon Study", 
              #align = "justify"))
         #),
      
      #hr(),
      
      # Seleccionamos el metodo 
      checkboxGroupInput("method", label = h4("1. Select one method"),
                         choices = list("Gene selection based on Conditional Mutual Information", 
                                        "Gene selection based on Regression based on Splines"), 
                         selected = "NULL"),
      
      #Si el metodo es CMI 
      conditionalPanel(
        condition = "input.method == 'Gene selection based on Conditional Mutual Information'",
      hr(),
      h4("2. Select computing conditions"),
      sliderInput("h", 
      label = "Bandwidth (h):",
      min = 0, max = 1, value = 0.3, step= 0.1),
      hr(),
      h4("3. Select filtering conditions"),
      sliderInput("r", 
      label = "Ratio (r) lower than:",
      min = 0, max = 1, value = 0.25, step= 0.01),
      sliderInput("cMI0", 
      label = "Unconditional cMI (cMI(0)) upper than:",
      min = 0, max = 1, value = 0.1, step= 0.01),
      
      actionButton("go1", "Show genes")
    ),
    
    #Si el metodo es B-splines
    conditionalPanel(
      condition = "input.method == 'Gene selection based on Regression based on Splines'",
      hr(),
      h4("2. Select filtering conditions"),
      sliderInput("nSI", 
                  label = "nSI - Number of observations with high expression and low methylation upper than:",
                  min = 0, max = 25, value = 3, step= 1),
      sliderInput("nID", 
                  label = "nID - Number of observations with low expression and high methylation upper than:",
                  min = 0, max = 25, value = 2, step= 1),
      sliderInput("nSD", 
                  label = "nSD - Number of observations with high expression and high methylation upper than:",
                  min = 0, max = 25, value = 1, step= 1),
      sliderInput("met_max", 
                  label = "met.max - Maximum methylation upper than:",
                  min = 0, max = 1, value = 0.05, step= 0.05),
      sliderInput("dif_met", 
                  label = "dif.met - Difference between minimum and maximum methylation upper than:",
                  min = 0, max = 1, value = 0.35, step= 0.05),
      actionButton("go2", "Show genes")
      )
    
    ),
    mainPanel(
      textOutput("method1"),
      tableOutput("table1"),
      tableOutput("table2")
    )
  )
))


  






