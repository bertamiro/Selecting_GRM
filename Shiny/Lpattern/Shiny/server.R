library(shiny)
library(lpattern)

data(MethData)
data(ExpData)

shinyServer(function(input, output) {
  
  output$method1 <- renderText({
    paste("Results for", input$method)
  })
  
  output$table1 <- renderTable({
  if (input$go1 == FALSE) {return()}
  isolate({
    genes <- if (input$method == "Gene selection based on Conditional Mutual Information") {
      cMI_genes(MethData, ExpData, input$h, input$r, input$cMI0)
  }
  print(genes)
  })
  })
  
  output$table2 <- renderTable({
    if (input$go2 == FALSE) {return()}
    isolate({
      #genes <- input$h
      genes <- if (input$method == "Gene selection based on Regression based on Splines") {
        splines_genes(MethData, ExpData, 
                      input$nSI, 
                      input$nID, 
                      input$nSD, 
                      input$met_max, 
                      input$dif_met)
      }
      print(genes)
    })
  })
  
  })




  
