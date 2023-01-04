## UI ####
aboutUI <- function(){
  tagList(
    htmlOutput("text1")
  )
}


## Serveur Functions ####
about <- function(){
  Text1 <- function(){
    text = tags$span(
      tags$h3("Explore cell states by MFIs", style = "color: #337ab7;"),
        tags$b("MFIs",style = "color: #337ab7;"),  "takes in scRNA-seq count matrix and estimate gene interactions. Here we show how to use these MFIs to explore cell states.",
        tags$br(),
      
      tags$h4("Data Visualization &  Analysis", style="color: #337ab7;"),
        tags$li("Table - A Summary statistics for deviating state", style="list-style-type: square;"),
        tags$li("Heatmaps - Over-representation test for MFIs and other cell annotations", style="list-style-type: square;"),
        tags$li("GO & KEGG for genes in each state", style="list-style-type: square;"),
        tags$li("rrvgo - Simplifying the redundance of GO sets", style="list-style-type: square;"),
        tags$li("Upset Plot", style="list-style-type: square;"),
        tags$li("DE analysis for mutually exclusive states", style="list-style-type: square;"),
        tags$br(),
      
      tags$h4("Tutorial", style="color: #337ab7;"),
        tags$br(),
        tags$br(),
      
      html = TRUE)
  }
  return(Text1)
}


## optionsTutorial UI:
InformationUI <- function(){
  text = tags$span(
    br("If you like MFIs and use it, please consider citing the related article:",
       tags$br(),
       a("Higher-order interactions in statistical physics and machine learning: A model-independent solution to the inverse problem at equilibrium", href = "https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.053314", target="_blank", style = "color: steelblue;"),
       a("Sjoerd Viktor Beentjes and Ava Khamseh, Phys. Rev. E. 2020 Nov 102, 053314", href= "https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.053314", target="_blank", style = "color: grey;")),
    
    tags$br(),
    tags$b("MFIs",style = "color: #337ab7;"), "application code is available through Github", a(icon("house"), href="https://github.com/AJnsm", target="_blank", style = "color: steelblue;"),
    tags$br(),
    br( "If you have any question, you can send an e-mail",a(icon("mail-bulk"), href="mailto: s1914230@ed.ac.uk",  style = "color: steelblue;")) 
  )
}





## Output to UI ####
aboutOutput <- function(output,Text1){
  output$text1 <- renderUI({  Text1() })
}
