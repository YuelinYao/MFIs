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
      tags$h3("Explore cell states by Stator", style = "color: #337ab7;"),
        tags$b("Stator",style = "color: #337ab7;"),  "takes in scRNA-seq count matrix, estimates higher-order gene interactions and defines cell states. Here we show how to use results from Stator to explore cell states.",
        tags$br(),
      
      tags$h4("Data Visualization &  Analysis", style="color: #337ab7;"),
        tags$li("Table - A summary statistics for deviating states.", style="list-style-type: square;"),
        tags$li("Heatmaps Cells - Over-representation test for MFIs and other external cell annotations.", style="list-style-type: square;"),  
        tags$li("Heatmaps Genes - Over-representation test for cell-state genes and external gene annotations.", style="list-style-type: square;"),
        tags$li("GO & KEGG for genes in each state.", style="list-style-type: square;"),
        tags$li("rrvgo - Simplifying the redundance of GO sets.", style="list-style-type: square;"),
        tags$li("Upset Plot - Visualisation of how many cells sharing among states.", style="list-style-type: square;"),
        tags$li("DE analysis - Differential expression analysis for two mutually exclusive states.", style="list-style-type: square;"),
        tags$li("Find Markers - Identify marker genes for a give cell state.", style="list-style-type: square;"),
        tags$li("Automatic annotations - Identify marker genes for all stator states and return DEGs in the provided gene list.", style="list-style-type: square;"),
        tags$li("Markov Blanket - Visualisation of Markov Blanket for given gene(s).", style="list-style-type: square;"),
        tags$li("UMAP Plot - Visualisation of UMAP for given a cell state.", style="list-style-type: square;"),
        tags$br(),
      
        ##Tutorial
        tags$h4("Tutorial", style="color: #337ab7;"),
        tags$b("Prepare and upload the dataset",style = "color: #337ab7;"),
        tags$br(),
        "The input files include:",
        tags$li("Count_matrix.csv: count matrix of scRNA-seq data, the same file indicated in the rawDataPath of Stator's JSON file", style="list-style-type: square;"),
        tags$li("Meta_data.csv: cell type annotations from other tools (e.g., clustering, NMFs, two coloum csv file, example in ./data folder), is only used in the overrepresentation test heatmap (heatmap-cell tabs), if you don't want to plot this heatmap, you can skip this file", style="list-style-type: square;"),
        tags$li("all_DTuples.csv: located in the output from stator pipeline (HOIsummaries folder)", style="list-style-type: square;"),
        tags$li("trainingData_.csv: which is also in the output from stator pipeline (output folder)", style="list-style-type: square;"),
        tags$li("GeneAnnotationSet.csv (optional): which is used in heatmap-genes, over-representation test between gene list. Example file format see ./data/CancerState.csv", style="list-style-type: square;"),
        tags$li("UMAP.csv (optional): which is used in UMAP plot. Example file format see ./data/UMAP_coords.csv", style="list-style-type: square;"),
        "The example of file can be found at: ",a(icon("house"), href="https://github.com/YuelinYao/MFIs/tree/main/data", target="_blank", style = "color: steelblue;"),
      
      tags$br(),
      tags$br(),
      
      #tags$b("Statistical analysis",style = "color: #337ab7;"),
      
      
      html = TRUE)
  }
  return(Text1)
}


## optionsTutorial UI:
InformationUI <- function(){
  text = tags$span(
    br("If you like Stator and use it, please consider citing the related article:",
       tags$br(),
       a("Higher-order interactions in statistical physics and machine learning: A model-independent solution to the inverse problem at equilibrium", href = "https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.053314", target="_blank", style = "color: steelblue;"),
       a("Sjoerd Viktor Beentjes and Ava Khamseh, Phys. Rev. E. 2020 Nov 102, 053314", href= "https://journals.aps.org/pre/abstract/10.1103/PhysRevE.102.053314", target="_blank", style = "color: grey;")),
    
    tags$br(),
    tags$b("Stator",style = "color: #337ab7;"), "application code is available through Github", a(icon("house"), href="https://github.com/AJnsm", target="_blank", style = "color: steelblue;"),
    tags$br(),
    br( "If you have any question, you can send an e-mail",a(icon("mail-bulk"), href="mailto: s1914230@ed.ac.uk",  style = "color: steelblue;")) 
  )
}





## Output to UI ####
aboutOutput <- function(output,Text1){
  output$text1 <- renderUI({  Text1() })
}
