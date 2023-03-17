#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#
# ====== import R libraries: ####
# Some initial setup:
#library(BiocManager)
#options(repos = BiocManager::repositories())
# Define any Python packages needed for the app in R:
#PYTHON_DEPENDENCIES = c('pip', 'numpy','pandas','igraph','argparse','scipy',
#'matplotlib','Pillow','seaborn')
#PYTHON_DEPENDENCIES = c('numpy','pandas','scipy')

# ------------------ App virtualenv setup (Do not edit) ------------------- #
# VIRTUALENV_NAME and PYTHON_PATH are definded in .Rprofile
#virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
#python_path = Sys.getenv('PYTHON_PATH')

# Create virtual env and install dependencies
#reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
#reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES, ignore_installed=TRUE)
#reticulate::use_virtualenv(virtualenv_dir, required = T)

library(shiny)            
library(shinythemes)      
import::from(shinycssloaders, withSpinner) 
reticulate::py_config()
library(shinyBS)        
library(shinyWidgets)     
library(gridExtra, verbose=FALSE)        
library(RColorBrewer, verbose=FALSE)
library(ComplexHeatmap)
library(Seurat)
library(stringr)
library(pheatmap)
library(dplyr)
library(rrvgo)
library(DT)
library(clusterProfiler)
library(data.table)
library(biomaRt)
library(ggplot2)
library(igraph)
library(ggrepel)
#library(org.Hs.eg.db)
#library(org.Mm.eg.db)
options(shiny.maxRequestSize = 5000*1024^2)

source("./app/general/general.R") # conditionalPanel
source("./app/tabs/about/about.R") # about
source("./app/tabs/heatmap/heatmap.R") # heatmap
source("./app/tabs/heatmap_gene/heatmap_gene.R") # heatmap-gene
source("./app/tabs/GO_plot/GO_plot.R") # GOplot
source("./app/tabs/upsetplot/upsetplot.R") 
source("./app/tabs/table/table.R") 
source("./app/tabs/DE/DE.R") 
source("./app/tabs/marker/marker.R") 
source("./app/tabs/rrvgo/rrvgo.R") 
source("./app/tabs/markovBlanket/markovBlanket.R") 

# Define UI for application
ui <- fluidPage(theme = shinytheme("spacelab"),
                # Application title
                titlePanel("",windowTitle = "Stator project"),
                tags$head(
                  tags$style(
                    HTML("@import url('//fonts.googleapis.com/css?family=Righteous|');"),
                    HTML(".shiny-output-error-validation {
                                        color: red;
                          }")
                  )
                ),
                headerPanel(                              ### add logos 
                 img(src= "logo.png",height  = 120, width = 700))
                      , ##class = "pull-left" 
                
  #Sidebar with a slider input for the cutoff
            sidebarLayout(
            sidebarPanel(
 
 ## Global conditional panel for read table 
            UploadFilesUI(),  
 ## About            
            conditionalPanel(
                      condition="input.tabs == 'about'",
                      InformationUI()),
 
 ## table
            conditionalPanel(condition= "input.tabs == 'table'",
                  TableInput()),
 
 ## heatmap
            conditionalPanel(condition= "input.tabs == 'heatmap'",
                  HeatmapInput()),
 
 ## heatmap_genes
            conditionalPanel(condition= "input.tabs == 'heatmap_genes'",
                  HeatmapGenesInput()), 
 

 ## GO enrichment        
            conditionalPanel(condition= "input.tabs == 'GO'", # GO & KEGG Tab
                             GOInput()),
                            
 ## rrvgo Tab   
            conditionalPanel(condition= "input.tabs == 'rrvgo'",
                             rrvgoInput()),
 ## Upset          
            conditionalPanel(condition= "input.tabs == 'upset'",  
                  UPsetInput()),                    
 ## DE Tab 
            conditionalPanel(condition= "input.tabs == 'DE'",
                             DEInput()),

 
 ## Find Marker Tab
            conditionalPanel(condition= "input.tabs == 'Marker'",
                             MarkerInput()),
          
 ## MB Tab
            conditionalPanel(condition= "input.tabs == 'mb'",
                             MBInput()),
 
        ),
        
        
        mainPanel(
            tabsetPanel(type = "pills", id = 'tabs',
                tabPanel("About", icon = icon("circle-info"), value = 'about', aboutUI()),
                tabPanel("Table", icon = icon("table"),value="table",TableUI()),
                tabPanel("Heatmap-Cells",icon = icon("map"), value="heatmap",
                         heatmapUI(),
                         NMF_UI(),
                         NMF_CelltypeUI(),
                         Test_CellstateUI()),
                tabPanel("Heatmap-Genes",icon = icon("map"), value="heatmap_genes",
                         heatmapGenesUI(),
                         stateGenesUI()),
                tabPanel("GO & KEGG", icon = icon("chart-line"),value="GO",
                         FunctionUI(),
                         Function_tableUI()),
                tabPanel("Using rrvgo", icon = icon("chart-line"),value="rrvgo",
                         rrvgoUI(),
                         rrvgo2UI(),
                         rrvgo3UI(),
                         rrvgo_tableUI()),
                tabPanel("Upset Plot",icon = icon("signal"),value="upset",
                         UPsetUI()),
                tabPanel("DE analysis",icon = icon("random"),value="DE",
                         DE_UI(),
                         DE_TableUI(),
                         DEGO_UI(),
                         DEGOtable_UI()),
                
                tabPanel("Find Markers",icon = icon("random"),value="Marker",
                         Marker_UI(),
                         Marker_TableUI(),
                         MarkerGO_UI(),
                         MarkerGOtable_UI()),
                
                
                tabPanel("Markov Blanket",icon = icon("circle-nodes"),value="mb",
                         MBUI()),
        ))   )

)



# Define server 
server <- function(input, output,session) {

  ## Load example data or read uploaded data: count matrix
  usedTable <- reactive({
    if(input$TestTable){
      load("./data/Example.Rdata")
      data<-list(count,srt)
      names(data)<-c("count","srt")
      return(data)
    }
    else{
      file_countmatrix<-input$Count_matrix
      data<-read_data(file_countmatrix$datapath)
      return(data)
    }
  })  
  
  ## Load example data or read uploaded data: topDeviatingHOIstates.csv and trainingData_.csv
  usedTable2 <- reactive({
    print(paste("Use HCC cancer data: ",input$TestTable))
    if(input$TestTable){
      devStates="./data/topDeviatingHOIstates.csv"
      trainDat="./data/trainingData_CL01_14698Cells_1000Genes.csv"
      data_path<-list(devStates,trainDat)
      names(data_path)<-c("devStates","trainDat")
      return(data_path)
    }
    
    else{
      file_devStates<-input$topDeviatingHOIstates
      file_trainDat<-input$trainingData_matrix
      data_path<-list(file_devStates$datapath,file_trainDat$datapath)
      names(data_path)<-c("devStates","trainDat")
      return(data_path)
    }
  })  
  
  ## Load example data or read uploaded data: Meta data
  usedMeta_data <- reactive({
    print(paste("Use HCC cancer data: ",input$TestTable))
    if(input$TestTable){
      Meta_data_path<-"./data/Meta_Data.csv"
      Meta_data<-read.csv(Meta_data_path,row.names = 1,colClasses = "character")
    } else{
      Meta_data_path<-input$Meta_Data 
      Meta_data<-read.csv(Meta_data_path$datapath,row.names = 1,colClasses = "character")
    }
      
      colnames(Meta_data)<-c("Cell_State","Cell_Types")
      return(Meta_data)
    
  })  
  
  
  ## Load example data or read uploaded data: Meta data
  usedMCMCGraph <- reactive({
    print(paste("Use HCC cancer data: ",input$TestTable))
    if(input$TestTable){
      MCMCGraph_path<-"./data/MCMCgraph_CL01_14698Cells_1000Genes.csv"
      MCMCGraph<-read.csv(MCMCGraph_path,header=T,row.names = 1,colClasses = "character")
   
    } else{
      MCMCGraph_path<-input$MCMCgraph
      MCMCGraph<-read.csv(MCMCGraph_path$datapath,header=T,row.names = 1,colClasses = "character")
    }
    
    MCMCGraph<-graph_from_adjacency_matrix(as(MCMCGraph, 'matrix'), mode="directed")
    return(MCMCGraph)
  })  
  
  
  
  
  ## Show example background genes for functional analysis
  # Background genes:
  BackgroundGenes<-reactive({
    background<-strsplit(input$background_genes, "\n")[[1]]
  }) 
  
  BackgroundGenes2<-reactive({
    background<-strsplit(input$background_genesrrgvo, "\n")[[1]]
  }) 
  
  BackgroundGenes3<-reactive({
    background<-strsplit(input$background_genesDE, "\n")[[1]]
  }) 
  
  
  BackgroundGenes4<-reactive({
    background<-strsplit(input$background_genesMarker, "\n")[[1]]
  }) 
  
  
  # Load background genes
  observe({
    if (input$bg_Liver0) {
      updateTextInput(session, "background_genes", value = background)
    }
  })
  
  observe({
    if (input$bg_Liver) {
      updateTextInput(session, "background_genesrrgvo", value = background)
    }
  })
  
  observe({
    if (input$bg_Liver2) {
      updateTextInput(session, "background_genesDE", value = background)
    }
  })
  
  
  observe({
    if (input$bg_Liver3) {
      updateTextInput(session, "background_genesMarker", value = background)
    }
  })
  

  ## Set up:
  # Summary Table
    summaryTable<-reactive({
      Table_cluster(input$cutoff,usedTable2(),input$minStateDeviation,input$minNoCells)
    }) %>% bindCache(input$cutoff,usedTable2(),input$minStateDeviation,input$minNoCells)
    
    # Get cell list, function in heatmap.R
    List<-reactive({
      GetCellList(usedTable()$count,summaryTable())#input$cutoff,usedTable()$count,
    })%>% bindCache(input$cutoff,usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) 
    
    
    # Get gene list: function in GO_plot.R
    GetGenesList<-reactive({
      GetGenes(summaryTable(),Remove=T)
    })%>% bindCache(input$cutoff,usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) 
    
    
    # Get gene list, this gene list consist all the genes, function in GO_plot.R
    GetGenesList_All<-reactive({
      GetGenes(summaryTable(),Remove=F)
    })%>% bindCache(input$cutoff,usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) 
    
    
    # About: include overview and tutorial
    Text1 <- about()
    aboutOutput(output,Text1)
    
    
    ## Table:
    showText<-eventReactive(input$action_table, {
      paste(dim(summaryTable())[1]," deviating MFIs in total, ",length(unique(summaryTable()$cluster)), "clusters." )
    })
    
    output$textsummary <- renderText({
      showText()
    })
    
    show_Table <- eventReactive(input$action_table, { 
      table<-summaryTable()
      table$cluster<-paste0("Cluster:",table$cluster)
      colnames(table)<-c("Genes","D-tuple","Enrichment","Pval_corrected","No.Cells","Cluster")
      table$`D-tuple`<-paste0("[",table$`D-tuple`,"]")
      table
    })
    
    output$table <- renderDT({
      show_Table()
    })
    
    output$downloadtable<-downloadHandler(
      filename=function(){
        paste0("State_Table-",Sys.Date(), '.csv')
      },
      content = function(file){
        write.csv(show_Table(),file)
      }
    )

    
    ## Heatmap-genes:
    usedCellStateGene <- reactive({
      print(paste("Use CellStateGene: ",input$CellStateGenes))
      if(!input$CellStateGenes=="UploadState"){
        CellStateGenes_path<-paste0("./data/",input$CellStateGenes,".csv")
      } else{
        CellStateGenes_path<-input$uploadCellStateGenes
        CellStateGenes_path<-CellStateGenes_path$datapath
      }
      return(CellStateGenes_path)
    })  
    
    # color by -log10FDR or enrichment fold
    colorHeatmapGene<- eventReactive(input$action_heatmap_genes, { 
      input$colorHeatmapGene
    })
    
    result_genes <- eventReactive(input$action_heatmap_genes,{
      result_genes<-heatmapGenes(usedCellStateGene(),GetGenesList_All(),N=as.integer(input$NO.background))
      #N=25678 the number of genes in whole human genome
      result_genes$cutoff=input$cutoff
      result_genes
    } )
    
    output$heatmap_GeneSet <- renderPlot({
      pheatmap::pheatmap(border_color = NA,result_genes()[[colorHeatmapGene()]],display_numbers = result_genes()[["Mydata_raw_m"]],
                         fontsize = 12,fontsize_number = 15,
                         main=result_genes()[["cutoff"]],
                         cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                         color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
    })
    
   ## Download csv
    output$downloadheatmap_genes<-downloadHandler(
      filename = function(){
        paste('HeatmapGenes-', colorHeatmapGene(),"-",Sys.Date(), '.csv', sep='')
      },
      content=function(heatmap1){
        write.csv(result_genes()[[colorHeatmapGene()]],heatmap1)
      }
    )
    
    ## Download heatmap pdf
    output$downloadheatmapgenes_plot<-downloadHandler(
      filename = function(){
        paste('HeatmapGenes-',colorHeatmapGene(),"-",Sys.Date(), '.pdf', sep='')
      },
      content=function(heatmap1){
        width=dim(result_genes()[["log10FDR"]])[2]*5/25.4+max(nchar(rownames(result_genes()[["log10FDR"]])))*5/25.4+2
        height=dim(result_genes()[["log10FDR"]])[1]*5/25.4+max(nchar(colnames(result_genes()[["log10FDR"]])))*5/25.4+5
        pdf(heatmap1,width =width,height = height)
        pheatmap::pheatmap(border_color = NA,result_genes()[[colorHeatmapGene()]],display_numbers = result_genes()[["Mydata_raw_m"]],
                           fontsize = 12,fontsize_number = 15,cellwidth = 15,cellheight = 15,
                           main=result_genes()[["cutoff"]],
                           cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
        dev.off()
        
      } )
    
    ### StateGenes:
    StateGenes <- eventReactive(input$action_heatmap_genes,{
      StateGenesResult<-ListTest(GetGenesList_All(),as.integer(input$NO.background))
      StateGenesResult$cutoff=input$cutoff
      StateGenesResult
    })
    
    output$heatmapStateGenes<- renderPlot({
      pheatmap::pheatmap(border_color = NA,   StateGenes()[[colorHeatmapGene()]],display_numbers = StateGenes()$Mydata_raw_m,
                         fontsize = 12,fontsize_number = 15,
                         cluster_cols = T,cluster_rows = T,main=StateGenes()[["cutoff"]],treeheight_row = 0, treeheight_col = 0,
                         color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
    })
    ## Download csv
    output$downloadheatmap_genes2<-downloadHandler(
      filename = function(){
        paste('HeatmapStateGenes-', colorHeatmapGene(),"-",Sys.Date(), '.csv', sep='')
      },
      content=function(heatmap1){
        write.csv(StateGenes()[[colorHeatmapGene()]],heatmap1)
      }
    )
    
    ## Download pdf
    output$downloadheatmapgenes_plot2<-downloadHandler(
      filename = function(){
        paste('HeatmapGenes-', colorHeatmapGene(),"-",Sys.Date(), '.pdf', sep='')
      },
      content=function(heatmap1){
        width=dim( StateGenes()[["log10FDR"]])[2]*5/25.4+max(nchar(rownames( StateGenes()[["log10FDR"]])))*5/25.4+2
        height=dim( StateGenes()[["log10FDR"]])[1]*5/25.4+max(nchar(colnames( StateGenes()[["log10FDR"]])))*5/25.4+5
        pdf(heatmap1,width =width,height = height)
        pheatmap::pheatmap(border_color = NA, StateGenes()[[colorHeatmapGene()]],display_numbers =  StateGenes()[["Mydata_raw_m"]],
                           fontsize = 12,fontsize_number = 15,cellwidth = 15,cellheight = 15,
                           main= StateGenes()[["cutoff"]],
                           cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
        dev.off()
        
      } )
    
  
    
    ## Heatmap for cells(Tab)
    colorHeatmapCells<- eventReactive(input$action_heatmap, { 
      input$colorHeatmapCells
    })
    
    result1 <- eventReactive(input$action_heatmap,{
      result1<-heatmap(usedMeta_data(),summaryTable(), List(),N=dim(usedTable()$count)[1])
      result1$cutoff=input$cutoff
      result1
    } )
    
    output$heatmap_celltypes <- renderPlot({
    pheatmap::pheatmap(border_color = NA,result1()[[colorHeatmapCells()]],display_numbers = result1()[["Mydata_raw_m"]],
                           fontsize = 12,fontsize_number = 15,
                            main=result1()[["cutoff"]],
                           cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
    })
    
    ## Download csv
    output$downloadheatmap1<-downloadHandler(
      filename = function(){
        paste('Heatmap1-', colorHeatmapCells(),"-",Sys.Date(), '.csv', sep='')
      },
      content=function(heatmap1){
        write.csv(result1()[[colorHeatmapCells()]],heatmap1)
      }
   )
    ## Download pdf
    output$downloadheatmap_plot1<-downloadHandler(
      filename = function(){
        paste('Heatmap1-',colorHeatmapCells(),"-",Sys.Date(), '.pdf', sep='')
      },
      content=function(heatmap1){
        width=dim(result1()[["log10FDR"]])[2]*5/25.4+max(nchar(rownames(result1()[["log10FDR"]])))*5/25.4+2
        height=dim(result1()[["log10FDR"]])[1]*5/25.4+max(nchar(colnames(result1()[["log10FDR"]])))*5/25.4+5
        pdf(heatmap1,width =width,height = height)
        pheatmap::pheatmap(border_color = NA,result1()[[colorHeatmapCells()]],display_numbers = result1()[["Mydata_raw_m"]],
                           fontsize = 12,fontsize_number = 15,cellwidth = 15,cellheight = 15,
                           main=result1()[["cutoff"]],
                           cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
        dev.off()

     } )
    
    
    # Heatmap2 (Tab)
    result2 <- eventReactive(input$action_heatmap,{  
      result2<-NMF_heatmap(usedMeta_data(),summaryTable(),List(),N=dim(usedTable()$count)[1])
      result2$cutoff=input$cutoff
      result2  } )
    
      output$heatmap_cellstates <- renderPlot({
      pheatmap::pheatmap(border_color = NA,result2()[[colorHeatmapCells()]],display_numbers = result2()[["Mydata_raw_m"]],
                         fontsize = 12,fontsize_number = 15,
                         main=result2()[["cutoff"]],
                         cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                         color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
    })
      ## download csv
      output$downloadheatmap2<-downloadHandler(
        filename = function(){
          paste('Heatmap2-', colorHeatmapCells(),"-",Sys.Date(), '.csv', sep='')
        },
        content=function(heatmap2){
          write.csv(result2()[[colorHeatmapCells()]],heatmap2)
        }
      )
      ## download pdf
      output$downloadheatmap_plot2<-downloadHandler(
        filename = function(){
          paste('Heatmap2-', colorHeatmapCells(),"-",Sys.Date(), '.pdf', sep='')
        },
        content=function(heatmap2){
          width=dim(result2()[["log10FDR"]])[2]*5/25.4+max(nchar(rownames(result2()[["log10FDR"]])))*5/25.4+2
          height=dim(result2()[["log10FDR"]])[1]*5/25.4+max(nchar(colnames(result2()[["log10FDR"]])))*5/25.4+5
          pdf(heatmap2,width =width,height = height)
          pheatmap::pheatmap(border_color = NA,result2()[[colorHeatmapCells()]],display_numbers = result2()[["Mydata_raw_m"]],
                             fontsize = 12,fontsize_number = 15,cellwidth = 15,cellheight = 15,
                             main=result2()[["cutoff"]],
                             cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                             color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
          dev.off()
        } )
      
      
      #Heatmap3(Tab):
      result <- eventReactive(input$action_heatmap,{
        StateVsType(usedMeta_data(),N=dim(usedTable()$count)[1])
      })
      output$cellstates_types <- renderPlot({
      pheatmap::pheatmap(border_color = NA,result()[[colorHeatmapCells()]],display_numbers = result()$Mydata_raw_m,
                           fontsize = 12,fontsize_number = 15,
                           cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
      })
      ## download csv
      output$downloadheatmap3<-downloadHandler(
        filename = function(){
          paste('Heatmap3-', colorHeatmapCells(),"-",Sys.Date(), '.csv', sep='')
        },
        content=function(heatmap3){
          write.csv(result()[[colorHeatmapCells()]],heatmap3)
        }
      )
      ## download pdf
      output$downloadheatmap_plot3<-downloadHandler(
        filename = function(){
          paste('Heatmap3-', colorHeatmapCells(),"-",Sys.Date(), '.pdf', sep='')
        },
        content=function(heatmap3){
          width=dim(result()[["log10FDR"]])[2]*5/25.4+max(nchar(rownames(result()[["log10FDR"]])))*5/25.4+2
          height=dim(result()[["log10FDR"]])[1]*5/25.4+max(nchar(colnames(result()[["log10FDR"]])))*5/25.4+5
          pdf(heatmap3,width =width,height = height)
          pheatmap::pheatmap(border_color = NA,result()[[colorHeatmapCells()]],display_numbers = result()[["Mydata_raw_m"]],
                             fontsize = 12,fontsize_number = 15,cellwidth = 15,cellheight = 15,
                             cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                             color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
          dev.off()
          
        } )
      
      
      ##Heatmap-Tab4: cellStates
      cellStates <- eventReactive(input$action_heatmap,{
        testResult<-ListTest(List(),N=dim(usedTable()$count)[1])
        testResult$cutoff=input$cutoff
        testResult
      })
      output$cellStates_cellStates<- renderPlot({
        pheatmap::pheatmap(border_color = NA,cellStates()[[colorHeatmapCells()]],display_numbers = cellStates()$Mydata_raw_m,
                           fontsize = 12,fontsize_number = 15,
                           cluster_cols = T,cluster_rows = T,main=cellStates()[["cutoff"]],treeheight_row = 0, treeheight_col = 0,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
      })
      
      ## download csv
      output$downloadheatmap4<-downloadHandler(
        filename = function(){
          paste('Heatmap4-', colorHeatmapCells(),"-",Sys.Date(), '.csv', sep='')
        },
        content=function(heatmap1){
          write.csv(cellStates()[[colorHeatmapCells()]],heatmap1)
        }
      )
      ## download pdf
      output$downloadheatmap_plot4<-downloadHandler(
        filename = function(){
          paste('Heatmap4-', colorHeatmapCells(),"-",Sys.Date(), '.pdf', sep='')
        },
        content=function(heatmap4){
          width=dim(cellStates()[["log10FDR"]])[2]*5/25.4+max(nchar(rownames(cellStates()[["log10FDR"]])))*5/25.4+2
          height=dim(cellStates()[["log10FDR"]])[1]*5/25.4+max(nchar(colnames(cellStates()[["log10FDR"]])))*5/25.4+5
          pdf(heatmap4,width =width,height = height)
          pheatmap::pheatmap(border_color = NA,cellStates()[[colorHeatmapCells()]],display_numbers = cellStates()[["Mydata_raw_m"]],
                             fontsize = 12,fontsize_number = 15,cellwidth = 15,cellheight = 15,
                             main=cellStates()[["cutoff"]],
                             cluster_cols = T,cluster_rows = T,treeheight_row = 0, treeheight_col = 0,
                             color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
          dev.off()
          
        } )
      
      
      ##GO Tab
      #http://www.genome.jp/kegg/catalog/org_list.html
      GO <- reactive({
        FunctionE(input$cutoff,input$selected_clusterGO,input$Mart,input$kegg_species,input$go_species,GetGenesList(),BackgroundGenes())
      }) %>%
        bindCache(input$cutoff,input$selected_clusterGO,input$Mart,input$kegg_species,input$go_species,BackgroundGenes(),usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) %>%
        bindEvent(input$action_GO)
      
      
      output$Functional_enrichment <- renderPlot(height=400,{
        p_go<-Plot_enrichment(GO())
        print(p_go)
        
      })
      
      output$Functional_table <- renderDT({
        go_table<-GO()[,-c(1,5,7)]
        go_table$cluster<-paste0("Cluster:",go_table$cluster)
        go_table
      })
      
      
      output$downloadGO<-downloadHandler(
        filename = function(){
          paste('GOtable-', Sys.Date(), '.csv', sep='')
        },
        content=function(gotable){
          write.csv(GO(),gotable)
        }
      )
      
      
      plot_size<-reactive({
        width=length(unique(GO()$cluster))*1.5/25.4+max(nchar(GO()$Description))*1.5/24.5+6
        all_function<-GO() %>%
          group_by(cluster,class) %>%
          slice_max(n = 4,order_by = Count)%>%
          slice_head(n=4)
        height=dim(all_function)[1]*5/25.4+2
        size_lists<-list(width,height)
        names(size_lists)<-c("width","height")
        return(size_lists)
      })
      
      
      output$downloadGOPlot<-downloadHandler(
        filename = function(){
          paste('GOPlot-', Sys.Date(), '.pdf', sep='')
        },
        content=function(gotable){
        ggsave(Plot_enrichment(GO()),filename = gotable,width = plot_size()[["width"]],height =plot_size()[["height"]] )
        }
      )
      
      
      ###rrvgo Tab:
      rrvGO<- reactive({
        rrvgo_FunctionE(input$cutoff,input$selected_clusterrrvgo,input$subOntology,input$go_species_rrgvo,input$Mart_rrgvo,GetGenesList(),BackgroundGenes2())
      }) %>%
        bindCache(input$cutoff,input$selected_clusterrrvgo,input$subOntology,input$go_species_rrgvo,input$Mart_rrgvo,BackgroundGenes2(),usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) %>%
        bindEvent(input$action_rrvgo)
      
      ###rrvgo wordcloudPlot:
      output$SimplifyingGO_wordcloud<- renderPlot(height=400,{
        par(mar = rep(0, 4))
        wordcloudPlot( rrvGO()$reducedTerms,scale=c(3.5,0.25),use.r.layout=T)

      })
      
      output$rrvgo_plot1<-downloadHandler(
        filename = function(){
          paste('wordcloudPlot-', Sys.Date(), '.pdf', sep='')
        },
        content=function(wordcloud){
          pdf(wordcloud)
          wordcloudPlot(rrvGO()$reducedTerms)
          dev.off()
        }
      )
      
      
      
      ###rrvgo treemapPlot
      output$SimplifyingGO_treemap<- renderPlot({
        treemapPlot( rrvGO()$reducedTerms)
      })
      
      output$rrvgo_plot2<-downloadHandler(
        filename = function(){
          paste('TreemapPlot-', Sys.Date(), '.pdf', sep='')
        },
        content=function(treemap){
          pdf(treemap)
          treemapPlot(rrvGO()$reducedTerms)
          dev.off()
        }
      )
      
      
      ###rrvgo ScatterPlot
      output$SimplifyingGO_Scatter<- renderPlot({
        scatterPlot(labelSize = 5,rrvGO()$simMatrix, rrvGO()$reducedTerms)
      })
      
      
      output$rrvgo_plot3<-downloadHandler(
        filename = function(){
          paste('ScatterPlot-', Sys.Date(), '.pdf', sep='')
        },
        content=function(scatter){
          scplot<-scatterPlot(labelSize = 5,rrvGO()$simMatrix, rrvGO()$reducedTerms)
          ggsave(scplot,filename = scatter,width = 7,height =7)
        }
      )
      
      
      ###rrvgo Table
      output$Reduced_Terms <- renderDT({
        rrvGO()$reducedTerms
        
      })
      
      output$rrvgo_table<-downloadHandler(
        filename = function(){
          paste('rrvgo_table-', Sys.Date(), '.csv', sep='')
        },
        content=function(rrvgotable){
       write.csv(rrvGO()$reducedTerms,rrvgotable)
        }
      )
      
      
      
      ### Upset Plot
      m <- reactive({
        Up_set(input$selected_cluster_upset,input$cutoff,List())
      }) %>%
        bindCache(input$selected_cluster_upset,input$cutoff,List(),usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) %>%
        bindEvent(input$action_upset)
      
      cutoff<- eventReactive(input$action_upset, { 
      input$cutoff})
      
      output$UPsetPlot<- renderPlot({
      ComplexHeatmap::UpSet(m(),top_annotation = upset_top_annotation(m(), add_numbers = TRUE),column_title =cutoff(),
           right_annotation = upset_right_annotation(m(), add_numbers = TRUE),comb_order = order(comb_size(m())))
    })
    
      output$upset_plot<-downloadHandler(
        filename = function(){
          paste('UPsetPlot-', Sys.Date(), '.pdf', sep='')
        },
        content=function(UPsetPlot){
          upsetplot<-ComplexHeatmap::UpSet(m(),top_annotation = upset_top_annotation(m(), add_numbers = TRUE),column_title =cutoff(),
                                right_annotation = upset_right_annotation(m(), add_numbers = TRUE),comb_order = order(comb_size(m())))
          pdf(UPsetPlot,width =8 ,height = 6)
          draw(upsetplot)
          dev.off()
        }
      )
  
  
    ####DEA heatmap
      srt<-reactive({
        if(input$TestTable){
          return(usedTable()$srt)
        } else{ 
       processe_srt(usedTable()$Count_matrix)
          }
      })
      
      
      p <- reactive({
        DE_set(input$selected_clusterDE,input$cutoff,usedTable()$count,srt(),input$logfc,input$Pvalue_DE,List())
      }) %>%
        bindCache(input$selected_clusterDE,input$cutoff,input$logfc,input$Pvalue_DE,List(),usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) %>%
        bindEvent(input$action_DE)
      
      
      p_heatmap <- eventReactive(input$action_DE, { 
        PlotDEheatmap(input$selected_clusterDE,input$cutoff,p(),srt(),List())
      })
  
      
      output$DE_heatmap<- renderPlot({
        doheatmap<-DoHeatmap(p_heatmap()[["sub"]],features = p_heatmap()[["top20"]],slot = "c",group.by = "state",
                     disp.min = 0,disp.max = 4) +scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
        
        print(doheatmap)
    })
      
      
      output$GeneHeatmap<-downloadHandler(
        filename = function(){
          paste('GeneHeatmap-', Sys.Date(), '.pdf', sep='')
        },
        content=function(genesheatmap){
          doheatmap<-DoHeatmap(p_heatmap()[["sub"]],features = p_heatmap()[["top20"]],slot = "c",group.by = "state",
                               disp.min = 0,disp.max = 4) +scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
          ggsave(doheatmap,filename = genesheatmap,width = 7,height =7)
        }
      )
    
      
      output$DEGtable <- renderDT({
        p()[,-1]
      })
      
      output$DegTable<-downloadHandler(
        filename = function(){
          paste('DegTable-', Sys.Date(), '.csv', sep='')
        },
        content=function(Degstable){
          write.csv(p(),Degstable)
        }
      )
      
      
      
 
    ####DE function enrichment
      p_plot<- reactive({
        DEGO(p(),input$selected_clusterDE,input$cutoff,input$Mart_DE,input$kegg_species_DE,input$go_species_DE,input$logfc,input$Pvalue_DE, BackgroundGenes3())
      }) %>%
        bindCache(p(),input$selected_clusterDE,input$cutoff,input$Mart_DE,input$kegg_species_DE,input$go_species_DE,input$logfc,input$Pvalue_DE, BackgroundGenes3(),List(),usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) %>%
        bindEvent(input$action_DE)
      
      output$DEG_enrichment<- renderPlot({
        Plot_DE_enrichment(p_plot())
      })
      
      
      output$DegAnnotationPlot<-downloadHandler(
        filename = function(){
          paste('DegAnnotationPlot-', Sys.Date(), '.pdf', sep='')
        },
        content=function(DegAnoPlot){
          p_DegAnoPlot<- Plot_DE_enrichment(p_plot())
          ggsave(p_DegAnoPlot,filename = DegAnoPlot)
          
        }
      )
      
      
      output$DEG_enrichmentTable <- renderDT({
        p_plot()[,-c(1,5,7)]
      })

      
      output$DegAnnotationTable<-downloadHandler(
        filename = function(){
          paste('DegAnnotationTable-', Sys.Date(), '.csv', sep='')
        },
        content=function(DegAnoTable){
        write.csv(p_plot(),DegAnoTable)
      
        }
      )
      
      ## Cell state marker
      Marker <- reactive({
        Marker_set(input$selected_clusterMarker,input$cutoff,usedTable()$count,srt(),input$logfcMarker,input$Pvalue_Marker,List())
      }) %>%  #selected_cluster,cutoff,count,srt,logfc,Pvalue,List
        bindCache(input$selected_clusterMarker,input$cutoff,input$logfcMarker,input$Pvalue_Marker,List(),usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) %>%
        bindEvent(input$action_Marker)
      
      #selected_cluster,cutoff,Marker,srt,List
      Marker_plot <- eventReactive(input$action_Marker, { 
       PlotVolcano(input$selected_clusterMarker,input$cutoff, Marker(),srt(),List())
      })
      
      
      output$VolcanoPlot<- renderPlot({
        print(Marker_plot())
      })
      
      output$VolcanoPlot_Download<-downloadHandler(
        filename = function(){
          paste('VolcanoPlot-', Sys.Date(), '.pdf', sep='')
        },
        content=function(VolcanoPlot){
          Volco<-Marker_plot()
          ggsave(Volco,filename = VolcanoPlot,width = 5,height =4)
        }
      )
      
      output$MarkerTable <- renderDT({
        Marker()[,-1]
      })
      
      output$MarkerTable_Downloaded<-downloadHandler(
        filename = function(){
          paste('MarkerTable-', Sys.Date(), '.csv', sep='')
        },
        content=function(MarkerTable){
          write.csv(Marker(),MarkerTable)
          
        }
      )
      
      marker_plot<- reactive({ #Marker,selected_cluster,cutoff,Mart,kegg_species,go_species,logfc,Pvalue,background_genes
        MarkerGO(Marker(),input$selected_clusterMarker,input$cutoff,input$Mart_Marker,input$kegg_species_Marker,input$go_species_Marker,input$logfcMarker,input$Pvalue_Marker, BackgroundGenes4())
      }) %>%
        bindCache(Marker(),input$selected_clusterMarker,input$cutoff,input$Mart_Marker,input$kegg_species_Marker,input$go_species_Marker,input$logfcMarker,input$Pvalue_Marker, BackgroundGenes4(),List(),usedTable(),usedTable2(),input$minStateDeviation,input$minNoCells) %>%
        bindEvent(input$action_Marker)
      
      output$Marker_enrichment<- renderPlot({
        Plot_Marker_enrichment(marker_plot())
      })
      
      
      
      output$MarkerAnnotationPlot<-downloadHandler(
        filename = function(){
          paste('MarkersAnnotationPlot-', Sys.Date(), '.pdf', sep='')
        },
        content=function(MarkersAnoPlot){
          p_Markers<-  Plot_Marker_enrichment(marker_plot())
          ggsave(p_Markers,filename = MarkersAnoPlot)
          
        }
      )
      
      
      output$Marker_enrichmentTable <- renderDT({
        marker_plot()[,-c(1,5,7)]
      })
      
      
      output$MarkerAnnotationTable<-downloadHandler(
        filename = function(){
          paste('MarkersAnnotationTable-', Sys.Date(), '.csv', sep='')
        },
        content=function(MarkersAnoTable){
          write.csv(marker_plot(),MarkersAnoTable)
          
        }
      )
      
      
      
      
      # Get markov blanket:
      NodeGene<-reactive({
        Gene<-strsplit(input$NodeGene, ",\\s*")[[1]]
        Gene
      })
      
      Markov_blanket<-reactive({
        MB<-mb(usedMCMCGraph(),NodeGene())
        MB
      })%>% bindCache(input$NodeGene,usedMCMCGraph()) 
      
      
      
      Subgraph <- eventReactive(input$action_mb, { 
        subgraph<-subgraph(usedMCMCGraph(), unique(c(unlist(Markov_blanket()),NodeGene())))
        V(subgraph)$color<-"Navy"
        V(subgraph)$color[V(subgraph)$name%in%NodeGene()]<-"red"
        subgraph
      })
      
      output$MarkovBlanket <- renderPlot({
        plot(Subgraph(),layout=layout_with_kk(Subgraph()),edge.arrow.size = 0.3,vertex.label.color=V(Subgraph())$color,vertex.shape="none",vertex.label.font=0.4)
        
      })
      
      
      
      
      output$MarkovBlanket_plot<-downloadHandler(
        filename = function(){
          paste('MarkovBlanket-',NodeGene(),"-",Sys.Date(), '.pdf', sep='')
        },
        content=function(mbPlot){
          pdf(mbPlot)
          plot(Subgraph(),layout=layout_with_kk(Subgraph()),edge.arrow.size = 0.3,vertex.label.color=V( Subgraph())$color,vertex.shape="none",vertex.label.font=0.4)
          dev.off()
        } )
      
      
      
      
      
      }

# Run the application 
shinyApp(ui = ui, server = server)
