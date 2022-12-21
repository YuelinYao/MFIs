#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# ====== import R libraries: ####
library(shiny)            #For shiny
library(shinythemes)      #For graphics of shiny interface
import::from(shinycssloaders, withSpinner) #For the spinner of load during the computing of the functions
reticulate::py_config()
library(shinyBS)          #For tooltips, popovers and alerts
library(shinyWidgets)     #For some shiny functions
library(gridExtra, verbose=FALSE)        #Grid display
library(RColorBrewer, verbose=FALSE)
library(reticulate)
library(ComplexHeatmap)
library(Seurat)
library(stringr)
library(pheatmap)
library(dplyr)
library(rrvgo)

# Some initial setup:
library(DT)
#library(org.Hs.eg.db)
library(clusterProfiler)
library(data.table)
library(biomaRt)
library(ggplot2)
#load("./data/Example.Rdata")
options(shiny.maxRequestSize = 1000*1024^2)
#source("./app/general/readData.R")
source("./app/general/general.R") # conditionalPanel
source("./app/tabs/about/about.R") # about
source("./app/tabs/heatmap/heatmap.R") # heatmap
source("./app/tabs/GO_plot/GO_plot.R") # GOplot
source("./app/tabs/upsetplot/upsetplot.R") 
source("./app/tabs/table/table.R") 
source("./app/tabs/DE/DE.R") 
source("./app/tabs/rrvgo/rrvgo.R") 



# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("spacelab"),
                
                # Application title
                titlePanel("",windowTitle = "MFIs project"),
                
                tags$head(
                  tags$style(
                    HTML("@import url('//fonts.googleapis.com/css?family=Righteous|');"),
                    HTML(".shiny-output-error-validation {
                                        color: red;
                          }")
                  )
                ),
                
                headerPanel(                              ### add logos 
                  
                  (img(src= "logo.png",height  = 120, width = 700)) ), ##class = "pull-left" 
                
                
                
                
                #Sidebar with a slider input for the cutoff
                sidebarLayout(
                  sidebarPanel(
            
                readTableUI(), 
            
            conditionalPanel(
                      condition="input.tabs == 'about'",
                      InformationUI()),
            
            conditionalPanel(condition= "input.tabs == 'upset'",
                             textInput("selected_cluster_upset", "Input cluster(s):",value = "5,18,6,11,19")),
            
            conditionalPanel(condition= "input.tabs == 'GO'",
                             textInput("selected_clusterGO", "Input cluster(s):",value = "5,18"),
                             selectInput("Mart", "Mart dataset:", choices=datasets_list, selected = "hsapiens_gene_ensembl", multiple = FALSE),
                             textInput("go_species", "GO OrgDb:",value = "org.Hs.eg.db"),
                             textInput("kegg_species", "KEGG organism:",value = "hsa"),
                             textAreaInput("background_genes", "Background genes (recommended): ",placeholder = "Just paste a list of genes (multiple-line gene list).",rows = 5),
                             
                             actionButton(inputId = "bg_Liver0",                                       #action button to display background genes
                                          label = NULL, icon = icon("tag"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                             
                             bsTooltip("bg_Liver0","Load background genes in HCC dataset.",placement = "bottom", trigger = "hover",
                                       options = NULL),
                             actionButton("action_GO","Submit",icon("paper-plane"), 
                                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                             
                             
                             ),
                            
            
            conditionalPanel(condition= "input.tabs == 'rrvgo'",
                             textInput("selected_clusterrrvgo", "Input a cluster:",value = "5"),
                             selectInput("Mart_rrgvo", "Mart dataset:", choices=datasets_list, selected = "hsapiens_gene_ensembl", multiple = FALSE),
                             textInput("go_species_rrgvo", "GO OrgDb:",value = "org.Hs.eg.db"),
                             radioButtons("subontology", "Select sub-ontology:",
                                                c("Biological Process" = "BP",
                                                  "Cellular Component" = "CC",
                                                  "Molecular Function" = "MF"),selected = "BP" ),

                             textAreaInput("background_genesrrgvo", "Background genes (recommended): ",placeholder = "Just paste a list of genes (multiple-line gene list).",rows = 5),

                             actionButton(inputId = "bg_Liver",                                       #action button to display background genes
                                          label = NULL, icon = icon("tag"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                             
                             bsTooltip("bg_Liver","Load background genes in HCC dataset.",placement = "bottom", trigger = "hover",
                                       options = NULL),
                              # submit button
                             actionButton("action_rrvgo","Submit",icon("paper-plane"), 
                                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                             ),
                    
            
            
            conditionalPanel(condition= "input.tabs == 'DE'",
                             textInput("selected_clusterDE", "Select cluster(s)",value = "6,19"),
                             selectInput("Mart_DE", "Mart dataset:", choices=datasets_list, selected = "hsapiens_gene_ensembl", multiple = FALSE),
                             textInput("go_species_DE", "GO OrgDb:",value = "org.Hs.eg.db"),
                             textInput("kegg_species_DE", "KEGG organism:",value = "hsa"),
                             textInput("logfc", "logFC:",value = 0.25),
                             sliderInput("Pvalue_DE",
                                         "Adjusted p value:",
                                         min = 0,
                                         max = 1,
                                         value = 0.05),
                             textAreaInput("background_genesDE", "Background genes (recommended): ",placeholder = "Just paste a list of genes (multiple-line gene list).",rows = 5),
                             
                             actionButton(inputId = "bg_Liver2",                                       #action button to display background genes
                                          label = NULL, icon = icon("tag"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                             
                             bsTooltip("bg_Liver2","Load background genes in HCC dataset.",placement = "bottom", trigger = "hover",
                                       options = NULL),
      
                             actionButton("action_DE","Submit",icon("paper-plane"), 
                                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                             ),
            
      
    

            ### submit button
            conditionalPanel(condition= "input.tabs == 'heatmap'",
                             actionButton("action_heatmap","Submit",icon("paper-plane"), 
                                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
            

            
            
            
            ### submit button
            conditionalPanel(condition= "input.tabs == 'upset'",
                             actionButton("action_upset","Submit",icon("paper-plane"), 
                                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
            
  
            ### table
            ### submit button
            conditionalPanel(condition= "input.tabs == 'table'",
                             actionButton("action_table","Submit",icon("paper-plane"), 
                                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
            
          
            
        ),
        
        
        mainPanel(
            tabsetPanel(type = "pills", id = 'tabs',
                
                        
                tabPanel("About", icon = icon("circle-info"), value = 'about', aboutUI()),
                
                #tabPanel("Tutorial", icon = icon("diamond"), value = 'about', aboutUI()), 
                
                tabPanel("Table", icon = icon("table"),value="table",TableUI()),

                tabPanel("Heatmap",icon = icon("map"), value="heatmap",
                         heatmapUI(), tags$br(),NMF_UI(),tags$br(),NMF_CelltypeUI()),
                
                tabPanel("GO & KEGG", icon = icon("chart-line"),value="GO",FunctionUI(),tags$br(),Function_tableUI()),
                
                
                tabPanel("Using rrvgo", icon = icon("chart-line"),value="rrvgo",
                         rrvgoUI(),tags$br(),rrvgo2UI(),tags$br(),
                         rrvgo3UI() ,tags$br(), rrvgo_tableUI()),
                        
                  
                tabPanel("Upset Plot",icon = icon("signal"),value="upset",UPsetUI()),
                
                tabPanel("DE analysis",icon = icon("random"),value="DE",DE_UI(),tags$br(),DE_TableUI(),
                         tags$br(),DEGO_UI(),
                         tags$br(),DEGOtable_UI())
          
                        
        
            
        ))   )

)






# Define server logic required to draw a histogram
server <- function(input, output,session) {
    
  # General:
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
  
  
  usedTable2 <- reactive({
    print(paste("Use HCC cancer data: ",input$TestTable))
    if(input$TestTable){
      devStates="./data/topDeviatingHOIstates.csv"
      trainDat="./data/trainingData_CL01_14698Cells_1000Genes.csv"
      #pcaCoords="./data/trainingData_CL01_14698Cells_1000Genes_PCAcoords.csv" 
      data_path<-list(devStates,trainDat)#, pcaCoords
      names(data_path)<-c("devStates","trainDat")#,"pcaCoords"
      return(data_path)
    }
    
    else{
      file_devStates<-input$topDeviatingHOIstates
      file_trainDat<-input$trainingData_matrix
      #file_pcaCoords<-input$trainingData_PCA
      data_path<-list(file_devStates$datapath,file_trainDat$datapath)#,file_pcaCoords$datapath
      names(data_path)<-c("devStates","trainDat")#,"pcaCoords"
      print(data_path)
      return(data_path)
    }
  })  
  
  
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
  
  # load background genes
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
  
  
  # About:
    Text1 <- about()
    aboutOutput(output,Text1)

    # setup:
    #summary Table
    summaryTable<-reactive({
      Table_cluster(input$cutoff,usedTable2())
    }) %>% bindCache(input$cutoff,usedTable2())
    
    # Get cell list
    List<-reactive({
      GetCellList(input$cutoff,usedTable()$count,summaryTable())
    })%>% bindCache(input$cutoff,usedTable(),usedTable2()) 
    
    
    # Get gene list:
    GetGenesList<-reactive({
      GetGenes(input$cutoff,summaryTable())
    })%>% bindCache(input$cutoff,usedTable(),usedTable2()) 
    
    
    
    ## Table:  
    showText<-eventReactive(input$action_table, {
      paste(dim(summaryTable())[1]," deviating MFIs in total, ",length(unique(summaryTable()$cluster)), "clusters." )
    })
    
    output$textsummary <- renderText({
      showText()
    })
    
    show_Table <- eventReactive(input$action_table, { 
      summaryTable()
    })
    
    output$table <- renderDT({
      show_Table()
    })
    
  

    # Heatmap1(Tab):
    result1 <- eventReactive(input$action_heatmap,{
      result1<-heatmap(input$cutoff,usedMeta_data(),summaryTable(), List())
      result1$cutoff=input$cutoff
      result1
    } )
    
    output$heatmap_celltypes <- renderPlot({
    pheatmap::pheatmap(border_color = NA,result1()[["Data_mtrix_log"]],display_numbers = result1()[["Mydata_raw_m"]],
                           fontsize = 12,fontsize_number = 15,fontsize_row = 15,main=result1()[["cutoff"]],
                           cluster_cols = T,cluster_rows = T,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
    })
    
    
  
    # Heatmap2(Tab):
    #main = input$cutoff,
    result2 <- eventReactive(input$action_heatmap,{  
      result2<-NMF_heatmap(input$cutoff,usedMeta_data(),summaryTable(),List())
      result2$cutoff=input$cutoff
      result2  } )
    
      output$heatmap_cellstates <- renderPlot({
      pheatmap::pheatmap(border_color = NA,result2()[["Data_mtrix_log"]],display_numbers = result2()[["Mydata_raw_m"]],
                         fontsize = 12,fontsize_number = 15,fontsize_row = 15,main=result2()[["cutoff"]],
                         cluster_cols = T,cluster_rows = T,
                         color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
    })
    
      #Heatmap3(Tab):
      result <- eventReactive(input$action_heatmap,{
        StateVsType(usedMeta_data())
      })
      output$cellstates_types <- renderPlot({
      pheatmap::pheatmap(border_color = NA,result()$Data_mtrix_log,display_numbers = result()$Mydata_raw_m,
                           fontsize = 12,fontsize_number = 15,fontsize_row = 15,
                           cluster_cols = T,cluster_rows = T,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
      })
    
      
      
      #  #http://www.genome.jp/kegg/catalog/org_list.html
      GO <- reactive({
        FunctionE(input$cutoff,input$selected_clusterGO,input$Mart,input$kegg_species,input$go_species,GetGenesList(),BackgroundGenes())
      }) %>%
        bindCache(input$cutoff,input$selected_clusterGO,input$Mart,input$kegg_species,input$go_species,BackgroundGenes(),usedTable(),usedTable2()) %>%
        bindEvent(input$action_GO)
      
      
      output$Functional_enrichment <- renderPlot(height=420,{
        p_go<-Plot_enrichment(GO())
        p_go
        
      })
      
      output$Functional_table <- renderDT({
        GO()[,-c(1,5,7)]
      })
      
      
      ###rrvgo function(cutoff,selected_cluster,subclass,go_species,Mart)
      # 
      rrvGO<- reactive({
        rrvgo_FunctionE(input$cutoff,input$selected_clusterrrvgo,input$subontology,input$go_species_rrgvo,input$Mart_rrgvo,GetGenesList(),BackgroundGenes2())
      }) %>%
        bindCache(input$cutoff,input$selected_clusterrrvgo,input$subontology,input$go_species_rrgvo,input$Mart_rrgvo,BackgroundGenes2(),usedTable(),usedTable2()) %>%
        bindEvent(input$action_rrvgo)
      
      output$SimplifyingGOBP_wordcloud<- renderPlot(height=440,{
       wordcloudPlot( rrvGO()$reducedTerms)

      })
      
      output$SimplifyingGOBP_treemap<- renderPlot({
        treemapPlot( rrvGO()$reducedTerms)
      })
      
      output$SimplifyingGOBP_Scatter<- renderPlot({
        scatterPlot(labelSize = 5,rrvGO()$simMatrix, rrvGO()$reducedTerms)
      })
      
      output$Reduced_Terms <- renderDT({
        rrvGO()$reducedTerms
        
      })
      
      
        
      ### upset
      m <- reactive({
        Up_set(input$selected_cluster_upset,input$cutoff,List())
      }) %>%
        bindCache(input$selected_cluster_upset,input$cutoff,List(),usedTable(),usedTable2()) %>%
        bindEvent(input$action_upset)
      
      cutoff<- eventReactive(input$action_upset, { 
      input$cutoff})
      
      output$UPsetPlot<- renderPlot({
      ComplexHeatmap::UpSet(m(),top_annotation = upset_top_annotation(m(), add_numbers = TRUE),column_title =cutoff(),
           right_annotation = upset_right_annotation(m(), add_numbers = TRUE),comb_order = order(comb_size(m())))
    })
    
    
  
      
      
      
    ####DEA heatmap
    
      srt<-reactive({
        if(input$TestTable){
          return(usedTable()$srt)
        } else{ 
       processed_srt(usedTable()$Count_matrix)
          }
      })%>%bindEvent(input$action_DE)
      
      
      p <- reactive({
        DE_set(input$selected_clusterDE,input$cutoff,usedTable()$count,srt(),input$logfc,input$Pvalue_DE,List())
      }) %>%
        bindCache(input$selected_clusterDE,input$cutoff,input$logfc,input$Pvalue_DE,List(),usedTable(),usedTable2()) %>%
        bindEvent(input$action_DE)
      
      
      p_heatmap <- eventReactive(input$action_DE, { 
        PlotDEheatmap(input$selected_clusterDE,input$cutoff,p(),srt(),List())
      })
  
      output$DE_heatmap<- renderPlot({
        doheatmap<-DoHeatmap(p_heatmap()[["sub"]],features = p_heatmap()[["top20"]],slot = "c",group.by = "state",
                     disp.min = 0,disp.max = 4) +scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
        
        print(doheatmap)
    })
      
      output$DEGtable <- renderDT({
        p()[,-1]
      })
      
 
    ####DE function enrichment
      p_plot<- reactive({
        DEGO(p(),input$selected_clusterDE,input$cutoff,input$Mart_DE,input$kegg_species_DE,input$go_species_DE,input$logfc,input$Pvalue_DE, BackgroundGenes3())
      }) %>%
        bindCache(p(),input$selected_clusterDE,input$cutoff,input$logfc,input$Pvalue_DE, BackgroundGenes3(),List(),usedTable(),usedTable2()) %>%
        bindEvent(input$action_DE)
      
      output$DEG_enrichment<- renderPlot({
        Plot_DE_enrichment( p_plot())
      })
      
      output$DEG_enrichmentTable <- renderDT({
        p_plot()[,-c(1,5,7)]
      })

      }

# Run the application 
shinyApp(ui = ui, server = server)
