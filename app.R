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
PYTHON_DEPENDENCIES = c('numpy','pandas','scipy')

# ------------------ App virtualenv setup (Do not edit) ------------------- #
# VIRTUALENV_NAME and PYTHON_PATH are definded in .Rprofile
virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
python_path = Sys.getenv('PYTHON_PATH')


# Create virtual env and install dependencies
reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES, ignore_installed=TRUE)
reticulate::use_virtualenv(virtualenv_dir, required = T)

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
#library(org.Hs.eg.db)
#library(org.Mm.eg.db)
options(shiny.maxRequestSize = 5000*1024^2)

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
                             radioButtons("subOntology", "Select sub-ontology:",
                                                c("Biological Process" = "BP",
                                                  "Cellular Component" = "CC",
                                                  "Molecular Function" = "MF"),selected = "BP" ),
                             textAreaInput("background_genesrrgvo", "Background genes (recommended): ",placeholder = "Just paste a list of genes (multiple-line gene list).",rows = 5),
                             actionButton(inputId = "bg_Liver",                                       #action button to display background genes
                                          label = NULL, icon = icon("tag"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                             bsTooltip("bg_Liver","Load background genes in HCC dataset.",placement = "bottom", trigger = "hover",
                                       options = NULL),
                             actionButton("action_rrvgo","Submit",icon("paper-plane"),       # submit button
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
                tabPanel("Table", icon = icon("table"),value="table",TableUI()),
                tabPanel("Heatmap",icon = icon("map"), value="heatmap",
                         heatmapUI(),downloadButton("downloadheatmap1","Download as .csv"),downloadButton("downloadheatmap_plot1","Download as .pdf"),
                         NMF_UI(),downloadButton("downloadheatmap2","Download as .csv"), downloadButton("downloadheatmap_plot2","Download as .pdf"),tags$br(),
                         NMF_CelltypeUI(),downloadButton("downloadheatmap3","Download as .csv"), downloadButton("downloadheatmap_plot3","Download as .pdf")),
                tabPanel("GO & KEGG", icon = icon("chart-line"),value="GO",FunctionUI(),downloadButton("downloadGOPlot","Download as .pdf"),tags$br(),
                         Function_tableUI(),downloadButton("downloadGO","Download as .csv")),
                tabPanel("Using rrvgo", icon = icon("chart-line"),value="rrvgo",
                         rrvgoUI(),downloadButton("rrvgo_plot1","Download as .pdf"),tags$br(),
                         rrvgo2UI(),downloadButton("rrvgo_plot2","Download as .pdf"),tags$br(),
                         rrvgo3UI() ,downloadButton("rrvgo_plot3","Download as .pdf"),tags$br(), 
                         rrvgo_tableUI(),downloadButton("rrvgo_table","Download as .csv")),
                tabPanel("Upset Plot",icon = icon("signal"),value="upset",UPsetUI(),
                         downloadButton("upset_plot","Download as .pdf")),
                tabPanel("DE analysis",icon = icon("random"),value="DE",DE_UI(),downloadButton("GeneHeatmap","Download as .pdf"),tags$br(),
                         DE_TableUI(),downloadButton("DegTable","Download as .csv"),tags$br(),
                         DEGO_UI(),downloadButton("DegAnnotationPlot","Download as .pdf"),tags$br(),
                         DEGOtable_UI(),downloadButton("DegAnnotationTable","Download as .csv"),tags$br())
        ))   )

)



# Define server logic required to draw a histogram
server <- function(input, output,session) {

  # Load example data or read uploaded data
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
      data_path<-list(devStates,trainDat)
      names(data_path)<-c("devStates","trainDat")
      return(data_path)
    }
    
    else{
      file_devStates<-input$topDeviatingHOIstates
      file_trainDat<-input$trainingData_matrix
      data_path<-list(file_devStates$datapath,file_trainDat$datapath)
      names(data_path)<-c("devStates","trainDat")
      #print(data_path)
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
      table<-summaryTable()
      table$cluster<-paste0("Cluster:",table$cluster)
      colnames(table)<-c("Genes","D-tuple","Dev","Pval","Cluster")
      table
    })
    
    output$table <- renderDT({
      show_Table()
    })
    
    output$downloadtable<-downloadHandler(
      filename=function(){
        paste0("State_Table ",Sys.Date(), '.csv')
      },
      content = function(file){
        write.csv(show_Table(),file)
      }
    )

    # Heatmap1(Tab):
    result1 <- eventReactive(input$action_heatmap,{
      result1<-heatmap(input$cutoff,usedMeta_data(),summaryTable(), List())
      result1$cutoff=input$cutoff
      result1
    } )
    
    output$heatmap_celltypes <- renderPlot({
    pheatmap::pheatmap(border_color = NA,result1()[["Data_mtrix_log"]],display_numbers = result1()[["Mydata_raw_m"]],
                           fontsize = 12,fontsize_number = 15,
                            main=result1()[["cutoff"]],
                           cluster_cols = T,cluster_rows = T,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
    })
    
    
    output$downloadheatmap1<-downloadHandler(
      filename = function(){
        paste('Heatmap1-', Sys.Date(), '.csv', sep='')
      },
      content=function(heatmap1){
        write.csv(result1()[["Data_mtrix_log"]],heatmap1)
      }
   )

    output$downloadheatmap_plot1<-downloadHandler(
      filename = function(){
        paste('Heatmap1-', Sys.Date(), '.pdf', sep='')
      },
      content=function(heatmap1){
        width=dim(result1()[["Data_mtrix_log"]])[2]*5/25.4+max(nchar(rownames(result1()[["Data_mtrix_log"]])))*5/25.4+2
        height=dim(result1()[["Data_mtrix_log"]])[1]*5/25.4+max(nchar(colnames(result1()[["Data_mtrix_log"]])))*5/25.4+5
        pdf(heatmap1,width =width,height = height)
        pheatmap::pheatmap(border_color = NA,result1()[["Data_mtrix_log"]],display_numbers = result1()[["Mydata_raw_m"]],
                           fontsize = 12,fontsize_number = 15,cellwidth = 15,cellheight = 15,
                           main=result1()[["cutoff"]],
                           cluster_cols = T,cluster_rows = T,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
        dev.off()

     } )
    
    
  
    # Heatmap2(Tab):
    #main = input$cutoff,
    result2 <- eventReactive(input$action_heatmap,{  
      result2<-NMF_heatmap(input$cutoff,usedMeta_data(),summaryTable(),List())
      result2$cutoff=input$cutoff
      result2  } )
    
      output$heatmap_cellstates <- renderPlot({
      pheatmap::pheatmap(border_color = NA,result2()[["Data_mtrix_log"]],display_numbers = result2()[["Mydata_raw_m"]],
                         fontsize = 12,fontsize_number = 15,
                         main=result2()[["cutoff"]],
                         cluster_cols = T,cluster_rows = T,
                         color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
    })
    
      
      output$downloadheatmap2<-downloadHandler(
        filename = function(){
          paste('Heatmap2-', Sys.Date(), '.csv', sep='')
        },
        content=function(heatmap2){
          write.csv(result2()[["Data_mtrix_log"]],heatmap2)
        }
      )
      
      output$downloadheatmap_plot2<-downloadHandler(
        filename = function(){
          paste('Heatmap2-', Sys.Date(), '.pdf', sep='')
        },
        content=function(heatmap2){
          width=dim(result2()[["Data_mtrix_log"]])[2]*5/25.4+max(nchar(rownames(result2()[["Data_mtrix_log"]])))*5/25.4+2
          height=dim(result2()[["Data_mtrix_log"]])[1]*5/25.4+max(nchar(colnames(result2()[["Data_mtrix_log"]])))*5/25.4+5
          pdf(heatmap2,width =width,height = height)
          pheatmap::pheatmap(border_color = NA,result2()[["Data_mtrix_log"]],display_numbers = result2()[["Mydata_raw_m"]],
                             fontsize = 12,fontsize_number = 15,cellwidth = 15,cellheight = 15,
                             main=result2()[["cutoff"]],
                             cluster_cols = T,cluster_rows = T,
                             color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
          dev.off()
        } )
      
      
      #Heatmap3(Tab):
      result <- eventReactive(input$action_heatmap,{
        StateVsType(usedMeta_data())
      })
      output$cellstates_types <- renderPlot({
      pheatmap::pheatmap(border_color = NA,result()$Data_mtrix_log,display_numbers = result()$Mydata_raw_m,
                           fontsize = 12,fontsize_number = 15,
                           cluster_cols = T,cluster_rows = T,
                           color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
      })
    
      
      output$downloadheatmap3<-downloadHandler(
        filename = function(){
          paste('Heatmap3-', Sys.Date(), '.csv', sep='')
        },
        content=function(heatmap3){
          write.csv(result()[["Data_mtrix_log"]],heatmap3)
        }
      )
      
      output$downloadheatmap_plot3<-downloadHandler(
        filename = function(){
          paste('Heatmap3-', Sys.Date(), '.pdf', sep='')
        },
        content=function(heatmap3){
          width=dim(result()[["Data_mtrix_log"]])[2]*5/25.4+max(nchar(rownames(result()[["Data_mtrix_log"]])))*5/25.4+2
          height=dim(result()[["Data_mtrix_log"]])[1]*5/25.4+max(nchar(colnames(result()[["Data_mtrix_log"]])))*5/25.4+5
          pdf(heatmap3,width =width,height = height)
          pheatmap::pheatmap(border_color = NA,result()[["Data_mtrix_log"]],display_numbers = result()[["Mydata_raw_m"]],
                             fontsize = 12,fontsize_number = 15,cellwidth = 15,cellheight = 15,
                             cluster_cols = T,cluster_rows = T,
                             color = colorRampPalette(c("white","firebrick3"))(10),breaks = seq(0, 10, by = 1))
          dev.off()
          
        } )
      
      
      
      ##GO Tab
      #http://www.genome.jp/kegg/catalog/org_list.html
      GO <- reactive({
        FunctionE(input$cutoff,input$selected_clusterGO,input$Mart,input$kegg_species,input$go_species,GetGenesList(),BackgroundGenes())
      }) %>%
        bindCache(input$cutoff,input$selected_clusterGO,input$Mart,input$kegg_species,input$go_species,BackgroundGenes(),usedTable(),usedTable2()) %>%
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
        bindCache(input$cutoff,input$selected_clusterrrvgo,input$subOntology,input$go_species_rrgvo,input$Mart_rrgvo,BackgroundGenes2(),usedTable(),usedTable2()) %>%
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
        bindCache(input$selected_cluster_upset,input$cutoff,List(),usedTable(),usedTable2()) %>%
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
        bindCache(p(),input$selected_clusterDE,input$cutoff,input$Mart_DE,input$kegg_species_DE,input$go_species_DE,input$logfc,input$Pvalue_DE, BackgroundGenes3(),List(),usedTable(),usedTable2()) %>%
        bindEvent(input$action_DE)
      
      output$DEG_enrichment<- renderPlot({
        Plot_DE_enrichment( p_plot())
      })
      
      
      output$DegAnnotationPlot<-downloadHandler(
        filename = function(){
          paste('DegAnnotationPlot-', Sys.Date(), '.pdf', sep='')
        },
        content=function(DegAnoPlot){
          p_DegAnoPlot<- Plot_DE_enrichment( p_plot())
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
      }

# Run the application 
shinyApp(ui = ui, server = server)
