## general UI ####
UploadFilesUI <- function(){
  tagList(
    conditionalPanel(condition="input.tabs != 'about'",           #conditional panel: if the tab Tutorial is selected
                     br("Here we use scRNA-seq HCC dataset. \nTo upload your data, click the box:"),
                     switchInput(inputId = "TestTable",                            #switch button to upload your own data
                                 value = T, onStatus = "success"),
    conditionalPanel(condition = "!input.TestTable",              #conditional panel: if TestTable (example) is not selected, display the upload bar 
       tagList(
             tags$p(tags$strong("Upload your file:")),
              tags$div(
              br("Count_matrix.csv"),
               fileInput(inputId = "Count_matrix",label = NULL, multiple = FALSE,
                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
                                          ),style="display:inline-block;vertical-align:top;width:70%;"),
                tags$div(
                br("Meta_Data.csv"),
                 fileInput(inputId = "Meta_Data",label = NULL, multiple = FALSE,
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                          style="display:inline-block;vertical-align:top;width:70%;"
                                        ),
                tags$div(
                br("all_DTuples.csv"),
                fileInput(inputId = "topDeviatingHOIstates",label = NULL, multiple = FALSE,
                           accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv") ),
                                          style="display:inline-block;vertical-align:top;width:70%;"
                                        ),
                tags$div(
                br("trainingData_*.csv"),
                fileInput(inputId = "trainingData_matrix",label = NULL, multiple = FALSE,
                                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                          style="display:inline-block;vertical-align:top;width:70%;"
                                        ),
             
             tags$div(
               br("MCMCgraph_*.csv"),
               fileInput(inputId = "MCMCgraph",label = NULL, multiple = FALSE,
                         accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
               style="display:inline-block;vertical-align:top;width:70%;"
             ),
             
             
                                      )
                     )
    ),
    
    conditionalPanel(condition= "input.tabs != 'about' & input.tabs != 'mb'" ,
                     textInput("minStateDeviation", "Minimum enrichment factor:",value = 3),
                     textInput("minNoCells", "Minimum number of cells in each D-tuple:",value = 0),
                     textInput("stateDevAlpha", "Min. enrichment significance (corrected):",value = 0.05),
                     sliderInput("cutoff",
                                 "Dice distance:",
                                 min = 0,
                                 max = 1,
                                 value = 0.91))
    
    
  )
}

# Read count matrix:
read_data<-function(path_countmatrix){
  print("Read count matrix")
  ## Count matrix
  Count_matrix <- as.matrix(fread(path_countmatrix),rownames=1)
  ## Binary matrix
  count<-Count_matrix
  count[count>0]<-1
  count[count==0]<-0
  data<-list(count,Count_matrix)
  names(data)<-c("count","Count_matrix")
  return(data)
  
}

# Example background_genes:
background_genes<-read.table("./data/Background_genes.txt")
background<-paste0(background_genes$V1,collapse='\n')


## processed_seruat ftunction:
processe_srt<-function(Count_matrix){
  
  # Creat CreateSeuratObject
  Count_matrix<-t(Count_matrix) 
  
  ## preprocessing srt object
  print("Preprocessing srt object")
  srt<-CreateSeuratObject(counts =Count_matrix)
  srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
  srt <- ScaleData(srt)
  
  return(srt)
  
}






