


#DTOutput

TableUI <- function(){
  tagList(
    tags$h3(paste0("Table for deviating MFIs"), style = "color: steelblue;"),
    textOutput("textsummary"),tags$br(),
    DTOutput("table") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadtable","Download as .csv")
  )}




### Input function
TableInput<- function(){
  tagList( 
    actionButton("action_table","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}



Table_cluster<-function(cutoff,data_path,minStateDeviation,minNoCells,stateDevAlpha){
  
  command=paste("python ./Produce_devStates.py", cutoff,data_path$devStates,data_path$trainDat,minStateDeviation,minNoCells,stateDevAlpha,sep = " ")
  system(command)
  Devstates<-read.csv(paste0('./',cutoff,"_",minStateDeviation,"_",minNoCells,"_",stateDevAlpha,'_devStates.csv'),colClasses = c("character"))
  file.remove(paste0('./',cutoff,"_",minStateDeviation,"_",minNoCells,"_",stateDevAlpha,'_devStates.csv'))
  Devstates$enrichment<-as.numeric(Devstates$enrichment)
  Devstates$cluster<-as.numeric(Devstates$cluster)
  Devstates$pval_corrected<-as.numeric(Devstates$pval_corrected)
  Devstates$No.Cells<-as.numeric(Devstates$No.Cells)
  Devstates<-Devstates[order(Devstates$cluster),]
  Devstates<-Devstates[,-1]
  
  return(Devstates)
  
}


