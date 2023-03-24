### UI
UPsetUI <- function(){
  tagList(
    tags$h3(paste0("Upset Plot"), style = "color: steelblue;"),
    plotOutput(outputId ="UPsetPlot", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("upset_plot","Download as .pdf")
  )}


### Input function
UPsetInput<- function(){
  tagList( 
    textInput("selected_cluster_upset", "Input cluster(s):",value = "1,3,8"),
    actionButton("action_upset","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}




Up_set<-function(selected_cluster,cutoff,List){
  selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]
  selected_cluster<-paste0("cluster_C:",selected_cluster)
  m=make_comb_mat(List[selected_cluster],mode = "intersect")  
  return(m)

}
