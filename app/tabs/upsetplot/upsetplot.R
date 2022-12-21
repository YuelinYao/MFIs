### UI

UPsetUI <- function(){
  tagList(
    tags$h3(paste0("Upset Plot"), style = "color: steelblue;"),
    plotOutput(outputId ="UPsetPlot", width = "90%") %>% withSpinner(color="#4682B4")
  )}



Up_set<-function(selected_cluster,cutoff,List){
  
  selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]
  selected_cluster<-paste0("cluster_C:",selected_cluster)
  m=make_comb_mat(List[selected_cluster],mode = "intersect")  
  return(m)

}
