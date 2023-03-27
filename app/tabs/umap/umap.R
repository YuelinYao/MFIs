## UI ####
umapUI <- function(){
  tagList(
    tags$h3(paste0("UMAP plot"), style = "color: steelblue;"),
    plotOutput(outputId ="umap_plot", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadumap_plot","Download as .pdf")
  )}




### Input function
UMAPInput<- function(){
  tagList( 
    textInput("selected_umap", "Select cluster",value = "4"),
    actionButton("action_umap","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}





### function
cluster_umap<-function(selected_cluster,umap,srt,List){
  
  print("Plot umap")
  selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]
  selected_cluster<-paste0("cluster_C:",selected_cluster)
  selected_1<-List[[selected_cluster[1]]]
  print(length(selected_1)) 
  srt[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(srt))
  p<-DimPlot(object = srt,reduction = "umap",cells.highlight = selected_1, cols.highlight = "red", cols = "gray", order = TRUE)
  return(p)
  
}



