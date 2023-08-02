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
cluster_umap<-function(selected_cluster,umap,srt,List,Devstates){
  
  print("Plot umap")
  selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]
  names<-selected_cluster[1]
  selected_cluster<-paste0("cluster_C:",selected_cluster[1])
  selected_1<-List[[selected_cluster]]

  stats<-Devstates[Devstates$cluster==names,]

  Gene_pattern<-NULL
  for (r in 1:length(stats$genes)){
    
    Genes<-str_split(stats$genes[r], "_")[[1]]
    state<-stats$state[r]
    State<-as.numeric(strsplit(as.character(state),"")[[1]])
    State<-gsub("1","+",State)
    State<-gsub("0","-",State)
    Genes<-paste0(Genes,State,"\n")
    Gene_pattern<-c(Gene_pattern,Genes)
  }
  
  genes_cap<-paste0(head(Gene_pattern,5),collapse = "")
  sums<-paste0("(",length(selected_1)," cells)\n","(",dim(stats)[1]," d-tuples)")
  all<-paste0(genes_cap,sums)
  
  srt[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(srt))
  p<-DimPlot(object = srt,reduction = "umap",cells.highlight =selected_1,cols.highlight = "#3776ab",cols = "gray", order = TRUE)+NoAxes()+NoLegend()+ggtitle(names)+theme(plot.title = element_text(hjust = 0.5))
  p<-p+labs(caption =  all)+theme(plot.caption = element_text(hjust = 0.5))+ggtitle(names)
  return(p)
  
}



