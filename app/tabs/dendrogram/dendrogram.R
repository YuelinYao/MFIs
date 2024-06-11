## UI ####
dendrogramUI <- function(){
  tagList(
    tags$h3(paste0("Full dendrogram"), style = "color: steelblue;"),
    plotOutput(outputId ="Full_dendrogram", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadFull_dendrogram","Download as .pdf"),
    
    tags$h3(paste0("Sub dendrogram"), style = "color: steelblue;"),
    plotOutput(outputId ="Sub_dendrogram", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadSub_dendrogram","Download as .pdf"),
  )}




### Input function
DendrogramInput<- function(){
  tagList( 
    textInput("selected_dendrogram", "Select a Stator state:",value = "45"),
    actionButton("action_dendrogram","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}


### function

dice = function(mat) {
  present = !is.na(mat) & mat > 0
  as.matrix(arules::dissimilarity(t(present), method = 'dice'))
}


GetSubend<-function(dend,cluster){
  dend_list <- lapply(unique(cluster), function(cluster.id) {
    find_dendrogram(dend, names(which(cluster == cluster.id)))
  })
  class(dend_list) <- "dendlist"
  dend_list
}



hcl_construct<-function(binRep,Devstates){
  print("construct tree")
  Devstates$mapping<-paste0(Devstates$genes,":",Devstates$state)
  binRep<-binRep[Devstates$mapping,]
  d<-dice(t(binRep))
  hcl <- hclust(as.dist(d),method='average')
  return(hcl)
  
}

hcl_cluster<-function(hcl,cutoff,Devstates){
  print("perform hcl")
  print(cutoff)
  Devstates$mapping<-paste0(Devstates$genes,":",Devstates$state)
  cluster<-cutree(hcl,h = cutoff)
  print(table(names(cluster)==Devstates$mapping))
  return(cluster)  
  
}

plot_full<-function(hcl,cutoff,cluster){
  
  print("plot full")
  avg_dend_obj <- as.dendrogram(hcl)
  avg_col_dend <- color_branches(avg_dend_obj, h = cutoff, groupLabels = unique(cluster[order.dendrogram(avg_dend_obj)]))
  return(avg_col_dend)  
}


plot_sub<-function(hcl,cluster,selected_cluster){
  
  print("plot sub")
  selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]
  dend<-as.dendrogram(hcl)
  dend_list<-GetSubend(dend,cluster)
  names(dend_list)<-c(1:length(dend_list))
  sub<-dend_list[[selected_cluster]]
  return(sub)  
}


