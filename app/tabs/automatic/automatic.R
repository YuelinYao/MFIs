## UI
AllMarker_TableUI <- function(){
  tagList(
    tags$h3(paste0("Find markers for all cell states"), style = "color: steelblue;"),
    DTOutput(outputId ="AllMarkerTable", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("AllMarkerTable_Downloaded","Download as .csv")
  )}




AllMarkerAnnotationTable_UI<- function(){
  tagList(
    tags$h3(paste0("Automatic Annotation"), style = "color: steelblue;"),
    DTOutput(outputId ="AllMarker_AutomaticAnnotation", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("AllMarker_AutomaticAnnotationTable","Download as .csv")
  )}


Automatic_Heatmap_UI <- function(){
  tagList(
    tags$h3("Automatic Annotation Heatmap", style = "color: steelblue;"),
    plotlyOutput(outputId ="Automatic_Annotation_Heatmap", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadAutomatic_Annotation_Heatmap","Download as .csv")#, downloadButton("downloadAutomatic_Annotation_Heatmap_plot","Download as .pdf")
  )
}


### Input function
AllMarkerInput<- function(){
  tagList( 

    textInput("logfcAllMarker", "logFC:",value = 0.25),
    sliderInput("Pvalue_AllMarker",
                "Adjusted p value:",
                min = 0,
                max = 1,
                value = 0.05),
    fileInput(inputId = "Auto_annot",label = "Upload annotation table:", multiple = FALSE,
              accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    radioButtons(inputId = "regulation", "The type of annotations:",
                 choices = c("Up_regulated genes" = "Up_regulated genes", "Down_regulated genes" = "Down_regulated genes"), 
                 selected = "Up_regulated genes", inline = TRUE),
    actionButton("action_AllMarkers","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}





MarkerAll_set<-function(srt,logfc,Pvalue,List){
  
  print("Find All marker genes")
  All_Markers<-NULL
  
  for (i in names(List)){
  print(i)
  selected_cells<-List[[i]]
  #The number of cells: print(length(selected_cells)) #The number of cells

  Marker<-FindMarkers(object = srt,ident.1 = selected_cells,logfc.threshold = logfc)
  print("done")
  Marker<-Marker[Marker$p_val_adj<Pvalue,]
  Marker$gene<-rownames(Marker)
  rownames(Marker)<-NULL
  Marker$cluster[Marker$avg_log2FC>0]<-paste0("Up_regulated genes")
  Marker$cluster[Marker$avg_log2FC<0]<-paste0("Down_regulated genes")
  print(table(Marker$cluster))
  Marker$State<-i
  All_Markers<-rbind(All_Markers,Marker)
  
  }
  return(All_Markers)
  
}


Markers_in_Annot<-function(MarkerAll,annotation,regulation){
  markers_in_annotation<-NULL
  for (c in colnames(annotation)){
    print(c)
    genes<-annotation[,c]
    genes<-genes[!genes%in%""]
    genes<-genes[!is.na(genes)]
    markers_subset<-MarkerAll[MarkerAll$gene%in%genes,]
    markers_subset<-markers_subset[markers_subset$cluster==regulation,]
    markers_subset<-markers_subset[order(markers_subset$p_val_adj),] 
    if(dim(markers_subset)[1]>0){
    
    markers_subset$Annotation<-NA
    markers_subset$Annotation<-c
    markers_in_annotation<-rbind(markers_in_annotation,markers_subset)}
    
  }
  return(markers_in_annotation)
}


number_all<-function(annotation){
  
  number_all<-NULL
  for (c in colnames(annotation)){
    genes<-annotation[,c]
    genes<-genes[!genes%in%""]
    genes<-genes[!is.na(genes)]
    number<-length(genes)
    names(number)<-c
    number_all<-c(number_all,number)

}
return(number_all)
}

summary_df<-function(Markers_in_Annot,number_all){
  
  df<-table(Markers_in_Annot$Annotation,Markers_in_Annot$State)
  colnames(df)<-gsub("cluster_C:","",colnames(df))
  number_all<-number_all[rownames(df)]
  rownames(df)<-paste0(rownames(df)," (",number_all,")")
  class(df)<-"matrix"
  return(df)
}

percent_summaryDF<-function(df,number_all){
  for (r in rownames(df)){
    df[r,]<-df[r,]/number_all[r]
  }
  class(df)<-"matrix"
 
  return(df)
}



get_genelist_annotation<-function(annotation,number_all){
  
  df<-table(annotation$Annotation,annotation$State)
  colnames(df)<-gsub("cluster_C:","",colnames(df))
  results<-array(data = NA,dim = dim(df),dimnames = dimnames(df))
  for (r in rownames(results)){
    
    for (c in colnames(results)){
      r1<-annotation[annotation$State==paste0("cluster_C:",c)&annotation$Annotation==r,]
      r1<-r1[order(r1$p_val_adj),]
      results[r,c]=""
      if( dim(r1)[1]>0){
        gene_list<-paste0(r1$gene,collapse = ", ")
        results[r,c]<-gene_list}
    }
  }
  
  number_all<-number_all[rownames(results)]
  rownames(results)<-paste0(rownames(results)," (",number_all,")")
  return(results)
}


