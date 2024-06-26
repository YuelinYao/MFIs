## UI ####
heatmapGenesUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap: Gene set"), style = "color: steelblue;"),
    plotOutput(outputId ="heatmap_GeneSet", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap_genes","Download as .csv"),
    downloadButton("downloadheatmapgenes_plot","Download as .pdf"),
    downloadButton("downloadheatmap_genes_pvalue","Download raw pvalue as .csv"),
  )}


stateGenesUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap: State genes"), style = "color: steelblue;"),
    plotOutput(outputId ="heatmapStateGenes", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap_genes2","Download as .csv"),
    downloadButton("downloadheatmapgenes_plot2","Download as .pdf"),
    downloadButton("downloadheatmap_genes2_pvalue","Download raw pvalue as .csv"),
  )}



iheatmapGenesUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap: Gene set"), style = "color: steelblue;"),
    plotlyOutput(outputId ="iheatmap_GeneSet", width = "90%") %>% withSpinner(color="#4682B4"),

    tags$h3(paste0("Heatmap: State genes"), style = "color: steelblue;"),
    plotlyOutput(outputId ="iheatmapStateGenes", width = "90%") %>% withSpinner(color="#4682B4"),
    
  )}


### Input function
HeatmapGenesInput<- function(){
  tagList( 
    radioButtons(inputId = "TestGenes", "Test:",
                choices = c("Enrichment" = "Over_representation", "Depletion" = "Under_representation", "Two.side" = "Fisher"), 
                 selected = "Over_representation", inline = TRUE),
    #radioButtons(inputId = "colorHeatmapGene", "Colored by:",
     #            choices = c("-log10FDR" = "log10FDR", "Fold" = "Fold"),
      #           selected = "log10FDR", inline = TRUE),
    selectInput("CellStateGenes", "Cell state gene set:", choices=c("Cancer cell state (Barkley, D.et.al., 2022)"="CancerState",
                                                                    "Cell cycle state (Tirosh et al, 2015)"="CellCycleState",
                                                                    "Upload state gene set" ="UploadState"), 
                selected = "CancerState", multiple = FALSE),
    textInput(inputId="NO.background", "The number of background genes", value = 25678),
    #N=25678
    conditionalPanel(
      condition = "input.CellStateGenes == 'UploadState'",
      fileInput(inputId = "uploadCellStateGenes",label = NULL, multiple = FALSE,
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
      )),
    actionButton("action_heatmap_genes","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}




makeGeneSetList<-function(csvfile){
  Cell_states<-read.csv(csvfile,row.names = 1)
  List=NULL
  for (i in 1:dim(Cell_states)[2]){
    
    cellstates<-Cell_states[,i]
    cellstates<-cellstates[!is.na(cellstates)]
    cellstates<-cellstates[!cellstates==""]
    cellstates<-list(cellstates)
    names(cellstates)<-colnames(Cell_states)[i]
    List<-c(List,cellstates)
    
  }
  return(List)
}

### function
#N=25678 the number of genes in whole human genome
heatmapGenes <- function(RefSet,Genelist,N=25678,test) { 

  RefSet<-makeGeneSetList(RefSet)
  #print(head(RefSet))
  r=length(names(RefSet))
  col=length(names(Genelist))
  
  Fold<-array(data=NA,dim = c(r,col))
  Data_mtrix<-array(data=NA,dim = c(r,col))
  
  colnames(Data_mtrix)<-names(Genelist)
  rownames(Data_mtrix)<-names(RefSet)
  
  dimnames(Fold)<-dimnames(Data_mtrix)
  
  Overlap<-array(data=NA,dim = c(r,col))
  dimnames(Overlap)<-dimnames(Data_mtrix)
  
  for (i in names(Genelist)){
    #print(i)
    Genes<-Genelist[[i]]
    P_set=rep(NA,length(names(RefSet)))
    names(P_set)<-names(RefSet)
    
    enrichment_set<-rep(NA,length(names(RefSet)))
    names(enrichment_set)<-names(RefSet)
    
    for (ct in names(RefSet)){
      #print(ct)
      Gene_types<-RefSet[[ct]]
      q=length(intersect(Genes,Gene_types))
      #print(q)
      m=length(Genes)
      n=N-m
      k=length(Gene_types)
      Overlap[ct,i]<-q
      if (test=="Over_representation") {
        p_value<-phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE) 
        fold_value=(q*N)/(m*k)
      } else if (test=="Under_representation") {
        #print("Under representation test")
        p_value<-phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE) 
        fold_value=(N*(m-q))/(m*(N-k))
      } else {
        fisherTest<-fisher.test(matrix(c(q, m-q, k-q, N-m-k+q), 2, 2), alternative='two.sided')
        p_value<-as.numeric(fisherTest$p.value)
        fold_value<-as.numeric(fisherTest$estimate)
        #as.numeric(fisherTest$estimate)
      } 
      
      #https://www.biostars.org/p/15548/
      #https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
      #print(p_value)
      P_set[ct]=p_value
      enrichment_set[ct]=fold_value
    }
    Data_mtrix[,i]<-P_set
    Fold[,i]<-enrichment_set
  }
  
  result=list()
  result$raw_pvalue<-Data_mtrix
  
  Data_mtrix<-Data_mtrix+2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  log10FDR<--log10(Mydata_raw_m)
  
  #print(log10FDR)
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  
  #breaksList = seq(0, 10, by = 1)
  dimnames(log10FDR)<-dimnames(Data_mtrix)
  
  df<-lengths(RefSet)
  percentage=df/sum(df)
  
  ##
  percentage_row<-percentage
  names(percentage_row)<-names(df)
  percentage_row<-percentage_row[names(RefSet)]
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(RefSet)]
  stats<-paste0(names(RefSet),stats)
  rownames(log10FDR)<-stats
  
  df<-lengths(Genelist)
  percentage=df/sum(df)
  
  ##
  percentage_col<-percentage
  names(percentage_col)<-names(df)
  percentage_col<-percentage_col[names(Genelist)]
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(Genelist)]
  stats<-paste0(names(Genelist),stats)
  colnames(log10FDR)<-stats
  colnames(log10FDR)<-gsub(pattern = "cluster_","",colnames(log10FDR))
  
  dimnames(Fold)<-dimnames(log10FDR)
  

  result$log10FDR=log10FDR
  result$Mydata_raw_m=Mydata_raw_m
  Fold[Fold==0]<-2.2e-16
  Fold[Fold==Inf]<-100
  result$Fold=log10(Fold)
  result$row=as.vector(percentage_row*100)
  result$col=as.vector(percentage_col*100)
  dimnames(Overlap)<-dimnames(log10FDR)
  result$Overlap<-as.matrix(Overlap)
  #print(result)
  print("done")
  
  return(result)
  
}



