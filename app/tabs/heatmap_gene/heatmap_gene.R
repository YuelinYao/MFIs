## UI ####
heatmapGenesUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap: Gene set"), style = "color: steelblue;"),
    plotOutput(outputId ="heatmap_GeneSet", width = "90%") %>% withSpinner(color="#4682B4")
  )}


stateGenesUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap: State genes"), style = "color: steelblue;"),
    plotOutput(outputId ="heatmapStateGenes", width = "90%") %>% withSpinner(color="#4682B4")
  )}


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
heatmapGenes <- function(cutoff,RefSet,Genelist,N=25678) { 

  RefSet<-makeGeneSetList(RefSet)
  #print(head(RefSet))
  r=length(names(RefSet))
  col=length(names(Genelist))
  
  EnrichmentM<-array(data=NA,dim = c(r,col))
  Data_mtrix<-array(data=NA,dim = c(r,col))
  
  colnames(Data_mtrix)<-names(Genelist)
  rownames(Data_mtrix)<-names(RefSet)
  
  dimnames(EnrichmentM)<-dimnames(Data_mtrix)
  
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
      p_value<-phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE) 
      #https://www.biostars.org/p/15548/
      #https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
      #print(p_value)
      P_set[ct]=p_value
      enrichment_set[ct]=(q*N)/(m*k)
    }
    Data_mtrix[,i]<-P_set
    EnrichmentM[,i]<-enrichment_set
  }
  
  Data_mtrix[Data_mtrix==0]<-2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  Data_mtrix_log<--log10(Mydata_raw_m)
  
  #print(Data_mtrix_log)
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  
  breaksList = seq(0, 10, by = 1)
  dimnames(Data_mtrix_log)<-dimnames(Data_mtrix)
  
  df<-lengths(RefSet)
  percentage=df/sum(df)
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(RefSet)]
  stats<-paste0(names(RefSet),stats)
  rownames(Data_mtrix_log)<-stats
  
  df<-lengths(Genelist)
  percentage=df/sum(df)
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(Genelist)]
  stats<-paste0(names(Genelist),stats)
  colnames(Data_mtrix_log)<-stats
  colnames(Data_mtrix_log)<-gsub(pattern = "cluster_","",colnames(Data_mtrix_log))
  
  dimnames(EnrichmentM)<-dimnames(Data_mtrix_log)
  
  result=list()
  result$Data_mtrix_log=Data_mtrix_log
  result$Mydata_raw_m=Mydata_raw_m
  result$Enrichment=EnrichmentM
  #print(result)
  print("done")
  
  return(result)
  
}


