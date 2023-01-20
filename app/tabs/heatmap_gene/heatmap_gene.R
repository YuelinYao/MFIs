## UI ####
heatmapGenesUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap: Gene set"), style = "color: steelblue;"),
    plotOutput(outputId ="heatmap_GeneSet", width = "90%") %>% withSpinner(color="#4682B4")
  )}


### function
heatmapGenes <- function(cutoff,RefSet,Genelist,N=10000) { 

  load(paste0("./data/",RefSet,".Rdata"))
  #print(head(RefSet))
  r=length(names(RefSet))
  col=length(names(Genelist))
  
  Data_mtrix<-array(data=NA,dim = c(r,col))
  colnames(Data_mtrix)<-names(Genelist)
  rownames(Data_mtrix)<-names(RefSet)
  
  
  for (i in names(Genelist)){
    #print(i)
    Genes<-Genelist[[i]]
    P_set=rep(NA,length(names(RefSet)))
    names(P_set)<-names(RefSet)
    
    for (ct in names(RefSet)){
      #print(ct)
      Gene_types<-RefSet[[ct]]
      q=length(intersect(Genes,Gene_types))-1
      #print(q)
      m=length(Genes)
      n=N-m
      k=length(Gene_types)
      p_value<-phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) 
      #https://www.biostars.org/p/15548/
      #https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
      print(p_value)
      P_set[ct]=p_value
    }
    Data_mtrix[,i]<-P_set
    
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
  
  result=list()
  result$Data_mtrix_log=Data_mtrix_log
  result$Mydata_raw_m=Mydata_raw_m
  #print(result)
  print("done")
  
  return(result)
  
}



