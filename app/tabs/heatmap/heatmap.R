## UI ####
heatmapUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap: Cell types (clustering + singleR)"), style = "color: steelblue;"),
    plotOutput(outputId ="heatmap_celltypes", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap1","Download as .csv"),downloadButton("downloadheatmap_plot1","Download as .pdf")
  )}

NMF_UI <- function(){
  tagList(
    tags$h3("Heatmap: Cell states (NMF)", style = "color: steelblue;"),
    plotOutput(outputId ="heatmap_cellstates", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap2","Download as .csv"), downloadButton("downloadheatmap_plot2","Download as .pdf")
  )
}


NMF_CelltypeUI<-function(){
  tagList(
    tags$h3("Heatmap: Cell states vs. Cell types", style = "color: steelblue;"),
    plotOutput(outputId ="cellstates_types", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap3","Download as .csv"), downloadButton("downloadheatmap_plot3","Download as .pdf")
  )
}

Test_CellstateUI<-function(){
  tagList(
    tags$h3("Heatmap: Cell states vs. Cell states", style = "color: steelblue;"),
    plotOutput(outputId ="cellStates_cellStates", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap4","Download as .csv"), downloadButton("downloadheatmap_plot4","Download as .pdf")
  )
}



### Input function
HeatmapInput<- function(){
  tagList( 
    radioButtons(inputId = "colorHeatmapCells", "Colored by:",
                 choices = c("-log10FDR" = "log10FDR", "Enrichment" = "Enrichment"),
                 selected = "log10FDR", inline = TRUE),
    actionButton("action_heatmap","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}





####function
### Get List of cell
GetCellList <- function(count,summaryTable) { 
  #cutoff,count
  print("Get cells in each state")
  Devstates<-summaryTable
  Devstates$cluster<-paste0("C:",(Devstates$cluster))
  
  List=NULL
  for (i in unique(Devstates$cluster)){

    cluster<-Devstates[Devstates$cluster==i,]
    Cell<-NULL
    for (c in 1:length(cluster$genes)){
      #print(paste0("Cell states: ",cluster$genes[c]))
      Genes<-cluster$genes[c]
      Genes<-str_split(Genes, "_")[[1]]
      
      state<-cluster$state[c]
      State<-as.numeric(strsplit(as.character(state),"")[[1]])

      
      if ( length(Genes)==3 ) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]),]
      } else if ( length(Genes)==4) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]),]
      } else if ( length(Genes)==5) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]),]
      } else if ( length(Genes)==6){
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]),]
      } else if ( length(Genes)==7) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]),]
      } else if ( length(Genes)==8) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]),]
      } else if ( length(Genes)==9) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]&count[,Genes[9]]==State[9]),]
      } else {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]&count[,Genes[9]]==State[9]&count[,Genes[10]]==State[10]),]
      }
      
      
      cellnames<-rownames(c1)
      Cell<-c(Cell,cellnames)
    }  
    
    Cell<-unique(Cell)
  
    list=list(Cell)
    names(list)<-paste0("cluster_",i)
    List<-c(List,list)

  }
  return(List)
}

##GetCell in each d-tuple:
### Get List of cell
GetCellList_d <- function(count,summaryTable) { 
  #cutoff,count
  print("Get cells in each d-tuple")
  Devstates<-summaryTable
  Devstates$cluster<-paste(Devstates$genes,Devstates$state)
  
  List=NULL
  for (i in unique(Devstates$cluster)){
    
    cluster<-Devstates[Devstates$cluster==i,]
    Cell<-NULL
    for (c in 1:length(cluster$genes)){
      #print(paste0("Cell states: ",cluster$genes[c]))
      Genes<-cluster$genes[c]
      Genes<-str_split(Genes, "_")[[1]]
      
      state<-cluster$state[c]
      State<-as.numeric(strsplit(as.character(state),"")[[1]])
      
      
      if ( length(Genes)==3 ) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]),]
      } else if ( length(Genes)==4) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]),]
      } else if ( length(Genes)==5) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]),]
      } else if ( length(Genes)==6){
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]),]
      } else if ( length(Genes)==7) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]),]
      } else if ( length(Genes)==8) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]),]
      } else if ( length(Genes)==9) {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]&count[,Genes[9]]==State[9]),]
      } else {
        c1<-count[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]&count[,Genes[9]]==State[9]&count[,Genes[10]]==State[10]),]
      }
      
      
      cellnames<-rownames(c1)
      Cell<-c(Cell,cellnames)
    }  
    
    Cell<-unique(Cell)
    list=list(Cell)
    #names(list)<-paste0("cluster_",i)
    List<-c(List,list)
    
  }
  return(List)
}





### function
heatmap <- function(Meta_data,summaryTable,List) { 
  
  print("Over-representation test with cell type")
  Devstates<-summaryTable
  #print(paste0("The number of clusters: ",length(unique(Devstates$cluster))))
  
  r=length(unique(Meta_data$Cell_Types))
  col=length(unique(Devstates$cluster))

  Data_mtrix<-array(data=NA,dim = c(r,col))
  EnrichmentM<-array(data=NA,dim = c(r,col))
  
  Devstates$cluster<-paste0("C:",(Devstates$cluster))
  
  colnames(Data_mtrix)<-unique(Devstates$cluster)
  rownames(Data_mtrix)<-unique(Meta_data$Cell_Types)
  dimnames(EnrichmentM)<-dimnames(Data_mtrix)

  
  for (i in unique(Devstates$cluster)){
  #print(i)
  cluster<-Devstates[Devstates$cluster==i,]
  #print(paste0("Cluster: ",i))
  names<-paste0("cluster_",i)
  Cell<-List[[names]]
  P_set=rep(NA,length(unique(Meta_data$Cell_Types)))
  names(P_set)<-unique(Meta_data$Cell_Types)
    
  enrichment_set<-rep(NA,length(unique(Meta_data$Cell_Types)))
  names(enrichment_set)<-unique(Meta_data$Cell_Types)
  
    for (ct in unique(Meta_data$Cell_Types)){
      cell_types<-rownames(Meta_data)[which(Meta_data$Cell_Types==ct)]
      q=length(intersect(Cell,cell_types))
      m=length(Cell)
      N=dim(Meta_data)[1]
      n=N-m
      k=length(cell_types)
      p_value<-phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE) 
      #https://www.biostars.org/p/15548/
      #https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
      P_set[ct]=p_value
      enrichment_set[ct]=(q*N)/(m*k)
      
    }
    Data_mtrix[,i]<-P_set
    EnrichmentM[,i]<-enrichment_set
    
  }
  
  Data_mtrix<-Data_mtrix+2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  log10FDR<--log10(Mydata_raw_m)
  
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  
  #breaksList = seq(0, 10, by = 1)
  dimnames(log10FDR)<-dimnames(Data_mtrix)
  

  
  df<-table(Meta_data$Cell_Types)
  percentage=df/sum(df)
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[unique(Meta_data$Cell_Types)]
  stats<-paste0(unique(Meta_data$Cell_Types),stats)
  rownames(log10FDR)<-stats

  
  df<-lengths(List)
  percentage=df/sum(df)
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-gsub("cluster_","",names(df))
  stats<-stats[unique(Devstates$cluster)]
  stats<-paste0(unique(Devstates$cluster),stats)
  colnames(log10FDR)<-stats
  colnames(log10FDR)<-gsub(pattern = "C:","",colnames(log10FDR))
  dimnames(EnrichmentM)<-dimnames(log10FDR)
  
  result=list()
  result$log10FDR=log10FDR
  result$Mydata_raw_m=Mydata_raw_m
  result$Enrichment=EnrichmentM
  
  print("done")
  
  return(result)
  
}



NMF_heatmap<-function(Meta_data,summaryTable,List){
  
  print("Over-representation test with NMFs")
  Devstates<-summaryTable
  #print(paste0("Cutoff: ",cutoff)) 
  #print(paste0("The number of clusters: ",length(unique(Devstates$cluster))))

  r=length(unique(Meta_data$Cell_State))
  col=length(unique(Devstates$cluster))

  EnrichmentM<-array(data=NA,dim = c(r,col))
  Data_mtrix<-array(data=NA,dim = c(r,col))
  
  Devstates$cluster<-paste0("C:",(Devstates$cluster))
  
  colnames(Data_mtrix)<-unique(Devstates$cluster)
  rownames(Data_mtrix)<-unique(Meta_data$Cell_State)
  
  dimnames(EnrichmentM)<-dimnames(Data_mtrix)
  
  for (i in unique(Devstates$cluster)){
    
    #print(i)
    cluster<-Devstates[Devstates$cluster==i,]
    #print(paste0("Cluster: ",i))
    names<-paste0("cluster_",i)
    Cell<-List[[names]]
    
    P_set=rep(NA,length(unique(Meta_data$Cell_State)))
    names(P_set)<-unique(Meta_data$Cell_State)
    
    enrichment_set<-rep(NA,length(unique(Meta_data$Cell_State)))
    names(enrichment_set)<-unique(Meta_data$Cell_State)
    
    
    for (ct in unique(Meta_data$Cell_State)){
      cell_types<-rownames(Meta_data)[which(Meta_data$Cell_State==ct)]
      q=length(intersect(Cell,cell_types))
      m=length(Cell)
      N=dim(Meta_data)[1]
      n=N-m
      k=length(cell_types)
      p_value<-phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE) #https://www.biostars.org/p/15548/
      P_set[ct]=p_value
      enrichment_set[ct]=(q*N)/(m*k)
    }
    
    Data_mtrix[,i]<-P_set
    EnrichmentM[,i]<-enrichment_set
  }
  
  Data_mtrix<-Data_mtrix+2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  log10FDR<--log10(Mydata_raw_m)
  
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  
  #breaksList = seq(0, 10, by = 1)
  dimnames(log10FDR)<-dimnames(Data_mtrix)
  
  
  df<-table(Meta_data$Cell_State)
  percentage=df/sum(df)
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[unique(Meta_data$Cell_State)]
  stats<-paste0(unique(Meta_data$Cell_State),stats)
  rownames(log10FDR)<-stats
  
  
  df<-lengths(List)
  percentage=df/sum(df)
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-gsub("cluster_","",names(df))
  stats<-stats[unique(Devstates$cluster)]
  stats<-paste0(unique(Devstates$cluster),stats)
  colnames(log10FDR)<-stats
  colnames(log10FDR)<-gsub(pattern = "C:","",colnames(log10FDR))
  dimnames(EnrichmentM)<-dimnames(log10FDR)
  
  result=list()
  result$log10FDR=log10FDR
  result$Mydata_raw_m=Mydata_raw_m
  result$Enrichment=EnrichmentM
  
  return(result)
  
}

StateVsType<-function(Meta_data){

    print("Over-representation test cell types with NMFs")
    r=length(unique(Meta_data$Cell_State))
    col=length(unique(Meta_data$Cell_Types))
    
    Data_mtrix<-array(data=NA,dim = c(r,col))
    EnrichmentM<-array(data=NA,dim = c(r,col))
    
    colnames(Data_mtrix)<-unique(Meta_data$Cell_Types)
    rownames(Data_mtrix)<-unique(Meta_data$Cell_State)

    dimnames(EnrichmentM)<-dimnames(Data_mtrix)
    
    for (i in unique(Meta_data$Cell_Types)){
      #print(i)
      Cell<-rownames(Meta_data)[which(Meta_data$Cell_Types==i)]
      Cell
      
      P_set=rep(NA,length(unique(Meta_data$Cell_State)))
      names(P_set)<-unique(Meta_data$Cell_State)
      
      
      enrichment_set<-rep(NA,length(unique(Meta_data$Cell_State)))
      names(enrichment_set)<-unique(Meta_data$Cell_State)
      
      for (ct in unique(Meta_data$Cell_State)){
        
        #print(ct)
        cell_types<-rownames(Meta_data)[which(Meta_data$Cell_State==ct)]
        q=length(intersect(Cell,cell_types))
        m=length(Cell)
        N=dim(Meta_data)[1]
        n=N-m
        k=length(cell_types)
        p_value<-phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE) #https://www.biostars.org/p/15548/
        P_set[ct]=p_value
        enrichment_set[ct]=(q*N)/(m*k)
        
      }
      Data_mtrix[,i]<-P_set
      EnrichmentM[,i]<-enrichment_set
    }
  
    Data_mtrix<-Data_mtrix+2.2e-16
    Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
    Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
    log10FDR<--log10(Mydata_raw_m)
    
    Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
    Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
    Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
    
    #breaksList = seq(0, 10, by = 1)
    dimnames(log10FDR)<-dimnames(Data_mtrix)
   
    df<-table(Meta_data$Cell_Types)
    percentage=df/sum(df)
    stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
    names(stats)<-names(df)
    stats<-stats[unique(Meta_data$Cell_Types)]
    stats<-paste0(unique(Meta_data$Cell_Types),stats)
    colnames(log10FDR)<-stats
    
    df<-table(Meta_data$Cell_State)
    percentage=df/sum(df)
    stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
    names(stats)<-names(df)
    stats<-stats[unique(Meta_data$Cell_State)]
    stats<-paste0(unique(Meta_data$Cell_State),stats)
    rownames(log10FDR)<-stats
    
    dimnames(EnrichmentM)<-dimnames(log10FDR)
    
    
    result=list()
    result$log10FDR=log10FDR
    result$Mydata_raw_m=Mydata_raw_m
    result$Enrichment=EnrichmentM
    
    return(result)
  
}


ListTest <- function(list,N) { 
  
  print("Over-representation test between states")
  #print(head(list))
  r=length(names(list))
  col=length(names(list))
  
  Data_mtrix<-array(data=NA,dim = c(r,col))
  colnames(Data_mtrix)<-names(list)
  rownames(Data_mtrix)<-names(list)
  
  EnrichmentM<-array(data=NA,dim = c(r,col))
  dimnames(EnrichmentM)<-dimnames(Data_mtrix)
  
  for (i in names(list)){
    #print(i)
    Lists<-list[[i]]
    P_set=rep(NA,length(names(list)))
    names(P_set)<-names(list)
   
    
    enrichment_set<-rep(NA,length(names(list)))
    names(enrichment_set)<-names(list)
    
    
    for (ct in names(list)){
      #print(ct)
      List_types<-list[[ct]]
      q=length(intersect(Lists,List_types))
      #print(q)
      m=length(Lists)
      n=N-m
      k=length(List_types)
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
  
  Data_mtrix<-Data_mtrix+2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  log10FDR<--log10(Mydata_raw_m)

  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  
  #breaksList = seq(0, 10, by = 1)
  dimnames(log10FDR)<-dimnames(Data_mtrix)
  
  df<-lengths(list)
  percentage=df/sum(df)
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(list)]
  stats<-paste0(names(list),stats)
  rownames(log10FDR)<-stats
  
  df<-lengths(list)
  percentage=df/sum(df)
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(list)]
  stats<-paste0(names(list),stats)
  colnames(log10FDR)<-stats
  colnames(log10FDR)<-gsub(pattern = "cluster_","",colnames(log10FDR))
  rownames(log10FDR)<-gsub(pattern = "cluster_","",colnames(log10FDR))
  
  colnames(log10FDR)<-gsub(pattern = "C:","",colnames(log10FDR))
  rownames(log10FDR)<-gsub(pattern = "C:","",colnames(log10FDR))
  
  
  dimnames(EnrichmentM)<-dimnames(log10FDR)
  
  result=list()
  result$log10FDR=log10FDR
  result$Mydata_raw_m=Mydata_raw_m
  result$Enrichment=EnrichmentM
  print("done")
  #print(result)
  return(result)
  
}






