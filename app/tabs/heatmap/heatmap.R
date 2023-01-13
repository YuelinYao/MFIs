## UI ####
heatmapUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap: Cell types (clustering + singleR)"), style = "color: steelblue;"),
    plotOutput(outputId ="heatmap_celltypes", width = "90%") %>% withSpinner(color="#4682B4")
  )}

NMF_UI <- function(){
  tagList(
    tags$h3("Heatmap: Cell states (NMF)", style = "color: steelblue;"),
    plotOutput(outputId ="heatmap_cellstates", width = "80%") %>% withSpinner(color="#4682B4")
  )
}


NMF_CelltypeUI<-function(){
  tagList(
    tags$h3("Heatmap: Cell states vs. Cell types", style = "color: steelblue;"),
    plotOutput(outputId ="cellstates_types", width = "80%") %>% withSpinner(color="#4682B4")
  )
}

####function
### Get List of cell
GetCellList <- function(cutoff,count,summaryTable) { 
  
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




### function
heatmap <- function(cutoff,Meta_data,summaryTable,List) { 
  
  print("Over-representation test with cell type")
  Devstates<-summaryTable
  #print(paste0("The number of clusters: ",length(unique(Devstates$cluster))))
  
  r=length(unique(Meta_data$Cell_Types))
  col=length(unique(Devstates$cluster))

  Data_mtrix<-array(data=NA,dim = c(r,col))
  Devstates$cluster<-paste0("C:",(Devstates$cluster))
  
  colnames(Data_mtrix)<-unique(Devstates$cluster)
  rownames(Data_mtrix)<-unique(Meta_data$Cell_Types)


  for (i in unique(Devstates$cluster)){
  #print(i)
  cluster<-Devstates[Devstates$cluster==i,]
  #print(paste0("Cluster: ",i))
  names<-paste0("cluster_",i)
  Cell<-List[[names]]
  P_set=rep(NA,length(unique(unique(Meta_data$Cell_Types))))
  names(P_set)<-unique(Meta_data$Cell_Types)
    
    for (ct in unique(Meta_data$Cell_Types)){
      cell_types<-rownames(Meta_data)[which(Meta_data$Cell_Types==ct)]
      q=length(intersect(Cell,cell_types))-1
      m=length(Cell)
      n=dim(Meta_data)[1]-m
      k=length(cell_types)
      p_value<-phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) 
      #https://www.biostars.org/p/15548/
      #https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
      
      P_set[ct]=p_value
    }
    Data_mtrix[,i]<-P_set
    
  }
  
  Data_mtrix[Data_mtrix==0]<-2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  Data_mtrix_log<--log10(Mydata_raw_m)
  
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  
  breaksList = seq(0, 10, by = 1)
  dimnames(Data_mtrix_log)<-dimnames(Data_mtrix)
  colnames(Data_mtrix_log)<-gsub(pattern = "C:","",colnames(Data_mtrix_log))

  
  result=list()
  result$Data_mtrix_log=Data_mtrix_log
  result$Mydata_raw_m=Mydata_raw_m

  print("done")
  
  return(result)
  
}



NMF_heatmap<-function(cutoff,Meta_data,summaryTable,List){
  
  print("Over-representation test with NMFs")
  Devstates<-summaryTable
  #print(paste0("Cutoff: ",cutoff)) 
  #print(paste0("The number of clusters: ",length(unique(Devstates$cluster))))

  r=length(unique(Meta_data$Cell_State))
  col=length(unique(Devstates$cluster))

  Data_mtrix<-array(data=NA,dim = c(r,col))
  Devstates$cluster<-paste0("C:",(Devstates$cluster))
  
  colnames(Data_mtrix)<-unique(Devstates$cluster)
  rownames(Data_mtrix)<-unique(Meta_data$Cell_State)
  
  for (i in unique(Devstates$cluster)){
    
    #print(i)
    cluster<-Devstates[Devstates$cluster==i,]
    #print(paste0("Cluster: ",i))
    names<-paste0("cluster_",i)
    Cell<-List[[names]]
    
    P_set=rep(NA,length(unique(unique(Meta_data$Cell_State))))
    names(P_set)<-unique(Meta_data$Cell_State)
    
    for (ct in unique(Meta_data$Cell_State)){
      cell_types<-rownames(Meta_data)[which(Meta_data$Cell_State==ct)]
      q=length(intersect(Cell,cell_types))-1
      m=length(Cell)
      n=dim(Meta_data)[1]-m
      k=length(cell_types)
      p_value<-phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) #https://www.biostars.org/p/15548/
      P_set[ct]=p_value
      
    }
    
    Data_mtrix[,i]<-P_set
    
  }
  
  Data_mtrix[Data_mtrix==0]<-2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  Data_mtrix_log<--log10(Mydata_raw_m)
  
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  
  breaksList = seq(0, 10, by = 1)
  dimnames(Data_mtrix_log)<-dimnames(Data_mtrix)
  colnames(Data_mtrix_log)<-gsub(pattern = "C:","",colnames(Data_mtrix_log))
  
  result=list()
  result$Data_mtrix_log=Data_mtrix_log
  result$Mydata_raw_m=Mydata_raw_m
  
  return(result)
  
}

StateVsType<-function(Meta_data){

    print("Over-representation test cell types with NMFs")
    r=length(unique(Meta_data$Cell_State))
    col=length(unique(Meta_data$Cell_Types))
    
    Data_mtrix<-array(data=NA,dim = c(r,col))
    colnames(Data_mtrix)<-unique(Meta_data$Cell_Types)
    rownames(Data_mtrix)<-unique(Meta_data$Cell_State)

  
    for (i in unique(Meta_data$Cell_Types)){
      #print(i)
      Cell<-rownames(Meta_data)[which(Meta_data$Cell_Types==i)]
      Cell
      
      P_set=rep(NA,length(unique(unique(Meta_data$Cell_State))))
      names(P_set)<-unique(Meta_data$Cell_State)
      
      for (ct in unique(Meta_data$Cell_State)){
        
        #print(ct)
        cell_types<-rownames(Meta_data)[which(Meta_data$Cell_State==ct)]
        q=length(intersect(Cell,cell_types))-1
        m=length(Cell)
        n=dim(Meta_data)[1]-m
        k=length(cell_types)
        p_value<-phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE) #https://www.biostars.org/p/15548/
        P_set[ct]=p_value
        
      }
      Data_mtrix[,i]<-P_set
    }
  
    Data_mtrix[Data_mtrix==0]<-2.2e-16
    Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
    Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
    Data_mtrix_log<--log10(Mydata_raw_m)
    
    Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
    Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
    Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
    
    breaksList = seq(0, 10, by = 1)
    dimnames(Data_mtrix_log)<-dimnames(Data_mtrix)
    colnames(Data_mtrix_log)<-gsub(pattern = "C:","",colnames(Data_mtrix_log))
    
    result=list()
    result$Data_mtrix_log=Data_mtrix_log
    result$Mydata_raw_m=Mydata_raw_m
    
    return(result)
  
}
