## UI ####
heatmapUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap 1"), style = "color: steelblue;"), # Cell types (clustering + singleR)
    tags$h5("Comparing Stator states vs. the annotations from the first column in meta data table."), 
    tags$h5("For liver cancer dataset: Stator states vs. Cell types (by clustering + singleR)."), 
    tags$br(),
    plotOutput(outputId ="heatmap_celltypes", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap1","Download as .csv"),
    downloadButton("downloadheatmap_plot1","Download as .pdf"),
    downloadButton("downloadheatmap1_pvalue","Download raw pvalue as .csv"),
  )}

NMF_UI <- function(){
  tagList(
    tags$h3("Heatmap 2", style = "color: steelblue;"), # Cell states (NMF)
    tags$h5("Comparing Stator states vs. the annotations from the second column in meta data table."), 
    tags$h5("For liver cancer dataset: Stator states vs. Cell states (by NMF)."), 
    tags$br(),
    plotOutput(outputId ="heatmap_cellstates", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap2","Download as .csv"), 
    downloadButton("downloadheatmap_plot2","Download as .pdf"),
    downloadButton("downloadheatmap2_pvalue","Download raw pvalue as .csv")
  )
}


NMF_CelltypeUI<-function(){
  tagList(
    tags$h3("Heatmap 3", style = "color: steelblue;"), #  Cell states vs. Cell types
    tags$h5("Comparing the annotations between the first and second column in meta data table."), 
    tags$h5("For liver cancer dataset: Cell states vs. Cell types."), 
    tags$br(),
    plotOutput(outputId ="cellstates_types", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap3","Download as .csv"), 
    downloadButton("downloadheatmap_plot3","Download as .pdf"),
    downloadButton("downloadheatmap3_pvalue","Download raw pvalue as .csv"), 
  )
}

Test_CellstateUI<-function(){
  tagList(
    tags$h3("Heatmap 4", style = "color: steelblue;"), # Cell states vs. Cell states
    tags$h5("Comparing between Stator states."), 
    tags$br(),
    plotOutput(outputId ="cellStates_cellStates", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("downloadheatmap4","Download as .csv"), 
    downloadButton("downloadheatmap_plot4","Download as .pdf"),
    downloadButton("downloadheatmap4_pvalue","Download raw pvalue as .csv"), 
  )
}


iheatmapUI <- function(){
  tagList(
    tags$h3(paste0("Heatmap 1"), style = "color: steelblue;"),
    plotlyOutput(outputId ="iheatmap_celltypes", width = "90%") %>% withSpinner(color="#4682B4"),
    
    tags$h3("Heatmap 2", style = "color: steelblue;"),
    plotlyOutput(outputId ="iheatmap_cellstates", width = "80%") %>% withSpinner(color="#4682B4"),
    
    tags$h3("Heatmap 3", style = "color: steelblue;"),
    plotlyOutput(outputId ="icellstates_types", width = "80%") %>% withSpinner(color="#4682B4"),
    
    tags$h3("Heatmap 4", style = "color: steelblue;"),
    plotlyOutput(outputId ="icellStates_cellStates", width = "80%") %>% withSpinner(color="#4682B4"),
    
  )}







### Input function
HeatmapInput<- function(){
  tagList( 
    radioButtons(inputId = "TestCells", "Test:",
                 choices = c("Enrichment" = "Over_representation", "Depletion" = "Under_representation", "Two.side" = "Fisher"), 
                 selected = "Over_representation", inline = TRUE),
    #radioButtons(inputId = "colorHeatmapCells", "Colored by:",
    #             choices = c("-log10FDR" = "log10FDR", "Fold" = "Fold"),
    #            selected = "log10FDR", inline = TRUE),
    actionButton("action_heatmap","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}





####function
### Get List of cell
GetCellList <- function(count,summaryTable) { 
  #cutoff,count
  # colnames(count)<-gsub("-",".",colnames(count))
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
      
      count_names<-rownames(count)
      
      if ( length(Genes)==3 ) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3])]
      } else if ( length(Genes)==4) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4])]
      } else if ( length(Genes)==5) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5])]
      } else if ( length(Genes)==6){
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6])]
      } else if ( length(Genes)==7) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7])]
      } else if ( length(Genes)==8) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8])]
      } else if ( length(Genes)==9) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]&count[,Genes[9]]==State[9])]
      } else {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]&count[,Genes[9]]==State[9]&count[,Genes[10]]==State[10])]
      }
      
      
      cellnames<-c1
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
  # colnames(count)<-gsub("-",".",colnames(count))
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
      
      
      count_names<-rownames(count)
      
      if ( length(Genes)==3 ) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3])]
      } else if ( length(Genes)==4) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4])]
      } else if ( length(Genes)==5) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5])]
      } else if ( length(Genes)==6){
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6])]
      } else if ( length(Genes)==7) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7])]
      } else if ( length(Genes)==8) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8])]
      } else if ( length(Genes)==9) {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]&count[,Genes[9]]==State[9])]
      } else {
        c1<- count_names[which(count[,Genes[1]]==State[1]&count[,Genes[2]]==State[2]&count[,Genes[3]]==State[3]&count[,Genes[4]]==State[4]&count[,Genes[5]]==State[5]&count[,Genes[6]]==State[6]&count[,Genes[7]]==State[7]&count[,Genes[8]]==State[8]&count[,Genes[9]]==State[9]&count[,Genes[10]]==State[10])]
      }
      
      
      cellnames<-c1
      Cell<-c(Cell,cellnames)
    }  
    
    Cell<-unique(Cell)
    list=list(Cell)
    names(list)<-paste0(i)
    List<-c(List,list)
    
  }
  return(List)
}





### function
heatmap <- function(Meta_data,summaryTable,List,N,test) { 
  
  print(test)
  Devstates<-summaryTable
  #print(paste0("The number of clusters: ",length(unique(Devstates$cluster))))
  
  r=length(unique(Meta_data$Cell_Types))
  col=length(unique(Devstates$cluster))
  
  Data_mtrix<-array(data=NA,dim = c(r,col))
  Fold<-array(data=NA,dim = c(r,col))
  
  
  
  Devstates$cluster<-paste0("C:",(Devstates$cluster))
  
  colnames(Data_mtrix)<-unique(Devstates$cluster)
  rownames(Data_mtrix)<-unique(Meta_data$Cell_Types)
  dimnames(Fold)<-dimnames(Data_mtrix)
  
  Overlap<-array(data=NA,dim = c(r,col))
  dimnames(Overlap)<-dimnames(Data_mtrix)
  
  for (i in unique(Devstates$cluster)){
    #print(i)
    cluster<-Devstates[Devstates$cluster==i,]
    #print(paste0("Cluster: ",i))
    names<-paste0("cluster_",i)
    Cell<-List[[names]]
    P_set=rep(NA,length(unique(Meta_data$Cell_Types)))
    names(P_set)<-unique(Meta_data$Cell_Types)
    
    fold_set<-rep(NA,length(unique(Meta_data$Cell_Types)))
    names(fold_set)<-unique(Meta_data$Cell_Types)
    
    for (ct in unique(Meta_data$Cell_Types)){
      cell_types<-rownames(Meta_data)[which(Meta_data$Cell_Types==ct)]
      q=length(intersect(Cell,cell_types))
      m=length(Cell)
      n=N-m #N=dim(Meta_data)[1]
      k=length(cell_types)
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
      #https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
      P_set[ct]=p_value
      fold_set[ct]=fold_value
    }
    Data_mtrix[,i]<-P_set
    Fold[,i]<-fold_set
    
  }
  result=list()
  result$raw_pvalue<-Data_mtrix
  
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
  
  ##
  percentage_row<-percentage
  names(percentage_row)<-names(df)
  percentage_row<-percentage_row[unique(Meta_data$Cell_Types)]
  ##
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[unique(Meta_data$Cell_Types)]
  stats<-paste0(unique(Meta_data$Cell_Types),stats)
  rownames(log10FDR)<-stats
  
  
  df<-lengths(List)
  percentage=df/sum(df)
  
  ##
  percentage_col<-percentage
  names(percentage_col)<-gsub("cluster_","",names(df))
  percentage_col<-percentage_col[unique(Devstates$cluster)]
  ##
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-gsub("cluster_","",names(df))
  stats<-stats[unique(Devstates$cluster)]
  stats<-paste0(unique(Devstates$cluster),stats)
  colnames(log10FDR)<-stats
  colnames(log10FDR)<-gsub(pattern = "C:","",colnames(log10FDR))
  dimnames(Fold)<-dimnames(log10FDR)
  #dimnames(Overlap)<-dimnames(log10FDR)
  
  
  result$log10FDR=log10FDR
  result$Mydata_raw_m=Mydata_raw_m
  Fold[Fold==0]<-2.2e-16
  Fold[Fold==Inf]<-100
  result$Fold=log10(Fold)
  
  result$row=as.vector(percentage_row*100)
  result$col=as.vector(percentage_col*100)
  #print(result$col)
  dimnames(Overlap)<-dimnames(log10FDR)
  result$Overlap<-as.matrix(Overlap)
  #print(Overlap)
  print("done")
  
  return(result)
  
}



NMF_heatmap<-function(Meta_data,summaryTable,List,N,test){
  
  #print("Over-representation test with NMFs")
  Devstates<-summaryTable
  #print(paste0("Cutoff: ",cutoff)) 
  #print(paste0("The number of clusters: ",length(unique(Devstates$cluster))))
  
  r=length(unique(Meta_data$Cell_State))
  col=length(unique(Devstates$cluster))
  
  Fold<-array(data=NA,dim = c(r,col))
  Data_mtrix<-array(data=NA,dim = c(r,col))
  
  Devstates$cluster<-paste0("C:",(Devstates$cluster))
  
  colnames(Data_mtrix)<-unique(Devstates$cluster)
  rownames(Data_mtrix)<-unique(Meta_data$Cell_State)
  
  dimnames(Fold)<-dimnames(Data_mtrix)
  
  Overlap<-array(data=NA,dim = c(r,col))
  dimnames(Overlap)<-dimnames(Data_mtrix)
  
  
  for (i in unique(Devstates$cluster)){
    
    #print(i)
    cluster<-Devstates[Devstates$cluster==i,]
    #print(paste0("Cluster: ",i))
    names<-paste0("cluster_",i)
    Cell<-List[[names]]
    
    P_set=rep(NA,length(unique(Meta_data$Cell_State)))
    names(P_set)<-unique(Meta_data$Cell_State)
    
    fold_set<-rep(NA,length(unique(Meta_data$Cell_State)))
    names(fold_set)<-unique(Meta_data$Cell_State)
    
    
    for (ct in unique(Meta_data$Cell_State)){
      cell_types<-rownames(Meta_data)[which(Meta_data$Cell_State==ct)]
      q=length(intersect(Cell,cell_types))
      m=length(Cell)
      #N=dim(Meta_data)[1]
      n=N-m
      k=length(cell_types)
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
      #https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
      
      P_set[ct]=p_value
      fold_set[ct]=fold_value
    }
    
    Data_mtrix[,i]<-P_set
    Fold[,i]<-fold_set
  }
  
  result=list()
  result$raw_pvalue<-Data_mtrix
  
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
  
  ##
  percentage_row<-percentage
  names(percentage_row)<-names(df)
  percentage_row<-percentage_row[unique(Meta_data$Cell_State)]
  ##
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[unique(Meta_data$Cell_State)]
  stats<-paste0(unique(Meta_data$Cell_State),stats)
  rownames(log10FDR)<-stats
  
  
  df<-lengths(List)
  percentage=df/sum(df)
  
  ##
  percentage_col<-percentage
  names(percentage_col)<-gsub("cluster_","",names(df))
  percentage_col<-percentage_col[unique(Devstates$cluster)]
  ##
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-gsub("cluster_","",names(df))
  stats<-stats[unique(Devstates$cluster)]
  stats<-paste0(unique(Devstates$cluster),stats)
  colnames(log10FDR)<-stats
  colnames(log10FDR)<-gsub(pattern = "C:","",colnames(log10FDR))
  dimnames(Fold)<-dimnames(log10FDR)
  
  
  result$log10FDR=log10FDR
  result$Mydata_raw_m=Mydata_raw_m
  Fold[Fold==0]<-2.2e-16
  Fold[Fold==Inf]<-100
  result$Fold=log10(Fold)
  
  dimnames(Overlap)<-dimnames(log10FDR)
  result$Overlap<-as.matrix(Overlap)
  
  result$row=as.vector(percentage_row*100)
  result$col=as.vector(percentage_col*100)
  
  return(result)
  
}

StateVsType<-function(Meta_data,N,test){
  
  #print("Over-representation test cell types with NMFs")
  r=length(unique(Meta_data$Cell_State))
  col=length(unique(Meta_data$Cell_Types))
  
  Data_mtrix<-array(data=NA,dim = c(r,col))
  Fold<-array(data=NA,dim = c(r,col))
  
  colnames(Data_mtrix)<-unique(Meta_data$Cell_Types)
  rownames(Data_mtrix)<-unique(Meta_data$Cell_State)
  
  dimnames(Fold)<-dimnames(Data_mtrix)
  
  Overlap<-array(data=NA,dim = c(r,col))
  dimnames(Overlap)<-dimnames(Data_mtrix)
  
  for (i in unique(Meta_data$Cell_Types)){
    #print(i)
    Cell<-rownames(Meta_data)[which(Meta_data$Cell_Types==i)]
    Cell
    
    P_set=rep(NA,length(unique(Meta_data$Cell_State)))
    names(P_set)<-unique(Meta_data$Cell_State)
    
    
    fold_set<-rep(NA,length(unique(Meta_data$Cell_State)))
    names(fold_set)<-unique(Meta_data$Cell_State)
    
    for (ct in unique(Meta_data$Cell_State)){
      
      #print(ct)
      cell_types<-rownames(Meta_data)[which(Meta_data$Cell_State==ct)]
      q=length(intersect(Cell,cell_types))
      m=length(Cell)
      #N=dim(Meta_data)[1]
      n=N-m
      k=length(cell_types)
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
      
      P_set[ct]=p_value
      fold_set[ct]=fold_value
      
    }
    Data_mtrix[,i]<-P_set
    Fold[,i]<-fold_set
  }
  
  
  result=list()
  result$raw_pvalue<-Data_mtrix
  
  
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
  
  percentage_col<-percentage
  names(percentage_col)<-names(df)
  percentage_col<-percentage_col[unique(Meta_data$Cell_Types)]
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[unique(Meta_data$Cell_Types)]
  stats<-paste0(unique(Meta_data$Cell_Types),stats)
  colnames(log10FDR)<-stats
  
  df<-table(Meta_data$Cell_State)
  percentage=df/sum(df)
  
  percentage_row<-percentage
  names(percentage_row)<-gsub("cluster_","",names(df))
  percentage_row<-percentage_row[unique(Meta_data$Cell_State)]
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[unique(Meta_data$Cell_State)]
  stats<-paste0(unique(Meta_data$Cell_State),stats)
  rownames(log10FDR)<-stats
  
  dimnames(Fold)<-dimnames(log10FDR)
  
  
  result$log10FDR=log10FDR
  result$Mydata_raw_m=Mydata_raw_m
  Fold[Fold==0]<-2.2e-16
  Fold[Fold==Inf]<-100
  result$Fold=log10(Fold)
  
  dimnames(Overlap)<-dimnames(log10FDR)
  result$Overlap<-as.matrix(Overlap)
  
  result$row=as.vector(percentage_row*100)
  result$col=as.vector(percentage_col*100)
  
  return(result)
  
}


ListTest <- function(list,N,test) { 
  
  r=length(names(list))
  col=length(names(list))
  
  Data_mtrix<-array(data=NA,dim = c(r,col))
  colnames(Data_mtrix)<-names(list)
  rownames(Data_mtrix)<-names(list)
  
  Fold<-array(data=NA,dim = c(r,col))
  dimnames(Fold)<-dimnames(Data_mtrix)
  
  Overlap<-array(data=NA,dim = c(r,col))
  dimnames(Overlap)<-dimnames(Data_mtrix)
  
  for (i in names(list)){
    #print(i)
    Lists<-list[[i]]
    P_set=rep(NA,length(names(list)))
    names(P_set)<-names(list)
    
    
    fold_set<-rep(NA,length(names(list)))
    names(fold_set)<-names(list)
    
    
    for (ct in names(list)){
      #print(ct)
      List_types<-list[[ct]]
      q=length(intersect(Lists,List_types))
      #print(q)
      m=length(Lists)
      n=N-m
      k=length(List_types)
      Overlap[ct,i]<-q
      if (test=="Over_representation") {
        p_value<-phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE) 
        fold_value=(q*N)/(m*k)
      } else if (test=="Under_representation") {
        # print("Under representation test")
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
      fold_set[ct]=fold_value
    }
    Data_mtrix[,i]<-P_set
    Fold[,i]<-fold_set
    
  }
  
  result=list()
  result$raw_pvalue<-Data_mtrix
  
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
  
  percentage_row<-percentage
  names(percentage_row)<-names(df)
  percentage_row<-percentage_row[names(list)]
  
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(list)]
  stats<-paste0(names(list),stats)
  rownames(log10FDR)<-stats
  
  df<-lengths(list)
  percentage=df/sum(df)
  
  percentage_col<-percentage
  names(percentage_col)<-names(df)
  percentage_col<-percentage_col[names(list)]
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(list)]
  stats<-paste0(names(list),stats)
  colnames(log10FDR)<-stats
  colnames(log10FDR)<-gsub(pattern = "cluster_","",colnames(log10FDR))
  rownames(log10FDR)<-gsub(pattern = "cluster_","",colnames(log10FDR))
  
  colnames(log10FDR)<-gsub(pattern = "C:","",colnames(log10FDR))
  rownames(log10FDR)<-gsub(pattern = "C:","",colnames(log10FDR))
  
  
  dimnames(Fold)<-dimnames(log10FDR)
  
  
  result$log10FDR=log10FDR
  result$Mydata_raw_m=Mydata_raw_m
  Fold[Fold==0]<-2.2e-16
  Fold[Fold==Inf]<-100
  result$Fold=log10(Fold)
  
  dimnames(Overlap)<-dimnames(log10FDR)
  result$Overlap<-as.matrix(Overlap)
  
  result$row=as.vector(percentage_row*100)
  result$col=as.vector(percentage_col*100)
  
  print("done")
  return(result)
  
}






