## UI ####
FunctionUI <- function(){
  tagList(
    tags$h3(paste0("Functional enrichment"), style = "color: steelblue;"),
    plotOutput(outputId ="Functional_enrichment", width = "90%") %>% withSpinner(color="#4682B4")
  )}


Function_tableUI<-function(){
  tagList(
    tags$h3(paste0("Functional enrichment table"), style = "color: steelblue;"),
    DTOutput(outputId ="Functional_table", width = "80%") %>% withSpinner(color="#4682B4")
  )}


## show all the reference set from ENSEMBL
mart <- useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(mart)
datasets_list<-as.list(datasets$dataset)
names(datasets_list)<-listDatasets(mart)$dataset


## Gene_list list in each states
GetGenes<-function(cutoff,Devstates){
  print("Get genes in each state")
  Geneset<-NULL
  for (i in unique(Devstates$cluster)){

    cluster<-Devstates[Devstates$cluster==i,]
    print(paste0("Cluster: ",i))
    Genes_set<-NULL
    
    for (c in 1:length(cluster$genes)){
      #print(paste0("Cell states: ",cluster$genes[c]))
      Genes<-cluster$genes[c]
      Genes<-str_split(Genes, "_")[[1]]
      
      #print(length(Genes))
      state<-cluster$state[c]
      State<-as.numeric(strsplit(as.character(state),"")[[1]])
      #print(State)
      
      Genes_state<-paste(Genes,State,sep = "_")
      Genes_set<-c(Genes_set,Genes_state)
      
    }
    
    Genes_set<-unique(Genes_set)
    Genes_set<- Genes_set[!grepl("_0", Genes_set)] # remove 0-state genes
    Genes_set<-gsub("_1","",Genes_set)
  
    list=list(Genes_set)
    names(list)<-paste0("cluster_",i)
    Geneset<-c(Geneset,list)
      
  }
  
  return(Geneset)
  
}



## Gene_list list in each states with 0-state genes
GetGenes2<-function(cutoff,Devstates){
  print("Get genes in each state")
  Geneset<-NULL
  for (i in unique(Devstates$cluster)){
    
    cluster<-Devstates[Devstates$cluster==i,]
    #print(paste0("Cluster: ",i))
    Genes_set<-NULL
    
    for (c in 1:length(cluster$genes)){
      #print(paste0("Cell states: ",cluster$genes[c]))
      Genes<-cluster$genes[c]
      Genes<-str_split(Genes, "_")[[1]]
      
      #print(length(Genes))
      state<-cluster$state[c]
      State<-as.numeric(strsplit(as.character(state),"")[[1]])
      #print(State)
      
      Genes_state<-paste(Genes,State,sep = "_")
      Genes_set<-c(Genes_set,Genes_state)
      
    }
    
    
    Genes_set<-gsub("_1","",Genes_set)
    Genes_set<-gsub("_0","",Genes_set)
    Genes_set<-unique(Genes_set)
    
    list=list(Genes_set)
    names(list)<-paste0("cluster_",i)
    Geneset<-c(Geneset,list)
    
  }
  
  return(Geneset)
  
}




FunctionE <- function(cutoff,selected_cluster,Mart,kegg_species,go_species,GenesList, background_genes){

  if (!any(rownames(installed.packages()) == go_species)){
    BiocManager::install(go_species,update = F)
  }
  library(go_species,character.only = T)
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset(Mart, mart) 
  
  selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]

  print(paste0("GO & KEGG for cluster(s): ",selected_cluster))
  print("GO & KEGG...")
  
  AllEnrichment<-NULL
  
  if (length(background_genes)==0){
    
    print("No user-defined background genes")
    
    for (i in selected_cluster){
      print(paste0("Cluster: ",i))
      names<-paste0("cluster_",i)
      Genes_set<-GenesList[[names]]
      print(length(Genes_set))
      if (length(Genes_set)>0) {
        
      Genes_set<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                       filter="external_gene_name", values=Genes_set, uniqueRows=TRUE)
      print(head(Genes_set))
      print("KEGG")
      cluster_kegg <- enrichKEGG(gene =  Genes_set$entrezgene_id,organism = kegg_species,
                                 pAdjustMethod = "BH",
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)
      
      
      if (is.null(cluster_kegg)) {
        print("no result")
      } else{
      
      cluster_kegg<-cluster_kegg@result
      cluster_kegg$cluster<-i
      cluster_kegg$class<-"KEGG"}
      
      print("GO")
      cluster_GOBP <- enrichGO(gene = Genes_set$entrezgene_id,
                               OrgDb= go_species,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               minGSSize = 1,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE)
      
      if (is.null(cluster_GOBP)) {
        print("no result")
      } else{
      cluster_GOBP<-cluster_GOBP@result
      cluster_GOBP$cluster<-i
      cluster_GOBP$class<-"GOBP"}
      
      
      cluster_GOCC <- enrichGO(gene = Genes_set$entrezgene_id,
                               OrgDb= go_species,
                               ont = "CC",
                               pAdjustMethod = "BH",
                               minGSSize = 1,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE)
      
      if (is.null(cluster_GOCC)) {
        print("no result")
      } else{
      cluster_GOCC<-cluster_GOCC@result
      cluster_GOCC$cluster<-i
      cluster_GOCC$class<-"GOCC"}
      
    
      cluster_GOMF <- enrichGO(gene = Genes_set$entrezgene_id,
                               OrgDb= go_species,
                               ont = "MF",
                               pAdjustMethod = "BH",
                               minGSSize = 1,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE)
      
      if (is.null(cluster_GOMF)) {
        print("no result")
      } else{
      
      cluster_GOMF<-cluster_GOMF@result
      cluster_GOMF$cluster<-i
      cluster_GOMF$class<-"GOMF"}
      
      
      ClusterEnrich<-rbind(cluster_GOBP,cluster_GOCC,cluster_kegg,cluster_GOMF)
      
      AllEnrichment<-rbind(AllEnrichment,ClusterEnrich)  }
      
    }
    
  }
  
  else{
  print("Use input background genes")
  background_genes<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                                               filter="external_gene_name", values=background_genes, uniqueRows=TRUE)
  
  
  for (i in selected_cluster){

  print(paste0("Cluster: ",i))
  names<-paste0("cluster_",i)
  Genes_set<-GenesList[[names]]
  print(length(Genes_set))
  if (length(Genes_set)>0) {
  Genes_set<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                   filter="external_gene_name", values=Genes_set, uniqueRows=TRUE)
  
  cluster_kegg <- enrichKEGG(gene =  Genes_set$entrezgene_id,organism = kegg_species,
                           pAdjustMethod = "BH",universe = as.character(background_genes$entrezgene_id),
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
  
  
  if (is.null(cluster_kegg)) {
    print("no result")
  } else{
  cluster_kegg<-cluster_kegg@result
  cluster_kegg$cluster<-i
  cluster_kegg$class<-"KEGG"}

  
  cluster_GOBP <- enrichGO(gene = Genes_set$entrezgene_id,universe = as.character(background_genes$entrezgene_id),
                       OrgDb= go_species,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
  
  if (is.null(cluster_GOBP)) {
    print("no result")
  } else{
  cluster_GOBP<-cluster_GOBP@result
  cluster_GOBP$cluster<-i
  cluster_GOBP$class<-"GOBP"}
  
  
  cluster_GOCC <- enrichGO(gene = Genes_set$entrezgene_id,universe = as.character(background_genes$entrezgene_id),
                       OrgDb= go_species,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
  
  if (is.null(cluster_GOCC)) {
    print("no result")
  } else{
  cluster_GOCC<-cluster_GOCC@result
  cluster_GOCC$cluster<-i
  cluster_GOCC$class<-"GOCC"}
  
  cluster_GOMF <- enrichGO(gene = Genes_set$entrezgene_id,universe = as.character(background_genes$entrezgene_id),
                       OrgDb= go_species,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       minGSSize = 1,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
  
  if (is.null(cluster_GOMF)) {
    print("no result")
  } else{
  
  cluster_GOMF<-cluster_GOMF@result
  cluster_GOMF$cluster<-i
  cluster_GOMF$class<-"GOMF"}
  
  ClusterEnrich<-rbind(cluster_GOBP,cluster_GOCC,cluster_kegg,cluster_GOMF)

  AllEnrichment<-rbind(AllEnrichment,ClusterEnrich)  }

  }

  
}
  
  AllEnrichment<-AllEnrichment[AllEnrichment$p.adjust<0.05,]
  return(AllEnrichment)
  
}






Plot_enrichment<-function(GO){
  print("Plot")
  all_function<-GO %>%
  group_by(cluster,class) %>%
  slice_max(n = 4,order_by = Count)%>%
  slice_head(n=4)#, Count

  all_function$Gene_number <- all_function$Count
  all_function$Description<-as.factor(all_function$Description)
  all_function$'|log10(FDR)|' <- -(log10(all_function$p.adjust))
  all_function<-all_function[order(all_function$cluster),]
  all_function$Description<-paste(all_function$class,all_function$Description,sep=":")
  all_function$Description<-factor(all_function$Description,levels=unique(all_function$Description))
  all_function<-all_function[!is.na(all_function$ID),]
  all_function<-all_function[all_function$p.adjust<0.05,]
  
  p<-ggplot(all_function, aes(x = Description, y = cluster)) +
  geom_point(data=all_function,aes(x=Description, y=cluster,size = Gene_number, colour = `|log10(FDR)|`), alpha=.7)+
  scale_color_gradient(low = "blue", high = "orange", limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt"),angle = 45, hjust = 1,size=13),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),size=13),
        axis.text = element_text(color = "black"),
        axis.title=element_blank())+
  #xlab("GO & KEGG functional enrichment")+
  labs(color="-log10(FDR)", size="Number\nof genes")
return(p)
}
