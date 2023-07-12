## UI
DE_UI <- function(){
  tagList(
    tags$h3(paste0("DE analysis for two states"), style = "color: steelblue;"),
    plotOutput(outputId ="DE_heatmap", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("GeneHeatmap","Download as .pdf")
  )}


DE_TableUI <- function(){
  tagList(
    tags$h3(paste0("DEGs result table"), style = "color: steelblue;"),
    DTOutput(outputId ="DEGtable", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("DegTable","Download as .csv")
  )}



DEGO_UI<- function(){
  tagList(
    tags$h3(paste0("GO & KEGG for DEGs"), style = "color: steelblue;"),
    plotOutput(outputId ="DEG_enrichment", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("DegAnnotationPlot","Download as .pdf")
  )}



DEGOtable_UI<- function(){
  tagList(
    tags$h3(paste0("GO & KEGG for DEGs table"), style = "color: steelblue;"),
    DTOutput(outputId ="DEG_enrichmentTable", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("DegAnnotationTable","Download as .csv")
  )}



### Input function
DEInput<- function(){
  tagList( 
    textInput("selected_clusterDE", "Select cluster(s)",value = "46,19"),
    selectInput("Mart_DE", "Mart dataset:", choices=datasets_list, selected = "hsapiens_gene_ensembl", multiple = FALSE),
    textInput("go_species_DE", "GO OrgDb:",value = "org.Hs.eg.db"),
    textInput("kegg_species_DE", "KEGG organism:",value = "hsa"),
    textInput("logfc", "logFC:",value = 0.25),
    sliderInput("Pvalue_DE",
                "Adjusted p value:",
                min = 0,
                max = 1,
                value = 0.05),
    textAreaInput("background_genesDE", "Background genes (recommended): ",placeholder = "Just paste a list of genes (multiple-line gene list).",rows = 5),
    actionButton(inputId = "bg_Liver2",                                       #action button to display background genes
                 label = NULL, icon = icon("tag"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    bsTooltip("bg_Liver2","Load background genes in HCC dataset.",placement = "bottom", trigger = "hover",
              options = NULL),
    actionButton("action_DE","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}






DE_set<-function(selected_cluster,cutoff,count,srt,logfc,Pvalue,List){
  
  print("Find DE genes")
  selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]
  selected_cluster<-paste0("cluster_C:",selected_cluster)
  
  selected_1<-List[[selected_cluster[1]]]
  print(length(selected_1))
  selected_2<-List[[selected_cluster[2]]]
  print(length(selected_2))
  
  overlap<-intersect(selected_1,selected_2)
  
  selected_1<-selected_1[!selected_1%in%overlap]
  selected_2<-selected_2[!selected_2%in%overlap]
  
  DE<-FindMarkers(object = srt,ident.1 = selected_1,ident.2 = selected_2,logfc.threshold = logfc)
  DE<-DE[DE$p_val_adj<Pvalue,]
  DE$gene<-rownames(DE)
  
  DE$cluster[DE$avg_log2FC>0]<-paste0("up_",selected_cluster[1])
  DE$cluster[DE$avg_log2FC<0]<-paste0("up_",selected_cluster[2])
  
  print(table(DE$cluster))
  return(DE)
  
  }
  


PlotDEheatmap<-function(selected_cluster,cutoff,DE,srt,List){
  
    print("Plot DEA Heatmap")
    selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]
    selected_cluster<-paste0("cluster_C:",selected_cluster)

    selected_1<-List[[selected_cluster[1]]]
    selected_2<-List[[selected_cluster[2]]]
    
    
    overlap<-intersect(selected_1,selected_2)
    
    selected_1<-selected_1[!selected_1%in%overlap]
    selected_2<-selected_2[!selected_2%in%overlap]    

    top20<-DE[order(DE$avg_log2FC),]

    top20<-c(head(top20$gene,20),tail(top20$gene,20))
    top20
  
    sub<-srt[,colnames(srt)%in%c(selected_1,selected_2)]
    sub$state<-NA
    sub$state[colnames(sub)%in%c(selected_1)]<-gsub("cluster_","",selected_cluster[1])
    sub$state[colnames(sub)%in%c(selected_2)]<-gsub("cluster_","",selected_cluster[2])
    
    DE_heatmap<-list(top20,sub)
    names(DE_heatmap)<-c("top20","sub")
    return(DE_heatmap)
  }






DEGO<-function(DE,selected_cluster,cutoff,Mart,kegg_species,go_species,logfc,Pvalue,background_genes){
  
  if (!any(rownames(installed.packages()) == go_species)){
    BiocManager::install(go_species,update = F)
  }
  library(go_species,character.only = T)
  selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]
  selected_cluster<-paste0("cluster_C:",selected_cluster)

  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset(Mart, mart) 
  
  AllEnrichment<-NULL
  
  print("DEG GO & KEGG")
  
  if (length(background_genes)==0){
    
    for (cluster in unique(DE$cluster)){

      marker<-DE$gene[DE$cluster==cluster]
      Genes_set<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                       filter="external_gene_name", values= marker, uniqueRows=TRUE)
      print(head(Genes_set))
      cluster_kegg <- enrichKEGG(gene =  Genes_set$entrezgene_id,organism = kegg_species,
                                 pAdjustMethod = "BH",
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)
      
      
      if (is.null(cluster_kegg)) {
        print("no result")
      } else{
      cluster_kegg<-cluster_kegg@result
      cluster_kegg$cluster<-cluster
      cluster_kegg$class<-"KEGG"}
      
      
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
      cluster_GOBP$cluster<-cluster
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
      cluster_GOCC$cluster<-cluster
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
      cluster_GOMF$cluster<-cluster
      cluster_GOMF$class<-"GOMF"}
      
      ClusterEnrich<-rbind(cluster_GOBP,cluster_GOCC,cluster_kegg,cluster_GOMF)
      
      AllEnrichment<-rbind(AllEnrichment,ClusterEnrich)  
      #print(AllEnrichment)
    }
    
    
    
  }
    
  else{
    
    background_genes<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                         filter="external_gene_name", values=background_genes, uniqueRows=TRUE)
  
    #print(head(background_genes))
  for (cluster in unique(DE$cluster)){
  
  marker<-DE$gene[DE$cluster==cluster]

  Genes_set<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                   filter="external_gene_name", values= marker, uniqueRows=TRUE)
  
  #print(head(Genes_set))
  
  cluster_kegg <- enrichKEGG(gene =  Genes_set$entrezgene_id,organism = kegg_species,
                             pAdjustMethod = "BH",universe = as.character(background_genes$entrezgene_id),
                             minGSSize = 1,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
  
  print(cluster_kegg)
  if (is.null(cluster_kegg)) {
    print("no result")
  } else{
  
  cluster_kegg<-cluster_kegg@result
  cluster_kegg$cluster<-cluster
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
  cluster_GOBP$cluster<-cluster
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
  cluster_GOCC$cluster<-cluster
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
  cluster_GOMF$cluster<-cluster
  cluster_GOMF$class<-"GOMF"}
  
  ClusterEnrich<-rbind(cluster_GOBP,cluster_GOCC,cluster_kegg,cluster_GOMF)
  AllEnrichment<-rbind(AllEnrichment,ClusterEnrich)  
  #print(AllEnrichment)
  }
  
  }
  
  DE_Enrichment<-AllEnrichment[AllEnrichment$p.adjust<0.05,]
  #print(DE_Enrichment)
  return(DE_Enrichment)
  
}


Plot_DE_enrichment<-function(DE_Enrichment){
  
  print("Plot DEG GO & KEGG")
  print(DE_Enrichment)
  all_function<-DE_Enrichment %>%
  group_by(cluster,class) %>%
  slice_max(n = 3, order_by = Count) %>%
  slice_head(n=4)

  all_function$Gene_number <- all_function$Count
  all_function$Description<-as.factor(all_function$Description)
  all_function$'|log10(FDR)|' <- -(log10(all_function$p.adjust))
  all_function<-all_function[order(all_function$cluster),]
  all_function$Description<-paste(all_function$class,all_function$Description,sep=":")
  all_function$Description<-factor(all_function$Description,levels=unique(all_function$Description))

  # Draw the plot with ggplot2
  # to represent -log10(FDR) and Number of genes 
  # of each GO biological process per cluster 
  #---------------------------------------------------
  all_function<-all_function[!is.na(all_function$ID),]
  all_function<-all_function[all_function$p.adjust<0.05,]
  
  p<-ggplot(all_function, aes(x = Description, y = cluster)) +
    geom_point(data=all_function,aes(x=Description, y=cluster,size = Gene_number, colour = `|log10(FDR)|`), alpha=.7)+
    scale_color_gradient(low = "blue", high = "orange", limits=c(0, NA))+
    coord_flip()+
    theme_bw()+
    theme(axis.ticks.length=unit(-0.1, "cm"),
          axis.text.x = element_text(margin=margin(5,5,0,5,"pt"),angle = 45, hjust = 1,size=10),
          axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),size=10),
          axis.text = element_text(color = "black"),
          axis.title=element_blank())+
    #xlab("GO & KEGG functional enrichment")+
    labs(color="-log10(FDR)", size="Number\nof genes")
  
  return(p)
}



