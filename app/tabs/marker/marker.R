## UI
Marker_UI <- function(){
  tagList(
    tags$h3(paste0("Find markers for a cell state"), style = "color: steelblue;"),
    plotOutput(outputId ="VolcanoPlot", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("VolcanoPlot_Download","Download as .pdf")
  )}


Marker_TableUI <- function(){
  tagList(
    tags$h3(paste0("Markers result table"), style = "color: steelblue;"),
    DTOutput(outputId ="MarkerTable", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("MarkerTable_Downloaded","Download as .csv")
  )}



MarkerGO_UI<- function(){
  tagList(
    tags$h3(paste0("GO & KEGG for Markers"), style = "color: steelblue;"),
    plotOutput(outputId ="Marker_enrichment", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("MarkerAnnotationPlot","Download as .pdf")
  )}



MarkerGOtable_UI<- function(){
  tagList(
    tags$h3(paste0("GO & KEGG for Marker table"), style = "color: steelblue;"),
    DTOutput(outputId ="Marker_enrichmentTable", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("MarkerAnnotationTable","Download as .csv")
  )}



### Input function
MarkerInput<- function(){
  tagList( 
    textInput("selected_clusterMarker", "Select cluster(s):",value = "6"),
    selectInput("Mart_Marker", "Mart dataset:", choices=datasets_list, selected = "hsapiens_gene_ensembl", multiple = FALSE),
    textInput("go_species_Marker", "GO OrgDb:",value = "org.Hs.eg.db"),
    textInput("kegg_species_Marker", "KEGG organism:",value = "hsa"),
    textInput("logfcMarker", "logFC:",value = 0.25),
    sliderInput("Pvalue_Marker",
                "Adjusted p value:",
                min = 0,
                max = 1,
                value = 0.05),
    textAreaInput("background_genesMarker", "Background genes (recommended): ",placeholder = "Just paste a list of genes (multiple-line gene list).",rows = 5),
    actionButton(inputId = "bg_Liver3",                                       #action button to display background genes
                 label = NULL, icon = icon("tag"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    bsTooltip("bg_Liver3","Load background genes as expressed genes in the current dataset.",placement = "bottom", trigger = "hover",
              options = NULL),
    actionButton("action_Marker","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}





Marker_set<-function(selected_cluster,srt,logfc,Pvalue,List){
  
  print("Find Marker genes")
  selected_cluster_list<-strsplit(selected_cluster, ",\\s*")[[1]]
  
  selected_1<-NULL
  for (i in selected_cluster_list){
  
  selected_cluster<-paste0("cluster_C:",i)
  selected_cells<-List[[selected_cluster]]
  selected_1<-unique(c(selected_1,selected_cells))
  #print(length(selected_1)) The number of cells
  }
  print(length(selected_1))
  Marker<-FindMarkers(object = srt,ident.1 = selected_1,logfc.threshold = logfc)
  print("done")
  Marker<-Marker[Marker$p_val_adj<Pvalue,]
  Marker$gene<-rownames(Marker)
  
  Marker$cluster[Marker$avg_log2FC>0]<-paste0("Up_regulated genes")
  Marker$cluster[Marker$avg_log2FC<0]<-paste0("Down_regulated genes")
  
  print(table(Marker$cluster))
  return(Marker)
  
}



PlotVolcano<-function(selected_cluster,Marker,List,GetGenesList_All){
  
  print("Plot VolcanoPlot")

  Marker[, "diff_pct"] <- Marker[, "pct.1"] - Marker[, "pct.2"]
  Marker[, "log10padj"] <- -log10(Marker[, "p_val_adj"])
  
  
  Marker<-Marker[order(Marker$avg_log2FC,decreasing = T),]
  up_genes<-Marker[Marker$cluster=="Up_regulated genes",]
  down_genes<-Marker[Marker$cluster=="Down_regulated genes",]
  top20<-c(head(up_genes$gene,10),tail(down_genes$gene,10))

  #shape = Ingroup
  Marker$Ingroup<-"Out"
  
  selected_cluster_list<-strsplit(selected_cluster, ",\\s*")[[1]]
  
  cluster_genes<-NULL
  
  for (i in selected_cluster_list){
  
  selected_cluster<-paste0("cluster_C:",i)
  names<-gsub("cluster_C:","cluster_",selected_cluster)
  cluster_genes_list<-GetGenesList_All[[names]]
  cluster_genes<-unique(c(cluster_genes,cluster_genes_list))
  }
  
  print(length(cluster_genes))
  Marker$Ingroup<-1
  Marker$Ingroup[Marker$gene%in%cluster_genes]<-"In"
  
  #label
  Marker$label<-NA
  Marker$label[Marker$gene%in%top20]<-Marker$gene[Marker$gene%in%top20]
  p<-ggplot(data=Marker, aes(x=avg_log2FC, y=log10padj, col=diff_pct,label=label,shape = Ingroup)) + 
    geom_point() + labs(x="log2FC",y="-log10(p-adjust)")+guides(shape = FALSE)+
    theme_bw() +scale_color_gradient2(name="Diff_pct",midpoint=0,low="#053061",mid = "white",high = "#67001f") +  geom_text_repel(max.overlaps = 100)
  #scale_color_gradientn(name="Diff_pct",colors = rev(brewer.pal(11 ,"RdBu")))
  #scale_color_gradient2(name="Diff_pct",midpoint=0,low="#053061",mid = "white",high = "#67001f")
  return(p)
}






MarkerGO<-function(Marker,Mart,kegg_species,go_species,logfc,Pvalue,background_genes){
  
  print("Markers GO & KEGG")
  if (!any(rownames(installed.packages()) == go_species)){
    BiocManager::install(go_species,update = F)
  }
  library(go_species,character.only = T)
  #print(class(selected_cluster))

  
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset(Mart, mart) 
  
  AllEnrichment<-NULL
  
  print("Marker GO & KEGG")
  
  if (length(background_genes)==0){
    
    for (cluster in unique(Marker$cluster)){
      
      marker<-Marker$gene[Marker$cluster==cluster]
      
      Genes_set<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                       filter="external_gene_name", values= marker, uniqueRows=TRUE)
      
      cluster_kegg <- enrichKEGG(gene =  Genes_set$entrezgene_id,organism = kegg_species,
                                 pAdjustMethod = "BH",
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)
      
      
      if (is.null(cluster_kegg)){
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
      } else {
      
      cluster_GOBP<-cluster_GOBP@result
      cluster_GOBP$cluster<-cluster
      cluster_GOBP$class<-"GOBP" }
      
      
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
      cluster_GOCC$class<-"GOCC" }
      
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
      cluster_GOMF$class<-"GOMF" }
      
      
      ClusterEnrich<-rbind(cluster_GOBP,cluster_GOCC,cluster_kegg,cluster_GOMF)
      AllEnrichment<-rbind(AllEnrichment,ClusterEnrich)  }
      
     
      
    }
    
  
  else{
    
    background_genes<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                            filter="external_gene_name", values=background_genes, uniqueRows=TRUE)
    
    for (cluster in unique(Marker$cluster)){
      
      marker<-Marker$gene[Marker$cluster==cluster]
      
      Genes_set<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                       filter="external_gene_name", values= marker, uniqueRows=TRUE)
      
      
      
      cluster_kegg <- enrichKEGG(gene =  Genes_set$entrezgene_id,organism = kegg_species,
                                 pAdjustMethod = "BH",universe = as.character(background_genes$entrezgene_id),
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)
      
      
      if (is.null(cluster_kegg)) {
        print("no result")
      }  else {
      cluster_kegg<-cluster_kegg@result
      cluster_kegg$cluster<-cluster
      cluster_kegg$class<-"KEGG" }
      
      
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
      } else {
      
      cluster_GOBP<-cluster_GOBP@result
      cluster_GOBP$cluster<-cluster
      cluster_GOBP$class<-"GOBP" }
      
      
      cluster_GOCC <- enrichGO(gene = Genes_set$entrezgene_id,universe = as.character(background_genes$entrezgene_id),
                               OrgDb= go_species,
                               ont = "CC",
                               pAdjustMethod = "BH",
                               minGSSize = 1,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE)
      
      if (is.null( cluster_GOCC)) {
        print("no result")
      } else {
        
      cluster_GOCC<-cluster_GOCC@result
      cluster_GOCC$cluster<-cluster
      cluster_GOCC$class<-"GOCC" }
      
      cluster_GOMF <- enrichGO(gene = Genes_set$entrezgene_id,universe = as.character(background_genes$entrezgene_id),
                               OrgDb= go_species,
                               ont = "MF",
                               pAdjustMethod = "BH",
                               minGSSize = 1,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE)
      
      if (is.null(cluster_GOMF)){
        print("no result")
      } else {
      cluster_GOMF<-cluster_GOMF@result
      cluster_GOMF$cluster<-cluster
      cluster_GOMF$class<-"GOMF" }
   
      ClusterEnrich<-rbind(cluster_GOBP,cluster_GOCC,cluster_kegg,cluster_GOMF)
      AllEnrichment<-rbind(AllEnrichment,ClusterEnrich)  
    }
   
  }
  
  Marker_Enrichment<-AllEnrichment[AllEnrichment$p.adjust<0.05,]
  return(Marker_Enrichment)
  
}


Plot_Marker_enrichment<-function(Marker_Enrichment){
  
  print("Plot Marker GO & KEGG")
  all_function<-Marker_Enrichment %>%
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



