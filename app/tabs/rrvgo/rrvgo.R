rrvgoUI <- function(){
  tagList(
    tags$h3(paste0("Word cloud"), style = "color: steelblue;"),
    plotOutput(outputId ="SimplifyingGOBP_wordcloud", width = "90%") %>% withSpinner(color="#4682B4")
  )}


rrvgo2UI <- function(){
  tagList(
    tags$h3(paste0("Treemap plot"), style = "color: steelblue;"),
    plotOutput(outputId ="SimplifyingGOBP_treemap", width = "90%") %>% withSpinner(color="#4682B4")
  )}


rrvgo3UI <- function(){
  tagList(
    tags$h3(paste0("Scatter plot"), style = "color: steelblue;"),
    plotOutput(outputId ="SimplifyingGOBP_Scatter", width = "90%") %>% withSpinner(color="#4682B4")
  )}




rrvgo_tableUI<-function(){
  
  tagList(
    tags$h3(paste0("Reduced Terms"), style = "color: steelblue;"),
    DTOutput(outputId ="Reduced_Terms", width = "80%") %>% withSpinner(color="#4682B4")
  )}


rrvgo_FunctionE <- function(cutoff,selected_cluster,subclass,go_species,Mart,GenesList,backgroud_genes){
  

  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset(Mart, mart)
  


  i=selected_cluster
  names<-paste0("cluster_",i)
  Genes_set<-GenesList[[names]]
  print(paste0("Cluster: ",i))
  
  Genes_set<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                   filter="external_gene_name", values=Genes_set, uniqueRows=TRUE)
  
  print("GO analysis")  
  
  if (length(backgroud_genes)==0){
    print("No user-defined background genes")
    cluster_GO <- enrichGO(gene = Genes_set$entrezgene_id,
                           OrgDb= go_species,
                           ont = subclass,
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
  } else{
  
  print("Use input background genes")
  
  backgroud_genes<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                           filter="external_gene_name", values=backgroud_genes, uniqueRows=TRUE)
  cluster_GO <- enrichGO(gene = Genes_set$entrezgene_id,universe = as.character(backgroud_genes$entrezgene_id),
                               OrgDb= go_species,
                               ont = subclass,
                               pAdjustMethod = "BH",
                               minGSSize = 1,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE)
      

      }

  
  if (dim( cluster_GO )[1]==0) {
    print("no result")
  }
        
  cluster_GO <- cluster_GO@result
  cluster_GO$cluster<-i
  cluster_GO$class<-subclass  
  cluster_GO<-cluster_GO[cluster_GO$p.adjust<0.05,]
  
      

    print("Simplify GO terms")  
      simMatrix <- calculateSimMatrix(cluster_GO$ID,
                                      orgdb=go_species,
                                      ont=subclass,
                                      method="Rel")
      
      scores <- setNames(-log10(cluster_GO$qvalue), cluster_GO$ID)
      
      reducedTerms <- reduceSimMatrix(simMatrix,
                                      scores,
                                      threshold=0.7, # this threhold might adjust later
                                      orgdb=go_species)
  
  
      
  
  reducedTerms<-list(reducedTerms,simMatrix)
  names(reducedTerms)<-c("reducedTerms","simMatrix")
  return(reducedTerms)
  
}





