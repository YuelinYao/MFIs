## UI
rrvgoUI <- function(){
  tagList(
    tags$h3(paste0("Word cloud"), style = "color: steelblue;"),
    plotOutput(outputId ="SimplifyingGO_wordcloud", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("rrvgo_plot1","Download as .pdf")
  )}


rrvgo2UI <- function(){
  tagList(
    tags$h3(paste0("Treemap plot"), style = "color: steelblue;"),
    plotOutput(outputId ="SimplifyingGO_treemap", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("rrvgo_plot2","Download as .pdf")
  )}


rrvgo3UI <- function(){
  tagList(
    tags$h3(paste0("Scatter plot"), style = "color: steelblue;"),
    plotOutput(outputId ="SimplifyingGO_Scatter", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("rrvgo_plot3","Download as .pdf")
  )}

rrvgo_tableUI<-function(){
  tagList(
    tags$h3(paste0("Reduced Terms"), style = "color: steelblue;"),
    DTOutput(outputId ="Reduced_Terms", width = "80%") %>% withSpinner(color="#4682B4"),
    downloadButton("rrvgo_table","Download as .csv")
  )}



## Input function
rrvgoInput<- function(){
  tagList( 
    textInput("selected_clusterrrvgo", "Input a cluster:",value = "5"),
    selectInput("Mart_rrgvo", "Mart dataset:", choices=datasets_list, selected = "hsapiens_gene_ensembl", multiple = FALSE),
    textInput("go_species_rrgvo", "GO OrgDb:",value = "org.Hs.eg.db"),
    radioButtons("subOntology", "Select sub-ontology:",
                 c("Biological Process" = "BP",
                   "Cellular Component" = "CC",
                   "Molecular Function" = "MF"),selected = "BP" ),
    textAreaInput("background_genesrrgvo", "Background genes (recommended): ",placeholder = "Just paste a list of genes (multiple-line gene list).",rows = 5),
    actionButton(inputId = "bg_Liver",                                       #action button to display background genes
                 label = NULL, icon = icon("tag"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    bsTooltip("bg_Liver","Load background genes as expressed genes in the current dataset.",placement = "bottom", trigger = "hover",
              options = NULL),
    actionButton("action_rrvgo","Submit",icon("paper-plane"),       # submit button
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}







rrvgo_FunctionE <- function(selected_cluster,subclass,go_species,Mart,GenesList,background_genes){
  
  if (!any(rownames(installed.packages()) == go_species)){
    BiocManager::install(go_species,update = F)
  }
  
  library(go_species,character.only = T)
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset(Mart, mart)
  
  i=selected_cluster
  names<-paste0("cluster_",i)
  Genes_set<-GenesList[[names]]
  print(paste0("Cluster: ",i))
  print(length(Genes_set))
  
  if (length(Genes_set)>0) {
  Genes_set<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                   filter="external_gene_name", values=Genes_set, uniqueRows=TRUE)
  
  print("GO analysis")  
  
  if (length(background_genes)==0){
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
  background_genes<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                           filter="external_gene_name", values=background_genes, uniqueRows=TRUE)
  cluster_GO <- enrichGO(gene = Genes_set$entrezgene_id,universe = as.character(background_genes$entrezgene_id),
                               OrgDb= go_species,
                               ont = subclass,
                               pAdjustMethod = "BH",
                               minGSSize = 1,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05,
                               readable = TRUE)
      

      }

  
  if (is.null(cluster_GO)) {
    print("no result")
  } else{
  cluster_GO <- cluster_GO@result
  cluster_GO$cluster<-i
  cluster_GO$class<-subclass  }
  
  
  cluster_GO<-cluster_GO[cluster_GO$p.adjust<0.05,]
  
  print("Simplify GO terms")  
  simMatrix <- calculateSimMatrix(cluster_GO$ID,
                                      orgdb=go_species,
                                      ont=subclass,
                                      method="Rel")
    
  
  if(any(is.na(cluster_GO$qvalue))){
    cluster_GO$qvalue<-cluster_GO$p.adjust
  }
  
    scores <- setNames(-log10(cluster_GO$qvalue), cluster_GO$ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                      scores,
                                      threshold=0.7, # this threhold might adjust later
                                      orgdb=go_species)
  
  reducedTerms<-list(reducedTerms,simMatrix)
  names(reducedTerms)<-c("reducedTerms","simMatrix")
  return(reducedTerms)
  }
}





