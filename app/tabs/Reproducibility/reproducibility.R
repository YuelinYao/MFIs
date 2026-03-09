## UI ####
reproducibilityUI <- function(){
  tagList(
    tags$h3("Reproducibility", style = "color: steelblue;"),
    tags$p(
      "Stator allows visualisation of reproducibility across two datasets using an overrepresentation test."
    ),
    tags$p(
      "Please upload the count matrix files for Dataset 1 and Dataset 2. ",
      "The required format is the same as the count matrix used in the other tabs."
    ),
    tags$p(
      "Please also upload D-tuple summary table 1 and D-tuple summary table 2. ",
      "The required format is the same as the table downloaded from the Table tab."
    ),
    tags$p(
      "If a D-tuple summary table contains an 'Annotation' column, the heatmap will include annotations. ",
      "Otherwise, the heatmap will be shown without annotation."
    ),
    br(),
    
    tags$h3("Reproducibility plot 1", style = "color: steelblue;"),
    tags$h5("Apply Stator states 1 and 2 to Dataset 1"),
    plotOutput(outputId = "Reproducibility_plot1", width = "90%") %>% 
      withSpinner(color = "#4682B4"),
    downloadButton("download_reproducibility_plot1", "Download as .pdf"),
    br(),
    br(),
    
    tags$h3("Reproducibility plot 2", style = "color: steelblue;"),
    tags$h5("Apply Stator states 1 and 2 to Dataset 2"),
    plotOutput(outputId = "Reproducibility_plot2", width = "90%") %>% 
      withSpinner(color = "#4682B4"),
    downloadButton("download_reproducibility_plot2", "Download as .pdf")
  )
}



reproducibilityInput <- function() {
  tagList(
    # Description
    div(
      p("Here, stator states are validated across different datasets."),
      style = "margin-bottom: 20px;"
    ),
    
    # Count Matrix 1
    div(
      fileInput(
        inputId = "reproduce_count_matrix_1",
        label = "Upload count matrix for dataset 1:",
        multiple = FALSE,
        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
      ),
      style = "display: inline-block; vertical-align: top; width: 70%;"
    ),
    
    # D-tuples Table 1
    div(
      fileInput(
        inputId = "reproduce_dtuples_table_1",
        label = "Upload dtuple table for dataset 1: ",
        multiple = FALSE,
        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
      ),
      style = "display: inline-block; vertical-align: top; width: 70%;"
    ),
    
    # Count Matrix 2
    div(
      fileInput(
        inputId = "reproduce_count_matrix_2",
        label = "Upload count matrix for dataset 2: ",
        multiple = FALSE,
        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
      ),
      style = "display: inline-block; vertical-align: top; width: 70%;"
    ),
    
    
    # D-tuples Table 2
    div(
      fileInput(
        inputId = "reproduce_dtuples_table_2",
        label = "Upload dtuple table for dataset 2: ",
        multiple = FALSE,
        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
      ),
      style = "display: inline-block; vertical-align: top; width: 70%;"
    ),
    
    # Submit Button
    br(),
    actionButton(
      inputId = "action_reproducibility",
      label = "Submit",
      icon = icon("paper-plane"),
      style = "color: #fff; background-color: #337ab7; border-color: #2e6da4;"
    )
  )
}

### function
## function
# N= the total number of background
heatmap_test <- function(RefSet,TestSet,N=25678,test="Over_representation") {
  
  
  r=length(names(RefSet))
  col=length(names(TestSet))
  
  Fold<-array(data=NA,dim = c(r,col))
  Data_mtrix<-array(data=NA,dim = c(r,col))
  
  colnames(Data_mtrix)<-names(TestSet)
  rownames(Data_mtrix)<-names(RefSet)
  
  dimnames(Fold)<-dimnames(Data_mtrix)
  
  Overlap<-array(data=NA,dim = c(r,col))
  dimnames(Overlap)<-dimnames(Data_mtrix)
  
  for (i in names(TestSet)){
    #print(i)
    Genes<-TestSet[[i]]
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
      #print(p_value)
      P_set[ct]=p_value
      enrichment_set[ct]=fold_value
    }
    Data_mtrix[,i]<-P_set
    Fold[,i]<-enrichment_set
  }
  
  result=list()
  result$raw_pvalue<-Data_mtrix
  
  #Data_mtrix[Data_mtrix==0]<-2.2e-16
  Data_mtrix<-Data_mtrix+2.2e-16
  Mydata_raw_FDR <- p.adjust(Data_mtrix,method = "BH")
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  log10FDR<--log10(Mydata_raw_m)
  
  #print(log10FDR)
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  
  Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Data_mtrix)[1],byrow = F)
  
  #breaksList = seq(0, 10, by = 1)
  dimnames(log10FDR)<-dimnames(Data_mtrix)
  
  df<-lengths(RefSet)
  percentage=df/sum(df)
  
  ##
  percentage_row<-percentage
  names(percentage_row)<-names(df)
  percentage_row<-percentage_row[names(RefSet)]
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(RefSet)]
  stats<-paste0(names(RefSet),stats)
  rownames(log10FDR)<-stats
  
  df<-lengths(TestSet)
  percentage=df/sum(df)
  
  ##
  percentage_col<-percentage
  names(percentage_col)<-names(df)
  percentage_col<-percentage_col[names(TestSet)]
  
  stats<-paste0(" (",df,", ",round(percentage*100,2),"%",")")
  names(stats)<-names(df)
  stats<-stats[names(TestSet)]
  stats<-paste0(names(TestSet),stats)
  colnames(log10FDR)<-stats
  colnames(log10FDR)<-gsub(pattern = "cluster_","",colnames(log10FDR))
  
  dimnames(Fold)<-dimnames(log10FDR)
  
  
  result$log10FDR=log10FDR
  result$Mydata_raw_m=Mydata_raw_m
  Fold[Fold==0]<-2.2e-16
  Fold[Fold==Inf]<-100
  result$Fold=log10(Fold)
  result$row=as.vector(percentage_row*100)
  result$col=as.vector(percentage_col*100)
  dimnames(Overlap)<-dimnames(log10FDR)
  result$Overlap<-as.matrix(Overlap)
  #print(result)
  print("done")
  
  return(result)
  
}


create_heatmap_plot <- function(results_data, dataset_name, row_dtuple, col_dtuple,
                                cellwidth = NULL, cellheight = NULL) {
  
  df_row <- row_dtuple[!duplicated(row_dtuple$Cluster), , drop = FALSE]
  df_col <- col_dtuple[!duplicated(col_dtuple$Cluster), , drop = FALSE]
  
  has_row_annotation <- ("Annotation" %in% colnames(df_row) &&
                           !all(is.na(df_row$Annotation)))
  has_col_annotation <- ("Annotation" %in% colnames(df_col) &&
                           !all(is.na(df_col$Annotation)))
  
  if (has_row_annotation) {
    df_row <- df_row[, c("Cluster", "Annotation"), drop = FALSE]
    df_row$Annotation <-factor(df_row$Annotation,levels = unique(df_row$Annotation))
  }
  
  if (has_col_annotation) {
    df_col <- df_col[, c("Cluster", "Annotation"), drop = FALSE]
    df_col$Annotation <-factor(df_col$Annotation,levels = unique(df_col$Annotation))
  }
  
  ha_row <- NULL
  if (has_row_annotation) {
    ha_row <- rowAnnotation(
      empty = anno_empty(border = TRUE, width = unit(0.01, "cm")),
      foo = anno_block(
        gp = gpar(col = "white"),
        labels = unique(df_row$Annotation),
        labels_gp = gpar(col = "black", fontsize = 12),
        width = unit(4, "cm"),
        height = unit(4, "cm"),
        labels_rot = 0
      )
    )
  }
  
  ha <- NULL
  if (has_col_annotation) {
    ha <- HeatmapAnnotation(
      empty = anno_empty(border = TRUE, height = unit(0.01, "cm")),
      foo = anno_block(
        gp = gpar(col = "white"),
        labels = unique(df_col$Annotation),
        labels_gp = gpar(col = "black", fontsize = 12),
        width = unit(4, "cm"),
        height = unit(4, "cm"),
        labels_rot = 90
      )
    )
  }
  
  pheatmap_args <- list(
    mat = results_data[["log10FDR"]],
    border_color = NA,
    display_numbers = results_data[["Mydata_raw_m"]],
    top_annotation = ComplexHeatmap::columnAnnotation(
      Dtuple1 = ComplexHeatmap::anno_barplot(
        border = FALSE,
        results_data[["col"]]
      )
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Dtuple2 = ComplexHeatmap::anno_barplot(
        border = FALSE,
        results_data[["row"]]
      )
    ),
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    fontsize = 12,
    fontsize_number = 15,
    name = paste0(dataset_name, ": -Log10FDR"),
    heatmap_legend_param = list(title_position = "lefttop-rot"),
    breaks = seq(0, 10, by = 0.1),
    
    # key logic here
    cluster_cols = !has_col_annotation,
    cluster_rows = !has_row_annotation,
    
    treeheight_row = 0,
    treeheight_col = 0,
    color = colorRampPalette(c("white", "firebrick3"))(100)
  )
  
  # column annotation goes to bottom
  if (!is.null(ha)) {
    pheatmap_args$bottom_annotation <- ha
  }
  
  # row annotation goes to right
  if (!is.null(ha_row)) {
    pheatmap_args$right_annotation <- ha_row
  }
  
  if (has_col_annotation) {
    pheatmap_args$column_split <- df_col$Annotation
    pheatmap_args$column_title_gp <- gpar(fontsize = 0)
  }
  
  if (has_row_annotation) {
    pheatmap_args$row_split <- df_row$Annotation
    pheatmap_args$row_title_gp <- gpar(fontsize = 0)
  }
  
  if (!is.null(cellwidth)) {
    pheatmap_args$cellwidth <- cellwidth
  }
  
  if (!is.null(cellheight)) {
    pheatmap_args$cellheight <- cellheight
  }
  
  do.call(ComplexHeatmap::pheatmap, pheatmap_args)
}

# Calculate heatmap dimensions for PDF export
calculate_heatmap_dimensions <- function(results_data) {
  log10fdr <- results_data[["log10FDR"]]
  
  width <- ncol(log10fdr) * 5 / 25.4 +
    max(nchar(rownames(log10fdr))) * 5 / 25.4 + 2
  
  height <- nrow(log10fdr) * 5 / 25.4 +
    max(nchar(colnames(log10fdr))) * 5 / 25.4 + 5
  
  return(list(width = width, height = height))
}


GetCellList_reproduce <- function(count, summaryTable) { 
  print("Get cells in each state")
  
  Devstates <- summaryTable
  Devstates$D.tuple <- gsub("\\[|\\]", "", Devstates$D.tuple)
  
  # build normalized gene-name lookup
  orig_colnames <- colnames(count)
  norm_colnames <- gsub("-", ".", orig_colnames, fixed = TRUE)
  gene_map <- setNames(orig_colnames, norm_colnames)
  
  count_names <- rownames(count)
  List <- list()
  
  for (i in unique(Devstates$Cluster)) {
    cluster <- Devstates[Devstates$Cluster == i, , drop = FALSE]
    Cell <- character(0)
    
    for (c in seq_len(nrow(cluster))) {
      # split gene tuple
      Genes <- strsplit(cluster$Genes[c], "_", fixed = TRUE)[[1]]
      
      # normalize tuple gene names so NR2F2-AS1 matches NR2F2.AS1
      Genes_norm <- gsub("-", ".", Genes, fixed = TRUE)
      Genes_mapped <- unname(gene_map[Genes_norm])
      
      # parse state string
      state <- as.character(cluster$D.tuple[c])
      State <- as.numeric(strsplit(state, "")[[1]])
      
      # check tuple length
      if (length(Genes_mapped) != length(State)) {
        warning(
          sprintf(
            "Skipping Cluster '%s', entry '%s / %s': number of genes (%d) != number of states (%d)",
            i, cluster$Genes[c], cluster$D.tuple[c], length(Genes_mapped), length(State)
          )
        )
        next
      }
      
      # check missing genes
      if (any(is.na(Genes_mapped))) {
        missing_genes <- Genes[is.na(Genes_mapped)]
        warning(
          sprintf(
            "Skipping Cluster '%s', entry '%s / %s': gene(s) not found in count: %s",
            i, cluster$Genes[c], cluster$D.tuple[c], paste(missing_genes, collapse = ", ")
          )
        )
        next
      }
      
      # match cells for any tuple length
      idx <- rep(TRUE, nrow(count))
      for (j in seq_along(Genes_mapped)) {
        idx <- idx & (count[, Genes_mapped[j]] == State[j])
      }
      
      Cell <- c(Cell, count_names[idx])
    }
    
    List[[as.character(i)]] <- unique(Cell)
  }
  
  return(List)
}


GetCellList_reproduce_v0 <- function(count,summaryTable) { 

  print("Get cells in each state")
  # print(head(summaryTable))
  Devstates <-summaryTable
  Devstates$D.tuple<-gsub('\\[',"",Devstates$D.tuple)
  Devstates$D.tuple<-gsub('\\]',"",Devstates$D.tuple)
  
  List=NULL
  for (i in unique(Devstates$Cluster)){
    
    cluster<-Devstates[Devstates$Cluster==i,]
    Cell<-NULL
    for (c in 1:length(cluster$Genes)){
      #print(paste0("Cell states: ",cluster$genes[c]))
      Genes<-cluster$Genes[c]
      Genes<-str_split(Genes, "_")[[1]]
      
      state<-cluster$D.tuple[c]
      
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
    names(list)<-i
    List<-c(List,list)
    
  }
  # print(head(List))
  return(List)
}

