# MFIs

## Files required:

1. **Count_matrix.csv** : count matrix of scRNA-seq data, it can be the same file indicated in the "rawDataPath" in Stator nextflow JSON file. It is used for finding cells that satisfy the d-tuples within a group and differential expression analysis. 
2. **Meta_data.csv** : cell's external annotations (e.g., from other tools clustering, NMFs, or experimental conditions), two coloum csv file, example in ./data folder), is only used in the overrepresentation test heatmap (heatmap-cell tabs), if you don't want to plot this heatmap, you can also skip this file.
3. **all_DTuples.csv** : located in the output from Stator nextflow pipeline (dtuples_output folder).
4. **trainingData_.csv** : which is also in the output from Stator nextflow pipeline (output folder).
5. **GeneAnnotationSet.csv** (optional): which is used in heatmap-genes, over-representation test between gene list. Example file format see ./data/CancerState.csv.
6. **UMAP.csv** (optional): which is used in UMAP Plot, Example file format see ./data/UMAP_coords.csv.
7. **MCMCgraph.csv** (optional): which is used in Markov Blanket Tab, located in the output folder from Stator nextflow pipeline. Example file format see ./data/MCMCgraph_14698Cells_1000Genes.csv.
8. **Binary_matrix.csv** (optional): upload your own binary matrix if needed. The same format as Count_matrix.csv.
## Libarary:
**Note:** Please use the following versions: R version 4.3.0, dbplyr_2.3.3 (by devtools::install_version("dbplyr", version = "2.3.3")), clusterProfiler_4.8.2 and Seurat_4.3.0.1 For faster version of DE analysis, update Seurat to V5 and install presto (see: https://github.com/immunogenomics/presto).

R:
```
#----- Install one by one...

if (!any(rownames(installed.packages()) == "shiny")){
  install.packages("shiny")
}
library(shiny)

if (!any(rownames(installed.packages()) == "shinythemes")){
  install.packages("shinythemes")
}
library(shinythemes)

if (!any(rownames(installed.packages()) == "arules")){
  install.packages("arules")
}
library(arules)

if (!any(rownames(installed.packages()) == "shinycssloaders")){
  install.packages("shinycssloaders")
}
library(shinycssloaders)

if (!any(rownames(installed.packages()) == "import")){
  install.packages("import")
}

import::from(shinycssloaders, withSpinner) 

if (!any(rownames(installed.packages()) == "reticulate")){
  install.packages("reticulate")
}
library(reticulate)
reticulate::py_config()

if (!any(rownames(installed.packages()) == "shinyBS")){
  install.packages("shinyBS")
}
library(shinyBS)

if (!any(rownames(installed.packages()) == "shinyWidgets")){
  install.packages("shinyWidgets")
}
library(shinyWidgets)

if (!any(rownames(installed.packages()) == "gridExtra")){
  install.packages("gridExtra")
}
library(gridExtra, verbose=FALSE) 

if (!any(rownames(installed.packages()) == "RColorBrewer")){
  install.packages("RColorBrewer")
}
library(RColorBrewer, verbose=FALSE)

if (!any(rownames(installed.packages()) == "ComplexHeatmap")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)

if (!any(rownames(installed.packages()) == "Seurat")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("Seurat")
}
library(Seurat)

if (!any(rownames(installed.packages()) == "stringr")){
  install.packages("stringr")
}
library(stringr)

if (!any(rownames(installed.packages()) == "pheatmap")){
  install.packages("pheatmap")
}
library(pheatmap)

if (!any(rownames(installed.packages()) == "dplyr")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("dplyr")
}
library(dplyr)

if (!any(rownames(installed.packages()) == "rrvgo")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("rrvgo")
}
library(rrvgo)

if (!any(rownames(installed.packages()) == "limma")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("limma")
}
library("limma")

if (!any(rownames(installed.packages()) == "org.Hs.eg.db")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}

if (!any(rownames(installed.packages()) == "org.Mm.eg.db")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("org.Mm.eg.db")
}

library(org.Hs.eg.db)
library(org.Mm.eg.db)


if (!any(rownames(installed.packages()) == "DT")){
  install.packages("DT")
}
library(DT)

if (!any(rownames(installed.packages()) == "clusterProfiler")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}
library(clusterProfiler)

if (!any(rownames(installed.packages()) == "data.table")){
  install.packages("data.table")
}
library(data.table)

if (!any(rownames(installed.packages()) == "biomaRt")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

if (!any(rownames(installed.packages()) == "ggplot2")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ggplot2")
}
library(ggplot2)


install.packages('heatmaply')
library(heatmaply)

install.packages("plotly")
library(plotly)


#------ alternatively: renv package: https://rstudio.github.io/renv/articles/renv.html
# uncomment "source("renv/activate.R")" (remove # ) in the .Rprofile file if you wish to use renv
# renv::restore() use the information from renv.lock file to retrieve and reinstall those packages in this project.
# just enter the following in the project directory:
renv::restore()



```


Install the Python packages by following R codes:

```
# Install python from https://www.python.org/downloads/
# Python 3.8
# Python packages: virtualenv and numpy
# pip install virtualenv
# pip install numpy
# https://github.com/ranikay/shiny-reticulate-app

# Define any Python packages needed for the app in R:
# PYTHON_DEPENDENCIES = c('pip', 'numpy','pandas','igraph','argparse','scipy',
#           'matplotlib','Pillow','seaborn')
PYTHON_DEPENDENCIES = c('numpy','pandas','scipy','scikit-learn') # The name of sklearn has been change to scikit-learn



# ------------------ App virtualenv setup (Do not edit) ------------------- #
# VIRTUALENV_NAME and PYTHON_PATH are definded in .Rprofile
virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
python_path = Sys.getenv('PYTHON_PATH')


# Create virtual env and install dependencies
reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES, ignore_installed=TRUE)
reticulate::use_virtualenv(virtualenv_dir, required = T)
```


## How to run:
**Share as R script with R studio:**

Once all the packages have been installed, 

Approach 1: Enter: ```shiny::runGitHub( "MFIs", "YuelinYao")```

Approach 2: Download and run locally: ```shiny::runApp()```


Or you can run with Docker:

i). Install docker on the local computer
  
ii). Git clone the repository, Local branch: 
   
    ```
    git clone https://github.com/YuelinYao/MFIs.git
    ```
iii). Enter the repository:
    
     ```
     docker build -t my-shiny-app .
     docker run --rm -p 3838:3838 my-shiny-app
     ```
For larger files go to Docker Desktop: Settings > Resources > Update the RAM to be more.

Dockerhub: https://hub.docker.com/r/yuelinyao120/stator-app

**Share as a web page:**


Shiny Server: https://shiny.igc.ed.ac.uk/MFIs/


## Notes
KEGG reference genome: https://www.genome.jp/kegg/catalog/org_list.html

GO reference genome: http://bioconductor.org/packages/release/BiocViews.html#___OrgDb

Ensembl dataset for more details: see ./data/ensembl_dataset.csv


 
