# MFIs

## Files required:

1. **Count_matrix.csv** : count matrix of scRNA-seq data, the same file indicated in the "rawDataPath" in Abel's JSON file.
2. **Meta_data.csv** : cell type annotations from other tools (e.g., clustering, NMFs, two coloum csv file, example in ./data folder), is only used in the overrepresentation test heatmap (heatmap-cell tabs), if you don't want to plot this heatmap, you can skip this file.
3. **topDeviatingHOIstates.csv** : located in the output from Abel's pipeline (HOIsummaries folder)
4. **trainingData_.csv** : which is also in the output from Abel's pipeline (output folder)
5. **GeneAnnotationSet.csv** (optional): which is used in heatmap-genes, over-representation test between gene list. Example file format see ./data/CancerState.csv 

## Libarary:
R:
```
#------renv package: https://rstudio.github.io/renv/articles/renv.html
# renv::restore() use the information from renv.lock file to retrieve and reinstall those packages in this project.
# just enter the following in the project directory:
renv::restore()


#----- alternatively, install one by one...

if (!any(rownames(installed.packages()) == "shiny")){
  install.packages("shiny")
}
library(shiny)

if (!any(rownames(installed.packages()) == "shinythemes")){
  install.packages("shinythemes")
}
library(shinythemes)


if (!any(rownames(installed.packages()) == "shinycssloaders")){
  install.packages("shinycssloaders")
}
library(shinycssloaders)


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
```


Python:
```
# Define any Python packages needed for the app in R:
#PYTHON_DEPENDENCIES = c('pip', 'numpy','pandas','igraph','argparse','scipy',
             #           'matplotlib','Pillow','seaborn')
PYTHON_DEPENDENCIES = c('numpy','pandas','scipy')



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

1. Enter: ```shiny::runGitHub( "MFIs", "YuelinYao")```

2. Download and run locally: ```shiny::runApp()```

3. Run with Docker:

    i). Install docker on the local computer
  
    ii). Git clone the repository: 
    ```
    git clone https://github.com/YuelinYao/MFIs.git
    ```
    iii). Enter the repository:
    
     ```
     docker build -t my-shiny-app .
     docker run --rm -p 3838:3838 my-shiny-app
     ```
     For larger files go to Docker Desktop: Settings > Resources > Update the RAM to be more.

**Share as a web page:**

1. shinyapps.io: https://yuelinyao120.shinyapps.io/MFIs-shinyapps/ (Server is provided by Rstudio, but the memory for free is not enough.)

2. Shiny Server (We need to provide a server to host this application, might consider put on: https://shiny.igc.ed.ac.uk/)


## Notes
KEGG reference genome: https://www.genome.jp/kegg/catalog/org_list.html

GO reference genome: http://bioconductor.org/packages/release/BiocViews.html#___OrgDb

Ensembl dataset for more details: see ./data/ensembl_dataset.csv


 
