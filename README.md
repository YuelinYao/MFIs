# MFIs

## Files required:

1. **Count_matrix.csv** : count matrix of scRNA-seq data, the same file indicated in the "rawDataPath" in Abel's JSON file.
2. **Meta_data.csv**: cell type annotations from other tools (e.g., clustering, NMFs, two coloum csv, file, example in ./data folder), is only used in the overrepresentation test heatmap (heatmap tabs), if you don't want to plot this heatmap, you can skip this file.
3. **topDeviatingHOIstates.csv**: located in the output from Abel's pipeline (HOIsummaries folder)
4. **trainingData_.csv**, which is also in the output from Abel's pipeline (output folder)

## Libarary:
R:
```
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
import sys
import pandas as pd
import numpy as np
import igraph as ig
import argparse
import sys
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cut_tree
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.offsetbox import OffsetImage,AnnotationBbox
import io
from PIL import Image
import seaborn as sns
sns.set(style='white')
sns.set_palette('colorblind')
```
## How to run:
**Share as R script with R studio:**

1. Enter runGitHub( "MFIs", "YuelinYao")

2. Download and run locally

**Share as a web page:**

1. shinyapps.io

2. Shiny Server
 
3. RStudio Connect
