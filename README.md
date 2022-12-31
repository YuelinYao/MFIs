# MFIs

## Files required:

1. **Count_matrix.csv** : count matrix of scRNA-seq data, the same file indicated in the "rawDataPath" in Abel's JSON file.
2. **Meta_data.csv**: cell type annotations from other tools (e.g., clustering, NMFs, two coloum csv, file, example in ./data folder), is only used in the overrepresentation test heatmap (heatmap tabs), if you don't want to plot this heatmap, you can skip this file.
3. **topDeviatingHOIstates.csv**: located in the output from Abel's pipeline (HOIsummaries folder)
4. **trainingData_.csv**, which is also in the output from Abel's pipeline (output folder)

## Libarary:
R:
```
library(shiny)            
library(shinythemes)      
import::from(shinycssloaders, withSpinner) 
library(reticulate)
reticulate::py_config()
library(shinyBS)         
library(shinyWidgets)     
library(gridExtra, verbose=FALSE)      
library(RColorBrewer, verbose=FALSE)
library(ComplexHeatmap)
library(Seurat)
library(stringr)
library(pheatmap)
library(dplyr)
library(rrvgo)
library(DT)
library(clusterProfiler)
library(data.table)
library(biomaRt)
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
**Share as R script:**

1. runGitHub( "MFIs", "YuelinYao")

2. Download and run locally

**Share as a web page:**

1. shinyapps.io

2. Shiny Server
 
3. RStudio Connect
