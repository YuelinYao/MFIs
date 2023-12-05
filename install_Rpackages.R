# install R packages for the shiny app
package_list_cran <- c(
    "shiny", "shinythemes", "shinycssloaders", "import",
    "reticulate", "shinyBS", "shinyWidgets", "gridExtra",
    "RColorBrewer", "stringr", "pheatmap","DT", 
    "data.table", "heatmaply", "plotly", "Seurat")
package_list_bioc <- c(
  "ComplexHeatmap", "dplyr", "rrvgo", "limma",
  "org.Hs.eg.db", "org.Mm.eg.db", "clusterProfiler", "biomaRt", "ggplot2")

for (package in package_list_cran) {
  if (!any(rownames(installed.packages()) == package)) {
    install.packages(package)
  }
}

for (package in package_list_bioc) {
  if (!any(rownames(installed.packages()) == package)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(package, update = TRUE, ask = FALSE)
  }
}

# load packages
for (package in c(package_list_cran, package_list_bioc)) {
  library(package, character.only = TRUE)
}

# Create virtual env and install dependencies
reticulate::virtualenv_create(envname = "example_env_name", python = "python3")
reticulate::virtualenv_install("example_env_name", packages = c("numpy","pandas","scipy","scikit-learn"), ignore_installed=TRUE)
reticulate::use_virtualenv("example_env_name", required = T)

# run the app
#shiny::runApp('app.R', host = '0.0.0.0', port = 3838)