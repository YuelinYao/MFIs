# MFIs

## Files required:

1. **Count_matrix.csv** : count matrix of scRNA-seq data, the same file indicated in the "rawDataPath" in Abel's JSON file.
2. **Meta_data.csv**: cell type annotations from other tools (e.g., clustering, NMFs, two coloum csv, file, example in ~/MFIs/data folder), is only used in the overrepresentation test heatmap (heatmap tabs), if you don't want to plot this heatmap, you can skip this file.
3. **topDeviatingHOIstates.csv**: located in the output from Abel's pipeline (HOIsummaries folder)
4. **trainingData_.csv**, which is also in the output from Abel's pipeline (output folder)
