# HCA (Hierarchical Cell Annotation)

![图片](https://github.com/user-attachments/assets/bdcd4e4b-7f14-4a06-a7c0-adb2615ff9ba)

## Note
This package has several dependencies with version constraints: <br>
1. Seurat: version < 5. <br>
2. scGate: version = 1.2.0. <br>
3. future: version = 1.31.0 <br>
```r
# scGate
library(remotes)
remotes::install_github("carmonalab/scGate", ref="v1.2.0")

# future
# download future_1.31.0.tar.gz from https://cran.r-project.org/src/contrib/Archive/future/
install.packages("future_1.31.0.tar.gz", repos = NULL, type = "source")

```

## Abstract 
This tutorial demonstrates how to use R packages for single-cell RNA data analysis, including cell type annotation, tumor cell annotation, and data integration.
## Loading Required R Packages

```r
library(HCA)
library(viridis)
library(Seurat)
library(iCNA)
library(stringr)
library(scGate)
library(future)
```

load the data:

```r
scGate_DB <- readRDS("data/scGate_DB.rds")
datafilt <- readRDS("data/sc_datafilt.rds")
```

## Celltype annotation 

For Human Single-Cell RNA Data:
```r
non_epi <- c("EPCAM-", "CDH1-", "KRT7-", "KRT18-", "KRT19-", "ALB-", "AFP-")

#Annotating Immune Cells
dataimmu <- anno_immune(datafilt, scGate_DB = scGate_DB, organism = 'human', non_epi = non_epi, min_cell = 100, ncore = 1)

#Annotating Tumor Cells
datacanc <- anno_tumor(datafilt, scGate_DB = scGate_DB, 
                       organism = 'human', 
                       thres_sig = 0.005, # Adjust this threshold based on scatter_plot.png
                       thres_cor = 0.5, # Adjust this threshold based on scatter_plot.png
                       ncore = 1, # Multi-core functionality is not available on Windows
                       isFilter = TRUE)
```
Note: If the code runs successfully, an image (inferCNV/scatter_plot.png) will be generated in the current path. You can select the threshold range based on the scatter plot positions in the image.

For Mouse Single-Cell RNA Data:
```r
non_epi <- c("Krt5-", "Krt14-", "Krt6a-", "Dsp-", "Krt17-", "Lgals7-")

#Annotating Immune Cells
dataimmu <- anno_immune(datafilt, scGate_DB = scGate_DB, organism = 'mouse', non_epi = non_epi, min_cell = 100, ncore = 1)

#Annotating Tumor Cells
datacanc <- anno_tumor(datafilt, scGate_DB = scGate_DB, 
                       organism = 'mouse', 
                       thres_sig = 0.005, 
                       thres_cor = 0.5, 
                       ncore = 1, 
                       isFilter = TRUE)
```

## Data Integration
```r

dataintg <- integrate(dataimmu, datacanc,
                      min_tumor = 50,
                      rm_doublet = FALSE,
                      prop_doublet = 0.075)

saveRDS(dataintg, 'data/sc_datafilt_anno.rds')

```

After running, the Seurat object will include celltype_sig2, representing the annotation results.

## Further filtering
```r

dataintg <- autocluster(dataintg, nfeatures = 2000,
                        ndim = 15, neigh = 20,
                        dist = 1, res = 3) # Set a higher resolution (res) to capture more clusters

# celltype visualization
dimplot_new(dataintg,
            reduction = "umap",
            pt.size = 0.2, label = T,
            group.by = c("celltype_sig2"))

# High-resolution clustering visualization
dimplot_new(dataintg,
            reduction = "umap",
            pt.size = 0.5, label = T,
            group.by = c("seurat_clusters"))
```

## Visualization

This visualization specifically delineates the comparison between the control and experimental groups
```r

# UMAP density plot
prop_density(datafilt = datafilt,
             group = "group",
             coord = "umap_harmony")

# Back-to-back plot
prop_back2back(datafilt = datafilt,
               group = "group",
               cluster = "seurat_clusters",
               order = TRUE)

prop_back2back_lollipop(datafilt = dataimmu,
                        group = "group",
                        group1 = "H",
                        group2 = "L",
                        cluster = "celltype_sig2")

# Sample-level proportional distribution difference
input <- data.frame(table(dataimmu$sample, dataimmu$celltype_sig2))
prop_plot_hca(input, rotate = 45, decreasing = T, species = "human")

```

Note: The iCNA package is essentially a more installable version of the infercna package (see https://github.com/jlaffy/infercna), created to address the challenges often encountered with installing infercna across different environments. If you use our package, please cite both our study (https://doi.org/10.1016/j.ccell.2024.10.008) and the related article for the infercna package (https://doi.org/10.1016/j.cell.2019.06.024).

