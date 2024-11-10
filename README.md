# HCA (Hierarchical Cell Annotation)

![图片](https://github.com/user-attachments/assets/bdcd4e4b-7f14-4a06-a7c0-adb2615ff9ba)

This tutorial shows how to use HCA R package for single-cell RNA data analysis, including immune cell annotation, tumor cell annotation, and data integration.

 <br>

## Note
### Dependency issues
This package has several dependencies with version constraints: <br>
1. Seurat: version <5 (excluding 4.9); 4.3 recommended <br>
2. scGate: version = 1.2.0 <br>
3. future: version = 1.31.0 <br>
4. Matrix:version=1.5-3 <br>
```r
# scGate
library(remotes)
remotes::install_github("carmonalab/scGate", ref = "v1.2.0")

# future
# download future_1.31.0.tar.gz from https://cran.r-project.org/src/contrib/Archive/future/
install.packages("future_1.31.0.tar.gz", repos = NULL, type = "source") # Install from a local directory

#Matrix
# download Matrix_1.5-3 .tar.gz from https://cran.r-project.org/src/contrib/Archive/Matrix/
install.packages("Matrix_1.5-3 .tar.gz", repos = NULL, type = "source")


# install our package
devtools::install_github('Liuzhicheng048/iCNA')
devtools::install_github('YangJAT/HCA')

```
### Memory issues
We recommend utilizing a high-memory server for optimal performance with this package. Running it on a personal computer is possible, though it may result in substantially slower processing speeds.

 <br>

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

"datafilt" is a Seurat object that requires only basic cell filtering and includes a column labeled "sample" to define the cell-to-sample correspondence, with no need for additional processing. 

 <br>

## Celltype annotation 

### Annotating Immune Cells
```r
non_epi <- c("EPCAM-", "CDH1-", "KRT7-", "KRT18-", "KRT19-", "ALB-", "AFP-") # for human
non_epi <- c("Krt5-", "Krt14-", "Krt6a-", "Dsp-", "Krt17-", "Lgals7-") # for mouse

dataimmu <- anno_immune(datafilt,
                        scGate_DB = scGate_DB,
                        organism = "human", # or mouse
                        non_epi = non_epi,
                        min_cell = 100,
                        ncore = 1) # Multi-core functionality is not available on Windows
```

### Annotating Tumor Cells
Note 1: This step is optional. If your data has undergone CD45 sorting, then you only need to run immune cell annotation, and data integration can also be skipped. <br>
Note 2: The input Seurat object must include a column labeled "sample" to define the cell-to-sample correspondence.
```r
datacanc <- anno_tumor(datafilt,
                       scGate_DB = scGate_DB, 
                       organism = "human", # or mouse
                       thres_sig = 0.005, # Adjust this threshold based on scatter_plot.png
                       thres_cor = 0.5, # Adjust this threshold based on scatter_plot.png
                       ncore = 1, # Multi-core functionality is not available on Windows
                       isFilter = TRUE)
```

If the code runs successfully, an image (inferCNV/scatter_plot.png) will be generated in the current path. You can select the threshold range based on the scatter plot positions in the image.

![图片](https://github.com/user-attachments/assets/1f66831c-f017-4c73-82c8-136e53c79f85)

 <br>

## Data Integration
```r

dataintg <- integrate(dataimmu, datacanc,
                      min_tumor = 50,
                      rm_doublet = FALSE,
                      prop_doublet = 0.075)

# If you skipped the annotation of tumor cells, please run
# dataintg <- dataimmu
```

After running, the Seurat object will contain a column labeled "celltype_sig2," representing the annotation results. 

 <br>

## Further filtering
```r

source("data/Other functions.R")

# clustering

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

Simultaneously review the cell type annotations and Seurat clustering results, removing clusters that encompass cells from divergent lineages (e.g., myeloid and lymphoid lineages within a single cluster) or clusters with atypical spatial positioning on the UMAP plot (e.g., T cell subsets positioned in close proximity to myeloid cells).

![图片](https://github.com/user-attachments/assets/6cb03d4b-832e-4563-b8e3-13a8ebe285f3)

```r
# exclude any problematic clusters

select = c("31","35","39","40","51")
dataintg <- dataintg[,!(dataintg$seurat_clusters %in% select)]

# re-analyze

dataintg <- autocluster(dataintg, nfeatures = 2000,
                        ndim = 15, neigh = 20,
                        dist = 1, res = 3)

dimplot_new(dataintg,
            reduction = "umap",
            pt.size = 0.2, label = T,
            group.by = c("celltype_sig2"))
```

![图片](https://github.com/user-attachments/assets/58d3954d-0928-4415-a78e-94c692456517)

 <br>

## Visualization
```r

source("data/Other functions.R")

# dotplot visualization (with default gene set)

name = "marker_dotplot.pdf"
dotplot_marker(dataintg,
               group.by = "celltype_sig2",
               marker = NULL,
               species = "human", # or mouse
               output = name,
               height = 6)

# dotplot visualization (manually selected gene set)

Tcell = c("Cd3d", "Cd3e")
CD8T = c("Cd8a", "Cd8b1")
gene_list <- list(name1 = Tcell,
                  name2 = CD8T)

name = "marker_dotplot.pdf"
dotplot_marker(dataintg,
               group.by = "celltype_sig2",
               marker = gene_list,
               species = NULL,
               output = name,
               height = 6)
```
![图片](https://github.com/user-attachments/assets/af7ef1fa-193f-4715-acdc-5765d7ca6e90)


This visualization specifically delineates the comparison between the control and experimental groups
```r

source("data/Other functions.R")

# UMAP density plot

prop_density(datafilt = datafilt,
             group = "group", # grouping information
             coord = "umap")

# Back-to-back plot

prop_back2back(datafilt = datafilt,
               group = "group", # grouping information
               cluster = "seurat_clusters",
               order = TRUE)

# Sample-level proportional distribution difference

input <- data.frame(table(dataimmu$sample, dataimmu$celltype_sig2))
prop_plot_hca(input, rotate = 45, decreasing = T, species = "human")
```
![图片](https://github.com/user-attachments/assets/35e35a93-5075-4392-adfd-0348b746436f)

More analysis and visualization capabilities will be introduced in upcoming updates.

 <br>

## How to cite
The iCNA package is essentially a more installable version of the infercna package (see https://github.com/jlaffy/infercna), created to address the challenges often encountered with installing infercna across different environments. If you use our package, please cite both our study (https://doi.org/10.1016/j.ccell.2024.10.008) and the related article for the infercna package (https://doi.org/10.1016/j.cell.2019.06.024).

