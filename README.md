# HCA (Hierarchical Cell Annotation)
## Abstract 
This tutorial shows how to use HCA to identify cells with usual ‘gene set’ (i.e. gene signatures) in single-cell RNA-seq data. 
## How to install this R package:

>##You need to install iCNA R package first: \
>devtools::install_github('Liuzhicheng048/iCNA') \
>##Then you can install this package:\
>devtools::install_github('YangJAT/HCA')

## How to use this R package:

## Preprocess data:

>library(HCA) \
>library(viridis) \
>library(UCell) \
>library(scGate) \
>library(tibble) \
>library(stringr) \
>library(Seurat) \
>library(iCNA) \
>library(ggplot2) \
>library(future.apply) 
>
>
>mycol <- c("#0070B2", "#5CB3DA", "#B8E3EA", "#DA1735", "#F15E4C", "#FF9F99", \
>           "#A231A1", "#A37CB7", "#F2D7EE", "#B91372", "#E93B8C", "#ECB2C8", \
>           "#FF7149", "#F7AE24", "#FBDD7E", "#679436", "#8BBE53", "#CDE391", \
>           "#067D69", "#00A385", "#98D4C6", "#114B5F", "#028090", "#B2DBBF", \
>          "#A23E48", "#CD6981", "#FBD0C0", "#788585", "#9CAEA9", "#CCDAD1")
>
>
>scGate_DB <- readRDS("scgate/auto_anno/scGate_DB.rds")
>
>##load data
>
>datafilt <- readRDS("data/sc_datafilt.rds")
>
>datafilt=PercentageFeatureSet(datafilt, "^mt-", col.name = "percent_mito")
>
>par(mfrow=c(2,2))
>
>hist(datafilt$nFeature_RNA,breaks=300,prob=TRUE) \
>hist(datafilt$nCount_RNA,breaks=300,prob=TRUE) \
>hist(datafilt$percent_mito,breaks=300,prob=TRUE)
>
>dev.off()
>
>selected_c <- WhichCells(datafilt, expression = nFeature_RNA > 200 & nCount_RNA < 100000 & >nFeature_RNA < 10000) \
>selected_f <- rownames(datafilt)[Matrix::rowSums(datafilt@assays$RNA@counts > 0 ) > 3] \
>datafilt <- subset(datafilt, features = selected_f, cells = selected_c) \
>selected_mito <- WhichCells(datafilt , expression = percent_mito < 25)  \
>datafilt <- subset(datafilt, cells = selected_mito)
>
>dim(datafilt)

## Celltype annotation 
>##for human scRNA:
>
>non_epi <- c("EPCAM-", "CDH1-", "KRT7-", "KRT18-", "KRT19-", "ALB-", "AFP-")
>
>dataimmu <- anno_immune(datafilt, scGate_DB = scGate_DB,organism = 'human', non_epi = non_epi, min_cell = 100, ncore = 1)
>
>dataimmu<-autoumap(dataimmu)
>
>plot <- DimPlot(dataimmu, pt.size = 0.1, label = T, repel = T, cols = sample(mycol), 
>               raster = FALSE, label.size = 5, reduction = "umap", 
>                group.by = c("celltype_sig"))
>
>ggsave("result/umap_immu.png", plot, dpi = 300, width = 9, height = 7)
>
>datacanc <- anno_tumor(datafilt, scGate_DB = scGate_DB, 
>                       organism = 'human', 
>                      thres_sig = 0.005, 
>                     thres_cor = 0.5, 
>                    ncore = 10, 
>                     isFilter = TRUE)
>
>##for mouse scRNA:
>
>non_epi <- c("Krt5-", "Krt14-", "Krt6a-", "Dsp-", "Krt17-", "Lgals7-")
>
>dataimmu <- anno_immune(datafilt, scGate_DB = scGate_DB,organism = 'mouse', non_epi = non_epi, min_cell = 100, ncore = 1)
>
>dataimmu<-autoumap(dataimmu)
>
>plot <- DimPlot(dataimmu, pt.size = 0.1, label = T, repel = T, cols = sample(mycol), 
>               raster = FALSE, label.size = 5, reduction = "umap", 
>                group.by = c("celltype_sig"))
>
>ggsave("result/umap_immu.png", plot, dpi = 300, width = 9, height = 7)
>
>datacanc <- anno_tumor(datafilt, scGate_DB = scGate_DB, 
>                       organism = 'mouse', 
>                      thres_sig = 0.005, 
>                     thres_cor = 0.5, 
>                    ncore = 10, 
>                     isFilter = TRUE)

# Integration

>dataintg <- integrate(dataimmu, datacanc)
>
>plot <- DimPlot(dataintg, pt.size = 0.1, label = T, repel = T, cols = sample(mycol),
>                raster = FALSE, label.size = 5, reduction = "umap",
>               group.by = c("sample"))
>
>ggsave("result/umap_intg.png", plot, dpi = 300, width = 9, height = 7)
>
>
>
>info<-dataintg@meta.data[c('celltype_sig','celltype_sig2')] \
>datafilt<-AddMetaData(datafilt,info) \
>datafilt<-autoumap(datafilt) \
>plot <- DimPlot(datafilt, pt.size = 0.1, label = T, repel = T, cols = sample(mycol),
>               raster = FALSE, label.size = 5, reduction = "umap",
>                group.by = c("celltype_sig2"))
>
>ggsave("result/umap_clusters.png", plot, dpi = 300, width = 9, height = 7)
>
>saveRDS(datafilt,'data/sc_datafilt_anno.rds')


### Note: this package was built on seurat V4 and will be updated to support V5 in the future.
![图片](https://github.com/YangJAT/HCA/assets/70686083/d8fb4993-175e-453f-bff6-45bcd8c91ef3)
