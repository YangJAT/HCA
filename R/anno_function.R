

#' anno_function
#'
#' @param data your Seurat object(support V4)
#' @param sig_thres threshold for signature scores
#' @param min_cell  cell types contain a reasonable number of cells
#'
#' @return Seurat object with annotations
#' @export
#'
#' @import tibble
#' @import future.apply
#' @import Seurat
#' @examples

anno_function <- function(data = data,
                          sig_thres = 0.2, min_cell = 50){

  # signature score----------

  sigscore <- data@meta.data
  sigscore <- sigscore[,str_detect(colnames(sigscore), "UCell")]


  #Binary classification----------

  exp_dicho <- function(exp){

    exp <- as.matrix(exp)

    # loop

    dicho <- do.call(cbind, lapply(1:ncol(exp), function(q){

      # Calculate the density

      den <- density(exp[,q])
      den_diff <- diff(den$y, bw = 0.05)

      # Determine thresholds

      thres <- which(den_diff[-1] * den_diff[-length(den_diff)] < 0) + 1
      thres <- den$x[thres[2]]

      data.frame(ifelse(exp[,q] > thres, 1, 0))
    }))

    # result

    colnames(dicho) <- colnames(exp)
    dicho <- rownames_to_column(dicho, var = "id")

    # filter

    dicho[is.na(dicho)] <- 0

    return(dicho)
  }


  # Binary classification----------

  sigscore_dich <- exp_dicho(sigscore)
  sigscore_dich <- column_to_rownames(sigscore_dich, var = "id")


  # Annotate the cells----------

  plan(multisession, workers = ncol(sigscore))
  options(future.globals.maxSize = 3000*1024^2)

  cellanno <- do.call(rbind, future_lapply(colnames(sigscore), function(i){

    cell1 <- rownames(sigscore_dich)[sigscore_dich[,i] == 1]
    cell2 <- rownames(sigscore)[do.call(c, lapply(rownames(sigscore), function(j){
      sigscore[j,i] - max(sigscore[j,!(colnames(sigscore) == i)]) > sig_thres*sigscore[j,i]}))]

    cell <- intersect(cell1, cell2)
    if (length(cell) > 0) {
      data.frame(type = i, id = cell)
    } else {NULL}

  }, future.seed = TRUE))
  plan(sequential)


  #filter ----------

  number <- data.frame(table(cellanno$type))
  cellname <- as.character(number[,1][number[,2] < min_cell])
  cellanno <- cellanno[!(cellanno$type %in% cellname),]


  # Statistical results----------

  allcell <- data.frame(id = colnames(data))
  cellanno <- merge(cellanno, allcell, by = "id", all = T)


  # result ----------

  cellanno <- data.frame(id = cellanno$id, celltype_sig = cellanno$type)
  cellanno <- column_to_rownames(cellanno, var = "id")

  # output ----------

  data <- AddMetaData(data, cellanno)
  return(data)
}



#' anno_function2
#'
#' @param data your Seurat object(support V4)
#' @param sig_thres threshold for signature scores
#' @param min_cell  cell types contain a reasonable number of cells
#'
#' @return Seurat object with annotations
#' @export
#'
#' @import tibble
#' @import future.apply
#' @import Seurat
#' @examples

anno_function2 <- function(data = data,
                           sig_thres = 0.2, min_cell = 50){

  # signature score-------

  sigscore <- data@meta.data
  sigscore <- sigscore[,str_detect(colnames(sigscore), "UCell")]


  # Binary classification----------

  exp_dicho <- function(exp){

    exp <- as.matrix(exp)

    # loop

    dicho <- do.call(cbind, lapply(1:ncol(exp), function(q){
      data.frame(ifelse(exp[,q] > 0.15, 1, 0))}))

    # result

    colnames(dicho) <- colnames(exp)
    dicho <- rownames_to_column(dicho, var = "id")

    # filter

    dicho[is.na(dicho)] <- 0

    return(dicho)
  }


  # Binary classification ----------

  sigscore_dich <- exp_dicho(sigscore)
  sigscore_dich <- column_to_rownames(sigscore_dich, var = "id")


  # annotate the cells ----------

  plan(multisession, workers = ncol(sigscore))
  options(future.globals.maxSize = 3000*1024^2)

  cellanno <- do.call(rbind, future_lapply(colnames(sigscore), function(i){

    cell1 <- rownames(sigscore_dich)[sigscore_dich[,i] == 1]
    cell2 <- rownames(sigscore)[do.call(c, lapply(rownames(sigscore), function(j){
      sigscore[j,i] - max(sigscore[j,!(colnames(sigscore) == i)]) > sig_thres*sigscore[j,i]}))]

    cell <- intersect(cell1, cell2)
    if (length(cell) > 0) {
      data.frame(type = i, id = cell)
    } else {NULL}

  }, future.seed = TRUE))
  plan(sequential)


  # filter ----------

  number <- data.frame(table(cellanno$type))
  cellname <- as.character(number[,1][number[,2] < min_cell])
  cellanno <- cellanno[!(cellanno$type %in% cellname),]


  # result ----------

  allcell <- data.frame(id = colnames(data))
  cellanno <- merge(cellanno, allcell, by = "id", all = T)


  # Organize results  ----------

  cellanno <- data.frame(id = cellanno$id, celltype_sig = cellanno$type)
  cellanno <- column_to_rownames(cellanno, var = "id")

  # output ----------

  data <- AddMetaData(data, cellanno)
  return(data)
}



#' autoumap
#'
#' @param datafilt your Seurat object(support V4)
#' @param nfeatures 	Number of features to select as top variable features
#' @param ndim 	Which dimensions to use as input features
#' @param neigh This determines the number of neighboring points used in local approximations of manifold structure.
#' @param dist This controls how tightly the embedding is allowed compress points together.
#'
#' @return Seurat object with umap
#' @export
#'
#' @import Seurat
#' @examples

autoumap <- function(datafilt, nfeatures = 2000, ndim = 15,
                     neigh = 20, dist = 1){

  #  NormalizeData

  datafilt <- NormalizeData(datafilt, scale.factor = 10000,
                            normalization.method = "LogNormalize")

  # Find Variable Features

  datafilt <- FindVariableFeatures(datafilt, nfeatures = nfeatures,
                                   selection.method = "vst")

  # scale

  datafilt <- ScaleData(datafilt, features = VariableFeatures(datafilt))

  # PCA

  datafilt <- RunPCA(datafilt, assay = 'RNA', slot = 'scale.data')

  # umap
  
  datafilt <- RunUMAP(datafilt, dims = 1:ndim,
                      n.neighbors = neigh, min.dist = dist,
                      reduction = "pca", reduction.name = "umap")

  return(datafilt)
}


#' autocluster
#'
#' @param datafilt your Seurat object(support V4)
#' @param nfeatures 	Number of features to select as top variable features
#' @param ndim 	Which dimensions to use as input features
#' @param neigh This determines the number of neighboring points used in local approximations of manifold structure.
#' @param dist This controls how tightly the embedding is allowed compress points together.
#' @param res resolution for cluster
#'
#' @return Seurat object with cluster
#' @export
#'
#' @import Seurat
#' @examples

autocluster <- function(datafilt, nfeatures = 1000, ndim = 15,
                        neigh = 20, dist = 0.5, res = 0.5){

  # NormalizeData

  datafilt <- NormalizeData(datafilt, scale.factor = 10000,
                            normalization.method = "LogNormalize")

  # FindVariableFeatures

  datafilt <- FindVariableFeatures(datafilt, nfeatures = nfeatures,
                                   selection.method = "vst")

  # scale

  datafilt <- ScaleData(datafilt, features = VariableFeatures(datafilt))

  # PCA

  datafilt <- RunPCA(datafilt, assay = 'RNA', slot = 'scale.data')


  # cell cluster  ------------------------------

  datafilt <- FindNeighbors(datafilt, k.param = neigh,
                            dims = 1:ndim, reduction = "pca")

  datafilt <- FindClusters(datafilt, resolution = res, n.iter = 50)


  # run umap
  datafilt <- RunUMAP(datafilt, dims = 1:ndim,
                      n.neighbors = neigh, min.dist = dist,
                      reduction = "pca", reduction.name = "umap")

  return(datafilt)
}


# Install MAGIC-KNN
#install.packages("reticulate")
#reticulate::install_miniconda() # if no conda env in your PC or server
#reticulate::py_install("magic-impute") # this will install py packages: `magic` and `scprep`
#devtools::install_github("cran/Rmagic")# R interface for magic package

#' MAGIC_impute
#'
#' @param datafilt your Seurat object(support V4)
#' @param knn
#' @param t
#' @param npca
#' @param only_marker impute marker genes
#' @param ncore number of cores
#'
#' @return Seurat object
#' @export
#'
#' @import Seurat
#' @import Rmagic
#' @import Matrix
#' @examples
#'
MAGIC_impute<-function(datafilt, knn = 30, t = 3, npca = 50,
                       only_marker = TRUE,
                       ncore= 10){

  # impute

  magic_data<-magic(datafilt,knn = knn, t = t, npca = npca,
                    genes='all_genes',n.jobs=ncore)

  # filter ----------

  if (only_marker == TRUE) {

    data.magic<-magic_data@assays$MAGIC_RNA$data

    # marker gene

    sel_gene <- c("CD2", "CD3E", "CD3D", "CD3G",
                  "KLRD1", "NKG7", "NCR1", "NCAM1",
                  "CD79A",
                  "CD68", "CD163", "CD14", "CSF1R", "LYZ", "C1QC", "FCN1",
                  "FCER1A", "CD1C", "CLEC9A", "FLT3",
                  "KIT", "MS4A2", "GATA2", "CPA3", "TPSB2",
                  "FCGR3B", "CEACAM8", "CD177", "CSF3R", "ITGAM",
                  "COL1A1", "COL1A2", "FAP", "PDPN", "DCN", "COL3A1", "COL6A1", "ACTA2",
                  "VWF", "CLDN5",
                  "CD8A", "CD8B",
                  "CD4", "FOXP3", "IL2RA",
                  "SPP1", "APOE", "APOC1", "C1QC", "C1QA","C1QB",
                  "LYZ", "FCN1", "S100A8", "S100A9", "IL1B")

    sel_gene<-unique(sel_gene)
    gene<-intersect(sel_gene,rownames(magic_data))

    
    count<-magic_data@assays[['RNA']]$count
    count<-as.matrix(count)

    count[gene,]<-data.magic[gene,]
    count<-Matrix(data = count,sparse = TRUE)

    magic_data@assays[['MAGIC_RNA']]$data<-count

  }else {
    magic_data
  }
  return(magic_data)
}


#' immune cell annotation
#'
#' @param datafilt your Seurat object(support V4)
#' @param scGate_DB scGate gene list
#' @param non_epi epi marker genes
#' @param min_cell cell types contain a reasonable number of cells
#' @param ncore number of cores
#'
#' @return Seurat object with immune cell annotation
#' @export
#' @import Seurat
#' @import UCell
#' @import tibble
#' @import scGate
#' @import future.apply
#' @examples
#'
immune <- function(datafilt, scGate_DB = scGate_DB,
                   non_epi = non_epi, min_cell = 100,
                   ncore = 30){

  # Preliminary division of CD45+and CD45- cells ----------

  model <- scGate_DB$human$generic$CD8T
  model_CD45 <- model[model$levels == "level1",]

  # Calculate ratings

  datafilt <- scGate(datafilt,
                     pos.thr = 0.30, neg.thr = 0.2,
                     ncores = ncore, model = model_CD45)


  # annotate cells ==============================

  datapos <- datafilt[,datafilt$is.pure == "Pure"]

  # T cells

  Tcell <- c("CD2", "CD3E", "CD3D", "CD3G",
             "CD19-", "CD79A-", "MS4A1-", "TNFRSF17-")

  # NK cells

  NKcell <- c("KLRD1", "NKG7", "NCR1", "NCAM1",
              "CD8A-", "CD8B-", "CD4-", "FOXP3-",
              "CD19-", "CD79A-", "MS4A1-", "TNFRSF17-")

  # B cells

  Bcell <- c("CD79A",
             "CD3E-", "CD3D-", "CD8A-", "CD4-",
             "CD68-", "CD163-", "CD14-", "FCGR3A-",
             "FCER1A-", "CLEC9A-", "FLT3-",
             "CEACAM8-", "FCGR3B-")

  # Mono/Macro

  MoMac <- c("CD68", "CD163", "CD14", "CSF1R", "LYZ",
             "C1QC", "FCN1",
             "CD2-", "CD3E-", "CD3D-", "CD3G-",
             "CD19-", "CD79A-", "MS4A1-", "TNFRSF17-",
             "FCER1A-", "CLEC9A-", "FLT3-",
             "FCGR3B-", "CEACAM8-")

  # DCs

  DC <- c("FCER1A", "CD1C", "CLEC9A", "FLT3",
          "CD2-", "CD3E-", "CD3D-", "CD3G-",
          "CD19-", "CD79A-", "MS4A1-", "TNFRSF17-",
          "CD68-", "CD163-", "CD14-", "FCGR3A-",
          "FCGR3B-", "CEACAM8-")

  # MAST

  MAST <- c("KIT", "MS4A2", "GATA2", "CPA3", "TPSB2",
            "CD3E-", "CD3D-", "CD8A-", "CD4-",
            "CD19-", "CD79A-", "MS4A1-", "TNFRSF17-",
            "FCGR3B-", "CEACAM8-")

  # Neutrophils

  Neutrophil <- c("FCGR3B", "CEACAM8", "CD177", "CSF3R", "ITGAM",
                  "CD3E-", "CD3D-", "CD8A-", "CD4-",
                  "CD19-", "CD79A-", "MS4A1-", "TNFRSF17-",
                  "CD14-", "CD163-")

  # 得到signature

  genelist <- list(Tcell = Tcell,
                   NKcell = NKcell,
                   Bcell = Bcell,
                   MoMac = MoMac,
                   DC = DC,
                   MAST = MAST,
                   Neutrophil = Neutrophil)

  # calculate signature

  datapos@meta.data <- datapos@meta.data[,
                                         !(str_detect(colnames(datapos@meta.data), "UCell"))]

  datapos <- AddModuleScore_UCell(datapos, chunk.size = 5000,
                                  ncores = ncore, features = genelist)

  # annotate

  datapos <- anno_function(datapos, sig_thres = 0.15, min_cell = min_cell)


  dataneg <- datafilt[,datafilt$is.pure == "Impure"]


  # Further exclude CD45+cells from CD45- cells----------

  dataneg <- scGate(dataneg,
                    pos.thr = 0.2, neg.thr = 0.3,
                    ncores = ncore, model = model_CD45)

  dataneg <- dataneg[,dataneg$is.pure == "Impure"]


  # Fibroblasts

  Fibroblast <- c("COL1A1", "COL1A2", "FAP", "PDPN", "DCN",
                  "COL3A1", "COL6A1", "ACTA2",
                  "CD3E-", "CD3D-", "CD8A-", "CD4-",
                  "CD19-", "CD79A-", "MS4A1-", "IGKC-", "IGLC2-",
                  "CD68-", "CD163-", "CD14-", "FCGR3A-",
                  "FCER1A-", "CD207-", "CLEC9A-",
                  "FCGR3B-", "CEACAM8-", non_epi)

  # Endothelial cells

  Endothelial <- c("VWF", "CLDN5",
                   "CD3E-", "CD3D-", "CD8A-", "CD4-",
                   "CD19-", "CD79A-", "MS4A1-", "IGKC-", "IGLC2-",
                   "CD68-", "CD163-", "CD14-", "FCGR3A-",
                   "FCER1A-", "CD207-", "CLEC9A-",
                   "FCGR3B-", "CEACAM8-", non_epi)


  genelist <- list(Fibroblast = Fibroblast,
                   Endothelial = Endothelial)

  # calculate signature

  dataneg@meta.data <- dataneg@meta.data[,
                                         !(str_detect(colnames(dataneg@meta.data), "UCell"))]

  dataneg <- AddModuleScore_UCell(dataneg, chunk.size = 5000,
                                  ncores = ncore, features = genelist)

  # annotate cells

  dataneg <- anno_function2(dataneg, sig_thres = 0.3, min_cell = min_cell)


  # integrate==============================

  # filter cells

  datameta <- merge(datapos, dataneg)
  datameta <- datameta[,!is.na(datameta$celltype_sig)]


  datameta$celltype_sig <- do.call(rbind,
                                   strsplit(datameta$celltype_sig, '_'))[,1]


  # annotate T cells ==============================

  tcell <- datameta[,datameta$celltype_sig == "Tcell"]

  # CD8 Tcell

  CD8T <- c("CD8A", "CD8B",
            "CD4-", "FOXP3-", "IL2RA-", "CD40LG-")

  # CD4 Tcell

  CD4T <- c("CD4",
            "CD8A-", "CD8B-", "GZMA-", "KLRB1-",
            "FOXP3-", "IL2RA-", "TNFRSF4-", "TNFRSF18-")

  # Treg

  Treg <- c("CD4", "FOXP3", "IL2RA",
            "IL7R-", "CD8A-", "CD8B-",
            "GZMA-", "KLRB1-")

  # signature

  genelist <- list(CD8T = CD8T,
                   CD4T = CD4T,
                   Treg = Treg)

  # calculate signature

  tcell@meta.data <- tcell@meta.data[,
                                     !(str_detect(colnames(tcell@meta.data), "UCell"))]

  tcell <- AddModuleScore_UCell(tcell, chunk.size = 5000,
                                ncores = ncore, features = genelist)

  # annotate cells

  tcell <- anno_function(tcell, sig_thres = 0.1, min_cell = min_cell)


  # Further subdivision of MoMac cells ==============================

  momac <- datameta[,datameta$celltype_sig == "MoMac"]

  # Macrophage

  Macro <- c("SPP1", "APOE", "APOC1", "C1QC", "C1QA","C1QB",
             "S100A8-", "FCN1-")

  # Monocyte

  Mono <- c("LYZ", "FCN1", "S100A8", "S100A9", "IL1B",
            "APOE-", "SPP1-", "C1QC-")

  # signature

  genelist <- list(Macro = Macro,
                   Mono = Mono)

  # calculate signature

  momac@meta.data <- momac@meta.data[,
                                     !(str_detect(colnames(momac@meta.data), "UCell"))]

  momac <- AddModuleScore_UCell(momac, chunk.size = 5000,
                                ncores = ncore, features = genelist)

  # annotate

  momac <- anno_function(momac, sig_thres = 0.1, min_cell = min_cell)


  # integrate cells ==============================

  anno_tcell <- data.frame(id = colnames(tcell), celltype_sig2 = tcell$celltype_sig)
  anno_momac <- data.frame(id = colnames(momac), celltype_sig2 = momac$celltype_sig)
  anno <- rbind(anno_tcell, anno_momac)

  # preprocess

  anno <- na.omit(anno)
  rownames(anno) <- NULL

  anno$celltype_sig2 <- do.call(rbind,
                                strsplit(anno$celltype_sig2, '_'))[,1]

  # integrate data

  anno <- column_to_rownames(anno, var = "id")
  datameta <- AddMetaData(datameta, anno)


  cellname <- c("Bcell", "DC", "NKcell", "Fibroblast",
                "Endothelial", "MAST", "Neutrophil")

  select <- datameta$celltype_sig %in% cellname
  datameta$celltype_sig2[select] <- datameta$celltype_sig[select]
  datameta <- datameta[,!is.na(datameta$celltype_sig2)]

  return(datameta)
}


#' tumor annotation
#'
#' @param datafilt your Seurat object(support V4)
#' @param scGate_DB scGate gene list
#' @param ncore numbers of cores
#'
#' @return Seurat object
#' @export
#'
#' @import iCNA
#' @import scGate
#' @import tibble
#' @import Seurat
#' @import future.apply
#' @examples
tumor <- function(datafilt = datafilt,
                  scGate_DB = scGate_DB,
                  ncore = 30){

  useGenome('hg38')

  # ScGate screening for CD45 positive/negative cells------------------------------

  model <- scGate_DB$human$generic$CD8T
  model_CD45 <- model[model$levels == "level1",]


  datafilt <- scGate(datafilt,
                     pos.thr = 0.30, neg.thr = 0.05,
                     ncores = ncore, model = model_CD45)


  # Prepare for inferCNA analysis------------------------------

  # Generate reference

  normal <- datafilt[,datafilt$is.pure == "Pure"]


  if (ncol(normal) > 1000) {
    ref <- normal[,sample(colnames(normal), 1000)]
  } else {ref <- normal}
  ref$celltype <- "normal"


  refname <- list(ref = colnames(ref))


  # Generate inferCNA input file----------


  datafilt <- scGate(datafilt,
                     pos.thr = 0.2, neg.thr = 0.2,
                     ncores = ncore, model = model_CD45)

  # Extract CD45- cells

  input <- datafilt[,datafilt$is.pure == "Impure"]
  input$celltype <- "candidate"

  input <- merge(ref, input)

  # Generate annotation file

  info <- data.frame(id = colnames(input), type = input$celltype)
  info_p <- data.frame(id = colnames(input), sample = input$sample)


  #Generate input matrix

  if (ncol(input) > 50000){

    input$group <- 1:ncol(input)
    input$group <- ntile(input$group, n = 20)
    input_split <- SplitObject(input, split.by = "group")

    input <- do.call(cbind, lapply(input_split, function(i){
      as.matrix(GetAssayData(i, slot = "counts", assay = "RNA"))}))

    input <- t(t(input)/colSums(input) * 1000000)
    input <- log2(input/10 + 1)

  } else {

    input <- as.matrix(GetAssayData(input, slot = "counts", assay = "RNA"))
    input <- t(t(input)/colSums(input) * 1000000)
    input <- log2(input/10 + 1)

  }


  # inferCNA ----------

  cna <- iCNA::infercna(m = input, refCells = refname, window = 150,
                        n = 5000, noise = 0.1, isLog = TRUE, verbose = FALSE)

  # calculate cnaSignal

  signal <- data.frame(signal = cnaSignal(cna, refCells = refname))

  # calculate cnaCor

  cor <- do.call(rbind, lapply(unique(info_p$sample), function(i){

    name <- info_p$id[info_p$sample == i]
    common <- setdiff(name, unlist(refname))

    subdata <- cna[,name]
    correlation <- cnaCor(subdata, threshold = NULL, bySample = F,
                          samples = NULL, excludeFromAvg = common)

    data.frame(cor = correlation)
  }))

  # summary the result

  signal <- rownames_to_column(signal, var = "id")
  cor <- rownames_to_column(cor, var = "id")

  allscore <- merge(signal, cor, by = "id")
  allscore <- merge(allscore, info, by = "id")


  # output ----------

  list(cna = cna,
       refname = refname,
       allscore = allscore,
       info_p = info_p,
       normal = normal)
}




#' plot_thres
#'
#' @param cnadata
#' @param thres_sig
#' @param thres_cor
#'
#' @return
#' @export
#' @import iCNA
#' @import ggplot2
#' @examples
#'
plot_thres <- function(cnadata,
                       thres_sig = 0.005,
                       thres_cor = 0.5){

  allscore <- cnadata[["allscore"]]

  # color

  col_rb <- c("#343391","#0064af","#0090cc",
              "#00b6db","#01b7c2","#53c0a3",
              "#8dcb8a","#bbd967","#fbd324",
              "#f6bd25","#f4a02e","#ed6f32",
              "#ea5c2e","#d5452f","#c02e2f","#8b2a21")


  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }

  # Calculate the density of each point

  allscore <- do.call(rbind, lapply(unique(allscore$type), function(i){
    subdata <- allscore[allscore$type == i,]
    subdata$dens <- get_density(subdata$signal, subdata$cor, n = 1000)
    subdata$dens <- subdata$dens / max(subdata$dens)
    subdata
  }))

  # plot

  plot <- ggplot(allscore, aes(signal, cor)) +
    geom_point(aes(color = dens), size = 1) +
    scale_color_gradientn(colors = col_rb) +

    geom_vline(xintercept = thres_sig, color = "grey20", linetype = "longdash") +
    geom_hline(yintercept = thres_cor, color = "grey20", linetype = "longdash") +
    xlim(0, quantile(allscore$signal, probs=c(0.995))) +

    theme_bw() + theme(panel.grid = element_blank()) +
    facet_wrap(~type, nrow = 1)

  ggsave("inferCNV/scatter_plot.png", plot, dpi = 200,
         width = 11, height = 5, limitsize = FALSE)
}



#' Title
#'
#' @param datafilt your Seurat object(support V4)
#' @param cnadata
#' @param thres_sig
#' @param thres_cor
#' @param ncore number of core
#' @param isFilter filter Endothelial/Stomal cells
#'
#' @return
#' @export
#' @import iCNA
#' @import Seurat
#' @import tibble
#' @examples
#'
anno_tumor <- function(datafilt = datafilt,
                       cnadata = cnadata,
                       thres_sig = 0.005,
                       thres_cor = 0.5,
                       ncore = 30,
                       isFilter = FALSE){


  allscore <- cnadata[["allscore"]]
  refname <- cnadata[["refname"]]
  normal <- cnadata[["normal"]]
  cna <- cnadata[["cna"]]


  #Extracting tumor cells----------

  candidate <- allscore[allscore$type == "candidate",]


  cell_tumor <- candidate$id[candidate$signal > thres_sig & candidate$cor > thres_cor]


  # Integrate ----------

  # 由于前面只抽取了1000个ref，这里要把其他未被抽取的细胞也标注为normal

  anno <- data.frame(id = cell_tumor, sample = "tumor")
  colnames(anno)[2] <- "tumor_cl"


  rownames(anno) <- NULL
  anno <- column_to_rownames(anno, var = "id")


  datafilt <- AddMetaData(datafilt, anno)
  datafilt <- datafilt[,!is.na(datafilt$tumor_cl)]


  # filter ----------

  if (isFilter == TRUE) {


    Endothelial <- c("PECAM1", "VWF", "ENG", "CDH5")
    Endothelial <- gating_model(name = "Endothelial", signature = Endothelial)

    datafilt <- scGate(datafilt,
                       pos.thr = 0.15, neg.thr = 0.2,
                       ncores = ncore, model = Endothelial)

    datafilt <- datafilt[,datafilt$is.pure == "Impure"]



    Fibroblast <- c("COL1A1", "COL1A2", "FAP", "PDPN", "DCN", "COL3A1", "COL6A1", "ACTA2")
    Fibroblast <- gating_model(name = "Fibroblast", signature = Fibroblast)

    datafilt <- scGate(datafilt,
                       pos.thr = 0.15, neg.thr = 0.2,
                       ncores = ncore, model = Fibroblast)

    datafilt <- datafilt[,datafilt$is.pure == "Impure"]

  } else {
    datafilt
  }
}



#' Title
#'
#' @param cna_heatmap
#' @param thres_sig
#' @param thres_cor
#' @param max_cell
#'
#' @return
#' @export
#' @import iCNA
#' @import Seurat
#' @import ggplot2
#' @examples
cna_heatmap <- function(cnadata = cnadata,
                        thres_sig = 0.005,
                        thres_cor = 0.5,
                        max_cell = 15000){

  allscore <- cnadata[["allscore"]]
  refname <- cnadata[["refname"]]
  info_p <- cnadata[["info_p"]]
  cna <- cnadata[["cna"]]


  # Extracting tumor cells

  candidate <- allscore[allscore$type == "candidate",]


  cell_tumor <- candidate$id[candidate$signal > thres_sig & candidate$cor > thres_cor]
  cell_normal <- setdiff(candidate$id, cell_tumor)
  cell_ref <- unlist(refname)

  #Obtain various cell annotations

  anno_tumor <- info_p[info_p$id %in% cell_tumor,]
  anno_tumor <- anno_tumor[order(anno_tumor$sample),] %>% list()

  anno_normal <- data.frame(id = cell_normal, sample = "normal") %>% list()
  anno_ref <- data.frame(id = cell_ref, sample = "ref") %>% list()

  # summary

  anno <- do.call(rbind, c(anno_ref, anno_normal, anno_tumor))
  anno <- split(anno$id, anno$sample)

  #filter

  number <- do.call(c, lapply(anno, function(i){length(i)}))
  select <- names(number)[number > 50]

  # output

  anno <- anno[select]


  # plot ----------

  cna_ref <- cna[,cell_ref]
  cna_normal <- cna[,cell_normal]
  cna_tumor <- cna[,anno_tumor[[1]][["id"]]]


  if (length(cell_normal) < 5000) {
    cna_normal <- cna_normal
  } else {cna_normal <- cna_normal[,sort(sample(length(cell_normal), 5000))]}

  if (length(cell_tumor) < max_cell) {
    cna_tumor <- cna_tumor
  } else {cna_tumor <- cna_tumor[,sort(sample(length(cell_tumor), max_cell))]}


  # heatmap

  plot <- cnaPlot(cna = cna_ref)
  plot <- plot$p
  ggsave("inferCNV/cna_ref.png", plot, dpi = 100,
         width = 15, height = 20, limitsize = FALSE)

  # normal

  plot <- cnaPlot(cna = cna_normal)
  plot <- plot$p
  ggsave("inferCNV/cna_normal.png", plot, dpi = 100,
         width = 15, height = 20, limitsize = FALSE)

  # tumor

  plot <- cnaPlot(cna = cna_tumor)
  plot <- plot$p
  ggsave("inferCNV/cna_tumor.png", plot, dpi = 100,
         width = 15, height = 20, limitsize = FALSE)
}



#' auto_immune
#'
#' @param datafilt
#' @param scGate_DB
#' @param non_epi
#' @param min_cell
#' @param MAGIC_impute
#' @param ncore
#'
#' @return
#' @export
#' @import iCNA
#' @import Seurat
#' @import ggplot2
#' @import Rmagic
#' @import Matrix
#' @examples
auto_immune <- function(datafilt, scGate_DB = scGate_DB,
                        non_epi = non_epi, min_cell = 100,
                        MAGIC_impute = FALSE,
                        ncore = 30){


  if (ncol(datafilt)>100000){

    datafilt$split_group<-1:ncol(datafilt)
    datafilt$split_group<-dplyr::ntile(datafilt$split_group,n=ncol(datafilt)%/%100000+1)
    data_split<-SplitObject(datafilt,split.by = 'split_group')

    alldata<-lapply(data_split, function(i){

      ## MAGIC impute

      if (MAGIC_impute == TRUE) {
        data<-MAGIC_impute(i,ncore = ncore)
        data@active.assay <- 'MAGIC_RNA'
      } else {
        data<-i
      }

      # immune annotate

      datameta<-immune(data, scGate_DB = scGate_DB,
                       non_epi = non_epi, min_cell = min_cell,
                       ncore = ncore)

      return(datameta)

    })

    ##integrate

    dataintg<-alldata[1]

    for (i in 2:length(alldata)) {
      subdata <- alldata[[i]]
      dataintg <- merge(dataintg, subdata)}

    return(dataintg)

  } else {

    ## MAGIC impute

    if (MAGIC_impute == TRUE) {
      datafilt<-MAGIC_impute(datafilt,ncore = ncore)
      datafilt@active.assay <- 'MAGIC_RNA'
    } else {
      datafilt<-datafilt
    }

    datameta<-immune(datafilt, scGate_DB = scGate_DB,
                     non_epi = non_epi, min_cell = min_cell,
                     ncore = ncore)

    return(datameta)

  }
}


# ==============================================================================
# Tumor annotation
# ==============================================================================


#' infertumor
#'
#' @param datafilt
#' @param scGate_DB
#' @param thres_sig
#' @param thres_cor
#' @param isFilter
#' @param ncore
#'
#' @return
#' @export
#' @import iCNA
#' @import Seurat
#' @import Rmagic
#' @import Matrix
#' @examples
infertumor <- function(datafilt = datafilt,
                       scGate_DB = scGate_DB,
                       thres_sig = 0.005,
                       thres_cor = 0.5,
                       isFilter = FALSE,
                       ncore = 30){


  if (ncol(datafilt)>100000){

    datafilt$split_group<-1:ncol(datafilt)
    datafilt$split_group<-dplyr::ntile(datafilt$group,n=ncol(datafilt)%/%100000+1)
    data_split<-SplitObject(datafilt,split.by = 'group')

    alldata<-lapply(data_split, function(i){

      ## tumor_cna

      cnadata<-tumor(datafilt = i,
                     scGate_DB = scGate_DB,
                     ncore = ncore)



      datacanc <- anno_tumor(datafilt = i, cnadata = cnadata,
                             thres_sig = thres_sig,
                             thres_cor = thres_cor,
                             ncore = ncore,
                             isFilter = TRUE)

      return(datacanc)

    })

    ##integrate

    dataintg<-alldata[1]

    for (i in 2:length(alldata)) {
      subdata <- alldata[[i]]
      dataintg <- merge(dataintg, subdata)}

    return(dataintg)

  } else {

    ## tumor_cna

    cnadata<-tumor(datafilt = datafilt,
                   scGate_DB = scGate_DB,
                   ncore = ncore)

    # annotate

    datacanc <- anno_tumor(datafilt, cnadata = cnadata,
                           thres_sig = thres_sig,
                           thres_cor = thres_cor,
                           ncore = ncore,
                           isFilter = TRUE)
    return(datacanc)

  }
}


#' integrate
#'
#' @param dataimmu
#' @param datacanc
#' @param min_tumor
#' @param rm_doublet
#' @param prop_doublet
#'
#' @return
#' @export
#' @import Seurat
#' @import DoubletFinder
#' @import dplyr
#' @import future.apply
#' @examples
integrate <- function(dataimmu = dataimmu,
                      datacanc = datacanc,
                      min_tumor = 50,
                      rm_doublet = FALSE,
                      prop_doublet = 0.075){


  select <- intersect(colnames(dataimmu), colnames(datacanc))
  dataimmu <- dataimmu[,!(colnames(dataimmu) %in% select)]
  datacanc <- datacanc[,!(colnames(datacanc) %in% select)]
  dataintg <- merge(dataimmu, datacanc)


  dataintg$celltype_sig2[is.na(dataintg$celltype_sig2)] <- "tumor"


  # filter doublets

  if (rm_doublet == TRUE){

    dataintg <- autocluster(dataintg, nfeatures = 2000, ndim = 15,
                            neigh = 20, dist = 0.5, res = 0.5)

    # seek  k numbers

    sweep.res.list <- paramSweep_v3(dataintg, PCs = 1:20)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)

    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>%
      as.character() %>% as.numeric()


    homotypic.prop <- modelHomotypic(dataintg$Seurat_clusters)

    nExp_poi <- round(prop_doublet*ncol(dataintg))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


    dataintg <- doubletFinder_v3(dataintg, PCs = 1:20, pN = 0.25,
                                 pK = pK_bcmvn, reuse.pANN = F,
                                 nExp = nExp_poi.adj)

    names(dataintg@meta.data) <- gsub("DF.classifications.*",
                                      "DF.classifications",
                                      colnames(dataintg@meta.data))

    infodata <- data.frame(id = rownames(dataintg@meta.data),
                           DFclass = dataintg@meta.data[["DF.classifications"]])

    # integrate date

    rownames(infodata) <- NULL
    infodata <- column_to_rownames(infodata, var = "id")
    dataintg <- AddMetaData(dataintg, infodata)

    # exclude doublet

    dataintg <- dataintg[,dataintg$DFclass == "Singlet"]
  }


  # filter tumor

  select <- data.frame(cell = colnames(dataintg),
                       id = dataintg$sample,
                       type = dataintg$celltype_sig2)

  select <- select[select$type %in% "tumor",]
  number <- data.frame(table(select$id))

  number$Var1 <- as.character(number$Var1)
  number <- number$Var1[number$Freq < min_tumor]

  select <- select$cell[select$id %in% number]
  dataintg <- dataintg[,!(colnames(dataintg) %in% select)]

  # umap

  dataintg <- autoumap(dataintg, nfeatures = 2000, ndim = 15,
                       neigh = 20, dist = 0.75, res = 0.5)
}




