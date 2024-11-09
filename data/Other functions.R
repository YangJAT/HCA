
# dimplot_new ====================

dimplot_new <- function(data = datafilt,
                        reduction = "umap",
                        pt.size = 1, label = T,
                        group.by = c("seurat_clusters")) {
  
  library(Seurat)
  
  info <- data@meta.data
  number <- length(unique(info[,group.by]))
  
  if (number < 17) {
    col <- c("#A6D719", "#176EBF", "#00A8DE", "#AEE0E8",
             "#00A9A3", "#FBD324", "#F28A24", "#A52828",
             "#A37CB7", "#F2D7EE", "#CD6981", "#FBD0C0",
             "#F15E4C", "#ECB2C8", "#B2DBBF", "#CCDAD1")
  } else {
    col <- c("#B8E3EA", "#5CB3DA", "#0070B2", "#FBDD7E", "#F7AE24", "#FF7149", 
             "#F2D7EE", "#A37CB7", "#A231A1", "#ECB2C8", "#E93B8C", "#B91372", 
             "#FF9F99", "#F15E4C", "#DA1735", "#CDE391", "#8BBE53", "#679436", 
             "#98D4C6", "#00A385", "#067D69", "#B2DBBF", "#028090", "#114B5F", 
             "#FBD0C0", "#CD6981", "#A23E48", "#CCDAD1", "#9CAEA9", "#788585")
    col <- colorRampPalette(col)(number)
  }
  
  DimPlot(data, pt.size = pt.size, label = label, repel = T, 
          raster = FALSE, label.size = 5, reduction = reduction,
          group.by = group.by) + 
    scale_color_manual(values = c(col)) + 
    
    theme_bw() +
    theme(legend.key.size = unit(0.5,'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          axis.title = element_text(colour = "black", size = 15),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1)) +   
    
    labs(x = 'UMAP1',y= 'UMAP2',title = '')
}


# dotplot marker ====================

dotplot_marker <- function(data = datafilt,
                           group.by = "cluster",
                           marker = select,
                           species = "human",
                           width = 15,
                           height = 6,
                           output = name){
  
  library(Seurat)
  library(ggpubr)
  
  if (is.null(species)) {
    
    gene_list <- marker
    
  } else if (species == "human") {
    
    Tcell = c("CD3D", "CD3E")
    CD8T = c("CD8A", "CD8B")
    CD4T = c("CD4")
    Treg = c("FOXP3", "IL2RA", "IL7R")
    NKcell = c("KLRD1", "NKG7", "NCAM1")
    Bcell = c("CD19", "CD79A", "MS4A1")
    Plasma = c("IGKC", "IGLC2", "TNFRSF17")
    Macro = c("APOC1", "SPP1", "C1QC")
    Mono = c("FCN1", "S100A8", "S100A9")
    cDC = c("FCER1A", "CD207", "XCR1")
    pDC = c("IL3RA", "LILRA4")
    MAST = c("KIT", "MS4A2")
    Neu = c("FCGR3B", "CEACAM8", "CSF3R")
    Fibro = c("COL1A1", "COL1A2")
    Endo = c("PECAM1", "VWF")
    hepatocyte = c("ALB", "CYP3A4", "HNF4A")
    biliary = c("EPCAM", "KRT7", "KRT19")
    
    gene_list <- list(Tcell = Tcell,
                      CD8T = CD8T,
                      CD4T = CD4T,
                      Treg = Treg,
                      NKcell = NKcell,
                      Bcell = Bcell,
                      Plasma = Plasma,
                      Macro = Macro,
                      Mono = Mono,
                      cDC = cDC,
                      pDC = pDC,
                      MAST = MAST,
                      Neu = Neu,
                      Fibro = Fibro,
                      Endo = Endo,
                      hepatocyte = hepatocyte,
                      biliary = biliary)
    
  } else if (species == "mouse") {
    
    Tcell = c("Cd3d", "Cd3e")
    CD8T = c("Cd8a", "Cd8b1")
    CD4T = c("Cd4")
    Treg = c("Foxp3", "Il2ra", "Il7r")
    NKcell = c("Klrd1", "Nkg7", "Ncam1")
    Bcell = c("Cd19", "Cd79a", "Ms4a1")
    Plasma = c("Igkc", "Tnfrsf17")
    Macro = c("Apoc1", "Spp1", "C1qc")
    Mono = c("S100a8", "S100a9")
    cDC = c("Fcer1a", "Cd207", "Xcr1")
    pDC = c("Il3ra", "Gm14548")
    MAST = c("Kit", "Ms4a2")
    Neu = c("Csf3r", "Fut4")
    Fibro = c("Col1a1", "Col1a2")
    Endo = c("Pecam1", "Vwf")
    hepatocyte = c("Alb", "Cyp3a11", "Hnf4a")
    biliary = c("Epcam", "Krt7", "Krt19")
    
    gene_list <- list(Tcell = Tcell,
                      CD8T = CD8T,
                      CD4T = CD4T,
                      Treg = Treg,
                      NKcell = NKcell,
                      Bcell = Bcell,
                      Plasma = Plasma,
                      Macro = Macro,
                      Mono = Mono,
                      cDC = cDC,
                      pDC = pDC,
                      MAST = MAST,
                      Neu = Neu,
                      Fibro = Fibro,
                      Endo = Endo,
                      hepatocyte = hepatocyte,
                      biliary = biliary)
    
  }
  
  col <- colorRampPalette(c("#0070b2","#009bc7","#b8e3ea",
                            "#f3f3f1","#fccdb9", "#f15e4c","#da1735"))(100)
  

  plot <- DotPlot(data, scale = T, col.min = -1,
                  group.by = group.by, col.max = 1, features = gene_list) + 
    scale_color_gradientn(colors = col) + 
    theme_bw()+
    theme(legend.position = "right", legend.box = "vertical",
          legend.margin = margin(t = 0, unit='cm'),
          panel.grid = element_blank(),
          axis.text = element_text(color = "black", size = 12),
          legend.text = element_text(size = 12,color = "black"),
          legend.title = element_text(size = 12,color = "black")) + 
    labs(x = '', y = '', title = '') + rotate_x_text(45)
  
  if (is.null(output)) {
    plot
  } else {
    ggsave(output, plot, width = width, height = height)
  }
  
}


# UMAP density plot ====================

prop_density <- function(datafilt = datafilt,
                         group = "group",
                         coord = "umap") {
  library(Seurat)
  library(ggplot2)
  library(viridis)
  
  info <- datafilt@meta.data
  input <- data.frame(Embeddings(datafilt, coord))
  input$group <- info[,group]
  
  ggplot(input, aes(x = input[,1],y = input[,2])) +
    stat_density_2d(aes(fill = ..density..),
                    geom = "raster", h = 2, contour = FALSE) + 
    geom_density2d(size = 0.1, colour = "#FDAF9199",
                   alpha = 0.3, bins = 15, h = 3.5) +
    scale_fill_viridis(option="magma") + 
    
    theme_bw() + 
    labs(x = 'UMAP1',y= 'UMAP2',title = '') + 
    theme(strip.background = element_rect(colour = NA, fill = 'grey90'),
          strip.text.x = element_text(size = 15), 
          axis.title = element_text(colour = "black", size = 15), 
          legend.key.size = unit(0.75, 'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_blank(), 
          panel.grid = element_blank(),
          strip.placement = 'outside',
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) + 
    facet_wrap(~group, nrow = 1)
}


# Back-to-back plot ====================

prop_back2back <- function(datafilt = datafilt,
                           group = "group",
                           cluster = "celltype",
                           order = TRUE) {
  library(Seurat)
  library(ggplot2)
  library(gtools)
  
  info <- datafilt@meta.data
  info <- data.frame(group = info[,group],
                     cluster = info[,cluster])
  
  input <- data.frame(table(info$group, info$cluster))
  names(input) <- c("group", "cluster", "prop")
  
  input <- do.call(rbind, lapply(unique(input$group), function(i){
    subinput <- input[input$group == i,]
    subinput$prop <- subinput$prop / sum(subinput$prop)
    subinput}))
  
  input$label <- paste0(round(input$prop * 100, 0), "%")
  
  if (order == TRUE) {
    name <- unique(mixedsort(input$cluster, decreasing = TRUE))
    input$cluster <- factor(input$cluster, levels = name)
  }
  
  ggplot(input, aes(x = cluster, fill = group)) + 
    scale_fill_manual(values = c("#A231A1", "#F2D7EE")) + 
    
    geom_bar(stat = "identity",
             data = subset(input, group == unique(input$group)[1]),
             aes(y = prop)) +
    geom_text(data = subset(input, group == unique(input$group)[1]), 
              aes(y = prop, label = label), size = 5, hjust = -0.1) +
    
    geom_bar(stat = "identity",
             data = subset(input, group == unique(input$group)[2]),
             aes(y = prop * (-1)) ) +
    geom_text(data = subset(input, group == unique(input$group)[2]), 
              aes(y = prop * (-1), label = label), size = 5, hjust = 1.1) +
    
    theme_bw() +
    theme(legend.key.size = unit(1,'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(colour = "black", size = 15),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1)) + 
    
    coord_flip() + ylab("Proportion") + xlab("") + 
    ylim(-max(input$prop) - 0.035, max(input$prop) + 0.035)
}


# Sample-level proportional distribution difference ====================

prop_plot_hca <- function(input = input,
                          rotate = 45,
                          decreasing = T,
                          species = "human") {
  
  library(ggplot2)
  library(ggpubr)
  library(gtools)
  
  # 变量排序
  
  input$Var2 <- as.character(input$Var2)
  name <- unique(mixedsort(input$Var2, decreasing = decreasing))
  input$Var2 <- factor(input$Var2, levels = name)
  
  # 设置颜色
  
  if (species == "human") {
    
    col <- c("CD8T" = "#B8E3EA",
             "CD4T" = "#0070B2",
             "Treg" = "#A6D719",
             "NKcell" = "#F15E4C",
             "Bcell" = "#A231A1",
             "DC" = "#D5B81B",
             "Macro" = "#028090",
             "Mono" = "#B2DBBF",
             "MAST" = "#B32226",
             "Neutrophil" = "#00A9A3",
             "Fibroblast" = "#FBD324",
             "Endothelial" = "#CD6981",
             "tumor" = "#9CAEA9",
             "None" = "#F6F5BD")
    
  } else if (species == "mouse") {
    
    col <- c("CD8T" = "#B8E3EA",
             "CD4T" = "#0070B2",
             "Treg" = "#A6D719",
             "NKcell" = "#F15E4C",
             "Bcell" = "#A231A1",
             "DC" = "#D5B81B",
             "MoMac" = "#B2DBBF", 
             "MAST" = "#B32226",
             "Neutrophil" = "#00A9A3",
             "Fibroblast" = "#FBD324",
             "Endothelial" = "#CD6981",
             "tumor" = "#9CAEA9",
             "None" = "#F6F5BD")
    
  }
  
  # 开始画图
  
  ggplot(input, aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity", position = "fill", width = 0.75) + 
    scale_fill_manual(values = col)+
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.key.size = unit(0.75,'cm'),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(colour = "black", size = 15)) + 
    scale_y_continuous(expand = c(0,0)) + 
    labs(x = '', y = 'Proportion', title = '') + 
    rotate_x_text(rotate)
}

