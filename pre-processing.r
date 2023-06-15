library(Seurat)
library(tidyverse)
library(glmGamPoi)

# read data
my_samples <- base::paste0('./data/raw/P', c(1:12), '/filtered_feature_bc_matrix')
names(my_samples) <- base::paste0('P', c(1:12))
parab.data <- Seurat::Read10X(
  data.dir = my_samples
)


# assemble Seurat object
parab.ssc <- Seurat::CreateSeuratObject(
  counts = parab.data,
  project = 'parabiosis',
  min.cells = 1,
  min.features = 1,
  names.delim = '_',
  names.field = 1
)
rm(parab.data)


# add metadata
parab.ssc[['exp.group']] <- dplyr::case_when(
  parab.ssc@meta.data$orig.ident == 'P1'  ~ 'Aged (Y-A)',
  parab.ssc@meta.data$orig.ident == 'P2'  ~ 'Young (Y-A)',
  parab.ssc@meta.data$orig.ident == 'P3'  ~ 'Aged (Y-A)',
  parab.ssc@meta.data$orig.ident == 'P4'  ~ 'Young (Y-A)',
  parab.ssc@meta.data$orig.ident == 'P5'  ~ 'Aged (Y-A)',
  parab.ssc@meta.data$orig.ident == 'P6'  ~ 'Young (Y-A)',
  parab.ssc@meta.data$orig.ident == 'P7'  ~ 'Aged (Y-A)',
  parab.ssc@meta.data$orig.ident == 'P8'  ~ 'Young (Y-A)',
  parab.ssc@meta.data$orig.ident == 'P9'  ~ 'Young (Y-Y)',
  parab.ssc@meta.data$orig.ident == 'P10' ~ 'Young (Y-Y)',
  parab.ssc@meta.data$orig.ident == 'P11' ~ 'Aged (A-A)',
  parab.ssc@meta.data$orig.ident == 'P12' ~ 'Aged (A-A)'
)
parab.ssc[['age']] <- base::ifelse(
  test = stringr::str_detect(string = parab.ssc@meta.data$exp.group, patter = 'Aged'),
  yes = 'aged',
  no = 'young'
)
parab.ssc[['parabiont']] <- dplyr::case_when(
  parab.ssc@meta.data$exp.group == 'Aged (A-A)'  ~ 'aged parabiont',
  parab.ssc@meta.data$exp.group == 'Young (Y-A)' ~ 'aged parabiont',
  parab.ssc@meta.data$exp.group == 'Young (Y-Y)' ~ 'young parabiont',
  parab.ssc@meta.data$exp.group == 'Aged (Y-A)'  ~ 'young parabiont'
)
parab.ssc[['batch']] <- dplyr::case_when(
  parab.ssc@meta.data$orig.ident %in% c('P1','P2','P3','P4') ~ 'batch 1',
  parab.ssc@meta.data$orig.ident %in% c('P5','P6','P7','P8') ~ 'batch 2',
  parab.ssc@meta.data$orig.ident %in% c('P9','P10') ~ 'batch 3',
  parab.ssc@meta.data$orig.ident %in% c('P11','P12') ~ 'batch 4'
)
parab.ssc[['percent.mt']] <- Seurat::PercentageFeatureSet(object = parab.ssc, pattern = '^mt-')


# QC
Seurat::VlnPlot(object = parab.ssc, features = 'percent.mt')
Seurat::VlnPlot(object = parab.ssc, features = 'nCount_RNA')
Seurat::VlnPlot(object = parab.ssc, features = 'nCount_RNA') + ggplot2::ylim(0, 2000)
Seurat::VlnPlot(object = parab.ssc, features = 'nFeature_RNA') 
Seurat::VlnPlot(object = parab.ssc, features = 'nFeature_RNA') + ggplot2::ylim(0, 2000)
parab.ssc <- base::subset(
  x = parab.ssc,
  subset = percent.mt < 15 & nCount_RNA < 19000 & nCount_RNA > 1000 & nFeature_RNA < 5000 & nFeature_RNA > 500  
)
parab.ssc


# normalize
parab.ssc <- Seurat::SCTransform(
  object = parab.ssc,
  method = 'glmGamPoi',
  vars.to.regress = 'percent.mt',
  conserve.memory	= T,
  return.only.var.genes = T,
  vst.flavor = 'v2'
)


# dimension reduction
parab.ssc <- Seurat::RunPCA(object = parab.ssc, features = Seurat::VariableFeatures(object = parab.ssc))
Seurat::ElbowPlot(object = parab.ssc, ndims = 30)
parab.ssc <- Seurat::RunTSNE(object = parab.ssc, dims = 1:20)
parab.ssc <- Seurat::RunUMAP(object = parab.ssc, dims = 1:20)


# clustering
parab.ssc <- Seurat::FindNeighbors(object = parab.ssc, dims = 1:20)
parab.ssc <- Seurat::FindClusters(object = parab.ssc, resolution = 0.05)
Seurat::DimPlot(object = parab.ssc, group.by = 'seurat_clusters')


# name clusters
ggpubr::ggarrange(
  plotlist = list(
    Seurat::DimPlot(object = parab.ssc),
    Seurat::FeaturePlot(object = parab.ssc, features = c(
      'Cldn5',    # endothelial cells
      'Hexb',     # microglias
      'Acta2',    # smooth muscle cells
      'Pdgfrb',   # pericytes
      'Cldn11',   # oligodendrocytes
      'Slc1a2'    # astrocytes
      )
    )
  )
)

EC  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = parab.ssc, reduction = 'umap'))
EC  <- data.frame(bc = EC, cluster = 'EC')
MG  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = parab.ssc, reduction = 'umap'))
MG  <- data.frame(bc = MG, cluster = 'MG')
SMC <- Seurat::CellSelector(plot = Seurat::DimPlot(object = parab.ssc, reduction = 'umap'))
SMC <- data.frame(bc = SMC, cluster = 'SMC')
PC  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = parab.ssc, reduction = 'umap'))
PC  <- data.frame(bc = PC, cluster = 'PC')
OC  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = parab.ssc, reduction = 'umap'))
OC  <- data.frame(bc = OC, cluster = 'OC')
AS  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = parab.ssc, reduction = 'umap'))
AS  <- data.frame(bc = AS, cluster = 'AS')

manual_clusters <- rbind(EC, MG, SMC, PC, OC, AS) %>% tibble::column_to_rownames('bc')
parab.ssc[['seurat_clusters']] <- manual_clusters
Seurat::Idents(parab.ssc) <- parab.ssc[['seurat_clusters']]
parab.ssc <- base::subset(x = parab.ssc, idents = c('EC', 'MG', 'SMC', 'PC', 'OC', 'AS'))


# save seurat object
base::saveRDS(
  object = parab.ssc,
  file = './data/single_cell_parab.rda'
)

