library(Seurat)
library(tidyverse)
library(glmGamPoi)

# read data
my_samples <- base::paste0('./data/raw/my_samples_', c(1:6), '/filtered_feature_bc_matrix')
names(my_samples) <- base::paste0('my_samples_', c(1:6))
my.data <- Seurat::Read10X(
  data.dir = my_samples
)


# assemble Seurat object
my.ssc <- Seurat::CreateSeuratObject(
  counts = my.data,
  project = 'test_project',
  min.cells = 1,
  min.features = 1,
  names.delim = '_',
  names.field = 1
)
rm(my.data)


# add metadata
my.ssc[['exp.group']] <- dplyr::case_when(
  my.ssc@meta.data$orig.ident == 'my_samples_1'  ~ 'my_group_1',
  my.ssc@meta.data$orig.ident == 'my_samples_2'  ~ 'my_group_1',
  my.ssc@meta.data$orig.ident == 'my_samples_3'  ~ 'my_group_2',
  my.ssc@meta.data$orig.ident == 'my_samples_4'  ~ 'my_group_2',
  my.ssc@meta.data$orig.ident == 'my_samples_3'  ~ 'my_group_3',
  my.ssc@meta.data$orig.ident == 'my_samples_3'  ~ 'my_group_3',
)
my.ssc[['age']] <- base::ifelse(
  test = stringr::str_detect(string = my.ssc@meta.data$exp.group, patter = 'my_group_1'),
  yes = 'aged',
  no = 'young'
)
my.ssc[['batch']] <- dplyr::case_when(
  my.ssc@meta.data$orig.ident %in% c('my_group_1', 'my_group_1') ~ 'batch 1',
  my.ssc@meta.data$orig.ident %in% c('my_group_3') ~ 'batch 2'
)
my.ssc[['percent.mt']] <- Seurat::PercentageFeatureSet(object = my.ssc, pattern = '^mt-')


# QC
Seurat::VlnPlot(object = my.ssc, features = 'percent.mt')
Seurat::VlnPlot(object = my.ssc, features = 'nCount_RNA')
Seurat::VlnPlot(object = my.ssc, features = 'nCount_RNA') + ggplot2::ylim(0, 2000)
Seurat::VlnPlot(object = my.ssc, features = 'nFeature_RNA') 
Seurat::VlnPlot(object = my.ssc, features = 'nFeature_RNA') + ggplot2::ylim(0, 2000)
my.ssc <- base::subset(
  x = my.ssc,
  subset = percent.mt < 15 & nCount_RNA < 19000 & nCount_RNA > 1000 & nFeature_RNA < 5000 & nFeature_RNA > 500  
)
my.ssc


# normalize
my.ssc <- Seurat::SCTransform(
  object = my.ssc,
  method = 'glmGamPoi',
  vars.to.regress = 'percent.mt',
  conserve.memory	= T,
  return.only.var.genes = T,
  vst.flavor = 'v2'
)


# dimension reduction
my.ssc <- Seurat::RunPCA(object = my.ssc, features = Seurat::VariableFeatures(object = my.ssc))
Seurat::ElbowPlot(object = my.ssc, ndims = 30)
my.ssc <- Seurat::RunTSNE(object = my.ssc, dims = 1:20)
my.ssc <- Seurat::RunUMAP(object = my.ssc, dims = 1:20)


# clustering
my.ssc <- Seurat::FindNeighbors(object = my.ssc, dims = 1:20)
my.ssc <- Seurat::FindClusters(object = my.ssc, resolution = 0.05)
Seurat::DimPlot(object = my.ssc, group.by = 'seurat_clusters')


# name clusters
ggpubr::ggarrange(
  plotlist = list(
    Seurat::DimPlot(object = my.ssc),
    Seurat::FeaturePlot(object = my.ssc, features = c(
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

EC  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = my.ssc, reduction = 'umap'))
EC  <- data.frame(bc = EC, cluster = 'EC')
MG  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = my.ssc, reduction = 'umap'))
MG  <- data.frame(bc = MG, cluster = 'MG')
SMC <- Seurat::CellSelector(plot = Seurat::DimPlot(object = my.ssc, reduction = 'umap'))
SMC <- data.frame(bc = SMC, cluster = 'SMC')
PC  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = my.ssc, reduction = 'umap'))
PC  <- data.frame(bc = PC, cluster = 'PC')
OC  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = my.ssc, reduction = 'umap'))
OC  <- data.frame(bc = OC, cluster = 'OC')
AS  <- Seurat::CellSelector(plot = Seurat::DimPlot(object = my.ssc, reduction = 'umap'))
AS  <- data.frame(bc = AS, cluster = 'AS')

manual_clusters <- rbind(EC, MG, SMC, PC, OC, AS) %>% tibble::column_to_rownames('bc')
my.ssc[['seurat_clusters']] <- manual_clusters
Seurat::Idents(my.ssc) <- my.ssc[['seurat_clusters']]
my.ssc <- base::subset(x = my.ssc, idents = c('EC', 'MG', 'SMC', 'PC', 'OC', 'AS'))


# save seurat object
base::saveRDS(
  object = my.ssc,
  file = './data/single_cell_my_project.rda'
)

