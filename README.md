# Single Cell 10X Workflow v1


### Project folders:
```bash
my_project
|-- metadata.md
|-- R
|   |-- pre-processing.R
|   |-- sub-cluster_analysis.R
|   |-- downstream_analysis_1.R
|   |-- downstream_analysis_2.R
|   |-- ...
|-- data
|   |-- raw_data
|   |-- gene_sets
|   |-- single_cell_seurat_object.rda
|   |-- single_cell_seurat_object_cell_type_1.rda
|   |-- single_cell_seurat_object_cell_type_2.rda
|   |-- ...
|-- plot
```


### Workflow steps:
1. raw_data -> pre-processing.R
2. pre-processing.R -> single_cell_seurat_object.rda
3. single_cell_seurat_object.rda -> sub-cluster_analysis.R
4. sub-cluster_analysis.R -> single_cell_seurat_object_cell_type_1.rda and single_cell_seurat_object_cell_type_2.rda ...
5. single_cell_seurat_object.rda -> downstream_analysis_1.R
6. single_cell_seurat_object.rda -> downstream_analysis_2.R
7. ...


### Styleguide
- naming convention:
  - project name: [costumer]_[animal]_[orgnan]_[treatment]_10X_single_cell
  - Seurat object names: 
  - plot names:
- Don’t use attach()!
- Librarys to the header!
- Explicitly qualify namespaces!
- Use explicit returns!
- https://google.github.io/styleguide/Rguide.html


### Pre-processing parameteres to set:
- quality control:
  - number of mitochondria genes per cell
  - number of features per cell
  - number of reads per cell
- normalizing and scaling:
  - method: LogNormalize vs. SCTransform
- dimension reduction:
  - method: t-SNE, UMAP etc...
  - number of relevant dimensions
- clustering cell types:
  - number of relevant dimensions
  - resolution


### Meta data:
- short narrative introduction
- experimental groups
- treatment descritpion
- cell isolation workflow
