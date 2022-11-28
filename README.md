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
|   |-- ...
|-- plot
```

### Naming convention:
- project name:
- Seurat object names:
- plot names:
- 

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

### Meat data:
