## Deconvolusion 

Users simply need to furnish the expression matrix and metadata of scRNA-seq, specifying cell types, to acquire the cellular composition results of SRT. These outcomes can be visualized on the interface. The provided example leverages single-cell data from Giotto, utilizing 10x public data for spatial transcriptomics.

``` r
library(Seurat)
library(shinySRT)

scmtx <- data.table::fread('https://github.com/drieslab/spatial-datasets/raw/master/data/2022_scRNAseq_mouse_brain/count_matrix/brain_sc_expression_matrix.txt.gz')
genes <- scmtx$V1
scmtx <- scmtx[, -1] %>% as.matrix()
cells <- colnames(scmtx)
# scmtx <- as(scmtx, 'dgCMatrix')
rownames(scmtx) <- genes

scmeta <- read.csv('https://github.com/drieslab/spatial-datasets/raw/master/data/2022_scRNAseq_mouse_brain/cell_metadata/brain_sc_metadata.csv')

dir.create('shinySRT')
setwd('shinySRT')

dat <- readRDS('merge.Rds')

shinySRT::CreateshinySRT(
  dat,
  title = 'Seurat visim',
  sp_normalize = F,
  gene.mapping = F,
  scmtx = scmtx,
  scmeta = scmeta, 
  normalize = T, #scmtx was normalized
  colcluster = 'Class', #identify the column which have informations of cell type with scmeta 
  sp_cols = 'seurat_clusters'
)


shiny::runApp('shinySRT/shinyspatial_app')

```

![](image/content7-2.png)