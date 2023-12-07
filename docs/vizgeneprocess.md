## vizgene process

Spatial transcriptome single-cell precision is a development trend, vizgene is a kind of single-cell precision spatial transcriptome platform that has been commercially available, the good point is that the data can be processed by the seurat, the following is the example code, the sample data is the vizgene public platform data [vizgene](https://console.cloud.google.com/storage/browser/public-datasets-vizgen-merfish;tab=objects?prefix=&forceOnObjectsSortingFiltering=false)

``` r
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(future)
plan("multisession", workers = 10)
library(ggplot2)
library(shinySRT)

vizgen.obj <- LoadVizgen(data.dir = "vizgen/", fov = "vizgene")

vizgen.obj <- SCTransform(vizgen.obj, assay = "Vizgen", clip.range = c(-10, 10))
vizgen.obj <- RunPCA(vizgen.obj, npcs = 30, features = rownames(vizgen.obj))
vizgen.obj <- RunUMAP(vizgen.obj, dims = 1:30)
vizgen.obj <- FindNeighbors(vizgen.obj, reduction = "pca", dims = 1:30)
vizgen.obj <- FindClusters(vizgen.obj, resolution = 0.3)

dir.create('shinySRT')
setwd('shinySRT')

CreateshinySRT(vizgen.obj,title = 'Seurat vizgene')
```