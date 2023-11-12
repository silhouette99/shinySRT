## vizgene process

空间转录组单细胞精度是一个发展趋势，vizgene是一种已经商用的单细胞精度空间转录组平台，好的一点是数据可以通过seurat处理，下面是实例代码，示例数据为vizgene公共平台数据[vizgene](https://console.cloud.google.com/storage/browser/public-datasets-vizgen-merfish;tab=objects?prefix=&forceOnObjectsSortingFiltering=false)

``` r
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(future)
plan("multisession", workers = 10)
library(ggplot2)

vizgen.obj <- LoadVizgen(data.dir = "vizgen/", fov = "vizgene")

vizgen.obj <- SCTransform(vizgen.obj, assay = "Vizgen", clip.range = c(-10, 10))
vizgen.obj <- RunPCA(vizgen.obj, npcs = 30, features = rownames(vizgen.obj))
vizgen.obj <- RunUMAP(vizgen.obj, dims = 1:30)
vizgen.obj <- FindNeighbors(vizgen.obj, reduction = "pca", dims = 1:30)
vizgen.obj <- FindClusters(vizgen.obj, resolution = 0.3)

dir.create('shinySRT')
setwd('shinySRT')

makespashiny(vizgen.obj,title = 'Seurat vizgene')
```