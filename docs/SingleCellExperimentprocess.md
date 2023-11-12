## SingleCellExperiment process

SingleCellExperiment (SCE) and SpatialExperiment (SPE) are the same S4 type of storage object, SCE is generally the storage object for single cell data, in the example data of spatialLIBD, we found that there is also the SCE object for spatial transcriptomes, the following is the example code:

``` r
library(spatialLIBD)

load(
  '/mnt/raid62/pzz/shiny_data/obj/Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata'
)

## 
sce2 <- sce[, which(sce$sample_name %in% c(151507, 151508))]
dat <- sce2
dat@colData <- dat@colData[,c("sample_name","tissue","imagerow","imagecol","Cluster","position","subject_position")]

dir.create('sce')
setwd('sce')

makespashiny(dat,title = 'single cell experiment',maxlevel = 20)

```