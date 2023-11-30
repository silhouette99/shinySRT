## SingleCellExperiment process

SingleCellExperiment (SCE) and SpatialExperiment (SPE) are the same S4 type of storage object, SCE is generally the storage object for single cell data, in the example data of spatialLIBD, we found that there is also the SCE object for spatial transcriptomes[SPE data](https://www.dropbox.com/s/f4wcvtdq428y73p/Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata?dl=1), the following is the example code:

``` r
library(spatialLIBD)

# spe <- fetch_data(type = "spe")
load(
  'Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata'
)

## 
sce2 <- sce[, which(sce$sample_name %in% c(151507, 151508))]
dat <- sce2
dat@colData <- dat@colData[,c("sample_name","tissue","imagerow","imagecol","Cluster","position","subject_position")]

dir.create('sce')
setwd('sce')

CreatshinySRT(dat,title = 'single cell experiment',maxlevel = 20)

```