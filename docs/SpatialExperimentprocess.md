## Spatial Experiment process

Spatial Experiment object is a S4 type of data storage, the following is for the use of SPE data, but it is worth mentioning that SpatialExperiment provides a very small number of sample data spots, the display is not good.The following also includes the use of 10x public data [Mouse Brain Serial Section 2 (Sagittal-Anterior)](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard), read by SpatialExperiment to generate SPE objects:

``` r
## Spatial Experiment sample data
library(SpatialExperiment)
example(read10xVisium, echo = FALSE)


## 10x Visium load in SPE
spe2 <- SpatialExperiment::read10xVisium(
  samples = 'anterior/',
  sample_id = "anterior",
  type = "sparse",
  data = "filtered",
  images = "lowres",
  load = TRUE
)

dir.create('spe')
setwd('spe')
CreatshinySRT(dat = spe2,title = 'spatial experiment',gex.assay = 'counts')
```