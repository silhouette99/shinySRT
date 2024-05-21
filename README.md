# shinySRT

---
![](image/shinysrt.png)

## Description


`shinySRT` is a web application developed utilizing the Shiny framework, explicitly designed for the sharable and interactive visualization of spatially resolved transcriptomics data. This application is adept at processing various formats of spatial transcriptome data, facilitating the development of an interactive interface conducive to thorough data analysis. The interactive interface of shinySRT is open-source, offering significant customization potential to align with the unique requirements of users.




---

## Features

- Developed under R, utilizing a Shiny application that generates an interactive and shareable interface via the web.
- Provides compatibility with various prominent formats of spatial transcriptome data.
- Facilitates the import of multiple ST datasets into the Shiny web application.
- Allows for the customization of spatial spot selection.
- Supports multivariate comparisons, enabling the analysis of how the dependent variable changes across different groups of independent variables.
- Offers features for visualizing images and downloading data sheets.
- Simplifies the creation of the Shiny interface in a single step, and the Shiny app is entirely open source and customizable.


| Data      | Object Type | source |
| ----------- | ----------- | ----------- |
| 10x Visim     | Seurat       | [seurat obj](https://www.10xgenomics.com/resources/datasets?menu%5Bproducts.name%5D=Spatial%20Gene%20Expression&query=&page=1&configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000) |
| SingleCellExperiment   | SCE        | [SCE obj](docs/SingleCellExperimentprocess.md) |
| SpatialExperiment   | SPE        | [SPE obj](docs/SpatialExperimentprocess.md) |
| Vizgene  | Seurat        | [viz_seurat obj](docs/vizgeneprocess.md) |
| 10x Visim   | h5ad        | [h5ad class](docs/scanpyprocess.md) |
| Customizable list   | list        | [lists](docs/customlistprocess.md) |


Users could deploy their application utilizing a customized list that includes an expression matrix, metadata, imagery, as well as coordinate information.

---

## Installation

To begin, it's important to verify whether the necessary installation packages for `shinySRT` have already been installed (The best version of R is 4.2 or above, to prevent some R packages are not good to install or incompatibility):

``` r

# If you employ Conda, it is imperative to install the essential packages requisite for the analysis of single-cell and spatial genomic data.
# mamba install conda-forge::r-cairo conda-forge::r-hdf5r conda-forge::r-curl conda-forge::r-devtools conda-forge::r-stringi conda-forge::r-biocmanager conda-forge::r-rfast conda-forge::quadprog -y

if (!require('pacman')) install.packages('pacman')

BiocManager::install(c('scran','Seurat'))


devtools::install_github(c('YuLab-SMU/ggtree', 'silhouette99/shinySRT'))
# remotes::install_github(c('YuLab-SMU/ggtree', 'silhouette99/shinySRT'))
# install.package('Seurat')

pacman::p_load(
  'Cairo',
  'hdf5r',
  'data.table',
  'magrittr',
  'dplyr',
  'glue',
  'readr',
  'reticulate',
  'ggplot2',
  'graphics',
  'gridExtra',
  'patchwork',
  'RColorBrewer',
  'maps',
  'grid',
  'ggtree',
  'ggrepel',
  'ggdendro',
  'Matrix',
  'scales',
  'aplot',
  'keys',
  'ggiraph',
  'ggpubr',
  'shiny',
  'shinyhelper',
  'DT',
  'shinydashboard',
  'scran',
  'scatterpie',
  'quadprog',
  'Rfast',
  'SingleCellExperiment',
  'SpatialExperiment'
)

```



---

## Quick Start

The primary procedure entails `shinySRT` producing the necessary configuration files and web application file by employing SRT data objects. The specific operational code is delineated as follows:

### For a Seurat object

``` r
# 10x Visim
suppressMessages(library(shinySRT))
library(Seurat)
library(SeuratData)

InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

CreateshinySRT(brain,title = 'shinySRT exmaple', sp_normalize = F, gene.mapping = F)

# SpatialExperiment
library(SpatialExperiment)
example(read10xVisium, echo = FALSE)

CreateshinySRT(spe,title = 'shinySRT exmaple',maxlevel = 20)

# run shiny app
shiny::runApp('shinyspatial_app/')
```

### For a h5ad object generating by SCANPY

The SRT was processed using scanpy to obtain the h5ad file, while the following URL was used to access the source data from [10X](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard).


``` r
# h5ad
CreateshinySRT(dat = 'Anterior.h5ad',title = 'spatial experiment', sp_normalize = F)

# run shiny app
shiny::runApp('shinyspatial_app/')

```

Upon running a single line of code, a new directory named `/shinyspatial_app` will be generated in the current directory, where the Shiny app is located. Users can utilize `shiny::runApp` to locally run the app within R env. Furthermore, the app can be deployed remotely by placing the shiny app's directory into the `/srv/shiny-server` directory of a server that has a proxy. It's worth noting that Shiny apps can also be deployed on various web platforms using alternative methods. For more comprehensive information, please refer to [shinyapps.io](https://www.shinyapps.io/).

---

For detailed tutorial, please refer to [Tutorial](https://github.com/silhouette99/shinySRT-guide/blob/main/README.md). We also compared shinySRT with other SRT visualization tools in [this](https://github.com/silhouette99/shinySRT-guide/blob/main/shiny-display.md).

## Q&A

Visit [issues](https://github.com/silhouette99/ShinySRT/issues) or contact [Pan](https://github.com/silhouette99)


If you use the tool in your publication, please cite by




