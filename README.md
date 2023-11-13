# ShinySRT

---

### Description

`ShinySRT` is a Shiny-based web application developed for the analysis of spatially resolved transcriptomics data. This application is designed to handle multiple formats of spatial transcriptome data and allows for the creation of an interactive interface that supports comprehensive data analysis. The interactive interface provided by ShinySRT is open-source and can be highly customized to meet the specific needs of users.




---

### Features

- Developed under R, utilizing a Shiny application that generates an interactive interface deployable on a server or shareable via the web.
- Provides compatibility with various prominent formats of spatial transcriptome data.
- Facilitates the import of multiple ST datasets into the Shiny web application.
- Allows for the customization of spatial spot annotations and spot selection criteria.
- Supports multivariate comparisons, enabling the analysis of how the dependent variable changes across different groups of independent variables.
- Offers features for visualizing images and downloading data sheets.
- Simplifies the creation of the Shiny interface in a single step, and the Shiny app is entirely open source and customizable.


| Data      | Object Type | source |
| ----------- | ----------- | ----------- |
| 10x Visim (Seurat)     | Seurat       | [seurat obj](https://www.10xgenomics.com/resources/datasets?menu%5Bproducts.name%5D=Spatial%20Gene%20Expression&query=&page=1&configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000) |
| SingleCellExperiment   | SCE        | [SCE obj](docs/SingleCellExperimentprocess.md) |
| SpatialExperiment   | SPE        | [SPE obj](docs/SpatialExperimentprocess.md) |
| Vizgene (Seurat)  | Seurat        | [viz_seurat obj](docs/vizgeneprocess.md) |
| 10x Visim (scanpy)   | h5ad        | [h5ad class](docs/scanpyprocess.md) |
| Customizable list   | list        | [lists](docs/customlistprocess.md) |


Users could build a customized list containing expression matrix, metadata, image as well as coordinate information.

---

### Installation
To begin, it's important to verify whether the necessary installation packages for `ShinySRT` have already been installed:

``` r
packages <- c(
  'SingleCellExperiment',
  'SpatialExperiment',
  'data.table',
  'dplyr',
  'glue',
  'hdf5r',
  'readr',
  'reticulate',
  'ggplot2',
  'graphics',
  'gridExtra',
  'patchwork',
  'RColorBrewer',
  'maps',
  'Cairo',
  'grid',
  'ggtree',
  'aplot',
  'magrittr',
  'ggrepel',
  'ggdendro',
  'Matrix',
  'scales',
  'aplot',
  'keys'
)
packages = packages[!(packages %in% installed.packages()[, "Package"])]

if (length(packages)) {
  install.packages(newPkg)
}

```

Then check that all required shiny-related packages are installed:

``` r
packages <- c('shiny','shinyhelper','DT','shinydashboard')
packages = packages[!(packages %in% installed.packages()[,"Package"])]
if (length(packages)) {
  install.packages(newPkg)
}
```


---

### Content and Guide

The fundamental process involves `ShinySRT` generating the essential configuration files and web application file using spatial transcriptome data objects. The precise operational code is outlined as follows:

``` r
# 10x Visim
library(ShinySRT)
library(Seurat)
library(SeuratData)

InstallData("stxBrain")
brain <- LoadData("stxBrain"， type = "anterior1")

makespashiny(brain,title = 'ShinySRT exmaple')

# SpatialExperiment
library(SpatialExperiment)
example(read10xVisium, echo = FALSE)

makespashiny(spe,title = 'ShinySRT exmaple',gex.assay = 'counts')

## run shiny app
shiny::runApp('shinyspatial_app/')
```

The ST was processed using scanpy to obtain the h5ad file, while the following URL was used to access the source data from [10X](https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-2-sagittal-anterior-1-standard).


``` r
# h5ad
makespashiny(dat = 'Anterior.h5ad',title = 'spatial experiment')
## run shiny app
shiny::runApp('shinyspatial_app/')
```

Upon running a single line of code, a new directory named `/shinyspatial_app` will be generated in the current directory, where the Shiny app is located. Users can utilize `shiny::runApp` to locally run the app within R env. Furthermore, the app can be deployed remotely by placing the shiny app's directory into the `/srv/shiny-server` directory of a server that has a proxy. It's worth noting that Shiny apps can also be deployed on various web platforms using alternative methods. For more comprehensive information, please refer to [shinyapps.io](https://www.shinyapps.io/).

We use a 10x spatial transcriptome data within the lab as an example to demonstrate content of `ShinySRT`.

The Shiny app created by `ShinySRT` comprises six primary modules, as indicated by the module names in the leftmost menu bar in the figure. At the bottom of the menu bar (highlighted in an orange box), there is an input field for importing annotations for new spots (highlighted in a red box).


The current view represents the initial module titled "SpotInfo vs GeneExpr," primarily illustrating the connection between spatial spot annotations and gene expression. You can switch the spot annotations using the dropdown menu labeled "Spot information" and choose the displayed genes from the dropdown menu labeled "Gene expression." On the right-hand side, there is a section for personalized spot selection, allowing you to pick an area either by drawing a circle or by clicking. Remember to provide a name for the selected annotation area before proceeding with the selection.


![](image/content1.png)


For storing new annotation information, you need to first create a column name for this group of annotation information, and then stored in the form of columns in the following figure in the metadata list, the list can be downloaded and exported, and then use this annotation can be copied from the column and pasted in the left input area. metadata on the right side of the gene in the display of the expression of the spots information statistics.


![](image/content2.png)


The second module "GeneExpr vs GeneExpr" focuses on the relationship between the spatial expression of two genes, in the same way as the selection of spots in the first module.


![](image/content3.png)


In addition to this, Module 4 demonstrates the relationship between two genes by the spatial co-determination of two genes in the same slice sample. Thresholds for gene expression can be set using the sliders on the right, and the table below shows the statistics for the co-determination of the genes.


![](image/content4.png)


The third module, "Viol / Box data chart", plots traditional statistical box and violin plots based on the annotation information and gene expression or scoring of the spot, showing the relationship between the annotation information and genes.


![](image/content5.png)


The fourth module, "Portion data chart", shows the proportion of one grouping of information over another, e.g., the proportion or number of seurat subgroups in each anatomical region of the brain, with most of the spots in the HT region belonging to the four subgroups.


![](image/content6.png)


The sixth module shows the expression of genes in different regions by means of "bubble charts or heatmaps", which can also be clustered in rows and columns respectively.


![](image/content7.png)


---
