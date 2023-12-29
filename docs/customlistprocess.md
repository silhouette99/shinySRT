## Customized Lists process
There are many types of spatial transcriptome data and it is difficult to cover all of them, for this you can construct a list on your own and use it to use `shinySRT`, the following is the code for the way and format of constructing the list, the data used is the public data from SPATA2 [SPATA 313_T](https://www.dropbox.com/s/iptcb9bqo5oo5xj/313_T.RDS?dl=1)

``` r
library(png)
library(SPATA2)
library(shinySRT)

## eg. SPATA2 obj
dat <- readRDS('313_T.RDS')
dat <- SPATA2::updateSpataObject(dat)

## background image
images <- list(dat@images$`313_T`@image)
names(images) <- unique(dat@fdata$`313_T`$sample)

## coordination
coordinates <- list(dat@coordinates$`313_T`[, c("barcodes", "x", "y")])
names(coordinates) <- unique(dat@fdata$`313_T`$sample)

coordinates <- lapply(names(coordinates), function(x) {
  cos <- coordinates[[x]]%>%as.data.frame()
  rownames(cos) <- paste(names(coordinates), cos$barcodes, sep = '_')
  cos <- cos[,-1]
  
  colnames(cos) <- c('y','x')
  cos <- cos[c(2,1)]
  ## convert x direction
  cos$x <- ncol(image[[x]]) - cos$x
  cos
})

coordinates <- do.call(rbind,coordinates)

## meta 
meta <- dat@fdata$`313_T` %>% as.data.frame()
rownames(meta) <- paste(meta$sample, meta$barcodes, sep = '_')
meta$slice_sample <- meta$sample

## matrix
data = dat@data$`313_T`$counts
colnames(data) <- paste(meta$sample, meta$barcodes, sep = '_')

## names
obj <- list(image = images,data = data,coordination = coordinates,metadata = meta,reduction = NULL)

dir.create('lists')
setwd('lists')

CreateshinySRT(obj,title = 'shiny list')
```