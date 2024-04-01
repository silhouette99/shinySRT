#' Input the SRT data
#'
#' iuput SRT matrix to generate a seurat obj or SRT list
#'
#'
#' @param matx gene-barcode matrix file path, formation:.txt, .csv, .tsv, .xlsx
#' @param meta meta data file path, formation:.txt, .csv, .tsv, .xlsx
#' @param coordi file path, Spatial location information, formation:.txt, .csv, .tsv, .xlsx
#' @param image image file path, formation: .png, .jepg, .tiff. 
#'
#'
#' @import data.table hdf5r reticulate dplyr Seurat SummarizedExperiment scran dplyr scatterpie
#' 
#' 
#' @return object for process
#'
#'
#' @export
load_spatial <-
  function(matx = NULL,
           meta = NULL,
           coordi = NULL,
           image = NULL,
           species = 'hg',
           x_reverse = F,
           min_cells = 3,
           scale_factors_file = NULL) {
    if (length(grep(pattern = '.xlsx', matx)) > 0) {
      data <- read.xlsx(matx, rowNames = T) %>% as.data.frame()
    } else{
      data <- data.table::fread(matx) %>% as.data.frame()
    }
    tmp = read.table(system.file("extdata", paste0(species, '_map.txt.gz'), package = 'shinySRT'),
                     header = T)
    
    if (length(intersect(rownames(data), tmp$EnsemblID)) == 0 |
        length(intersect(rownames(data), tmp$GeneName)) == 0) {
      if (length(intersect(data[, 1], tmp$EnsemblID)) > 0 |
          length(intersect(data[, 1], tmp$GeneName)) > 0) {
        rownames(data) <- data[, 1] %>% make.unique()
        data  <- data[, -1]
      } else if (length(intersect(data[, ncol(data)], tmp$EnsemblID)) > 0 |
                 length(intersect(data[, ncol(data)], tmp$GeneName)) > 0) {
        rownames(data) <- data[, ncol(data)] %>% make.unique()
        data  <- data[, -ncol(data)]
      } else{
        stop('Can not find the geneID!')
      }
    }
    
    object <-
      Seurat::CreateSeuratObject(data, assay = 'Spatial', min.cells = min_cells)
    
    if (length(meta) > 0) {
      if (length(grep(pattern = '.xlsx', meta)) > 0) {
        coordination <- read.xlsx(meta, rowNames = T)
      } else {
        coordination <- data.table::fread(meta) %>% as.data.frame()
        rownames(coordination)  <-
          coordination[, grep(pattern = 'barcode|Barcode',
                              colnames(coordination),
                              value = T)]
      }
      
      object@meta.data <-
        cbind(object@meta.data, meta[match(rownames(object@meta.data), rownames(meta)), ])
    }
    
    
    if (length(coordi) == 0) {
      stop('Do not find coordination infotmation!')
    }
    
    if (length(grep(pattern = '.xlsx', coordi)) > 0) {
      coordination <- read.xlsx(coordi, rowNames = T)
    } else {
      coordination <- data.table::fread(coordi) %>% as.data.frame()
      rownames(coordination)  <-
        coordination[, grep(pattern = 'barcode|Barcode',
                            colnames(coordination),
                            value = T)]
    }
    
    colnames(coordination)[grep(pattern = 'x|X|row|Row', colnames(coordination))] <-
      'imagerow'
    colnames(coordination)[grep(pattern = 'y|Y|col|Col', colnames(coordination))] <-
      'imagecol'
    
    if (length(grep(pattern = 'tissue', colnames(coordination))) > 0) {
      colnames(coordination)[grep(pattern = 'tissue', colnames(coordination))] <-
        'tissue'
      coordination <-
        coordination[which(coordination$tissue == 1), ]
    }
    ## image object
    if (length(image) > 0) {
      img_file <- image
      
      if (length(grep(pattern = '.png$', img_file)) > 0) {
        img <- png::readPNG(img_file)
      } else if (length(grep(pattern = '.jpeg$', img_file)) > 0) {
        img <- jpeg::readJPEG(img_file)
      } else if (length(grep(pattern = 'tiff$', img_file)) > 0) {
        img <- tiff::readTIFF(img_file)
      }
      
      # coordination <- coordination[,c('imagerow','imagecol')]
      if (x_reverse == TRUE) {
        coordination$imagecol <- ncol(img) - coordination$imagecol
      }
      
      # scale_factors_file <- Sys.glob(paths = file.path(x, "scalefactors_json.json"))
      
      if (length(scale_factors_file) > 0) {
        scale_factors <- jsonlite::fromJSON(txt = scale_factors_file)
        
        unnormalized.radius <-
          scale_factors$fiducial_diameter_fullres *
          scale_factors$tissue_lowres_scalef
        spot.radius <- unnormalized.radius / max(dim(x = image))
      } else {
        scale_factors <-  list(
          tissue_hires_scalef = 1,
          tissue_lowres_scalef = 1,
          fiducial_diameter_fullres = 1,
          spot_diameter_fullres = 1
        )
        unnormalized.radius <- 1
        spot.radius <- 1
      }
      img_obj <- new(
        Class = "VisiumV1",
        image = img,
        scale.factors = Seurat::scalefactors(
          spot = scale_factors$spot_diameter_fullres,
          fiducial = scale_factors$fiducial_diameter_fullres,
          hires = scale_factors$tissue_hires_scalef,
          scale_factors$tissue_lowres_scalef
        ),
        coordinates = coordination,
        spot.radius = spot.radius
      )
      
      img_obj <- img_obj[Seurat::Cells(x = object)]
      Seurat::DefaultAssay(object = img_obj) <- 'Spatial'
      object[['slice']] <- img_obj
    } else{
      object  <-
        list(seurat_obj = object, coordination = coordination)
      # names(object) <- 'spatial'
    }
    return(object)
  }




#' Seurat process
#'
#' Traditional seurat process with spatial data (using SCT)
#'
#'
#' @param obj seurat object (None-processing)
#' @param npcs Dimensions used of pca
#' @param resolution precision of clustering
#'
#'
#' @import data.table hdf5r reticulate dplyr Seurat SummarizedExperiment scran dplyr scatterpie
#' 
#' 
#' @return a object for shinySRT
#'
#'
#' @export
seurat_sp_process <- function(obj,
                              species = 'hg', npcs = 20, resolution = 1) {
  if (species == 'hg') {
    mt <- toupper('^mt-')
  } else{
    mt <- '^mt-'
  }
  
  obj@meta.data$percent.mt <-
    Seurat::PercentageFeatureSet(obj, pattern = mt)
  ccgenes <- Seurat::cc.genes
  
  if (species == 'mm') {
    ccgenes$s.genes <- HgMM_ex(ccgenes$s.genes)
    ccgenes$g2m.genes <- HgMM_ex(ccgenes$g2m.genes)
  }
  
  obj <- Seurat::SCTransform(obj, assay = Seurat::Assays(obj)[1])
  
  # if(length(VariableFeatures(obj)) < 3000){
  #   if(3000 > length(rownames(obj))){
  #     obj <- FindVariableFeatures(obj,selection.method = 'vst',nfeatures = rownames(obj))
  #   }else{
  #     obj <- FindVariableFeatures(obj,selection.method = 'vst',nfeatures = 3000)
  #   }
  # }
  
  # obj <-
  #   Seurat::CellCycleScoring(
  #     object = obj,
  #     s.features = ccgenes$s.genes,
  #     g2m.features = ccgenes$g2m.genes
  #   )
  
  obj <-
    RunPCA(obj,
           features = VariableFeatures(obj),
           npcs = 50,
           verbose = F)
  
  obj <- FindNeighbors(obj, dims = 1:npcs, verbose = F)
  obj <- FindClusters(obj, verbose = F,resolution = resolution)
  obj <- RunUMAP(obj, dims = 1:npcs, verbose = F)
  # obj <- RunTSNE(obj, dims = 1:npcs, verbose = F)
  return(obj)
}

## process
HgMM_ex <- function(gene) {
  gene <- tolower(gene)
  substring(gene, first = 1, last = 1) <-
    toupper(substring(gene, first = 1, last = 1))
  gene
}




#' Seurat process for multi-object(object in list)
#'
#' Traditional seurat process with spatial data (using SCT)
#'
#'
#' @param obj_list seurat object list (None-processing)
#'
#'
#' @import data.table hdf5r reticulate dplyr Seurat SummarizedExperiment scran dplyr scatterpie
#' 
#' 
#' @return a object for shinySRT
#'
#'
#' @export
obj_list_process <- function(obj_list,
                             species = 'hg',
                             npcs = 20,resolution = 1) {
  multi_types <- lapply(1:length(obj_list), function(x) {
    class(obj_list[[x]])
  })
  
  if (length(which(unique(multi_types) == 'list')) > 0) {
    extr_coordi = T
  } else{
    extr_coordi = F
  }
  
  if (extr_coordi) {
    obj_list_t <- lapply(1:length(obj_list), function(x) {
      obj <- obj_list[[x]]
      if (class(obj) == 'list') {
        object <- obj[['seurat_obj']]
        object@meta.data$slice_sample <- paste0('slice_',x)
      } else{
        object <- obj
        object@meta.data$slice_sample <- paste0('slice_',x)
      }
      object <-
        RenameCells(object,
                    new.names = paste(
                      object@meta.data$slice_sample,
                      rownames(object@meta.data),
                      sep = '_'
                    ))
      return(object)
    })
    
    obj_t <- Reduce(merge, obj_list_t)
    obj_t <-
      seurat_sp_process(obj = obj_t, species = species, npcs = npcs,resolution = resolution)
    
    
    coordi_list_t <- lapply(1:length(obj_list), function(x) {
      obj <- obj_list[[x]]
      if (class(obj) == 'list') {
        coordi <- obj[['coordination']]
        coordi <- coordi[, c('imagerow',  'imagecol')]
        coordi$slice_sample <- paste0('slice_',x)
        
      } else{
        coordi <- GetTissueCoordinates(obj)
        coordi <- coordi[, c('imagerow',  'imagecol')]
        coordi$slice_sample <- paste0('slice_',x)
        
      }
      rownames(coordi) <-
        paste(coordi$slice_sample, rownames(coordi), sep = '_')
      return(coordi)
    })
    
    mat <- GetAssayData(obj_t)
    meta <- obj_t@meta.data
    images <- lapply(1:length(obj_list), function(x) {
      NULL
    })
    names(images) <- paste0('slice_', 1:length(obj_list))
    names(coordi_list_t) <- paste0('slice_', 1:length(obj_list))
    SRT_object <-
      list(
        image = images,
        data = mat,
        coordination = coordi_list_t,
        metadata = meta,
        reduction = NULL,
        normalize = T
      )
    
  } else{
    obj_list_t <- lapply(1:length(obj_list), function(x) {
      obj <- obj_list[[x]]
      object <- obj
      object@meta.data$slice_sample <- paste0('slice_',x)
      
      object <-
        Seurat::RenameCells(object,
                            new.names = paste(
                              object@meta.data$slice_sample,
                              rownames(object@meta.data),
                              sep = '_'
                            ))
      return(object)
    })
    
    obj_t <- Reduce(merge, obj_list_t)
    obj_t <-
      seurat_sp_process(obj = obj_t, species = species, npcs = npcs,resolution = resolution)
    
    SRT_object <- obj_t
  }
  return(SRT_object)
}


#' Load the SRT data from dir
#'
#' 
#'
#'
#' @param x directory of SRT data, matrix, position information, image, et.al.
#'
#'
#' @import data.table hdf5r reticulate dplyr Seurat SummarizedExperiment scran dplyr scatterpie
#' 
#' 
#' @return a object for process
#'
#'
#' @export
## load SRT data(dir)
spatial_load_dir <- function(x,
                             species = 'hg',
                             x_reverse = F,
                             min_cells = 3,
                             ...) {
  if (length(grep(pattern = 'filtered_feature_bc_matrix.h5$', dir(x))) == 1) {
    ## Raw data processing requires seurat
    h5 <-
      grep(pattern = 'filtered_feature_bc_matrix.h5$', dir(x), value = T)
    data <- Seurat::Read10X_h5(filename = file.path(x, h5))
    object <-
      Seurat::CreateSeuratObject(counts = data,
                                 assay = 'Spatial',
                                 min.cells = min_cells)
    img_obj <-
      Seurat::Read10X_Image(image.dir = file.path(x, "spatial"),
                            filter.matrix = T)
    img_obj <- img_obj[Seurat::Cells(x = object)]
    Seurat::DefaultAssay(object = img_obj) <- 'Spatial'
    object[['slice']] <- img_obj
    
  } else if (length(which(dir(x) == "filtered_feature_bc_matrix")) == 1) {
    data <- Seurat::Read10X(file.path(x, "filtered_feature_bc_matrix"))
    object <-
      Seurat::CreateSeuratObject(data, assay = 'Spatial', min.cells = min_cells)
    
    img_obj <-
      Seurat::Read10X_Image(image.dir = file.path(x, "spatial"),
                            filter.matrix = T)
    img_obj <- img_obj[Seurat::Cells(x = object)]
    
    Seurat::DefaultAssay(object = img_obj) <- 'Spatial'
    object[['slice']] <- img_obj
    
  } else if (length(grep(pattern = 'matrix', dir(x))) > 0) {
    mtx_file <- Sys.glob(paths = file.path(x, "*matrix*"))
    meta <- Sys.glob(paths = file.path(x, "*meta*"))
    scale_factors_file <-
      Sys.glob(paths = file.path(x, "scalefactors_json.json"))
    tissue.positions.path <-
      Sys.glob(paths = file.path(x, "*position*"))
    img_file <- Sys.glob(paths = file.path(x, "*image*"))
    object <- load_spatial(
      matx = mtx_file,
      meta = meta,
      coordi = tissue.positions.path,
      image = img_file,
      x_reverse = x_reverse,
      species = species,
      scale_factors_file = scale_factors_file
    )
  }
  return(object)
}



#' multi-dir generate obj for shiny app
#'
#' 
#'
#'
#' @param multi_dir directory of multi SRT data directory
#'
#'
#' @import data.table hdf5r reticulate dplyr Seurat SummarizedExperiment scran dplyr scatterpie
#' 
#' 
#' @return a object for process
#'
#'
#' @export
multi_dir_spatial <-
  function(multi_dir,
           species = 'hg',
           x_reverse = F,min_cells = 3,resolution =1,
           npcs = 20) {
    if (length(dir(multi_dir)) == 0) {
      stop('Do not find SRT dir!')
    }
    
    obj_list <- lapply(dir(multi_dir, full.names = T), function(x) {
      spatial_load_dir(x = x,
                       species = species,
                       x_reverse = x_reverse,min_cells = min_cells)
    })
    SRT_object <-
      obj_list_process(obj_list, species = species, npcs = 20,resolution = resolution)
    return(SRT_object)
  }



#' Load the SRT data from dir and process
#'
#' 
#'
#'
#' @param dir directory of SRT data, matrix, position information, image, et.al.
#'
#'
#' @import data.table hdf5r reticulate dplyr Seurat SummarizedExperiment scran dplyr scatterpie
#' 
#' 
#' @return a object for shinySRT
#'
#'
#' @export
## load SRT data(dir)
single_op_dir <- function(dir,
                          species = 'hg',
                          x_reverse = F,
                          min_cells = 3,resolution = 1,
                          npcs = 20) {
  object <- spatial_load_dir(
    x = dir,
    species = species,
    x_reverse = x_reverse,
    min_cells = min_cells
  )
  SRT_object <- obj_list_process(obj_list = list(object),resolution = resolution,
                                 species = species,
                                 npcs = npcs)
  return(SRT_object)
}




#' Load the SRT data  and process
#'
#' 
#'
#'
#' @param matx gene-barcode matrix file path, formation:.txt, .csv, .tsv, .xlsx
#' @param meta meta data file path, formation:.txt, .csv, .tsv, .xlsx
#' @param coordi file path, Spatial location information, formation:.txt, .csv, .tsv, .xlsx
#' @param image image file path, formation: .png, .jepg, .tiff. 
#'
#'
#' @import data.table hdf5r reticulate dplyr Seurat SummarizedExperiment scran dplyr scatterpie
#' 
#' 
#' @return a object for shinySRT
#'
#'
#' @export
## load SRT data(dir)
single_op_file <- function(matx = NULL,
                           meta = NULL,
                           coordi = NULL,
                           image = NULL,
                           species = 'hg',
                           x_reverse = F,
                           min_cells = 3,resolution = 1,
                           scale_factors_file = NULL,
                           npcs = 20) {
  object <-
    load_spatial(
      matx = matx,
      meta = meta,
      coordi = coordi,
      image = image,
      species = species,
      x_reverse = x_reverse,
      min_cells = min_cells,
      scale_factors_file = scale_factors_file
    )
  SRT_object <- obj_list_process(obj_list = list(object),resolution = resolution,
                                 species = species,
                                 npcs = npcs)
  return(SRT_object)
}


