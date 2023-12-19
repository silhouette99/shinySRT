#' Prepare data files required for shiny app
#'
#' Generate data files required for shiny app.Six files will be generated,
#' (1) shinySRT config \code{meta_group.Rds}
#' (2) the gene mapping object config \code{genesets.Rds}
#' (3) the gene expression \code{data.h5}
#' (4) the spatial metadata \code{meta.Rds}
#' (5) the defaults for the shiny app \code{df_select.Rds}
#' (6) the spatial image \code{image.Rds}
#' Note that both \code{preparedata_shinyspatial} and \code{prepare_code} functions are ran when
#' running the wrapper function \code{CreateshinySRT}.
#'
#' @param dat imported data for shinySRT
#' @param meta.to.include display the meta.data colnames
#' @param maxlevel maximum number of levels allowed for categorical metadata.
#'  maximum number of levels allowed for categorical metadata.
#' @param shiny.dir specify directory to create the shiny app in. Default is
#'   to create a new directory named "shinyspatial_app"
#' @param chunkSize number of genes written to h5file at any one time. Lower
#'   this number to reduce memory consumption. Should not be less than 10
#' @param gex.assay assay in spatially resolved transcriptomic data object to use for plotting
#'   gene expression, which must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: "RNA" or "integrated" assay,
#'       default is "RNA"
#'     \item{SCE objects}: "logcounts" or "normcounts" or "counts",
#'       default is "logcounts"
#'     \item{h5ad files}: "X" or any assay in "layers",
#'       default is "X"
#'   }
#' @param gene.mapping specifies whether to convert human / mouse Ensembl gene
#'   IDs (e.g. ENSG000xxx / ENSMUSG000xxx) into "user-friendly" gene symbols.
#'   Set this to \code{TRUE} if you are using Ensembl gene IDs. Default is
#'   \code{FALSE} which is not to perform any conversion. Alternatively, users
#'   can supply a named vector where \code{names(gene.mapping)} correspond
#'   to the actual gene identifiers in the gene expression matrix and
#'   \code{gene.mapping} correspond to new identifiers to map to
#' @param shiny.prefix specify file prefix
#' @param shiny.dir specify directory to create the shiny app in
#' @param default.gene1 specify primary default gene to show
#' @param default.gene2 specify secondary default gene to show
#' @param default.multigene character vector specifying default genes to
#'   show in bubbleplot / heatmap
#'
#'
#'
#'
#'
#'
#'
#'
#' @import data.table hdf5r reticulate dplyr SpatialExperiment SingleCellExperiment Seurat
#'
#'
#' @examples
#' preparedata_shinyspatial(dat,
#' meta.to.include = NA, #spot information
#' maxlevel = 50, # spot information selected
#' shiny.dir = 'shinyspatial_app', # dir
#' chunkSize = 500,
#' gex.assay = NA, # assay spesific
#' gex.slot = c("data", "scale.data","counts"),
#' gene.mapping = TRUE,
#' default.gene1 = NA,
#' default.gene2 = NA,default.multigene = NA)
#'
#'
#' @export
preparedata_shinyspatial <- function(dat,
                                     meta.to.include = NA,
                                     maxlevel = 50,
                                     shiny.dir = 'shinyspatial_app',
                                     chunkSize = 500,
                                     gex.assay = NA,
                                     gex.slot = c("data", "scale.data","counts"),
                                     gene.mapping = TRUE,
                                     default.gene1 = NA,
                                     default.gene2 = NA,
                                     default.multigene = NA
){
  drExist = TRUE # image array
  dmExist = TRUE # dim matrix
  
  if (class(dat)[1] == "Seurat") {
    ## get the image array and coordination matrix
    
    
    
    if (length(names(dat@images)) == 0) {
      drExist = FALSE
    }
    if (!drExist) {
      stop(paste0("shinSRT did not detect any coordination data"))
    }
    
    
    if (class(dat@images[[1]])[[1]] == 'FOV') { ## vizgene
      image <- dat@images
      #
      # for (i in names(image)) {
      #   png(paste0(i, '.png'))
      #   vizgene_bgplot(image[[i]])
      #   dev.off()
      # }
      #
      # ima <- lapply(paste0(names(image), '.png'), function(x) {
      #   imgs <- png::readPNG(x)
      # })
      
      ima <- list(rep(NULL,length(image)))
      
      
      names(ima) <- names(image)
      coordi <- lapply(names(ima), function(x) {
        coordi <- Seurat::GetTissueCoordinates(image[[x]])
        coordi$slice_sample <- x
        rownames(coordi) <- coordi$cell
        coordi <- coordi[,c(2,1,3,4)]
        
        colnames(coordi) <- c('imagerow','imagecol','cell', 'slice_sample')
        fs <- (max(coordi$imagerow) - min(coordi$imagerow))/480
        coordi$imagerow <- coordi$imagerow/fs
        coordi$imagecol <- coordi$imagecol/fs
        return(coordi)
      })
      
      names(coordi) <- NULL
      coordi <- do.call(rbind, coordi)
    } else if (class(dat@images[[1]])[[1]] == 'VisiumV1') { ## visium
      image <- dat@images
      ### seurat spot coordination are needed to flip
      coordi <- lapply(names(image), function(x) {
        coordi <- Seurat::GetTissueCoordinates(image[[x]])
        coordi$slice_sample <- x
        xrange <- c(0, dim(image[[x]])[2])
        yrange <- c(0, dim(image[[x]])[1])
        
        ## Image Orientation
        coordi$imagerow <-
          (yrange[2] - yrange[1]) / 2 - (coordi$imagerow - (yrange[2] - yrange[1]) / 2)
        return(coordi)
      })
      names(coordi) <- NULL
      coordi <- do.call(rbind, coordi)
      
      ### get image data
      ima <- lapply(names(image), function(x) {
        image[[x]]@image
      })
      names(ima) <- names(image)
    }
    
    
  }else if (class(dat)[1] == "SingleCellExperiment") {
    if (length(SingleCellExperiment::reducedDimNames(dat)) ==
        0) {
      dmExist = FALSE
    } else{
      dmExist = TRUE
    }
    
    if (!dmExist) {
      warning(
        paste0(
          "shinySRT did not detect any dimension reduction data \n",
          "       e.g. umap / tsne. Has any analysis been performed?"
        )
      )
    }
    
    if (length(dat@metadata$image$sample) == 0) {
      drExist = FALSE
    }
    if (!drExist) {
      stop(paste0("shinySRT did not detect any coordination data"))
    }
    
    ima <- lapply(dat@metadata$image$grob, function(x) {
      x[[1]]
    })
    names(ima) <- dat@metadata$image$sample
    
    slice_samples <- lapply(colnames(dat@colData), function(x) {
      if (length(intersect(unique(dat@colData[, x]), names(ima))) > 0) {
        samples <- x
      } else{
        samples <- NULL
      }
      samples
    }) %>% unlist()
    
    ima <- ima[which(names(ima) %in% dat@colData[, slice_samples[1]])]
    
    coordi <- lapply(1:length(names(ima)), function(x) {
      coordi <- dat@colData@listData %>% as.data.frame()
      
      rownames(coordi) <-
        paste(dat@colData[[slice_samples[1]]], dat@colData@rownames, sep = '_')
      
      coordi <- coordi[, c('imagerow', 'imagecol', slice_samples[1])]
      
      coordi_ <- subset(coordi, coordi[, slice_samples[1]] == names(ima)[x])
      
      coordi_$slice_sample <- coordi_[, slice_samples[1]]
      
      
      xrange <- c(0, dim(ima[[names(ima)[x]]])[2])
      yrange <- c(0, dim(ima[[names(ima)[x]]])[1])
      
      coordi_$imagerow <-
        (yrange[2] - yrange[1]) / 2 - (coordi_$imagerow - (yrange[2] - yrange[1]) / 2)
      return(coordi_)
    })
    names(coordi) <- NULL
    coordi <- do.call(rbind, coordi)
  } else if (class(dat)[1] == "SpatialExperiment") {
    if (length(SingleCellExperiment::reducedDimNames(dat)) ==
        0) {
      dmExist = FALSE
    } else{
      dmExist = TRUE
    }
    
    if (!dmExist) {
      warning(
        paste0(
          "shinySRT did not detect any dimension reduction data \n",
          "       e.g. umap / tsne. Has any analysis been performed?"
        )
      )
    }
    
    if (length(dat@int_metadata$imgData@listData$sample_id) == 0) {
      drExist = FALSE
    }
    if (!drExist) {
      stop(paste0("shinySRT did not detect any coordination data"))
    }
    
    ima <- dat@int_metadata$imgData$data
    
    ima <- lapply(1:length(ima), function(x){
      # if(class(ima[[x]])[1] == 'StoredSpatialImage'){
      #
      #   image <- SpatialImage(ima[[x]])
      #   image <- as(image, "LoadedSpatialImage")
      #   image <- imgRaster(dat,
      #                      sample_id = dat@int_metadata$imgData@listData$sample_id[[x]],
      #                      image_id = dat@int_metadata$imgData@listData$image_id[[x]])
      # }else{
      #   image <- ima[[x]]
      # }
      
      
      
      image <- SpatialImage(ima[[x]])
      image <- as(image, "LoadedSpatialImage")
      image <- SpatialExperiment::rotateImg(image, +90)
      image <- SpatialExperiment::mirrorImg(image, "v")
      as.raster(image)
    })
    
    names(ima) <- dat@int_metadata$imgData$sample_id
    
    slice_samples <- lapply(colnames(dat@colData), function(x) {
      if (length(intersect(unique(dat@colData[, x]), names(ima))) > 0) {
        samples <- x
      } else{
        samples <- NULL
      }
      samples
    }) %>% unlist()
    
    ima <- ima[which(names(ima) %in% dat@colData[, slice_samples[1]])]
    
    ### seurat spot coordination are needed to flip
    coordi <- lapply(1:length(names(ima)), function(x) {
      coordi <- dat@int_colData$spatialCoords %>% as.data.frame()
      
      coordi[,slice_samples[1]] <- dat@colData[[slice_samples[1]]]
      
      rownames(coordi) <-
        paste(dat@colData[[slice_samples[1]]], dat@colData@rownames, sep = '_')
      
      coordi_ <- coordi[,c(1,2,3)]
      
      colnames(coordi) <- c('imagerow', 'imagecol', 'slice_sample')
      
      coordi_ <- subset(coordi,coordi$slice_sample == names(ima)[x])
      
      coordi_$imagerow <- coordi_$imagerow*dat@int_metadata$imgData@listData$scaleFactor[x]
      coordi_$imagecol <- coordi_$imagecol*dat@int_metadata$imgData@listData$scaleFactor[x]
      
      xrange <- c(0, dim(ima[[x]])[2])
      yrange <- c(0, dim(ima[[x]])[1])
      
      coordi_$imagerow <-
        (yrange[2] - yrange[1]) / 2 - (coordi_$imagerow - (yrange[2] - yrange[1]) / 2)
      return(coordi_)
    })
    names(coordi) <- NULL
    coordi <- do.call(rbind, coordi)
  } else if (class(dat)[1] == "list") {
    if (length(dat[['redutcion']]) ==
        0) {
      dmExist = FALSE
    } else{
      dmExist = TRUE
    }
    
    if (!dmExist) {
      warning(
        paste0(
          "shinySRT did not detect any dimension reduction data \n",
          "       e.g. umap / tsne. Has any analysis been performed?"
        )
      )
    }
    
    if (length(dat[['metadata']]$slice_sample) == 0) {
      drExist = FALSE
    }
    if (!drExist) {
      stop(paste0("shinySRT did not detect any coordination data"))
    }
    
    ima <- dat[['image']]
    
    slice_samples <- lapply(colnames(dat[['metadata']]), function(x) {
      if (length(intersect(unique(dat[['metadata']][, x]), names(ima))) > 0) {
        samples <- x
      } else{
        samples <- NULL
      }
      samples
    }) %>% unlist()
    
    ima <- ima[which(names(ima) %in% dat[['metadata']][, slice_samples[1]])]
    coordi <- dat[['coordination']]
    
    coordi <- coordi[,c(2,1)]
    
    colnames(coordi) <- c('imagerow','imagecol')
  }else if(tolower(tools::file_ext(dat)) == "h5ad"){
    obj <- hdf5r::H5File$new(dat, mode = "r")
    
    if (length(grep(pattern = '^X_',obj[['obsm']]$names)) ==
        0) {
      dmExist = FALSE
    } else{
      dmExist = TRUE
    }
    
    if (!dmExist) {
      warning(
        paste0(
          "shinySRT did not detect any dimension reduction data \n",
          "       e.g. umap / tsne. Has any analysis been performed?"
        )
      )
    }
    
    if (length(obj[['uns']][['spatial']]$names) == 0) {
      drExist = FALSE
    }
    if (!drExist) {
      stop(paste0("shinySRT did not detect any coordination data"))
    }
    
    meta <- lapply(obj[['obs']]$names, function(i){
      if(class(obj[['obs']][[i]])[1] == "H5Group"){
        meta_inf <- obj[['obs']][[i]][['codes']]$read()
        meta_inf <- factor(meta_inf, levels = unique(meta_inf))
        levels(meta_inf) <- obj[['obs']][[i]][['categories']]$read()
      }else if(class(obj[['obs']][[i]])[1] == "H5D"){
        meta_inf <- obj[['obs']][[i]]$read()
      }
      meta_inf <- as.data.frame(meta_inf)
      if(i == '_index'){
        colnames(meta_inf) <- 'spots'
      }else{
        colnames(meta_inf) <- i
      }
      
      return(meta_inf)
    })
    
    meta <- do.call(cbind, meta)
    
    ima <- lapply(obj[['uns']][['spatial']]$names, function(m){
      image <- obj[['uns']][['spatial']][[m]][['images']][['lowres']]$read()
      array(c(image[1,,],image[2,,],image[3,,]),c(rev(dim(image[1,,])),3))
    })
    
    names(ima) <- obj[['uns']][['spatial']]$names
    
    ### seurat spot coordination are needed to flip
    coordi <- lapply(names(ima), function(x) {
      coordi<- obj[['obsm']][['spatial']]$read()%>%t()%>%as.data.frame()
      if(length(names(ima)) > 1){
        coordi$sample_id <- meta$library_id
      }else{
        coordi$sample_id <- x
      }
      
      rownames(coordi) <- meta$spots
      coordi_ <- subset(coordi,coordi$sample_id == x)
      
      coordi_$V1 <- coordi_$V1*obj[['uns']][['spatial']][[x]][['scalefactors']][['tissue_lowres_scalef']]$read()
      coordi_$V2 <- coordi_$V2*obj[['uns']][['spatial']][[x]][['scalefactors']][['tissue_lowres_scalef']]$read()
      
      colnames(coordi_) <- c('imagerow', 'imagecol','sample_id')
      
      coordi_$slice_sample <- x
      
      xrange <- c(0, dim(ima[[x]])[2])
      yrange <- c(0, dim(ima[[x]])[1])
      
      coordi_$imagerow <-
        (yrange[2] - yrange[1]) / 2 - (coordi_$imagerow - (yrange[2] - yrange[1]) / 2)
      return(coordi_)
    })
    names(coordi) <- NULL
    coordi <- do.call(rbind, coordi)
  }
  
  
  
  
  ## meta
  if (class(dat)[1] == "Seurat") {
    meta = dat@meta.data
    
    ### reduction matrix
    if (length(names(dat@reductions)) == 0) {
      dmExist <- FALSE
    } else{
      dmExist <- TRUE
    }
    if (!dmExist) {
      warning(
        paste0(
          "shinySRT did not detect any dimension reduction data \n",
          "       e.g. umap / tsne. Has any analysis been performed?"
        )
      )
    }
    if(dmExist){
      embedding <- lapply(dat@reductions, function(x) {
        embedding <- x@cell.embeddings %>% as.data.frame()
        if (ncol(embedding) > 5) {
          embedding <- embedding[, 1:5]
        } else{
          embedding <- embedding
        }
        return(embedding)
      })
      embedding <- do.call(cbind, embedding)
      meta <- cbind(meta, embedding)
    }
    
    cells <- intersect(rownames(coordi),rownames(meta))
    coordi <- coordi[which(rownames(coordi) %in% cells),]
    meta <- meta[which(rownames(meta) %in% cells),]
    dat <- dat[,cells]
    
    
    coordi <- cbind(coordi[rownames(x = meta),],meta)
    
    if (is.na(meta.to.include[1])) {
      meta.to.include = colnames(meta)
    }
    if (length(meta.to.include) < 2) {
      stop("At least 2 metadata is required!")
    }
    if (length(which(is.na(coordi$orig.ident))) > 0) {
      stop('some cells can not be identified !')
    }
  }else if (class(dat)[1] == "SingleCellExperiment") {
    
    meta <- dat@colData@listData %>% as.data.frame()
    if (dmExist) {
      embedding <-
        lapply(SingleCellExperiment::reducedDims(dat), function(x) {
          embedding <- x %>% as.data.frame()
          if (ncol(embedding) > 5) {
            embedding <- embedding[, 1:5]
          } else{
            embedding <- embedding
          }
          rownames(embedding) <-
            paste(dat@colData[[slice_samples]], dat@colData@rownames, sep = '_')
          return(embedding)
        })
      embedding <- do.call(cbind, embedding)
      meta <- cbind(meta, embedding)
    }
    rownames(meta) <-
      paste(dat@colData[[slice_samples]], dat@colData@rownames, sep = '_')
    
    
    cells <- intersect(rownames(coordi),rownames(meta))
    coordi <- coordi[which(rownames(coordi) %in% cells),]
    meta <- meta[which(rownames(meta) %in% cells),]
    
    coordi <- cbind(coordi[rownames(x = meta),],meta)
    
    coordi <- cbind(coordi[rownames(x = meta), ], meta)
    if (is.na(meta.to.include[1])) {
      meta.to.include = colnames(meta)
    }
    if (length(meta.to.include) < 2) {
      stop("At least 2 metadata is required!")
    }
    if (length(which(is.na(coordi$orig.ident))) > 0) {
      stop('some cells can not be identified !')
    }
  } else if (class(dat)[1] == "SpatialExperiment") {
    meta <- dat@colData@listData %>% as.data.frame()
    if(dmExist){
      embedding <- lapply(SingleCellExperiment::reducedDims(dat), function(x) {
        embedding <- x %>% as.data.frame()
        if (ncol(embedding) > 5) {
          embedding <- embedding[, 1:5]
        } else{
          embedding <- embedding
        }
        
        
        return(embedding)
      })
      embedding <- do.call(cbind, embedding)
      
      meta <- cbind(meta, embedding)
    }
    
    rownames(meta) <-
      paste(dat@colData$sample_id, dat@colData@rownames, sep = '_')
    
    cells <- intersect(rownames(coordi),rownames(meta))
    coordi <- coordi[which(rownames(coordi) %in% cells),]
    meta <- meta[which(rownames(meta) %in% cells),]
    
    coordi <- cbind(coordi[rownames(x = meta),],meta)
    
    if (is.na(meta.to.include[1])) {
      meta.to.include = colnames(meta)
    }
    if (length(meta.to.include) < 2) {
      stop("At least 2 metadata is required!")
    }
    if (length(which(is.na(coordi$orig.ident))) > 0) {
      stop('some cells can not be identified !')
    }
  }else if (class(dat)[1] == "list") {
    
    meta <- dat[['metadata']] %>% as.data.frame()
    if (dmExist) {
      embedding <-
        lapply(dat[['reduction']], function(x) {
          embedding <- x %>% as.data.frame()
          if (ncol(embedding) > 5) {
            embedding <- embedding[, 1:5]
          } else{
            embedding <- embedding
          }
          # rownames(embedding) <-
          #   paste(dat@colData[[slice_samples]], dat@colData@rownames, sep = '_')
          return(embedding)
        })
      embedding <- do.call(cbind, embedding)
      meta <- cbind(meta, embedding)
    }
    
    # rownames(meta) <-
    #   paste(dat@colData[[slice_samples]], dat@colData@rownames, sep = '_')
    cells <- intersect(rownames(coordi),rownames(meta))
    coordi <- coordi[which(rownames(coordi) %in% cells),]
    meta <- meta[which(rownames(meta) %in% cells),]
    
    
    coordi <- coordi[match(rownames(meta),rownames(coordi)),]
    coordi <- cbind(coordi,meta)
    
    
    if (is.na(meta.to.include[1])) {
      meta.to.include = colnames(meta)
    }
    if (length(meta.to.include) < 2) {
      stop("At least 2 metadata is required!")
    }
    if (length(which(is.na(coordi$orig.ident))) > 0) {
      stop('some cells can not be identified !')
    }
  }else if(tolower(tools::file_ext(dat)) == "h5ad"){
    obj <- hdf5r::H5File$new(dat, mode = "r")
    
    meta <- lapply(obj[['obs']]$names, function(i){
      if(class(obj[['obs']][[i]])[1] == "H5Group"){
        meta_inf <- obj[['obs']][[i]][['codes']]$read()
        meta_inf <- factor(meta_inf, levels = unique(meta_inf))
        levels(meta_inf) <- obj[['obs']][[i]][['categories']]$read()
      }else if(class(obj[['obs']][[i]])[1] == "H5D"){
        meta_inf <- obj[['obs']][[i]]$read()
      }
      meta_inf <- as.data.frame(meta_inf)
      if(i == '_index'){
        colnames(meta_inf) <- 'spots'
      }else{
        colnames(meta_inf) <- i
      }
      return(meta_inf)
    })
    
    meta <- do.call(cbind, meta)
    if(dmExist){
      embedding <- lapply(grep('^X_',obj[['obsm']]$names,value = T), function(x) {
        embedding <- obj[['obsm']][[x]]$read()%>%t()%>%as.data.frame()
        if (ncol(embedding) > 5) {
          embedding <- embedding[, 1:5]
        } else{
          embedding <- embedding
        }
        colnames(embedding) <- paste(x,1:ncol(embedding),sep = '_')
        
        rownames(embedding) <-
          meta$spots
        
        return(embedding)
      })
      
      embedding <- do.call(cbind, embedding)
      
      meta <- cbind(meta, embedding)
    }
    
    
    rownames(meta) <- meta$spots
    
    cells <- intersect(rownames(coordi),rownames(meta))
    coordi <- coordi[which(rownames(coordi) %in% cells),]
    meta <- meta[which(rownames(meta) %in% cells),]
    
    coordi <- cbind(coordi[rownames(x = meta),],meta)
    
    if (is.na(meta.to.include[1])) {
      meta.to.include = colnames(meta)
    }
    if (length(meta.to.include) < 2) {
      stop("At least 2 metadata is required!")
    }
    if (length(which(is.na(coordi$orig.ident))) > 0) {
      stop('some cells can not be identified !')
    }
    
    obj$close_all()
  }
  
  coordi <- coordi %>% distinct(.keep_all = T)
  
  colnames(coordi)[which(duplicated(colnames(coordi)))] <-
    paste(colnames(coordi)[which(duplicated(colnames(coordi)))], sample(length(colnames(coordi)[which(duplicated(colnames(coordi)))])))
  
  
  ### config
  meta_group <- lapply(colnames(coordi), function(x) {
    if (length(grep(pattern = 'factor',class(coordi[[x]]))) > 0 |
        length(grep(pattern = 'character',class(coordi[[x]]))) > 0) {
      if (length(grep(pattern = 'character',class(coordi[[x]]))) > 0 & length(unique(coordi[, x])) < maxlevel) {
        coordi[, x] <- factor(coordi[, x], levels = unique(coordi[, x]))
      }
      Unit <- paste(levels(coordi[, x]), collapse = '|')
      color <-
        paste(col_box()[1:length(levels(coordi[, x]))], collapse = '|')
      metas <- data.frame(matrix(nrow = 1, ncol = 5,))
      colnames(metas) <-
        c('group', 'unit', 'color', 'info', 'default')
      if(length(unique(coordi[, x])) < maxlevel){
        metas[1,] <- c(x, Unit, color, TRUE, 0)
      }else{
        metas[1,] <- c(x, Unit, color, FALSE, 0)
      }
      
      metas
    } else if (length(grep(pattern = '\\.', unique(coordi[, x]))) == 0 &
               length(unique(coordi[, x])) < maxlevel) {
      coordi[, x] <- factor(coordi[, x], levels = unique(coordi[, x]))
      Unit <- paste(levels(coordi[, x]), collapse = '|')
      color <-
        paste(col_box()[1:length(levels(coordi[, x]))], collapse = '|')
      metas <- data.frame(matrix(nrow = 1, ncol = 5,))
      colnames(metas) <-
        c('group', 'unit', 'color', 'info', 'default')
      metas[1,] <- c(x, Unit, color, TRUE, 0)
      metas
    } else{
      metas <- data.frame(matrix(nrow = 1, ncol = 5,))
      colnames(metas) <-
        c('group', 'unit', 'color', 'info', 'default')
      metas[1,] <- c(x, NA, NA, FALSE, 0)
      metas
    }
  })
  
  meta_group <- do.call(rbind, meta_group)
  meta_group <- data.table::data.table(meta_group)
  
  def1 = grep("ident|library|Ident|Librar", meta_group$group[meta_group$info == TRUE], ignore.case = TRUE)[1]
  def2 = grep("clust|Clust", meta_group$group[meta_group$info == TRUE], ignore.case = TRUE)[1]
  def2 = setdiff(def2, def1)[1]
  if (is.na(def1)) {
    def1 = setdiff(c(1, 2), def2)[1]
  }
  if (is.na(def2)) {
    def2 = setdiff(c(1, 2), def1)[1]
  }
  meta_group[meta_group$info == TRUE][def1]$default = 1
  meta_group[meta_group$info == TRUE][def2]$default = 2
  
  for (i in as.character(meta_group[!is.na(unit)]$group)) {
    coordi[[i]] = factor(coordi[[i]], levels = strsplit(meta_group[group == i]$unit, "\\|")[[1]])
    levels(coordi[[i]]) = strsplit(meta_group[group == i]$unit,
                                   "\\|")[[1]]
    meta_group[group == i]$unit = meta_group[group == i]$unit
  }
  meta_group$group = as.character(meta_group$group)
  coordi = data.table(sampleID = rownames(coordi),
                      coordi)
  
  
  
  
  ## default information
  if (class(dat)[1] == "Seurat") {
    if (is.na(gex.assay[1])) {
      gex.assay = Seurat::DefaultAssay(dat)
    }
    
    ### data dimation
    
    gex.matdim = dim(slot(dat@assays[[gex.assay[1]]], gex.slot[1]))
    gex.rownm = rownames(slot(dat@assays[[gex.assay[1]]],
                              gex.slot[1]))
    gex.colnm = colnames(slot(dat@assays[[gex.assay[1]]],
                              gex.slot[1]))
    defGenes = Seurat::VariableFeatures(dat)[1:10]
    if (is.na(defGenes[1])) {
      warning(
        paste0(
          "Variable genes for seurat object not found! Have you ",
          "ran `FindVariableFeatures` or `SCTransform`?"
        )
      )
      defGenes = gex.rownm[1:10]
    }
  } else if (class(dat)[1] == "SpatialExperiment") {
    if (is.null(colnames(dat)[1])) {
      colnames(dat) = paste0("cell_", seq(ncol(dat)))
    }
    if (is.na(gex.assay[1])) {
      if(length(grep(pattern = 'logcounts', assayNames(dat))) > 0){
        gex.assay = "logcounts"
      }else{
        gex.assay = assayNames(dat)[1]
      }

    }
    gex.matdim = dim(SummarizedExperiment::assay(dat, gex.assay[1]))
    gex.rownm = rownames(SummarizedExperiment::assay(dat,
                                                     gex.assay[1]))
    gex.colnm = paste(dat@colData[[slice_samples]], colnames(SummarizedExperiment::assay(dat,
                                                                                         gex.assay[1])),sep = '_')
    defGenes = gex.rownm[1:10]
  }  else if (class(dat)[1] == "SingleCellExperiment") {
    if (is.null(colnames(dat)[1])) {
      colnames(dat) = paste0("cell_", seq(ncol(dat)))
    }

    if (is.na(gex.assay[1])) {
      if(length(grep(pattern = 'logcounts', assayNames(dat))) > 0){
        gex.assay = "logcounts"
      }else{
        gex.assay = assayNames(dat)[1]
      }

    }
    gex.matdim = dim(SummarizedExperiment::assay(dat, gex.assay[1]))
    gex.rownm = rownames(SummarizedExperiment::assay(dat,
                                                     gex.assay[1]))
    gex.colnm = paste(dat@colData[[slice_samples]], colnames(SummarizedExperiment::assay(dat,
                                                                                         gex.assay[1])),sep = '_')
    defGenes = gex.rownm[1:10]
    
  } else if (class(dat)[1] == "list") {
    gex.matdim = dim(dat[['data']])
    gex.rownm = rownames(dat[['data']])
    gex.colnm = colnames(dat[['data']])
    defGenes = gex.rownm[1:10]
  } else if(tolower(tools::file_ext(dat)) == "h5ad"){
    obj <- hdf5r::H5File$new(dat, mode = "r")
    if (is.na(gex.assay[1])) {
      gex.assay = "X"
    }
    gex.rownm = obj[['var']][['gene_ids']]$read()
    gex.colnm = meta$spots
    gex.matdim = c(length(gex.rownm),length(gex.colnm))
    defGenes = gex.rownm[1:10]
    
    obj$close_all()
  }
  
  ## gene mapping(ensembl id to symbol id)
  if (gene.mapping[1] == TRUE) {
    if (sum(grepl("^ENSG000", gex.rownm)) >= sum(grepl("^ENSMUSG000",gex.rownm))) {
      tmp1 = fread(system.file("extdata",'hg_map.txt.gz',package = 'shinySRT'))
    } else {
      tmp1 = fread(system.file("extdata",'mm_map.txt.gz',package = 'shinySRT'))
    }
    gene.mapping = tmp1$GeneName
    names(gene.mapping) = tmp1$EnsemblID
  }
  
  if (gene.mapping[1] == FALSE) {
    gene.mapping = gex.rownm
    names(gene.mapping) = gex.rownm
  } else {
    if (!all(gex.rownm %in% names(gene.mapping))) {
      tmp1 = gex.rownm[gex.rownm %in% names(gene.mapping)]
      tmp1 = gene.mapping[tmp1]
      tmp2 = gex.rownm[!gex.rownm %in% names(gene.mapping)]
      names(tmp2) = tmp2
      gene.mapping = c(tmp1, tmp2)
    }
    gene.mapping = gene.mapping[gex.rownm]
  }
  
  
  defGenes = gene.mapping[defGenes]
  default.gene1 = default.gene1[1]
  default.gene2 = default.gene2[1]
  if (is.na(default.gene1)) {
    default.gene1 = defGenes[1]
  }
  if (is.na(default.gene2)) {
    default.gene2 = defGenes[2]
  }
  if (is.na(default.multigene[1])) {
    default.multigene = defGenes
  }
  if (default.gene1 %in% gene.mapping) {
    default.gene1 = default.gene1
  } else if (default.gene1 %in% names(gene.mapping)) {
    default.gene1 = gene.mapping[default.gene1]
  } else {
    warning("default.gene1 doesn't exist in gene expression, using defaults...")
    default.gene1 = defGenes[1]
  }
  if (default.gene2 %in% gene.mapping) {
    default.gene2 = default.gene2
  } else if (default.gene2 %in% names(gene.mapping)) {
    default.gene2 = gene.mapping[default.gene2]
  } else {
    warning("default.gene2 doesn't exist in gene expression, using defaults...")
    default.gene2 = defGenes[2]
  }
  if (all(default.multigene %in% gene.mapping)) {
    default.multigene = default.multigene
  } else if (all(default.multigene %in% names(gene.mapping))) {
    default.multigene = gene.mapping[default.multigene]
  } else {
    warning(paste0(
      "default.multigene doesn't exist in gene expression, ",
      "using defaults..."
    ))
    default.multigene = defGenes
  }
  
  names(default.multigene) <- default.multigene
  names(default.gene1) <- default.gene1
  names(default.gene2) <- default.gene2
  
  
  if (!isTRUE(all.equal(coordi$sampleID, gex.colnm))) {
    coordi$sampleID = factor(coordi$sampleID, levels = gex.colnm)
    coordi = coordi[order(sampleID)]
    coordi$sampleID = as.character(coordi$sampleID)
  }
  
  genesets = seq(gex.matdim[1])
  names(gene.mapping) = NULL
  names(genesets) = gene.mapping
  genesets = genesets[order(names(genesets))]
  genesets = genesets[order(nchar(names(genesets)))]
  
  if(!is.null(ima[[1]])){
    boxsize <- lapply(ima, function(x) {
      dim(x)[2]
    }) %>% unlist()
  }else{
    boxsize <- c(max(coordi$imagerow)-min(coordi$imagerow),max(coordi$imagecol)-min(coordi$imagecol))
  }
  
  
  
  
  ## expression gene matrix
  
  ### Empty file
  if (!dir.exists(shiny.dir)) {
    dir.create(shiny.dir)
  }
  filename = paste0(shiny.dir, '/', "data.h5")
  mat_exp <- H5File$new(filename, mode = "w")
  mat_exp.grp <- mat_exp$create_group("grp")
  mat_exp.grp.data <- mat_exp.grp$create_dataset(
    "data",
    dtype = h5types$H5T_NATIVE_FLOAT,
    space = H5S$new("simple",
                    dims = gex.matdim, maxdims = gex.matdim),
    chunk_dims = c(1,
                   gex.matdim[2])
  )
  chk = chunkSize
  while (chk > (gex.matdim[1] - 8)) {
    chk = floor(chk / 2)
  }
  
  
  if (class(dat)[1] == "Seurat") {
    for (i in 1:floor((gex.matdim[1] - 8) / chk)) {
      mat_exp.grp.data[((i - 1) * chk + 1):(i * chk),] <-
        as.matrix(slot(dat@assays[[gex.assay[1]]],
                       gex.slot[1])[((i - 1) * chk + 1):(i * chk),])
    }
    mat_exp.grp.data[(i * chk + 1):gex.matdim[1],] <-
      as.matrix(slot(dat@assays[[gex.assay[1]]],
                     gex.slot[1])[(i * chk + 1):gex.matdim[1],])
    mat_exp$close_all()
  } else if (class(dat)[1] == "SpatialExperiment" |
             class(dat)[1] == "SingleCellExperiment") {
    for (i in 1:floor((gex.matdim[1] - 8) / chk)) {
      mat_exp.grp.data[((i - 1) * chk + 1):(i * chk),] <-
        as.matrix(SummarizedExperiment::assay(dat,
                                              gex.assay[1])[((i - 1) * chk + 1):(i * chk),])
    }
    mat_exp.grp.data[(i * chk + 1):gex.matdim[1],] <-
      as.matrix(SummarizedExperiment::assay(dat,
                                            gex.assay[1])[(i * chk + 1):gex.matdim[1],])
    mat_exp$close_all()
  }else if (class(dat)[1] == 'list') {
    for (i in 1:floor((gex.matdim[1] - 8) / chk)) {
      mat_exp.grp.data[((i - 1) * chk + 1):(i * chk),] <-
        as.matrix(dat[['data']][((i - 1) * chk + 1):(i * chk),])
    }
    mat_exp.grp.data[(i * chk + 1):gex.matdim[1],] <-
      as.matrix(dat[['data']][(i * chk + 1):gex.matdim[1],])
    mat_exp$close_all()
  }else if(tolower(tools::file_ext(dat)) == "h5ad"){
    obj <- hdf5r::H5File$new(dat, mode = "r")
    da <- rep(NA, length(gex.colnm)*length(gex.rownm))
    da[obj[['X']][['indices']]$read()] <- obj[['X']][['data']]$read()
    da[which(is.na(da))] <- 0
    da <- matrix(da,nrow = length(gex.colnm)) %>% t()
    
    for (i in 1:floor((gex.matdim[1] - 8) / chk)) {
      mat_exp.grp.data[((i - 1) * chk + 1):(i * chk),] <-
        da[((i - 1) * chk + 1):(i * chk),]
    }
    mat_exp.grp.data[(i * chk + 1):gex.matdim[1],] <-
      da[(i * chk + 1):gex.matdim[1],]
    mat_exp$close_all()
    obj$close_all()
  }
  
  
  if(!is.null(ima[[1]])){
    ranges <- lapply(ima, function(x){
      list(
        xrange = c(0, dim(x)[2]),
        yrange = c(0, dim(x)[1])
      )
    })
  }else{
    ranges <- lapply(ima, function(x){
      list(
        xrange = c(0,(max(coordi$imagecol) - min(coordi$imagecol)) + 0.01*(max(coordi$imagecol) - min(coordi$imagecol))),
        yrange = c(0,(max(coordi$imagerow) - min(coordi$imagerow)) + 0.01*(max(coordi$imagerow) - min(coordi$imagerow)))
      )
    })
  }
  names(ranges) <- names(ima)
  
  
  coordi$tooltip <- paste(
    'sampleID: ',
    coordi$sampleID,
    '\n',
    'slice :',coordi$slice_sample,
    '\n',
    meta_group[default == 1 &
                 info == TRUE]$group,
    ': ',
    coordi[[meta_group[default == 1 & info == TRUE]$group]],
    meta_group[default == 2 &
                 info == TRUE]$group,'\n',
    ': ',
    coordi[[meta_group[default == 2 & info == TRUE]$group]]
  )
  
  
  df_select = list()
  df_select$meta1 = meta_group[default == 1 & info == TRUE]$group
  df_select$meta2 = meta_group[default == 2 & info == TRUE]$group
  df_select$gene1 = default.gene1
  df_select$gene2 = default.gene2
  df_select$genes = default.multigene
  df_select$ranges = ranges
  df_select$boxsize = max(boxsize) + 0.1 * max(boxsize)
  df_select$maxarea <-
    ceiling(abs(max(coordi$imagecol) - min(coordi$imagecol)) * 10)
  df_select$slice <- levels(coordi$slice_sample)
  df_select$image_size <- max(boxsize)
  
  saveRDS(coordi, paste0(shiny.dir, '/meta.Rds'))
  saveRDS(meta_group, paste0(shiny.dir, '/meta_group.Rds'))
  saveRDS(genesets, paste0(shiny.dir, '/genesets.Rds'))
  saveRDS(df_select, paste0(shiny.dir, '/df_select.Rds'))
  saveRDS(ima, paste0(shiny.dir, '/image.Rds'))
}

## clusters colors displayed on shiny web
col_box <- function(){
  c(
    '#ff3030','#1e90ff','#ffd700','#ff6eb4','#bf3eff',
    '#c1ffc1','#00ffff','#ff1493','#7fff00','#ff0000',
    '#0000ff','#00ff00','#ffc125','#ffb5c5','#9b30ff',
    '#54ff9f','#87ceff','#ffaeb9','#c0ff3e','#ff4040',
    '#836fff','#00ff7f','#ffff00','#ff3e96','#ff83fa',
    '#9aff9a','#4876ff','#ff3e96','#00f5ff','#8b1a1a',
    '#104e8b','#008b00','#8b7500','#8b3a62','#68228b',
    '#698b69','#008b8b','#8b0a50','#458b00','#eedfcc',
    '#76eec6','#eec591','#dbe9e9','#8ee5ee','#ee6a50',
    '#eee8cd','#eead0e','#bcee68','#ee7600','#8deeee',
    '#00b2ee','#ee6363','#eeeee0','#eee685','#b2dfee',
    '#eedc82','#ee00ee','#ee30a7','#006400'
  )
}
## gradual change of color
heat_col <- function(){
  c("#0000FF", "#0013FF", "#0028FF", "#003CFF", "#0050FF", "#0065FF", "#0078FF", "#008DFF", "#00A1FF", "#00B5FF",
    "#00CAFF", "#00DEFF", "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", "#FFC900", "#FFB700",
    "#FFA500", "#FF9200", "#FF8000", "#FF6D00", "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000")
}
## get vizgene
# vizgene_bgplot <- function(viz_im) {
#   coordinations <- Seurat::GetTissueCoordinates(viz_im)
#   graphics::plot.new()
#   graphics::par(
#     pty = 's',
#     bg = 'white',
#     mar = c(0, 0, 1, 0),
#     oma = c(0, 0, 1, 0)
#   )
#   graphics::plot(
#     x = coordinations$x,
#     y = coordinations$y,
#     col = ggplot2::alpha("white",
#                          0),
#     xlab = NA_character_,
#     ylab = NA_character_,
#     axes = F,
#     xlim = unit(c(
#       0, max(coordinations$x) - min(coordinations$x)
#     ), 'npc'),
#     ylim = unit(c(
#       0, max(coordinations$y) - min(coordinations$y)
#     ), 'npc')
#   )
#   graphics::points(
#     x = coordinations$x,
#     y = coordinations$y,
#     pch = 19,
#     # cex = pt_size + pt_size * 0.1,
#     cex = unit(1, "npc") * 0.7,
#     col = ggplot2::alpha(colour = 'grey90',
#                          alpha = 1),
#     asp = 1
#   )
# }
## prepared shiny data
