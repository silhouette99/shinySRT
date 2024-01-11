#' Add new group for spatial spots
#'
#' add some information into meta.Rds and meta_group.Rds
#'
#'
#' @param add_meta new group of spots, characters or data frame/data table
#' @param dirs directory shinyapp
#' @param colname group names, if add_meta is the data frame, colname could be used to identify the column of add_meta
#'
#'
#' @import data.table hdf5r reticulate dplyr SpatialExperiment SingleCellExperiment Seurat SummarizedExperiment scran dplyr
#' 
#' 
#'
#'
#'
#' @export


add_meta <- function(add_meta,dirs = 'shinyspatial_app',colname = 'new_group'){
  if(length(which(dir(dirs) == "meta.Rds")) == 0 | length(which(dir(dirs) == "meta_group.Rds")) == 0){
      stop('Do not find the meta for shinySRT! ')
  }
 
  meta <- readRDS(paste0(dirs,'/meta.Rds')) 
  meta_group <- readRDS(paste0(dirs,'/meta_group.Rds'))
  
  if(class(add_meta) == 'character'){
    if(length(add_meta) != nrow(meta)){
      stop('length of new group do not equal to rows of meta!')
    }
    if(!is.null(names(add_meta)) & length(intersect(names(add_meta), meta$sampleID)) == nrow(meta)){
      add_meta <- add_meta[match(meta$sampleID, names(add_meta))]
    }
    meta[,c(colname) := add_meta]
      
  }else if(class(add_meta) == 'data.frame'|class(add_meta) == 'data.table'){
    if(nrow(add_meta) != nrow(meta)){
      stop('length of new group do not equal to rows of meta!')
    }
    if(!is.null(rownames(add_meta)) & length(intersect(rownames(add_meta), meta$sampleID)) == nrow(meta)){
      add_meta <- add_meta[match(meta$sampleID, rownames(add_meta)),]
    }
    meta[,c(colname) := add_meta[[colname]]]
  }
  
  meta[[colname]] <- as.factor(meta[[colname]])
  unit <- paste(levels(meta[[colname]]),collapse = '|')
  color <- paste(col_box()[1:length(levels(meta[[colname]]))],collapse = '|')
  news <- data.table(colname, unit, color, TRUE, 0)
  colnames(news) <- colnames(meta_group)
  
  meta_group <- rbind(meta_group,news)
  
  saveRDS(meta,paste0(dirs,'/meta.Rds'))
  saveRDS(meta_group, paste0(dirs,'/meta_group.Rds'))
  
}