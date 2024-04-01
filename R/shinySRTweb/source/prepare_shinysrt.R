library(shiny)
library(shinyhelper)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(hdf5r)
library(ggdendro)
library(gridExtra)
library(dplyr)
library(scales)
library(patchwork)
library(RColorBrewer)
library(Seurat)
library(shinydashboard)
library(maps)
library(Cairo)
library(dplyr)
library(grid)
library(ggtree)
library(aplot)
library(ggiraph)
library(scatterpie)
library(ggpubr)




preparedata_shinyspatial <- function(
  dat,  
  meta.to.include = NA,
         tooltip_col = NULL,
         maxlevel = 50,
         shiny.dir = 'shinyspatial_app',
         chunkSize = 500,
         gex.assay = NA,
         gex.slot = c("data", "scale.data", "counts"),
         gene.mapping = F,
         colcluster = NULL,
         sp_cols = NULL,
         default.gene1 = NA,
         default.gene2 = NA,
         default.multigene = NA) {
  
  col_box <- c(
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
  
  drExist = TRUE # image array
  dmExist = TRUE # dim matrix
  
  if (length(names(dat@images)) == 0) {
    drExist = FALSE
  }
  if (!drExist) {
    stop(paste0("shinSRT did not detect any coordination data"))
  }
  
  
  ## visium
  image <- dat@images
  ### seurat spot coordination are needed to flip
  coordi <- lapply(names(image), function(x) {
    coordi <- Seurat::GetTissueCoordinates(image[[x]])
    coordi$slice_sample <- x
    xrange <- c(0, dim(image[[x]])[2])
    yrange <- c(0, dim(image[[x]])[1])
    
    ## Image Orientation
    # coordi$imagerow <-
    #   (yrange[2] - yrange[1]) / 2 - (coordi$imagerow - (yrange[2] - yrange[1]) / 2)
    return(coordi)
  })
  names(coordi) <- NULL
  coordi <- do.call(rbind, coordi)
  
  ### get image data
  ima <- lapply(names(image), function(x) {
    image[[x]]@image
  })
  names(ima) <- names(image)
  
  
  
  meta = dat@meta.data
  #
  # ### reduction matrix
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
  if (dmExist) {
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
  #
  cells <- intersect(rownames(coordi), rownames(meta))
  coordi <- coordi[which(rownames(coordi) %in% cells), ]
  meta <- meta[which(rownames(meta) %in% cells), ]
  dat <- dat[, cells]
  #
  #
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
  #
  #
  coordi <- coordi %>% dplyr::distinct(.keep_all = T)
  
  colnames(coordi)[which(duplicated(colnames(coordi)))] <-
    paste(colnames(coordi)[which(duplicated(colnames(coordi)))], sample(length(colnames(coordi)[which(duplicated(colnames(coordi)))])))
  
  
  
  ### config
  meta_group <- lapply(colnames(coordi), function(x) {
    if (length(grep(pattern = 'factor', class(coordi[[x]]))) > 0 |
        length(grep(pattern = 'character', class(coordi[[x]]))) > 0) {
      if (length(grep(pattern = 'character', class(coordi[[x]]))) > 0 &
          length(unique(coordi[[x]])) < maxlevel) {
        coordi[[x]] <- factor(coordi[[x]], levels = unique(coordi[[x]]))
      }
      Unit <- paste(levels(coordi[[x]]), collapse = '|')
      color <-
        paste(col_box[1:length(levels(coordi[[x]]))], collapse = '|')
      metas <- data.frame(matrix(nrow = 1, ncol = 5))
      colnames(metas) <-
        c('group', 'unit', 'color', 'info', 'default')
      if (length(unique(coordi[[x]])) < maxlevel) {
        metas[1, ] <- c(x, Unit, color, TRUE, 0)
      } else{
        metas[1, ] <- c(x, Unit, color, FALSE, 0)
      }
      
      metas
    } else if (length(grep(pattern = '\\.', unique(coordi[[x]]))) == 0 &
               length(unique(coordi[[x]])) < maxlevel) {
      coordi[[x]] <- factor(coordi[[x]], levels = unique(coordi[[x]]))
      Unit <- paste(levels(coordi[[x]]), collapse = '|')
      color <-
        paste(col_box[1:length(levels(coordi[[x]]))], collapse = '|')
      metas <- data.frame(matrix(nrow = 1, ncol = 5))
      colnames(metas) <-
        c('group', 'unit', 'color', 'info', 'default')
      metas[1, ] <- c(x, Unit, color, TRUE, 0)
      metas
    } else{
      metas <- data.frame(matrix(nrow = 1, ncol = 5))
      colnames(metas) <-
        c('group', 'unit', 'color', 'info', 'default')
      metas[1, ] <- c(x, NA, NA, FALSE, 0)
      metas
    }
  })
  #
  meta_group <- do.call(rbind, meta_group)
  meta_group <- data.table::data.table(meta_group)
  
  def1 = grep("ident|library|Ident|Librar",
              meta_group$group[meta_group$info == TRUE],
              ignore.case = TRUE)[1]
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
  
  
  if (is.na(gex.assay[1])) {
    gex.assay = Seurat::DefaultAssay(dat)
  }
  
  ### data dimation
  if(class(dat@assays[[gex.assay[1]]]) == 'Assay5'){
    gex.matdim <- dim(dat@assays[[gex.assay[1]]])
    gex.rownm  <- rownames(dat@assays[[gex.assay[1]]])
    gex.colnm <- colnames(dat@assays[[gex.assay[1]]])
  }else{
    gex.matdim = dim(slot(dat@assays[[gex.assay[1]]], gex.slot[1]))
    gex.rownm = rownames(slot(dat@assays[[gex.assay[1]]],
                              gex.slot[1]))
    gex.colnm = colnames(slot(dat@assays[[gex.assay[1]]],
                              gex.slot[1]))
    
  }
  
  if(is.null(Seurat::VariableFeatures(dat))){
    if(length(intersect(Layers(dat),'data')) == 0){
      dat <- NormalizeData(dat,normalization.method = 'LogNormalize')
    }
    dat <- FindVariableFeatures(dat,selection.method = 'vst')
  }
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
  
  
  if (gene.mapping[1] == TRUE) {
    if (sum(grepl("^ENSG000", gex.rownm)) >= sum(grepl("^ENSMUSG000",gex.rownm))) {
      tmp1 = read.table(system.file("extdata",'hg_map.txt.gz',package = 'shinySRT'),header = T)
      # tmp1 = data.table::fread(system.file("extdata",'hg_map.txt.gz',package = 'shinySRT'))
    } else {
      tmp1 <- read.table(system.file("extdata",'mm_map.txt.gz',package = 'shinySRT'),header = T)
      # tmp1 = data.table::fread(system.file("extdata",'mm_map.txt.gz',package = 'shinySRT'))
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
  
  
  
    da <- Seurat::GetAssayData(dat,slot = gex.slot[1])
  
  for (i in 1:floor((gex.matdim[1] - 8) / chk)) {
    mat_exp.grp.data[((i - 1) * chk + 1):(i * chk),] <-
      as.matrix(da[((i - 1) * chk + 1):(i * chk),])
  }
  mat_exp.grp.data[(i * chk + 1):gex.matdim[1],] <-
    as.matrix(da[(i * chk + 1):gex.matdim[1],])
  mat_exp$close_all()
  
  
  
  
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
  
  
  ## spot lab
  if (length(intersect(tooltip_col, colnames(coordi))) == 0 |
      length(intersect(tooltip_col, meta_group[default == 1 &
                                               info == TRUE]$group)) > 0 |
      length(intersect(tooltip_col, meta_group[default == 2 &
                                               info == TRUE]$group) > 0)) {
    coordi$tooltip <- paste(
      'sampleID: ',
      coordi$sampleID,
      '\n',
      'slice :',
      coordi$slice_sample,
      '\n',
      meta_group[default == 1 &
                   info == TRUE]$group,
      ': ',
      coordi[[meta_group[default == 1 & info == TRUE]$group]],'\n',
      meta_group[default == 2 &
                   info == TRUE]$group,
      ': ',
      coordi[[meta_group[default == 2 & info == TRUE]$group]]
    )
    
  } else{
    coordi$tooltip <- paste(
      'sampleID: ',
      coordi$sampleID,
      '\n',
      'slice :',
      coordi$slice_sample,
      '\n',
      meta_group[default == 1 &
                   info == TRUE]$group,
      ': ',
      coordi[[meta_group[default == 1 & info == TRUE]$group]],
      '\n',
      meta_group[default == 2 &
                   info == TRUE]$group,
      ': ',
      coordi[[meta_group[default == 2 & info == TRUE]$group]],
      '\n',
      tooltip_col,
      ' :',
      coordi[[tooltip_col]]
    )
  }
  
  

  
  
  
  
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




prepare_code <- function(shiny.dir = 'shinyspatial_app',title = 'spatial_example'){
  
  filename = paste0(shiny.dir, "/server.R")
  df_select <- readRDS(paste0(shiny.dir,"/df_select.Rds"))
  # library
  readr::write_file(lib_server(),file = filename)
  # color
  readr::write_file(col_server(),file = filename,append = T)
  # load data
  readr::write_file(data_server(df_select = df_select),file = filename,append = T)
  # function
  readr::write_file(fun_server(),file = filename,append = T)
  # page
  readr::write_file(server_heads(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p1(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p2(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p3(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p4(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p5(df_select = df_select),file = filename,append = T)
  
  if(length(grep(pattern = 'deconvolusion',dir(shiny.dir),value = T)) > 0){
    readr::write_file(sp_server_p6(df_select = df_select),file = filename,append = T)
  }else{
    readr::write_file('}',file = filename,append = T)
  }
  
  ##
  filename = paste0(shiny.dir, "/ui.R")
  readr::write_file(ui_load(), file = filename)
  readr::write_file(
    glue::glue(
      'container <- function(...) {{\n',
      '  shiny::fluidRow(shiny::column(...))\n',
      '}}\n\n\n'
    ),
    file = filename,
    append = T
  )
  
  readr::write_file(ui_head(title,df_select = df_select), file = filename, append = T)
  ## page1
  readr::write_file(ui_p1(df_select), file = filename, append = T)
  ## page2
  readr::write_file(ui_p2(df_select), file = filename, append = T)
  ## page3
  readr::write_file(ui_p3(df_select), file = filename,append = T)
  ## page4
  readr::write_file(ui_p4(df_select = df_select), file = filename,append = T)
  ## page5
  readr::write_file(ui_p5(df_select = df_select), file = filename,append = T)
  ## page6
  if(!is.null(df_select[['cell']])){
    readr::write_file(ui_p6(df_select = df_select), file = filename,append = T)
  }
  
  readr::write_file(glue::glue(')))'), file = filename,append = T)
}









CreateshinySRT <- function(dat,
                           meta.to.include = NA,
                           tooltip_col = NULL,
                           maxlevel = 50,
                           shiny.dir = 'shinyspatial_app',
                           title = 'spatial_example',
                           chunkSize = 500,
                           gex.assay = NA,
                           gex.slot = c("data", "scale.data", "counts"),
                           gene.mapping = TRUE,
                           colcluster = NULL,
                           sp_cols = NULL,
                           default.gene1 = NA,
                           default.gene2 = NA,
                           default.multigene = NA) {
  preparedata_shinyspatial(
    dat,
    meta.to.include = meta.to.include,
    tooltip_col = tooltip_col,
    maxlevel = maxlevel,
    shiny.dir = shiny.dir,
    chunkSize = chunkSize,
    gex.assay = gex.assay,
    gex.slot = gex.slot,
    gene.mapping = gene.mapping,
    colcluster = colcluster,
    sp_cols = sp_cols,
    default.gene1 = default.gene1,
    default.gene2 = default.gene2,
    default.multigene = default.multigene
  )
  prepare_code(shiny.dir,title)
}
