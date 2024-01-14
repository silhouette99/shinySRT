#' Deconvolusion of Giotto in shinySRT
#'
#' 
#' This function could be used to predict the cellular component of each spots, refering to Giotto's DWLS deconvolusion function. 
#' 
#' 
#' @param norm_mtx normalized data matrix (or sparse)
#' @param logbase base used for log normalization
#' @param cell_metadata spatial meta
#' @param cluster_column name of cluster column
#' @param sign_matrix signature gene-cell matrix for deconvolusion
#' @param  number of cells per spot
#'
#'
#'
#' @return data.table of cellular component 
#'
#'
#' @import data.table hdf5r reticulate dplyr SpatialExperiment SingleCellExperiment Seurat SummarizedExperiment quadprog Rfast scran
#'
#'
#' @export


## Giotto DWLS deconvolusion function

runDWLSDeconv <-
  function (norm_mtx,
            logbase = 2,
            cell_metadata,
            cluster_column,
            sign_matrix,
            n_cell = 50,
            cutoff = 2,
            name = NULL) {
    # package_check(pkg_name = "quadprog", repository = "CRAN")
    # package_check(pkg_name = "Rfast", repository = "CRAN")
    
    
    
    nolog_expr = logbase ^ (norm_mtx) - 1
    
    if (!cluster_column %in% colnames(cell_metadata)) {
      stop("\n cluster column not found \n")
    }
    cluster = cell_metadata[[cluster_column]]
    sign_matrix = as.matrix(sign_matrix)
    intersect_gene = intersect(rownames(sign_matrix), rownames(nolog_expr))
    filter_Sig = sign_matrix[intersect_gene,]
    filter_expr = nolog_expr[intersect_gene,]
    filter_log_expr = norm_mtx[intersect_gene,]
    enrich_spot_proportion = enrich_deconvolution(
      expr = filter_expr,
      log_expr = filter_log_expr,
      cluster_info = cluster,
      ct_exp = filter_Sig,
      cutoff = cutoff
    )
    
    resolution = (1 / n_cell)
    binarize_proportion = ifelse(enrich_spot_proportion >= resolution,
                                 1, 0)
    spot_proportion <- spot_deconvolution(
      expr = filter_expr,
      cluster_info = cluster,
      ct_exp = filter_Sig,
      binary_matrix = binarize_proportion
    )
    deconvolutionDT = data.table::data.table(sampleID = colnames(spot_proportion))
    deconvolutionDT = cbind(deconvolutionDT, data.table::as.data.table(t(spot_proportion)))
    
    rownames(enrich_spot_proportion) <-
      paste0('score_', rownames(enrich_spot_proportion))
    enrich_spot_proportion <- t(enrich_spot_proportion)
    deconvolutionDT <-
      cbind(deconvolutionDT, enrich_spot_proportion)
    
    for(i in sapply(strsplit(grep(pattern = '^score_',colnames(deconvolutionDT), value = T),split = '_'),'[[',2)){
      deconvolutionDT[[i]] <- round(deconvolutionDT[[i]],digits = 5)
    }
    
    toolstip <- lapply(grep(pattern = 'score_',colnames(deconvolutionDT),value = T), function(i){
      cell <- strsplit(i,split = '_')[[1]][2]
      paste0(cell,' : ',signif(deconvolutionDT[[i]]), '\n')
    })
    
    toolstip <- do.call(paste0, toolstip)
    
    deconvolutionDT$tooltip <- paste0('sampleID', ' : ', deconvolutionDT$sampleID, '\n', toolstip)
    
    return(deconvolutionDT)
  }

#' 
#' Scran normalization
#' 
#' Normalize the count data present in a given assay.
#' 
#' 
#' @param mymatrix transcriptomic data matrix
#' @param log_offset offset
#' @param logbase base used for log normalization
#' 
#' @return a normalized data matrix 
#' 
#' @export

scran_norm <-
  function(mymatrix,
           log_offset = 1,
           logbase = 2,
           scalefactor = 6000,
           library_size_norm = TRUE,
           log_norm = TRUE) {
    cells <- colnames(mymatrix)
    genes <- rownames(mymatrix)
    
    if (class(mymatrix)[[1]] == 'matrix') {
      mymatrix <- as(mymatrix, 'dgCMatrix')
    } else if (class(mymatrix)[[1]] == 'data') {
      mymatrix <- as.matrix(mymatrix)
      mymatrix <- as(mymatrix, 'dgCMatrix')
    }
    
    if (library_size_norm == TRUE) {
      libsizes <- Matrix::colSums(mymatrix)
      
      if (any(libsizes == 0)) {
        warning(
          
            'Total library size or counts for individual spat units are 0.
                     This will likely result in normalization problems.
                     filter (filterGiotto) or impute (imputeGiotto) spatial units.'
          
        )
      }
      
      mymatrix = Matrix::t(Matrix::t(mymatrix) / libsizes) * scalefactor
    }
    
    if (log_norm == TRUE) {
      mymatrix@x = log(mymatrix@x + log_offset) / log(logbase)
    }
    
    rownames(mymatrix) <- genes
    colnames(mymatrix) <- cells
    
    return(mymatrix)
  }


#'
#' Find highly variable gene
#' 
#' 
#' Identifies features
#' 
#' @param mtx normalized data matrix
#' 
#'
#'
#'
#'
#'
#'
#'
#'@export


scran_hvg <- function(mtx,
                      meta,
                      colcluster,
                      method = "scran",
                      pval = 0.01,
                      logFC = 0.5,
                      min_feats = 3) {
  if (class(meta)[1] != "data.table") {
    metas <- data.table::data.table(meta)
    metas[, `:=`(sampleID, rownames(meta))]
  }
  
  if(class(metas[[colcluster]]) != 'factor'){
    metas[[colcluster]] <- factor(metas[[colcluster]], levels = unique(metas[[colcluster]]))
  }
  # if (method == "scran") {
  deg <- lapply(levels(metas[[colcluster]]), function(x) {
    g1 = x
    g2 = levels(metas[[colcluster]])[levels(metas[[colcluster]]) != x]
    
    group_1_name = paste0(g1, collapse = "_")
    group_2_name = paste0(g2, collapse = "_")
    
    pairwise_select_comp = NULL
    metas[, `:=`(pairwise_select_comp,
                 ifelse(get(colcluster) %in%
                          g1, group_1_name, group_2_name))]
    
    deg <-
      scran::findMarkers(mtx, groups = metas$pairwise_select_comp)
    genes = clusters = NULL
    
    deg = lapply(names(deg), function(i) {
      dfr = deg[[i]]
      DT = data.table::as.data.table(dfr)
      DT[, `:=`(genes, rownames(dfr))]
      DT[, `:=`(clusters, i)]
    })
    
    select_bool = unlist(lapply(deg, function(m) {
      unique(m$clusters) == x
    }))
    
    selected_table = data.table::as.data.table(deg[select_bool])
    col_ind_keep = !grepl("summary", colnames(selected_table))
    selected_table = selected_table[, col_ind_keep, with = F]
    
    data.table::setnames(selected_table, colnames(selected_table)[4],
                         "logFC")
    data.table::setnames(selected_table, colnames(selected_table)[5],
                         "genes")
    
    filtered_table = selected_table[logFC > 0]
    filtered_table[, `:=`("ranking", rank(-logFC))]
    p.value = ranking = NULL
    
    filtered_table = filtered_table[(p.value <= pval &
                                       logFC >= logFC) |
                                      (ranking <= min_feats)]
    # pb(message = c("cluster ", clus_i, "/", length(uniq_clusters)))
    return(filtered_table)
    
  })
  deg <- do.call(rbind, deg)
  # }
  return(deg)
}



#' 
#' 
#' 
#' Get DWLS signature gene-cell matrix for deconvolusion
#' 
#' This function uses the hvg of each cell type as cell type signatures
#' 
#' 
#' @param matrix normalized scRNA-seq matrix
#' @param sign_gene signature gene
#' @param cell_cols colums of cell type
#' @param meta scRNA-seq meta
#'
#'@export

DWLSmatrix <- function (matrix, sign_gene, cell_cols, meta)
{
  if (ncol(matrix) != length(meta[[cell_cols]])) {
    stop("ncol(matrix) needs to be the same as length(cell_type_vector)")
  }
  if (!is.character(sign_gene)) {
    stop("\n sign_gene needs to be a character vector of cell type specific genes \n")
  }
  intersect_sign_gene = intersect(rownames(matrix), sign_gene)
  matrix_subset = matrix[intersect_sign_gene,]
  signMatrix = matrix(
    data = NA,
    nrow = nrow(matrix_subset),
    ncol = length(unique(meta[[cell_cols]]))
  )
  for (cell_type_i in seq_along(unique(meta[[cell_cols]]))) {
    cell_type = unique(meta[[cell_cols]])[cell_type_i]
    selected_cells = colnames(matrix_subset)[meta[[cell_cols]] == cell_type]
    mean_expr_in_selected_cells = Matrix::rowMeans(matrix_subset[, selected_cells])
    signMatrix[, cell_type_i] = mean_expr_in_selected_cells
  }
  rownames(signMatrix) = rownames(matrix_subset)
  colnames(signMatrix) = unique(meta[[cell_cols]])
  return(signMatrix)
}