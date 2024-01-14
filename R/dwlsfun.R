#' function of Giotto DWLS deconvolusion
#'
#'
#'
#'
#'
#'
#'
#' @export

enrich_deconvolution <- function(expr,
                                 log_expr,
                                 cluster_info,
                                 ct_exp,
                                 cutoff) {
  #####generate enrich 0/1 matrix based on expression matrix
  ct_exp <- ct_exp[rowSums(ct_exp) > 0, ]
  enrich_matrix <- matrix(0, nrow = dim(ct_exp)[1], ncol = dim(ct_exp)[2])
  rowmax_col <- Rfast::rowMaxs(ct_exp)
  for (i in seq_along(rowmax_col)) {
    enrich_matrix[i, rowmax_col[i]] = 1
  }
  colsum_ct_binary <- colSums(enrich_matrix)
  for (i in seq_along(colsum_ct_binary)) {
    if (colsum_ct_binary[i] <= 2) {
      rank <- rank(-ct_exp[, i])
      enrich_matrix[rank <= 2, i] = 1
    }
  }
  rownames(enrich_matrix) <- rownames(ct_exp)
  colnames(enrich_matrix) <- colnames(ct_exp)
  
  #####page enrich
  enrich_result <- enrich_analysis(log_expr, enrich_matrix)
  #####initialize dwls matrix
  dwls_results <-
    matrix(0, nrow = dim(enrich_matrix)[2], ncol = dim(expr)[2])
  rownames(dwls_results) <- colnames(enrich_matrix)
  colnames(dwls_results) <- colnames(expr)
  cluster_sort <- sort(unique(cluster_info))
  cluster_info <- cluster_info
  for (i in seq_along(cluster_sort)) {
    cluster_i_enrich <-
      enrich_result[, which(cluster_info == cluster_sort[i])]
    row_i_max <- Rfast::rowMaxs(cluster_i_enrich, value = TRUE)
    ct <- rownames(enrich_result)[which(row_i_max > cutoff)]
    if (length(ct) < 2) {
      sort_rank <- sort(row_i_max, decreasing = T)
      ct <- rownames(enrich_result)[which(row_i_max >= sort_rank[2])]
    }
    ct_gene <- c()
    for (j in seq_along(ct)) {
      sig_gene_j <- rownames(enrich_matrix)[which(enrich_matrix[, ct[j]] == 1)]
      ct_gene <- c(ct_gene, sig_gene_j)
    }
    uniq_ct_gene <- intersect(rownames(expr), unique(ct_gene))
    select_sig_exp <- ct_exp[uniq_ct_gene, ct]
    cluster_i_cell <- which(cluster_info == cluster_sort[i])
    cluster_cell_exp <- expr[uniq_ct_gene, cluster_i_cell]
    
    cluster_i_dwls <-
      optimize_deconvolute_dwls(cluster_cell_exp, select_sig_exp)
    dwls_results[ct, cluster_i_cell] <- cluster_i_dwls
  }
  #####remove negative values
  for (i in dim(dwls_results)[1]) {
    negtive_index <- which(dwls_results[i, ] < 0)
    dwls_results[i, negtive_index] == 0
  }
  return(dwls_results)
}


enrich_analysis <- function(expr_values,
                            sign_matrix) {
  # output enrichment
  # only continue with genes present in both datasets
  interGene = intersect(rownames(sign_matrix), rownames(expr_values))
  filterSig = sign_matrix[interGene, ]
  signames = rownames(filterSig)[which(filterSig[, 1] == 1)]
  # calculate mean gene expression
  #mean_gene_expr = rowMeans(expr_values)
  mean_gene_expr = log2(Matrix::rowMeans(2 ^ expr_values - 1, dims = 1) +
                          1)
  geneFold = expr_values - mean_gene_expr
  # calculate sample/spot mean and sd
  cellColMean = apply(geneFold, 2, mean)
  cellColSd = apply(geneFold, 2, stats::sd)
  
  # get enrichment scores
  enrichment = matrix(
    data = NA,
    nrow = dim(filterSig)[2],
    ncol = length(cellColMean)
  )
  for (i in (1:dim(filterSig)[2])) {
    signames = rownames(filterSig)[which(filterSig[, i] == 1)]
    sigColMean = apply(geneFold[signames, ], 2, mean)
    m = length(signames)
    vectorX = NULL
    for (j in(seq_along(cellColMean))) {
      Sm = sigColMean[j]
      u = cellColMean[j]
      sigma = cellColSd[j]
      zscore = (Sm - u) * m ^ (1 / 2) / sigma
      vectorX = append(vectorX, zscore)
    }
    enrichment[i, ] = vectorX
  }
  rownames(enrichment) = colnames(filterSig)
  colnames(enrichment) = names(cellColMean)
  return(enrichment)
}


optimize_deconvolute_dwls <- function(exp,
                                      Signature) {
  ######overlap signature with spatial genes
  Genes <- intersect(rownames(Signature), rownames(exp))
  S <- Signature[Genes, ]
  S <- Matrix::as.matrix(S)
  Bulk <- Matrix::as.matrix(exp)
  subBulk = Bulk[Genes, ]
  allCounts_DWLS <- NULL
  all_exp <- Matrix::rowMeans(exp)
  
  solution_all_exp <- solve_OLS_internal(S, all_exp[Genes])
  
  constant_J <-
    find_dampening_constant(S, all_exp[Genes], solution_all_exp)
  for (j in 1:(dim(subBulk)[2])) {
    B <- subBulk[, j]
    if (sum(B) > 0) {
      solDWLS <- optimize_solveDampenedWLS(S, B, constant_J)
    } else{
      solDWLS <- rep(0, length(B))
      names(solDWLS) <- names(B)
    }
    allCounts_DWLS <- cbind(allCounts_DWLS, solDWLS)
  }
  colnames(allCounts_DWLS) <- colnames(exp)
  return(allCounts_DWLS)
}



solve_OLS_internal <- function(S,
                               B) {
  D = t(S) %*% S
  d = t(S) %*% B
  A = cbind(diag(dim(S)[2]))
  bzero = c(rep(0, dim(S)[2]))
  
  
  out = tryCatch(
    expr = {
      quadprog::solve.QP(
        Dmat = D,
        dvec = d,
        Amat = A,
        bvec = bzero
      )$solution
    },
    
    error = function(cond) {
      message('Original error message: \n')
      message(cond)
      
      message('\n Will try to fix error with Matrix::nearPD()')
      return(NULL)
    }
  )
  
  
  if (is.null(out)) {
    D = Matrix::nearPD(D)
    D = as.matrix(D$mat)
    
    out = tryCatch(
      expr = {
        quadprog::solve.QP(
          Dmat = D,
          dvec = d,
          Amat = A,
          bvec = bzero
        )$solution
      },
      
      error = function(cond) {
        message('Original error message: \n')
        message(cond)
        
        message('\n nearPD() did not fix the error')
        return(NULL)
      }
    )
    
    if (is.null(out)) {
      stop('Errors could not be fixed')
    } else {
      names(out) = colnames(S)
    }
    
  } else {
    names(out) = colnames(S)
  }
  
  return(out)
}




find_dampening_constant <- function(S,
                                    B,
                                    goldStandard) {
  solutionsSd = NULL
  
  #goldStandard is used to define the weights
  sol = goldStandard
  ws = as.vector((1 / (S %*% sol)) ^ 2)
  wsScaled = ws / min(ws)
  wsScaledMinusInf = wsScaled
  
  if (is.na(max(wsScaled))) {
    wsScaled = wsScaled[-which(is.na(wsScaled))]
    wsScaledMinusInf = wsScaled
  }
  
  #ignore infinite weights
  if (max(wsScaled) == "Inf") {
    wsScaledMinusInf = wsScaled[-which(wsScaled == "Inf")]
  }
  
  
  
  #try multiple values of the dampening constant (multiplier)
  #for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))) {
    multiplier = 1 * 2 ^ (j - 1)
    wsDampened = wsScaled
    wsDampened[which(wsScaled > multiplier)] = multiplier
    solutions = NULL
    seeds = c(1:100)
    for (i in 1:100) {
      set.seed(seeds[i]) #make nondeterministic
      subset = sample(length(ws), size = length(ws) * 0.5) #randomly select half of gene set
      
      #solve dampened weighted least squares for subset
      fit = stats::lm (B[subset] ~ -1 + S[subset, ], weights = wsDampened[subset])
      sol = fit$coef * sum(goldStandard) / sum(fit$coef)
      solutions = cbind(solutions, sol)
    }
    solutionsSd = cbind(solutionsSd, apply(solutions, 1, stats::sd))
  }
  
  #choose dampening constant that results in least cross-validation variance
  j = which.min(colMeans(solutionsSd ^ 2))
  return(j)
}




optimize_solveDampenedWLS <- function(S,
                                      B,
                                      constant_J) {
  #first solve OLS, use this solution to find a starting point for the weights
  solution = solve_OLS_internal(S, B)
  #now use dampened WLS, iterate weights until convergence
  iterations = 0
  changes = c()
  #find dampening constant for weights using cross-validation
  j = constant_J
  change = 1
  
  while (change > .01 & iterations < 1000) {
    newsolution = solve_dampened_WLSj(S, B, solution, j)
    #decrease step size for convergence
    solutionAverage = rowMeans(cbind(newsolution,
                                     matrix(
                                       solution,
                                       nrow = length(solution),
                                       ncol = 4
                                     )))
    change = norm(Matrix::as.matrix(solutionAverage - solution))
    solution = solutionAverage
    iterations = iterations + 1
    changes = c(changes, change)
  }
  
  
  return(solution / sum(solution))
}




solve_dampened_WLSj <- function(S,
                                B,
                                goldStandard,
                                j) {
  multiplier = 1 * 2 ^ (j - 1)
  sol = goldStandard
  ws = as.vector((1 / (S %*% sol)) ^ 2)
  wsScaled = ws / min(ws)
  wsDampened = wsScaled
  wsDampened[which(wsScaled > multiplier)] = multiplier
  W = diag(wsDampened)
  D = t(S) %*% W %*% S
  d = t(S) %*% W %*% B
  A = cbind(diag(dim(S)[2]))
  bzero = c(rep(0, dim(S)[2]))
  sc = norm(D, "2")
  
  D_positive_definite <- Matrix::nearPD(x = D / sc)
  
  solution <-
    quadprog::solve.QP(
      Dmat = as.matrix(D_positive_definite$mat),
      dvec = d / sc,
      Amat = A,
      bvec = bzero
    )$solution
  
  names(solution) = colnames(S)
  return(solution)
}


#######################

spot_deconvolution <- function(expr,
                               cluster_info,
                               ct_exp,
                               binary_matrix) {
  #####generate enrich 0/1 matrix based on expression matrix
  enrich_matrix <- matrix(0, nrow = dim(ct_exp)[1], ncol = dim(ct_exp)[2])
  rowmax_col <- Rfast::rowMaxs(ct_exp)
  for (i in seq_along(rowmax_col)) {
    enrich_matrix[i, rowmax_col[i]] = 1
  }
  rownames(enrich_matrix) <- rownames(ct_exp)
  colnames(enrich_matrix) <- colnames(ct_exp)
  
  cluster_sort <- sort(unique(cluster_info))
  ####initialize dwls matrix
  dwls_results <- matrix(0, nrow = dim(ct_exp)[2], ncol = dim(expr)[2])
  rownames(dwls_results) <- colnames(ct_exp)
  colnames(dwls_results) <- colnames(expr)
  
  for (i in seq_along(cluster_sort)) {
    cluster_i_matrix <-
      binary_matrix[, which(cluster_info == cluster_sort[i])]
    row_i_max <- Rfast::rowMaxs(cluster_i_matrix, value = TRUE)
    ct_i <- rownames(cluster_i_matrix)[which(row_i_max == 1)]
    ########calculate proportion based on binarized deconvolution results at first step
    if (length(ct_i) == 1) {
      dwls_results[ct_i[1], which(cluster_info == cluster_sort[i])] == 1
    } else {
      ct_gene <- c()
      for (j in seq_along(ct_i)) {
        sig_gene_j <-
          rownames(enrich_matrix)[which(enrich_matrix[, ct_i[j]] == 1)]
        ct_gene <- c(ct_gene, sig_gene_j)
      }
      uniq_ct_gene <- intersect(rownames(expr), unique(ct_gene))
      select_sig_exp <- ct_exp[uniq_ct_gene, ct_i]
      cluster_i_cell <- which(cluster_info == cluster_sort[i])
      cluster_cell_exp <- expr[uniq_ct_gene, cluster_i_cell]
      ######calculate
      ######overlap signature with spatial genes
      all_exp <- Matrix::rowMeans(cluster_cell_exp)
      solution_all_exp <- solve_OLS_internal(select_sig_exp, all_exp)
      constant_J <-
        find_dampening_constant(select_sig_exp, all_exp, solution_all_exp)
      ######deconvolution for each spot
      for (k in 1:(dim(cluster_cell_exp)[2])) {
        B <- Matrix::as.matrix(cluster_cell_exp[, k])
        ct_spot_k <-
          rownames(cluster_i_matrix)[which(cluster_i_matrix[, k] == 1)]
        if (length(ct_spot_k) == 1) {
          dwls_results[ct_spot_k[1], colnames(cluster_cell_exp)[k]] <- 1
        } else {
          ct_k_gene <- c()
          for (m in seq_along(ct_spot_k)) {
            sig_gene_k <-
              rownames(enrich_matrix)[which(enrich_matrix[, ct_spot_k[m]] == 1)]
            ct_k_gene <- c(ct_k_gene, sig_gene_k)
          }
          uniq_ct_k_gene <-
            intersect(rownames(ct_exp), unique(ct_k_gene))
          S_k <- Matrix::as.matrix(ct_exp[uniq_ct_k_gene, ct_spot_k])
          solDWLS <-
            optimize_solveDampenedWLS(S_k, B[uniq_ct_k_gene, ], constant_J)
          dwls_results[names(solDWLS), colnames(cluster_cell_exp)[k]] <-
            solDWLS
        }
      }
    }
  }
  #####remove negative values
  for (i in dim(dwls_results)[1]) {
    negtive_index <- which(dwls_results[i, ] < 0)
    dwls_results[i, negtive_index] == 0
  }
  return(dwls_results)
}


#' enrich
#'
#'
#'
#'
#'
#'
#'
#'
#'


page_dt_method = function(sign_matrix,
                          expr_values,
                          min_overlap_genes = 5,
                          logbase = 2,
                          reverse_log_scale = TRUE,
                          output_enrichment = c('original', 'zscore'),
                          p_value = FALSE,
                          include_depletion = FALSE,
                          n_times = 1000,
                          max_block = 20e6,
                          verbose = TRUE) {
  
  
  # data.table variables
  Var1 = value = Var2 = V1 = marker = nr_markers = fc = cell_ID = zscore = colmean = colSd = pval = NULL
  mean_zscore = sd_zscore = pval_score = NULL
  
  # output enrichment
  output_enrichment = match.arg(output_enrichment, choices = c('original', 'zscore'))
  
  ## identify available cell types
  all_genes = rownames(expr_values)
  sign_matrix = as.matrix(sign_matrix)
  sign_matrix_DT = data.table::as.data.table(reshape2::melt(sign_matrix))
  sign_matrix_DT = sign_matrix_DT[Var1 %in% all_genes]
  detected_DT = sign_matrix_DT[, sum(value), by = Var2]
  
  lost_cell_types_DT = detected_DT[V1 <= min_overlap_genes]
  if(nrow(lost_cell_types_DT) > 0) {
    for(row in 1:nrow(lost_cell_types_DT)) {
      output = paste0("Warning, ",lost_cell_types_DT[row][['Var2']]," only has ",lost_cell_types_DT[row][['V1']]," overlapping genes. Will be removed.")
      if(verbose) print(output)
    }
  }
  available_ct = as.character(detected_DT[V1 > min_overlap_genes][['Var2']])
  
  if (length(available_ct) == 1){
    stop("Only one cell type available.")
  }
  
  # create subset of sinature matrix
  interGene = intersect(rownames(sign_matrix), rownames(expr_values))
  filterSig = sign_matrix[interGene, available_ct]
  
  # create fold expression for each gene in each spot
  # calculate mean gene expression
  if(reverse_log_scale == TRUE) {
    mean_gene_expr = log(rowMeans(logbase^expr_values-1, dims = 1)+1)
  } else {
    mean_gene_expr = rowMeans(expr_values)
  }
  geneFold = expr_values - mean_gene_expr
  
  # calculate sample/spot mean and sd
  cellColMean = colMeans(geneFold)
  cellColSd = apply(geneFold, 2, stats::sd)
  cellColMeanSd =  data.table::data.table(cell_ID = names(cellColMean),
                                          colmean = cellColMean,
                                          colSd = cellColSd)
  
  filterSig_DT = data.table::as.data.table(reshape2::melt(filterSig))
  colnames(filterSig_DT) = c('gene', 'cell_type', 'marker')
  sub_ct_DT = filterSig_DT[marker == 1]
  sub_ct_DT[, nr_markers := .N, by = cell_type]
  
  ## reshape gene fold-expression
  geneFold_DT = data.table::as.data.table(reshape2::melt(geneFold))
  colnames(geneFold_DT) = c('gene', 'cell_ID', 'fc')
  
  mergetest = data.table::merge.data.table(sub_ct_DT, geneFold_DT, by = 'gene')
  mergetest = mergetest[, mean(fc), by = .(cell_type, cell_ID, nr_markers)]
  if (is.integer(mergetest$cell_ID) && is.character(cellColMeanSd$cell_ID)){
    mergetest$cell_ID = as.character(mergetest$cell_ID)
  }
  mergetest = data.table::merge.data.table(mergetest, cellColMeanSd, by = 'cell_ID')
  mergetest[, zscore := ((V1 - colmean)* nr_markers^(1/2)) / colSd]
  
  if(output_enrichment == 'zscore') {
    mergetest[, zscore := scale(zscore), by = 'cell_type']
  }
  
  
  
  
  ## return p-values based on permutations ##
  if(p_value == TRUE) {
    
    ## 1. get number of markers instructions ##
    sample_intrs = unique(sub_ct_DT[,.(cell_type, nr_markers)])
    
    
    ## 2. first create the random samples all together ##
    cell_type_list = list()
    perm_type_list = list()
    for(row in 1:nrow(sample_intrs)) {
      
      cell_type = sample_intrs[row][['cell_type']]
      nr_genes = as.numeric(sample_intrs[row][['nr_markers']])
      
      gene_list = list()
      perm_list = list()
      for(i in 1:n_times) {
        sampled_genes = sample(rownames(expr_values), size = nr_genes)
        gene_list[[i]] = sampled_genes
        perm_list[[i]] = rep(paste0('p_',i), nr_genes)
      }
      
      gene_res = unlist(gene_list)
      names(gene_res) = rep(cell_type, length(gene_res))
      cell_type_list[[row]] = gene_res
      
      perm_res = unlist(perm_list)
      perm_type_list[[row]] = perm_res
      
    }
    
    cell_type_perm = unlist(cell_type_list)
    perm_round = unlist(perm_type_list)
    
    cell_type_perm_DT = data.table::data.table(cell_type = names(cell_type_perm),
                                               gene = cell_type_perm,
                                               round = perm_round)
    
    sample_intrs_vec = sample_intrs$nr_markers
    names(sample_intrs_vec) = sample_intrs$cell_type
    cell_type_perm_DT[, nr_markers := sample_intrs_vec[cell_type]]
    
    
    ## 3. decide on number of blocks to process ##
    nr_perm_lines = as.numeric(nrow(cell_type_perm_DT))
    nr_spots = as.numeric(ncol(expr_values))
    total_lines = nr_spots * nr_perm_lines
    nr_groups = round(total_lines / max_block)
    
    ## 4. create groups
    all_perms = unique(perm_round)
    all_perms_num = seq_along(all_perms)
    names(all_perms_num) = all_perms
    group_labels = paste0('group_',1:nr_groups)
    groups_vec = cut(all_perms_num, breaks = nr_groups, labels = group_labels)
    names(all_perms) = groups_vec
    
    
    ## 5. do random enrichment per block
    res_list = list()
    for(group_i in seq_along(group_labels)) {
      
      group = group_labels[group_i]
      sub_perms = all_perms[names(all_perms) == group]
      cell_type_perm_DT_sub = cell_type_perm_DT[round %in% sub_perms]
      
      mergetest_perm_sub = data.table::merge.data.table(cell_type_perm_DT_sub, geneFold_DT, allow.cartesian = TRUE)
      mergetest_perm_sub = mergetest_perm_sub[, mean(fc), by = .(cell_type, cell_ID, nr_markers, round)]
      if (is.integer(mergetest_perm_sub$cell_ID) && is.character(cellColMeanSd$cell_ID)){
        mergetest_perm_sub$cell_ID = as.character(mergetest_perm_sub$cell_ID)
      }
      mergetest_perm_sub = data.table::merge.data.table(mergetest_perm_sub, cellColMeanSd, by = 'cell_ID')
      mergetest_perm_sub[, zscore := ((V1 - colmean)* nr_markers^(1/2)) / colSd]
      
      res_list[[group_i]] = mergetest_perm_sub
      
    }
    
    res_list_comb = do.call('rbind', res_list)
    res_list_comb_average = res_list_comb[, .(mean_zscore = mean(zscore), sd_zscore = stats::sd(zscore)), by = c('cell_ID', 'cell_type')]
    mergetest_final = data.table::merge.data.table(mergetest, res_list_comb_average, by = c('cell_ID', 'cell_type'))
    
    ## calculate p.values based on normal distribution
    if(include_depletion == TRUE) {
      mergetest_final[, pval := stats::pnorm(abs(zscore), mean = mean_zscore, sd = sd_zscore, lower.tail = FALSE, log.p = FALSE)]
    } else {
      mergetest_final[, pval := stats::pnorm(zscore, mean = mean_zscore, sd = sd_zscore, lower.tail = FALSE, log.p = FALSE)]
    }
    
    data.table::setorder(mergetest_final, pval)
    
    ## calculate pval_score
    if(include_depletion == TRUE) {
      mergetest_final[, pval_score := sign(zscore)*-log10(pval)]
    } else {
      mergetest_final[, pval_score := -log10(pval)]
    }
    
    
    resultmatrix = data.table::dcast(mergetest_final, formula = cell_ID~cell_type, value.var = 'pval_score')
    return(list(DT = mergetest_final, matrix = resultmatrix))
    
    
  } else {
    
    resultmatrix = data.table::dcast(mergetest, formula = cell_ID~cell_type, value.var = 'zscore')
    return(list(DT = mergetest, matrix = resultmatrix))
    
  }
  
}


