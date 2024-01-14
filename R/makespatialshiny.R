#' Make a shiny app
#'
#'
#' make a shiny app based on the spatial transcriptome
#'
#'
#' @param dat inported data for shinySRT
#' @param meta.to.include display the meta.data colnames
#' @param maxlevel maximum number of levels allowed for categorical metadata.
#'  maximum number of levels allowed for categorical metadata.
#' @param shiny.dir specify directory to create the shiny app in. Default is
#'   to create a new directory named "shinyspatial_app"
#' @param chunkSize number of genes written to h5file at any one time. Lower
#'   this number to reduce memory consumption. Should not be less than 10
#' @param gex.assay assay in single-cell data object to use for plotting
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
#' @param normalize normlize the scRNA expression data
#' @param sp_normalize normalize the ST RNA expression data
#' @param scmtx scRNA expression matrix
#' @param scmeta scRNA metadata
#' @param shiny.dir specify directory to create the shiny app in
#' @param default.gene1 specify primary default gene to show
#' @param default.gene2 specify secondary default gene to show
#' @param default.multigene character vector specifying default genes to
#'   show in bubbleplot / heatmap
#'
#'
#'
#' @return data files and codes required for shiny app
#'
#'
#'
#' @examples
#' CreateshinySRT(
#' dat,
#' meta.to.include = NA,
#' maxlevel = 50,
#' shiny.dir = 'shinyspatial_app',
#' chunkSize = 500,
#' gex.assay = NA,
#' gex.slot = c("data", "scale.data", "counts"),
#' gene.mapping = TRUE,
#' default.gene1 = NA,
#' default.gene2 = NA,
#' default.multigene = NA)
#'
#'
#'
#' @export
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
                          scmtx = NULL,
                          scmeta = NULL,
                          normalize = T,
                          sp_normalize = T,
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
    scmtx = scmtx,
    scmeta = scmeta,
    normalize = normalize,
    sp_normalize = sp_normalize,
    colcluster = colcluster,
    sp_cols = sp_cols,
    default.gene1 = default.gene1,
    default.gene2 = default.gene2,
    default.multigene = default.multigene
  )
  prepare_code(shiny.dir,title)
}
