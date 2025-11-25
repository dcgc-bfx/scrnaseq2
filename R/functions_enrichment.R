#' Perform Over-Representation Analysis (ORA) for a Single Cluster
#'
#' This function performs Over-Representation Analysis (ORA) for a given set of 
#' genes within a cluster. It uses the [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) 
#' package to identify enriched biological processes or pathways.
#'
#' ORA determines whether a set of genes (e.g., marker genes for a cluster) is 
#' significantly over-represented in a predefined set of gene sets or pathways.
#' The analysis uses the hypergeometric test with Benjamini-Hochberg correction
#' for multiple testing.
#'
#' @param genes Character vector. Gene symbols to perform ORA on. These are typically
#'   marker genes or differentially expressed genes for a specific cluster.
#' @param universe Character vector. All genes in the background/universe for the 
#'   enrichment analysis. This is typically all genes detected in the dataset.
#' @param cluster Character or numeric. The cluster identifier for which ORA is being 
#'   performed. If \code{NA} (default), the cluster column will not be added to the 
#'   results. Useful when analyzing multiple clusters separately.
#' @param term2gene A data frame with two columns mapping gene sets (terms) to their 
#'   corresponding genes. First column should be the term/pathway ID, second column 
#'   should be the gene symbol.
#'
#' @return A data frame containing the ORA results with the following columns:
#'   \itemize{
#'     \item \code{ID} – the term/pathway identifier
#'     \item \code{Description} – description of the term (if available)
#'     \item \code{GeneRatio} – ratio of input genes found in the term
#'     \item \code{BgRatio} – ratio of background genes found in the term
#'     \item \code{pvalue} – raw p-value from hypergeometric test
#'     \item \code{p.adjust} – BH-adjusted p-value
#'     \item \code{qvalue} – q-value for FDR control
#'     \item \code{geneID} – genes from input found in the term
#'     \item \code{Count} – number of input genes found in the term
#'     \item \code{Cluster} – cluster identifier (if \code{cluster} was provided)
#'   }
#'   Returns \code{NULL} if no genes are provided or no significant enrichments are found.
#'
#' @import clusterProfiler
#' @export
#'
#' @examples
#' \dontrun{
#' # Define gene sets (term to gene mapping)
#' term2gene <- data.frame(
#'   term = c("GO:0001", "GO:0001", "GO:0002", "GO:0002"),
#'   gene = c("GeneA", "GeneB", "GeneC", "GeneD")
#' )
#'
#' # Marker genes for cluster 1
#' cluster_markers <- c("GeneA", "GeneB", "GeneE")
#'
#' # Background universe
#' all_genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF")
#'
#' # Run ORA for cluster 1
#' ora_result <- perform_ora(
#'   genes = cluster_markers,
#'   universe = all_genes,
#'   cluster = 1,
#'   term2gene = term2gene
#' )
#' }
perform_ora <- function(genes, universe, cluster=NA, term2gene) {
  if (length(genes) == 0) {
    return(NULL)
  }
  
  ora_result <- clusterProfiler::enricher(gene = genes,
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.2,
                                        universe = universe,
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        TERM2GENE = term2gene)
  
  if (is.null(ora_result)) { return(NULL) }

  if (nrow(ora_result@result) > 0) {
    if (is.na(cluster)) {
      return(ora_result@result)
    } else {
      ora_result@result$Cluster <- cluster
      return(ora_result@result)
    }
  } else {
    return(NULL)
  }
}