#' Perform Over-Representation Analysis
#'
#' Run clusterProfiler over-representation analysis for a gene set and optionally
#' annotate the resulting enrichment table with a cluster identifier.
#'
#' @param genes Character vector of gene symbols to test for enrichment.
#' @param universe Character vector of background genes used as the enrichment universe.
#' @param cluster Cluster identifier to append to the result table. Default: NA, which leaves the cluster column absent.
#' @param term2gene Data frame mapping gene-set terms to their member genes.
#' 
#' @import clusterProfiler
#' 
#' @return A data frame of enrichment results, optionally with a `Cluster` column, or NULL when no input genes or enrichments are available.
#' @note AI-assisted documentation
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