#' Run over-representation analysis.
#'
#' @description Performs clusterProfiler over-representation analysis for a gene set and optional cluster label.
#'
#' @param genes Gene identifiers or symbols to analyze.
#' @param universe Background gene universe.
#' @param cluster Optional cluster label added to the output. Default is NA.
#' @param term2gene Term-to-gene mapping table.
#' @return A data frame with ORA results, optionally annotated with the cluster identifier.
#'
#' @note AI-assisted documentation.
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
