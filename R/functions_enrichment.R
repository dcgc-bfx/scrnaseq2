#' Perform Over-Representation Analysis (ORA) for a single cluster
#'
#' This function performs Over-Representation Analysis (ORA) for a given set of 
#' genes within a cluster. It uses the clusterProfiler package to identify 
#' enriched biological processes or pathways.
#'
#' @param genes A vector of gene symbols to perform ORA on.
#' @param universe A vector of all genes in the universe of interest.
#' @param cluster The cluster identifier for which ORA is being performed. Default NA; cluster column not added
#' @param term2gene A data frame mapping gene sets (terms) to their corresponding genes.
#' 
#' @import clusterProfiler
#' 
#' @return A data frame containing the results of the ORA, including the cluster identifier.
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