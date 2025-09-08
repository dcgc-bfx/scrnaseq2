#' Sum of Top-N Entries per Row or Column in a Sparse Matrix
#'
#' This function computes the sum of the top `n` entries for each row (features)
#' or column (barcodes) of a sparse [`dgCMatrix`][Matrix::dgCMatrix-class] or
#' iterable (`IterableMatrix`) object. If an iterable is provided, it is processed
#' in chunks to avoid storing the entire matrix in memory.
#'
#' @param matrix A sparse [`dgCMatrix`][Matrix::dgCMatrix-class] or
#'   iterable (`IterableMatrix`) matrix.
#' @param top_n Integer vector. The number(s) of top entries to sum. Default is `50`.
#' @param margin Integer. The margin to operate over:
#'   * `1` = rows (top-N barcodes per feature)
#'   * `2` = columns (top-N features per barcode).  
#'   Default is `1`.
#' @param chunk_size Integer or `NULL`. If not `NULL`, the matrix will be processed
#'   in chunks of this size to reduce memory usage. Default is `NULL` (no chunking).
#'
#' @return A list of numeric vectors, one for each value of `top_n`. Each vector
#'   contains the top-N sums for all rows (`margin=1`) or columns (`margin=2`).
#'
#' @import Matrix
#' @import assertthat
#' @importFrom purrr map flatten_int
#' @importFrom furrr future_map furrr_options
#' @importFrom progressr progressor
#' @export
#' 
#' @examples
#' library(Matrix)
#'
#' # Create a random sparse matrix with 100 rows (features) and 50 columns (barcodes)
#' set.seed(123)
#' mat <- rsparsematrix(nrow = 100, ncol = 50, density = 0.1, rand.x = function(n) rpois(n, 5))
#'
#' # Sum top-5 entries per row (features)
#' result_rows <- SumTopN(mat, top_n = 5, margin = 1)
#' head(result_rows[[1]])
#'
#' # Sum top-3 entries per column (barcodes)
#' result_cols <- SumTopN(mat, top_n = 3, margin = 2)
#' head(result_cols[[1]])
#'
#' # Multiple top_n values
#' result_multi <- SumTopN(mat, top_n = c(3, 10), margin = 1)
#' str(result_multi)
SumTopN = function(matrix, top_n=50, margin=1, chunk_size=NULL){
    # Checks
    assertthat::assert_that(margin %in% c("1", "2"),
                            msg="Margin can only be 1 - rows or 2 - columns.")
    
    # Split counts into chunks for processing (if requested)
    chunks = NULL
    if (!is.null(chunk_size)) {
        if (margin == 1) {
            indices = 1:nrow(matrix)
        } else {
            indices = 1:ncol(matrix)
        }
        
        if (chunk_size < length(indices)) {
            chunks = split(indices, ceiling(seq_along(indices)/chunk_size))
            chunks = purrr::map(chunks, function(c) {
                if (margin == 1) {
                    mt = matrix[c, ]
                } else {
                    mt = matrix[, c]
                }
                return(mt)
            })
        }
    }
    
    if (!is.null(chunks)) {
        # Calculate for each chunk the top n sum (n can be one or more values)
        msg = paste("Sum up top", n, ifelse(margin==1, "barcodes", "features"))
        progr = progressr::progressor(along=chunks, message=msg)
        
        top_n_counts = furrr::future_map(chunks, function(counts) {
            progr()
            
            if (margin == 2) {
                # Per barcode
                if (!is(counts, "dgCMatrix")) counts = as(counts, "dgCMatrix")
                totals = Matrix::colSums(counts)
                
                # See here: https://stackoverflow.com/questions/59545714/r-matrix-package-meaning-of-the-attributes-in-the-dgcmatrix-class-for-sparse-ma
                d = diff(counts@p) 
                col_lst = split(as.integer(counts@x)*(-1), rep.int(1:ncol(counts), d))  ## columns to list
                col_lst = lapply(col_lst, function(x) return(sort(x)*(-1)))
                
                top_n_cts = lapply(top_n, function(n) {
                    s = sapply(col_lst, function(x) return(sum(head(x, n))))
                    names(s) = names(totals[totals>0])
                    
                    if (length(s) == 0) s = c()
                    
                    if (sum(totals==0) > 0) {
                        s[names(totals[totals==0])] = 0
                        s = s[names(totals)]
                    }
                    
                    return(s)
                })
            } else {
                # Per feature
                if (!is(counts, "dgCMatrix")) counts = as(counts, "dgCMatrix")
                totals = Matrix::rowSums(counts)
                
                row_lst = split(as.er(counts@x)*(-1), counts@i) ## rows to list
                row_lst = lapply(row_lst, function(x) return(sort(x)*(-1)))
                
                top_n_cts = lapply(top_n, function(n) {
                    s = sapply(row_lst, function(x) return(sum(head(x, n))))
                    names(s) = names(totals[totals>0])
                    
                    if (length(s) == 0) s = c()
                    
                    if (sum(totals==0) > 0) {
                        s[names(totals[totals==0])] = 0
                        s = s[names(totals)]
                    }
                    
                    return(s)
                })
            }
            return(top_n_cts)
        }, .options = furrr::furrr_options(seed=getOption("random_seed"), globals=c("margin", "top_n")))
        progr(type='finish')
        
        # Now combine chunk results: each chunk has values for the top n sum where n can have multiple values
        top_n_counts = purrr::map(seq(top_n), function(i) {
            # Get top sums for this n value 
            top_n_cts = purrr::map(top_n_counts, i) %>% purrr::flatten_int()
            return(top_n_cts)
        })
    } else {
        # Convert to sparse matrix
        if (!is(matrix, "dgCMatrix")) matrix = as(matrix, "dgCMatrix")
        
        # Per barcode
        if (margin == 2) {
            # Per barcode
            totals = Matrix::colSums(matrix)
            
            d = diff(matrix@p) 
            col_lst = split(as.integer(matrix@x)*(-1), rep.int(1:ncol(matrix), d))  ## columns to list
            col_lst = lapply(col_lst, function(x) return(sort(x)*(-1)))
            
            top_n_counts = lapply(top_n, function(n) {
                s = sapply(col_lst, function(x) return(sum(head(x, n))))
                names(s) = names(totals[totals>0])
                
                if (length(s) == 0) s = c()
                
                if (sum(totals==0) > 0) {
                    s[names(totals[totals==0])] = 0
                    s = s[names(totals)]
                }
                
                return(s)
            })
        } else {
            # Per feature
            totals = Matrix::rowSums(matrix)
            
            row_lst = split(as.integer(matrix@x)*(-1), matrix@i) ## rows to list
            row_lst = lapply(row_lst, function(x) return(sort(x)*(-1)))
            
            top_n_counts = lapply(top_n, function(n) {
                s = sapply(row_lst, function(x) return(sum(head(x, n))))
                names(s) = names(totals[totals>0])
                
                if (length(s) == 0) s = c()
                
                if (sum(totals==0) > 0) {
                    s[names(totals[totals==0])] = 0
                    s = s[names(totals)]
                }
                
                return(s)
            })
        }
    }
    
    return(top_n_counts)
}

#' Column-wise Percentages of a Sparse Matrix
#'
#' This function calculates the percentage contribution of each entry
#' relative to the column sum in a sparse [`dgCMatrix`][Matrix::dgCMatrix-class].
#' Columns with zero total counts are set to denominator `1` to avoid division by zero.
#'
#' @param mat A sparse [`dgCMatrix`][Matrix::dgCMatrix-class].
#'
#' @return A sparse matrix of the same dimension as `mat`, where each entry
#'   represents the percentage (0–100) of the column total.
#'
#' @import Matrix
#' @export
#' 
#' @examples
#' library(Matrix)
#'
#' # Create a random sparse matrix with 10 rows (features) and 5 columns (barcodes)
#' set.seed(123)
#' mat <- rsparsematrix(nrow = 10, ncol = 5, density = 0.3, rand.x = function(n) rpois(n, 5))
#'
#' # Calculate column percentages
#' perc_mat <- CalculateColumnPerc(mat)
#'
#' # Inspect the first few rows
#' perc_mat[1:5, ]
CalculateColumnPerc = function(mat) {
  totals = colSums(mat)
  totals = ifelse(totals>0, totals, 1)
  return(t(t(mat) / totals) * 100)
}

#' Calculate Row or Column Medians of a Sparse Matrix
#'
#' This function computes the median of each row (features) or column (barcodes)
#' of a sparse [`dgCMatrix`][Matrix::dgCMatrix-class] or iterable (`IterableMatrix`) matrix.
#' For large iterable matrices, the computation can be done in chunks to reduce memory usage.
#' A user-defined function can optionally be applied to each chunk before calculating the median.
#'
#' @param matrix A sparse [`dgCMatrix`][Matrix::dgCMatrix-class] or iterable (`IterableMatrix`) matrix.
#' @param margin Integer. The margin to calculate the median over:
#'   * `1` = rows (median per feature)
#'   * `2` = columns (median per barcode).  
#'   Default is `1`.
#' @param chunk_size Integer or `NULL`. If not `NULL`, the matrix will be 
#'    processed in chunks of this size to reduce memory usage. Default is `NULL`.
#' @param fun Optional function to apply to each chunk before calculating the 
#'    median. The function should take a single matrix argument. Default is `NULL`.
#'
#' @return A numeric vector of medians, named by row or column names depending on `margin`.
#'
#' @import Matrix
#' @import assertthat
#' @importFrom purrr map flatten_dbl
#' @importFrom furrr future_map furrr_options
#' @importFrom progressr progressor
#' @importFrom sparseMatrixStats rowMedians colMedians
#' @export
#' 
#' @examples
#' library(Matrix)
#'
#' # Create a random sparse matrix with 20 rows (features) and 10 columns (barcodes)
#' set.seed(123)
#' mat <- rsparsematrix(nrow = 20, ncol = 10, density = 0.2, rand.x = function(n) rpois(n, 5))
#'
#' # Calculate medians per row (features)
#' medians_row <- CalculateMedians(mat, margin = 1)
#' head(medians_row)
#'
#' # Calculate medians per column (barcodes)
#' medians_col <- CalculateMedians(mat, margin = 2)
#' head(medians_col)
#'
#' # Apply a log transformation before calculating medians
#' medians_log <- CalculateMedians(mat, margin = 1, fun = function(x) log1p(x))
#' head(medians_log)
CalculateMedians = function(matrix, margin=1, chunk_size=NULL, fun=NULL){
  # Checks
  assertthat::assert_that(margin %in% c("1", "2"),
                          msg="Margin can only be 1 - rows or 2 - columns.")
  
  # Split data into chunks of size chunk_size for parallel processing
  chunks = NULL
  if (!is.null(chunk_size)) {
    if (margin == 1) {
      indices = 1:nrow(matrix)
    } else {
      indices = 1:ncol(matrix)
    }
    
    if (chunk_size < length(indices)) {
      chunks = split(indices, ceiling(seq_along(indices)/chunk_size))
      chunks = purrr::map(chunks, function(c) {
        if (margin == 1) {
          mt = matrix[c, ]
        } else {
          mt = matrix[, c]
        }
        return(mt)
      })
    }
  }
  
  if (!is.null(chunks)) {
    # Analyse chunks in parallel
    # The number of workers can be configured with future::plan (at top level)
    msg = paste("Calculate medians per ", ifelse(margin==1, "barcodes", "features"))
    progr = progressr::progressor(along=chunks, message=msg)
    medians = furrr::future_map(chunks, function(counts) {
      progr()
      if (margin == 1) {
        # Analyse columns (per barcode)
        
        # If data is still an on-disk IterableMatrix, convert to sparse matrix
        if (!is(counts, "dgCMatrix")) counts = as(counts, "dgCMatrix")
        
        # Apply user-defined function stored in variable 'fun'
        if (!is.null(fun)) counts = fun(counts)
        
        # Calculate medians with function for sparse matrices
        mds = sparseMatrixStats::rowMedians(counts)
      } else {
        # Analyse rows (per feature)
        
        # If data is still an on-disk IterableMatrix, convert to sparse matrix
        if (!is(counts, "dgCMatrix")) counts = as(counts, "dgCMatrix")
        
        # Apply user-defined function stored in variable 'fun'
        if (!is.null(fun)) counts = fun(counts)
        
        # 
        mds = sparseMatrixStats::colMedians(counts)
      }
      return(mds)
    }, .options = furrr::furrr_options(seed=getOption("random_seed"), globals=c("margin", "fun"))) %>% 
      purrr::flatten_dbl()
    progr(type='finish')
  } else {
    # If data is still an on-disk IterableMatrix, convert to sparse matrix
    if (!is(matrix, "dgCMatrix")) matrix = as(matrix, "dgCMatrix")
    
    # Apply user-defined function stored in variable 'fun'
    if (!is.null(fun)) matrix = fun(matrix)
    
    # Calculate medians
    if (margin == 1) {
      medians = sparseMatrixStats::rowMedians(matrix)
    } else {
      medians = sparseMatrixStats::colMedians(matrix)
    }
  }
  return(medians)
}

#' Calculate Boxplot Statistics for Rows or Columns of a Sparse Matrix
#'
#' This function computes the boxplot statistics for each row (features) or column (barcodes)
#' of a sparse [`dgCMatrix`][Matrix::dgCMatrix-class] or iterable (`IterableMatrix`) matrix. 
#' Statistics include minimum, 25th percentile (q25), median (q50), 75th percentile (q75),
#' maximum, interquartile range (IQR), and lower/upper whiskers (q50 ± 1.5*IQR, bounded by min/max).
#'
#' @param matrix A sparse [`dgCMatrix`][Matrix::dgCMatrix-class] or iterable (`IterableMatrix`) matrix.
#' @param margin Integer. The margin to calculate statistics over:
#'   * `1` = rows (statistics per feature)
#'   * `2` = columns (statistics per barcode).  
#'   Default is `1`.
#' @param chunk_size Integer or `NULL`. If not `NULL`, the matrix will be processed 
#'    in chunks of this size to reduce memory usage. Default is `NULL`.
#'
#' @return A data frame with columns:
#'   * `min`, `q25`, `q50`, `q75`, `max` – standard boxplot statistics
#'   * `IQR` – interquartile range (`q75 - q25`)
#'   * `lower_whisker` – `q50 - 1.5*IQR` bounded by `min`
#'   * `upper_whisker` – `q50 + 1.5*IQR` bounded by `max`
#'
#' @import Matrix
#' @import assertthat
#' @importFrom purrr map pmap_dbl
#' @importFrom furrr future_map_dfr furrr_options
#' @importFrom progressr progressor
#' @importFrom sparseMatrixStats rowQuantiles colQuantiles
#' @export
#'
#' @examples
#' library(Matrix)
#'
#' # Create a random sparse matrix with 20 rows (features) and 10 columns (barcodes)
#' set.seed(123)
#' mat <- rsparsematrix(nrow = 20, ncol = 10, density = 0.3, rand.x = function(n) rpois(n, 5))
#'
#' # Calculate boxplot statistics per row (features)
#' stats_row <- CalculateBoxplotStats(mat, margin = 1)
#' head(stats_row)
#'
#' # Calculate boxplot statistics per column (barcodes)
#' stats_col <- CalculateBoxplotStats(mat, margin = 2)
#' head(stats_col)
CalculateBoxplotStats = function(matrix, margin=1, chunk_size=NULL){
  # Checks
  assertthat::assert_that(margin %in% c("1", "2"),
                          msg="Margin can only be 1 - rows or 2 - columns.")
  
  # Define chunks
  chunks = NULL
  if (!is.null(chunk_size)) {
    if (margin == 1) {
      indices = 1:nrow(matrix)
    } else {
      indices = 1:ncol(matrix)
    }
    
    if (chunk_size < length(indices)) {
      chunks = split(indices, ceiling(seq_along(indices)/chunk_size))
      chunks = purrr::map(chunks, function(c) {
        if (margin == 1) {
          mt = matrix[c, ]
        } else {
          mt = matrix[, c]
        }
        return(mt)
      })
    }
  }
  
  if (!is.null(chunks)) {
    # Analyse chunks
    msg = paste("Calculate boxplot stats per ", ifelse(margin==1, "barcodes", "features"))
    progr = progressr::progressor(along=chunks, message=msg)
    boxplot_stats = furrr::future_map_dfr(chunks, function(counts) {
      progr()
      if (!is(counts, "dgCMatrix")) counts = as(counts, "dgCMatrix")
      
      if (margin == 1) {
        mds = sparseMatrixStats::rowQuantiles(mt)
      } else {
        mds = sparseMatrixStats::colQuantiles(mt)
      }
      return(as.data.frame(mds))
    }, .options = furrr::furrr_options(seed=getOption("random_seed"), globals=c("margin")))
    progr(type='finish')
  } else {
    # Convert to sparse matrix
    if (!is(matrix, "dgCMatrix")) matrix = as(matrix, "dgCMatrix")
    
    # Calculate medians
    if (margin == 1) {
      boxplot_stats = sparseMatrixStats::rowQuantiles(matrix)
    } else {
      boxplot_stats = sparseMatrixStats::colQuantiles(matrix)
    }
    boxplot_stats = as.data.frame(boxplot_stats)
  }
  
  colnames(boxplot_stats) = c("min", "q25", "q50", "q75", "max")
  iqr = boxplot_stats$q75 - boxplot_stats$q25
  boxplot_stats$IQR = ifelse(is.na(iqr), 0, iqr)
  
  boxplot_stats$lower_whisker = purrr::pmap_dbl(boxplot_stats, function(min, q50, IQR, ...) {
    if (IQR>0) {
      lower_whisker = q50 - 1.5*IQR
      if (lower_whisker < min) lower_whisker = min
    } else {
      lower_whisker = min
    }
    return(lower_whisker)
  })
  
  boxplot_stats$upper_whisker = purrr::pmap_dbl(boxplot_stats, function(q50, max, IQR, ...) {
    if (IQR > 0){
      upper_whisker = q50 + 1.5*IQR
      if (upper_whisker > max) upper_whisker = max
    } else {
      upper_whisker = max
    }
    return(upper_whisker)
  })
  
  return(boxplot_stats)
}

#' Calculate Module Scores per Cell Using UCell
#'
#' This function computes module scores for sets of features (gene sets) per cell
#' using the [UCell](https://github.com/carmonalab/UCell) scoring method.
#' Supports sparse [`dgCMatrix`][Matrix::dgCMatrix-class] or iterable (`IterableMatrix`) matrices.
#' Large matrices can be processed in chunks to reduce memory usage.
#'
#' @param matrix A sparse [`dgCMatrix`][Matrix::dgCMatrix-class] or iterable (`IterableMatrix`) matrix.
#' @param features A named list of feature vectors (gene sets) to score.
#' @param chunk_size Integer or `NULL`. If not `NULL`, the matrix will be processed 
#'    in chunks of this size. Default is `NULL`.
#'
#' @return A data frame of module scores, with rows corresponding to cells (columns of `matrix`)
#'   and columns corresponding to the named feature sets.
#'
#' @import Matrix
#' @import assertthat
#' @importFrom purrr map flatten
#' @importFrom furrr future_map_dfr furrr_options
#' @importFrom progressr progressor
#' @import BiocParallel
#' @import UCell
#' @export
#' 
#' @examples
#' library(Matrix)
#' library(UCell)
#'
#' # Create a random sparse matrix with 20 features and 10 cells
#' set.seed(123)
#' mat <- rsparsematrix(nrow = 20, ncol = 10, density = 0.3, rand.x = function(n) rpois(n, 5))
#'
#' # Define feature sets
#' features_list <- list(
#'   set1 = rownames(mat)[1:5],
#'   set2 = rownames(mat)[6:10]
#' )
#'
#' # Calculate UCell scores
#' scores <- CalculateModuleScoreUCell(mat, features_list)
#' head(scores)
CalculateModuleScoreUCell = function(matrix, features, chunk_size=NULL){
    #matrix = Seurat::GetAssayData(sc, layer="counts", assay=assay)
    #features = known_markers_list
    #chunk_size = 5000
    
    
    # Checks
    features_vect = purrr::flatten(features) %>% unlist() %>% unique()
    valid = features_vect %in% rownames(matrix)
    assertthat::assert_that(all(valid),
                            msg=FormatString("The following features were not found in the data: {features_vect[!valid]*}."))
    
    # Define chunks
    chunks = NULL
    if (!is.null(chunk_size)) {
        indices = 1:ncol(matrix)
        
        if (chunk_size < length(indices)) {
            chunks = split(indices, ceiling(seq_along(indices)/chunk_size))
            chunks = purrr::map(chunks, function(c) {
                mt = matrix[, c]
                return(mt)
            })
        }
    }
    
    if (!is.null(chunks)) {
        # Analyse chunks
        msg = paste("Calculate UCell scores for feature sets")
        progr = progressr::progressor(along=chunks, message=msg)
        ucell_scores = furrr::future_map_dfr(chunks, function(counts) {
            progr()
            if (!is(counts, "dgCMatrix")) counts = as(counts, "dgCMatrix")
            
            # Calculate scores with UCell
            ucell = UCell:::calculate_Uscore(counts,
                                    features=features,
                                    maxRank=1500,
                                    chunk.size=ncol(counts),
                                    BPPARAM=BiocParallel::SerialParam(),
                                    ncores=1,
                                    w_neg=1,
                                    ties.method="average",
                                    storeRanks=FALSE,
                                    force.gc=FALSE,
                                    name = "")
            ucell = as.data.frame(ucell[[1]])
            rownames(ucell) = colnames(counts)
            return(ucell)
        }, .options = furrr::furrr_options(seed=getOption("random_seed"), globals=c("features")))
        progr(type='finish')
    } else {
        # Convert to sparse matrix
        if (!is(matrix, "dgCMatrix")) matrix = as(matrix, "dgCMatrix")
        
        # Calculate scores with UCell
        ucell_scores = UCell:::calculate_Uscore(matrix,
                                         features=features,
                                         maxRank=1500,
                                         chunk.size=ncol(matrix),
                                         BPPARAM=BiocParallel::SerialParam(),
                                         ncores=1,
                                         w_neg=1,
                                         ties.method="average",
                                         storeRanks=FALSE,
                                         force.gc=FALSE,
                                         name = "")
        # check and use correct name of UCell output
        if ("cells_AUC" %in% names(ucell_scores[[1]])) {
          ucell_scores = as.data.frame(ucell_scores[[1]][["cells_AUC"]])
        } else if ("cells_U" %in% names(ucell_scores[[1]])) {
          ucell_scores = as.data.frame(ucell_scores[[1]][["cells_U"]])
        } else {
          stop("The name of the UCell:::calculate_Uscore function (in functions_analysis.R/CalculateModuleScoreUCell) has changed - please check!")
        }
        rownames(ucell_scores) = colnames(matrix)
        
    }
    
    return(ucell_scores)
}



#' Calculate Cell Cycle Scores for Seurat v5 Objects
#'
#' This function computes cell cycle scores for each cell in a Seurat v5 object
#' based on provided S-phase and G2M-phase gene sets. The results include
#' the cell cycle phase as well as S.Score, G2M.Score, and CC.Difference 
#' (S.Score - G2M.Score) and are added to the cell metadata of the Seurat object.
#'
#' @param sc A [Seurat](https://satijalab.org/seurat/) v5 object.
#' @param genes_s Character vector of gene names characteristic for S-phase.
#' @param genes_g2m Character vector of gene names characteristic for G2M-phase.
#' @param assay Character string specifying the assay to use. Default is the default 
#'    assay of the Seurat object.
#' @param verbose Logical. Whether to print messages during scoring. Default is `TRUE`.
#'
#' @return Updated Seurat object with cell cycle scores added to the metadata. 
#'    The following columns are added:
#'   * `Phase` – factor with levels `G1`, `G2M`, `S`
#'   * `S.Score` – numeric S-phase score
#'   * `G2M.Score` – numeric G2M-phase score
#'   * `CC.Difference` – numeric difference `S.Score - G2M.Score`
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom purrr map
#' @importFrom furrr future_map_dfr furrr_options
#' @export
#' 
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(furrr)
#'
#' # Example Seurat object with RNA assay
#' sc <- CreateSeuratObject(matrix(rpois(200, lambda=5), nrow=20, ncol=10))
#' genes_s <- rownames(sc)[1:10]
#' genes_g2m <- rownames(sc)[11:20]
#'
#' # Calculate cell cycle scores
#' sc <- CCScoring(sc, genes_s, genes_g2m)
#' head(sc[[]])
#' }
CCScoring = function(sc, genes_s, genes_g2m, assay=NULL, verbose=TRUE){
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  if (length(genes_s) >= 20 & length(genes_g2m) >= 20) {
    # For each layer (dataset)
    # since CellCycleScoring and AddModuleScore still cannot work with layers
    # In this case, we just keep the data layers
    layers = SeuratObject::Layers(sc, assay=assay, search="^data\\.")
    sc_split = purrr::map(layers, function(l) {
      # do not calculate nCount and nFeature
      op = options(Seurat.object.assay.calcn = FALSE)
      on.exit(expr = options(op), add = TRUE)
      
      data = SeuratObject::LayerData(sc, assay=assay, layer=l)
      s = CreateAssay5Object(data=data)
      s = CreateSeuratObject(s, assay=assay)
      return(s)
    })
    
    # Calculate cell cycle scores
    cell_cycle_scores = furrr::future_map_dfr(sc_split, function(s) {
        # Check that the genes exist
        genes_s_exists = genes_s %in% rownames(s[[assay]])
        genes_g2m_exists = genes_g2m %in% rownames(s[[assay]])
        
        
        if (sum(genes_s_exists) >= 20 & sum(genes_g2m_exists) >= 20) {
            s = Seurat::CellCycleScoring(s,
                                         s.features=genes_s[genes_s_exists],
                                         g2m.features=genes_g2m[genes_g2m_exists],
                                         assay=assay,
                                         verbose=verbose)
            cc_scores = s[[c("Phase", "S.Score", "G2M.Score")]]
            cc_scores[["CC.Difference"]] = cc_scores[["S.Score"]] - cc_scores[["G2M.Score"]]
        } else {
            barcodes = Cells(s)
            cc_scores = data.frame(Phase=rep(NA, length(barcodes)) %>% as.character(), 
                                   S.Score=rep(0, length(barcodes)) %>% as.numeric(), 
                                   G2M.Score=rep(0, length(barcodes)) %>% as.numeric(), 
                                   CC.Difference=rep(0, length(barcodes)) %>% as.numeric(),
                                   row.names=barcodes)
        }
        return(cc_scores)
    }, .options = furrr::furrr_options(seed=getOption("random_seed")))
  } else {
      barcodes = Cells(sc)
      cell_cycle_scores = data.frame(Phase=rep(NA, length(barcodes)) %>% as.character(), 
                                     S.Score=rep(0, length(barcodes)) %>% as.numeric(), 
                                     G2M.Score=rep(0, length(barcodes)) %>% as.numeric(), 
                                     CC.Difference=rep(0, length(barcodes)) %>% as.numeric(),
                                     row.names=barcodes)
  }
  
  # Add Phase factor levels
  cell_cycle_scores[["Phase"]] = factor(cell_cycle_scores[["Phase"]], levels=c("G1", "G2M", "S"))

  # Add to barcode metadata
  sc = Seurat::AddMetaData(sc, cell_cycle_scores)
  
  # Log command
  sc = Seurat::LogSeuratCommand(sc)
  
  return(sc)
}

#' Transform Data and Save as New Layers in a Seurat v5 Object
#'
#' This function applies either an identity or log transformation to the data
#' in one or more layers of a Seurat v5 object and saves the result as new layers.
#' By default, raw counts layers are transformed and stored as `data.X`, `data.Y`, etc.
#'
#' @param sc A [Seurat](https://satijalab.org/seurat/) v5 object.
#' @param assay Character string specifying the assay to transform. Default is 
#'    the default assay of the Seurat object.
#' @param layer Character vector specifying which layers to transform. Default is `"counts"`.
#' @param save Character vector specifying the names of the new layers. Default 
#'    is `"data"`. If multiple layers are transformed, names are made unique automatically.
#' @param log Logical. If `TRUE`, apply `log1p` transformation. If `FALSE`, the 
#'    identity transformation is applied. Default is `FALSE`.
#'
#' @return Updated Seurat v5 object with new layers added containing the transformed data.
#' 
#' @import Seurat
#' @import SeuratObject
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Create a small Seurat object
#' sc <- CreateSeuratObject(matrix(rpois(100, lambda=5), nrow=10, ncol=10))
#'
#' # Apply log transformation to counts and save as new layers
#' sc <- TransformData(sc, layer="counts", save="data", log=TRUE)
#' }
TransformData = function(sc, assay=NULL, layer="counts", save="data", log=FALSE) {
    if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
    
    # Iterate over layers (datasets) and get size factors
    olayers = layers = unique(layer)
    layers = SeuratObject::Layers(sc[[assay]], layer)
    if (length(save) != length(layers)) {
        save = make.unique(names=gsub(pattern=olayers, replacement=save, x=layers))
    }
    
    for (i in seq_along(layers)) {
        l = layers[i]
        
        # Get counts
        counts = SeuratObject::LayerData(sc[[assay]], layer=l, fast=NA)
        
        # If requested apply log transformation, else it is just identity
        if (log) counts = log1p(counts)
        
        LayerData(sc[[assay]], 
                  layer=save[i], 
                  features=SeuratObject::Features(sc[[assay]], layer=l),
                  cells=SeuratObject::Cells(sc[[assay]], layer=l)) = counts
        
    }
    
    # Log command
    sc = Seurat::LogSeuratCommand(sc)
    
    return(sc)
}

#' Apply Scran Normalization to Counts Data in a Seurat v5 Object
#'
#' This function normalizes counts data using pooled size factors computed by
#' the [scran](https://bioconductor.org/packages/release/bioc/html/scran.html) package.
#' Large counts matrices are automatically split into chunks to reduce memory usage.
#' Normalized data are log-transformed and stored as new layers.
#'
#' @param sc A [Seurat](https://satijalab.org/seurat/) v5 object.
#' @param assay Character string specifying the assay to normalize. Default is 
#'    the default assay of the Seurat object.
#' @param layer Character vector specifying which layers to normalize. Default is `"counts"`.
#' @param save Character vector specifying names of the new layers. Default is 
#'    `"data"`. If multiple layers are transformed, names are made unique automatically.
#' @param chunk_size Integer. Maximum number of cells per chunk when computing 
#'    size factors. Large matrices are split into chunks to save memory. Default is 50,000.
#'
#' @return Updated Seurat v5 object with new normalized layers (`log1p` transformed) added.
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom purrr map flatten
#' @importFrom furrr future_map furrr_options
#' @importFrom progressr progressor
#' @import SingleCellExperiment
#' @import scran
#' @import scuttle
#' @import Matrix
#' @export
#' 
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scran)
#' library(scuttle)
#'
#' # Create a small Seurat object
#' sc <- CreateSeuratObject(matrix(rpois(100, lambda=5), nrow=10, ncol=10))
#'
#' # Apply scran normalization to counts layer
#' sc <- NormalizeDataScran(sc, layer="counts", save="data")
#' }
NormalizeDataScran = function(sc, assay=NULL, layer="counts", save="data", chunk_size=50000) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  # Iterate over layers (datasets) and get size factors
  olayers = layers = unique(layer)
  layers = SeuratObject::Layers(sc[[assay]], layer)
  if (length(save) != length(layers)) {
    save = make.unique(names=gsub(pattern=olayers, replacement=save, x=layers))
  }
  
  for (i in seq_along(layers)) {
    l = layers[i]
    
    # Define chunks
    chunks = NULL
    if (!is.null(chunk_size)) {
      barcodes = SeuratObject::Cells(sc[[assay]], layer=l)
      
      if (chunk_size < length(barcodes)) {
        chunks = split(barcodes, ceiling(seq_along(barcodes)/chunk_size))
        chunks = purrr::map(chunks, function(c) {
          return(SeuratObject::LayerData(sc[[assay]], layer=l, fast=NA, cells=c))
        })
      }
    }
    
    # Calculate size factors per chunk
    if (!is.null(chunks)) {
      # Analyse chunks
      msg = paste("Calculate size factors for layer", l)
      progr = progressr::progressor(along=chunks, message=msg)
      size_factors = furrr::future_map(chunks, function(counts) {
        progr()
        # Convert to dgCMatrix and create
        if (!is(data, "dgCMatrix")) counts = as(counts, "dgCMatrix")
      
        # Convert to SingleCellExperiment, cluster and compute size factors per cluster
        # Note: These are already centered
        counts = SingleCellExperiment::SingleCellExperiment(list(counts=counts))
        clusters = scran::quickCluster(counts, min.size=100)
        sf = scuttle::pooledSizeFactors(counts, clusters=clusters)
        sf = setNames(sf, colnames(counts))
        return(sf)
      }, .options = furrr::furrr_options(seed=getOption("random_seed"), globals=c()))
      size_factors = purrr::flatten(size_factors) %>% unlist()
      progr(type='finish')
    } else {
      # Convert to dgCMatrix and create
      counts = SeuratObject::LayerData(sc[[assay]], layer=l, fast=NA)
      if (!is(counts, "dgCMatrix")) counts = as(counts, "dgCMatrix")
      
      # Convert to SingleCellExperiment, cluster and compute size factors per cluster
      # Note: These are already centered
      counts = SingleCellExperiment::SingleCellExperiment(list(counts=counts))
      clusters = scran::quickCluster(counts, min.size=100)
      size_factors = scuttle::pooledSizeFactors(counts, clusters=clusters)
      size_factors = setNames(size_factors, colnames(counts))
    }
    
    # Now apply size factors and normalise
    counts = SeuratObject::LayerData(sc[[assay]], layer=l, fast=NA)
    counts = log1p(Matrix::t(Matrix::t(counts) / size_factors))
    
    LayerData(sc[[assay]], 
              layer=save[i], 
              features=SeuratObject::Features(sc[[assay]], layer=l),
              cells=SeuratObject::Cells(sc[[assay]], layer=l)) = counts
    
  }
  
  # Log command
  sc = Seurat::LogSeuratCommand(sc)
  
  return(sc)
}

#' Identify Highly Variable Features Using Scran
#'
#' This function identifies highly variable features (genes) in a Seurat v5 object
#' using the [scran](https://bioconductor.org/packages/release/bioc/html/scran.html) package.
#' Variable features are determined using mean-variance modeling. Data can be analyzed
#' either combined across all layers or separately per layer.
#'
#' @param sc A [Seurat](https://satijalab.org/seurat/) v5 object.
#' @param assay Character string specifying the assay to analyze. Default is the 
#'    default assay of the Seurat object.
#' @param nfeatures Integer. Number of highly variable features to identify. Default is 2000.
#' @param combined Logical. If `TRUE`, all layers are analyzed together (more memory-intensive). 
#'                 If `FALSE`, features are identified per layer and then combined. Default is `TRUE`.
#'
#' @return Updated Seurat v5 object with highly variable features stored in:
#'   * `VariableFeatures(sc[[assay]])` – character vector of variable features
#'   * Metadata columns for each feature containing variance statistics, ranks, and variable status.
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom purrr map reduce flatten_chr
#' @importFrom dplyr filter arrange bind_cols
#' @import scran
#' @export
#' 
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scran)
#'
#' # Create a small Seurat object
#' sc <- CreateSeuratObject(matrix(rpois(100, lambda=5), nrow=10, ncol=10))
#'
#' # Identify highly variable features using scran
#' sc <- FindVariableFeaturesScran(sc, nfeatures=5)
#' VariableFeatures(sc)
#' }
FindVariableFeaturesScran = function(sc, assay=NULL, nfeatures=2000, combined=TRUE) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  # Checks
  layers = SeuratObject::Layers(sc[[assay]], "^data")
  assertthat::assert_that(length(layers) > 0,
                          msg=FormatString("Could not find normalized data for assay {assay}."))
  
  # Get normalized data
  if (combined) {
    # Find features using all layers (dataset) together
    data = purrr::map(layers, function(l) {
      dt = SeuratObject::LayerData(sc, assay=assay, layer=l)
      
      # Convert to dgCMatrix matrix
      if (!is(dt, "dgCMatrix")) {
        dt = as(dt, "dgCMatrix")
      }
      return(dt)
    })
    
    # Record to which dataset a cell belongs
    block = purrr::map(seq_along(layers), function(i) {
      b = rep(layers[i], ncol(data[[i]]))
      return(b)
    }) %>% purrr::flatten_chr()
    block = factor(block, levels=unique(block))
    
    # Find variable features using the entire dataset
    data = purrr::reduce(data, cbind)
    hvf_info = scran::modelGeneVar(data, block=block)
    hvf = scran::getTopHVGs(hvf_info, n=nfeatures)
    
    # Combine hvf_info tables
    hvf_info = purrr::map(layers, function(l) {
      hvfi = hvf_info[["per.block"]][[l]]
      h = scran::getTopHVGs(hvfi, n=nfeatures)
      
      hvfi = as.data.frame(hvfi)
      hvfi$variable = rownames(hvfi) %in% h
      hvfi$rank = match(rownames(hvfi), h)
      colnames(hvfi) = paste0("vf_scran_", l, "_", colnames(hvfi))
      return(hvfi)
    })
    hvf_info = dplyr::bind_cols(hvf_info_lst)
  } else {
    # Find for each layer (dataset) separately
    
    # Find highly variable features per layer (dataset)
    hvf_info_lst = purrr::map(layers, function(l) {
      data = SeuratObject::LayerData(sc, assay=assay, layer=l)
      
      # Convert to dgCMatrix matrix
      if (!is(data, "dgCMatrix")) {
        data = as(data, "dgCMatrix")
      }
      
      # Find highly variable genes with scran
      hvf_info = scran::modelGeneVar(data)
      hvf = scran::getTopHVGs(hvf_info, n=nfeatures)
      
      hvf_info = as.data.frame(hvf_info)
      hvf_info$variable = rownames(hvf_info) %in% hvf
      hvf_info$rank = match(rownames(hvf_info), hvf)
      
      colnames(hvf_info) = paste0("vf_scran_", l, "_", colnames(hvf_info))
      return(hvf_info)
    })
    names(hvf_info_lst) = layers
    
    # Get variable features per layer (dataset) and rank
    hvf_lst =purrr::map(layers, function(l) {
      var_col = paste0("vf_scran_", l, "_variable")
      rank_col = paste0("vf_scran_", l, "_rank")
      hvf = hvf_info_lst[[l]] %>% 
        dplyr::filter(!!sym(var_col)) %>%
        dplyr::arrange(!!sym(rank_col)) %>%
        rownames()
      return(hvf)
    })
    
    # Find a good overall set (from SeuratObject:::VariableFeatures.StdAssay)
    hvf = SeuratObject:::.SelectFeatures(hvf_lst, 
                                         all.features=intersect(x = slot(sc[[assay]], name = "features")[, layers]), 
                                         nfeatures = nfeatures)
    
    # Combine hvf_info tables
    hvf_info = dplyr::bind_cols(hvf_info_lst)
  }
  
  # Add to feature metadata
  hvf_info$var.features = rownames(hvf_info) %in% hvf
  hvf_info$var.features.rank = match(rownames(hvf_info), hvf)
  sc[[assay]] = SeuratObject::AddMetaData(sc[[assay]], metadata=hvf_info)
  
  # Set variable features
  SeuratObject::VariableFeatures(sc[[assay]]) = hvf
  
  # Log command
  sc = Seurat::LogSeuratCommand(sc)
  
  return(sc)
}

#' Wrapper for Highly Variable Feature Selection
#'
#' This function provides a unified interface to identify highly variable
#' features (genes) in a Seurat v5 object using either:
#' * Seurat’s built-in `vst` method, or
#' * the scran-based mean–variance modeling method (via `FindVariableFeaturesScran`).
#'
#' @param sc A [Seurat](https://satijalab.org/seurat/) v5 object.
#' @param feature_selection_method Character string specifying the method to use.
#'   Options are:
#'   * `"vst"` – variance stabilizing transformation (Seurat default),
#'   * `"scran"` – mean–variance modeling with the [scran](https://bioconductor.org/packages/scran) package.
#' @param num_variable_features Integer. Number of features to identify. Default is `2000`.
#' @param assay Character string. Assay to analyze. If `NULL`, defaults to the active assay of the Seurat object.
#' @param verbose Logical. Whether to print messages during processing. Default is `TRUE`.
#'
#' @return A Seurat v5 object with variable features set for the specified assay.
#'   Results are stored in:
#'   * `VariableFeatures(sc[[assay]])` – character vector of variable features.
#'   * Feature-level metadata (if `scran` is used).
#'
#' @import Seurat
#' @importFrom assertthat assert_that
#' @export
#' 
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Create a toy Seurat object
#' sc <- CreateSeuratObject(matrix(rpois(100, lambda=5), nrow=10, ncol=10))
#'
#' # Identify HVFs using vst
#' sc <- FindVariableFeaturesWrapper(sc, feature_selection_method="vst", num_variable_features=5)
#' VariableFeatures(sc)
#'
#' # Identify HVFs using scran
#' sc <- FindVariableFeaturesWrapper(sc, feature_selection_method="scran", num_variable_features=5)
#' VariableFeatures(sc)
#' }
FindVariableFeaturesWrapper = function(sc, feature_selection_method, num_variable_features=2000, assay=NULL, verbose=TRUE) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  # Check
  valid_feature_selection_methods = c("vst", "scran")
  assertthat::assert_that(feature_selection_method %in% valid_feature_selection_methods,
                          msg=FormatString("Variable features method must be one of: {valid_feature_selection_methods*}."))
  
  # Find variable features
  if (feature_selection_method == "vst") {
    sc = Seurat::FindVariableFeatures(sc, 
                                      assay=assay, 
                                      selection.method="vst", 
                                      nfeatures=num_variable_features,
                                      verbose=verbose)
  } else if (feature_selection_method == "scran") {
    sc = FindVariableFeaturesScran(sc,
                                   assay=assay,
                                   nfeatures=num_variable_features,
                                   combined=TRUE)
  }
  
  return(sc)
}

#' Wrapper for Dimensionality Reduction
#'
#' Runs a dimensionality reduction on a Seurat v5 object.
#' Currently supports **PCA**, with results stored as a new reduction in the object.
#'
#' @param sc A [Seurat](https://satijalab.org/seurat/) v5 object.
#' @param method Character string. Dimensionality reduction method.
#'   Currently supported: `"pca"`.
#' @param name Character string. Name to assign to the reduction in the Seurat object.
#'   If `NULL`, defaults to the method name in lowercase (e.g., `"pca"`).
#' @param assay Character string. Assay to use. If `NULL`, defaults to the active assay of the Seurat object.
#' @param dim_n Integer. Number of dimensions to compute. Default is `50`.
#' @param verbose Logical. Whether to print messages during processing. Default is `TRUE`.
#'
#' @return A Seurat v5 object with a new dimensionality reduction stored in `sc@reductions`.
#'   Additionally:
#'   * The reduction is accessible via `sc[[<name>]]`.
#'   * Metadata (`title`, `method`) are stored in the `Misc` slot of the reduction.
#'   * The reduction is set as the default for downstream analyses.
#'
#' @import Seurat
#' @importFrom assertthat assert_that
#' @importFrom stringr str_to_title
#' @export
#' 
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Create a toy object
#' sc <- CreateSeuratObject(matrix(rpois(100, 5), nrow=10, ncol=10))
#'
#' # Run PCA with default settings
#' sc <- RunDimRedWrapper(sc, method="pca")
#'
#' # Access PCA results
#' Embeddings(sc[["pca"]])[1:5, 1:3]
#' }
RunDimRedWrapper = function(sc, method="pca", name=NULL, assay=NULL, dim_n=50, verbose=TRUE) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
    
  # Checks
  method = tolower(method)
  valid_methods = c("pca")
  assertthat::assert_that(method %in% valid_methods,
                          msg=FormatMessage("Method is {method} but must be one of: {valid_methods*}."))
  
  # Reduction name and key
  if (is.null(name)) {
    reduction_name = paste0(method, "_pca") %>% tolower()
  } else {
    reduction_name = name
  }
  reduction_key = gsub("[\\._]+", " ", reduction_name) %>%
    stringr::str_to_title() %>%
    gsub(" ", "", .)
  reduction_key = paste0(reduction_key, "_")

  # Run dimensionality reduction
  if (method == "pca") {
    sc = Seurat::RunPCA(sc,
                        assay=assay,
                        verbose=verbose, 
                        npcs=min(dim_n, ncol(sc)), 
                        seed.use=getOption("random_seed"),
                        reduction.name=reduction_name,
                        reduction.key=reduction_key)
    
    # Set title and method in misc slot
    SeuratObject::Misc(sc[[reduction_name]], slot="title") = paste0("PCA", " (", assay, ")")
    SeuratObject::Misc(sc[[reduction_name]], slot="method") = "PCA"
  }

  # Set default assay for dimensionality reduction
  SeuratObject::DefaultAssay(sc[[reduction_name]]) = assay
  
  # Set as active dimensionality reduction
  DefaultReduct(sc, assay=assay) = reduction_name

  return(sc)
}

#' Wrapper for Layer Integration
#'
#' Integrates multiple layers of a Seurat v5 object using a specified integration
#' method. This does **not** alter the raw data but instead converts an existing
#' dimensionality reduction into an integrated reduction suitable for downstream
#' analyses.  
#'
#' Supported integration methods include:
#' - `"CCAIntegration"`  
#' - `"RPCAIntegration"`  
#' - `"HarmonyIntegration"`  
#' - `"FastMNNIntegration"`  
#' - `"scVIIntegration"`  
#'
#' @param sc A [Seurat](https://satijalab.org/seurat/) v5 object.
#' @param integration_method Character string. Integration method to use.
#'   Must be one of: `"CCAIntegration"`, `"RPCAIntegration"`, `"HarmonyIntegration"`,
#'   `"FastMNNIntegration"`, `"scVIIntegration"`.
#' @param assay Character string. Assay to integrate. If `NULL`, defaults to the active assay.
#' @param orig_reduct Character string. Original dimensionality reduction to base
#'   the integration on. Default is `"pca"`. If `NULL`, uses the default reduction.
#' @param new_reduct Character string. Name of the new (integrated) reduction.
#'   If `NULL`, a name will be derived from the integration method.
#' @param new_reduct_suffix Character string or `NULL`. Optional suffix appended to the
#'   new reduction name (useful for distinguishing multiple runs).
#' @param additional_args Named list of additional arguments passed to the selected
#'   integration method (e.g., batch info, conda environment for scVI).
#' @param verbose Logical. Whether to print progress messages. Default is `TRUE`.
#'
#' @return A Seurat v5 object with a new integrated dimensionality reduction
#'   stored in `sc@reductions`. Specifically:
#'   * The reduction is accessible via `sc[[<new_reduct>]]`.
#'   * Metadata (`title`, `method`) are stored in the `Misc` slot of the reduction.
#'   * The reduction is set as the default for the assay.
#'   * For `"FastMNNIntegration"`, the intermediate corrected assay is dropped automatically.
#'
#' @import Seurat
#' @importFrom assertthat assert_that
#' @importFrom stringr str_to_title
#' @importFrom purrr exec
#' @export
#' 
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Toy Seurat object
#' sc <- CreateSeuratObject(matrix(rpois(100, 5), nrow=10, ncol=10))
#' sc <- NormalizeData(sc)
#' sc <- FindVariableFeatures(sc)
#' sc <- RunPCA(sc)
#'
#' # Run Harmony integration
#' sc <- IntegrateLayersWrapper(sc, integration_method="HarmonyIntegration")
#' Reductions(sc)
#' }
IntegrateLayersWrapper = function(sc, integration_method, assay=NULL, orig_reduct=NULL, new_reduct=NULL, new_reduct_suffix=NULL, additional_args=NULL, verbose=TRUE) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  if (is.null(orig_reduct)) orig_reduct = SeuratObject::DefaultDimReduc(sc)

  # Checks
  valid_integration_methods = c("CCAIntegration", "RPCAIntegration", "HarmonyIntegration", "FastMNNIntegration", "scVIIntegration")
  assertthat::assert_that(integration_method %in% valid_integration_methods,
                          msg=FormatString("Integration method method must be one of: {valid_integration_methods*}."))
  
  assertthat::assert_that(orig_reduct %in% SeuratObject::Reductions(sc),
                          msg=FormatString("Original reduction {orig_reduct} is not part of the Seurat object."))
  
  # Collect method-specific arguments that are always required (set only if they are not already set)
  integration_method_arg = integration_method
  if (integration_method == "CCAIntegration") {
    # Layers to use
    layers = SeuratObject::Layers(sc, assay=assay, search="^data\\.")
    
    # Name of new reduction
    new_reduct = paste0(assay, "_cca") %>% tolower()
    
    # Method call
    integration_method_arg = "CCAIntegration"
    
    # Reduction title
    new_reduct_title = "CCA"
  } else if (integration_method == "RPCAIntegration") {
    # Layers to use
    layers = SeuratObject::Layers(sc, assay=assay, search="^data\\.")
    
    # Name of new reduction
    new_reduct = paste0(assay, "_rpca") %>% tolower()
    
    # Method call
    integration_method_arg = "RPCAIntegration"
    
    # Reduction title
    new_reduct_title = "RPCA"
  } else if (integration_method == "HarmonyIntegration") {
    # Layers to use
    layers = SeuratObject::Layers(sc, assay=assay, search="^data\\.")
    
    # Name of new reduction
    new_reduct = paste0(assay, "_harmony") %>% tolower()
    
    # Method call
    integration_method_arg = "HarmonyIntegration"
    
    # Reduction title
    new_reduct_title = "Harmony"
  } else if (integration_method == "FastMNNIntegration") {
    # Layers to use
    layers = SeuratObject::Layers(sc, assay=assay, search="^data\\.")
    
    # Name of new reduction
    new_reduct = paste0(assay, "_mnn") %>% tolower()
    
    # Method call
    integration_method_arg = "FastMNNIntegration"
    
    # Name of batch-corrected assay
    if (!"reconstructed.assay" %in% names(additional_args)) additional_args[["reconstructed.assay"]] = paste0(assay, "mnn")
    
    # Add grouping information
    additional_args[["groups"]] = data.frame(group=SeuratObject::Idents(sc))
    
    # Reduction title
    new_reduct_title = "FastMNN"
  } else if (integration_method == "scVIIntegration") {
    # Layers to use
    layers = SeuratObject::Layers(sc, assay=assay, search="^counts\\.")
    
    # Name of new reduction
    new_reduct = paste0(assay, "_scvii") %>% tolower()
    
    # Method call
    integration_method_arg = scVIIntegration_Fixed

    # Conda environment for scVI
    if (!"conda_env" %in% names(additional_args)) additional_args[["conda_env"]] = reticulate::py_config()[["python"]]

    # Add grouping information
    additional_args[["groups"]] = data.frame(group=Idents(sc))
    
    # Reduction title
    new_reduct_title = "scVII"
  }
  
  # Normalization method argument
  if (!"normalization.method" %in% names(additional_args)) additional_args[["normalization.method"]] = ifelse(grepl(pattern="sct", x=assay), "SCT", "LogNormalize")
  
  # Additional suffix for name or new reduction
  if (!is.null(new_reduct_suffix)) {
    new_reduct = paste0(new_reduct, new_reduct_suffix)
  }
  
  # Call integration method
  sc = purrr::exec(Seurat::IntegrateLayers,
                   !!!c(list(sc,
                             method=integration_method_arg,
                             orig.reduction=orig_reduct,
                             assay=assay,
                             layers=layers,
                             new.reduction=new_reduct,
                             verbose=verbose),
                        additional_args)
                   )
  
  # Add title and method
  SeuratObject::Misc(sc[[new_reduct]], slot="title") = paste0(new_reduct_title, " (", assay,")")
  SeuratObject::Misc(sc[[new_reduct]], slot="method") = integration_method
  
  # Set as active dimensionality reduction
  DefaultReduct(sc, assay=assay) = new_reduct
  
  # Set default assay for dimensionality reduction
  SeuratObject::DefaultAssay(sc[[new_reduct]]) = assay
  
  # Fix key
  new_reduct_key = gsub("[\\._]+", " ", new_reduct) %>%
    stringr::str_to_title() %>%
    gsub(" ", "", .)
  new_reduct_key = paste0(new_reduct_key, "_")
  SeuratObject::Key(sc[[new_reduct]]) = new_reduct_key

  # Post-process
  if (integration_method == "FastMNNIntegration") {
    # Drop assay with corrected counts
    sc[[paste0(assay, "mnn")]] = NULL
  }
  
  return(sc)
}

# This function is a copy of the scVIIntegration (from the SeuratWrappers) with some bugs fixed.
# We keep the original code (e.g. <- instead of =) and only fix the bugs.

#' scVI Integration (Bug-Fixed Version)
#'
#' This function is a patched copy of `scVIIntegration()` from
#' [SeuratWrappers](https://github.com/satijalab/seurat-wrappers).  
#' The original code is preserved, except for two important bug fixes:
#' \enumerate{
#'   \item Only uses the specified `conda_env` if provided (does not force one).
#'   \item Ensures on-disk matrices (e.g. BPCells) are converted to
#'         `dgCMatrix` before creating the `AnnData` object.
#' }
#' 
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Load a Seurat object
#' sc <- CreateSeuratObject(matrix(rpois(2000, 5), nrow=100, ncol=20))
#'
#' # Run scVI integration
#' latent_list <- scVIIntegration_Fixed(
#'   object = sc,
#'   conda_env = "scvi-env",
#'   new.reduction = "scvi.latent",
#'   ndims = 20
#' )
#'
#' # Attach result back to Seurat object
#' sc[["scvi.latent"]] <- latent_list[["scvi.latent"]]
#' Reductions(sc)
#' }
#' 
#' @importFrom Seurat CreateDimReducObject
#' @importFrom SeuratObject LayerData
#' @importFrom future nbrOfWorkers
#' @importFrom Matrix t
#' @importFrom methods as
#' @importFrom reticulate import use_condaenv py_to_r r_to_py
#' @export
scVIIntegration_Fixed = function (object, groups = NULL, features = NULL, layers = "counts", 
                                  conda_env = NULL, new.reduction = "integrated.dr", ndims = 30, 
                                  nlayers = 2, gene_likelihood = "nb", max_epochs = NULL, ...) 
{
  # BUG/FIX - AP: Only use conda environment if specified
  if (!is.null(conda_env)) reticulate::use_condaenv(conda_env, required = TRUE)
  #
  
  sc <- reticulate::import("scanpy", convert = FALSE)
  scvi <- reticulate::import("scvi", convert = FALSE)
  anndata <- reticulate::import("anndata", convert = FALSE)
  scipy <- reticulate::import("scipy", convert = FALSE)
  if (is.null(max_epochs)) {
    max_epochs <- reticulate::r_to_py(max_epochs)
  }
  else {
    max_epochs <- as.integer(max_epochs)
  }
  batches <- SeuratWrappers:::.FindBatches(object, layers=layers)
  object <- JoinLayers(object = object, layers="counts")
  
  # BUG/FIX - AP: If on-disk matrices are used (with BPCells), convert to dgCMatrix first
  counts_matrix <- as(t(SeuratObject::LayerData(object, layer="counts")[features, ]), "dgCMatrix")
  adata <- sc$AnnData(X = scipy$sparse$csr_matrix(counts_matrix), obs = batches, var = object[[]][features,])
  
  # Set number of workers and batch size
  num_workers <- future::nbrOfWorkers()
  scvi$settings$dl_num_workers <- as.integer(num_workers)
  scvi$settings$num_threads <- as.integer(num_workers)
  scvi$settings$batch_size <- as.integer(512)
  
  scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")
  model <- scvi$model$SCVI(adata = adata, n_latent = as.integer(x = ndims), 
                           n_layers = as.integer(x = nlayers), gene_likelihood = gene_likelihood)
  model$train(max_epochs = max_epochs)
  latent <- model$get_latent_representation()
  latent <- as.matrix(latent)
  rownames(latent) <- reticulate::py_to_r(adata$obs$index$values)
  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))
  suppressWarnings(latent.dr <- CreateDimReducObject(embeddings = latent, 
                                                     key = new.reduction))
  output.list <- list(latent.dr)
  names(output.list) <- new.reduction
  return(output.list)
}


#' Enrichment of Cells per Sample per Cluster (Fisher's Exact Test)
#'
#' This function calculates enrichment of cells from each sample
#' (`orig.ident`) in each cluster (`seurat_clusters`) using Fisher's exact test
#' (`alternative = "greater"`). For each sample–cluster combination, a 2x2
#' contingency table is constructed and tested, returning odds ratios and p-values.
#'
#' @param sc Seurat object with `orig.ident` (sample) and `seurat_clusters` (cluster)
#'   in its metadata.
#'
#' @return A data frame with one row per sample, and columns for each cluster (`Cl_X`).
#'   Each cell contains a vector with:
#'   \itemize{
#'     \item \code{oddsRatio}: estimated odds ratio of enrichment.
#'     \item \code{p}: p-value from Fisher's exact test (formatted scientific notation).
#'   }
#' 
#' @importFrom dplyr pull filter count
#' @importFrom purrr reduce
#' @importFrom stats fisher.test
#' @export
#' 
#' @examples
#' library(Seurat)
#'
#' # Example Seurat object
#' sc <- CreateSeuratObject(matrix(rpois(2000, 5), nrow = 100, ncol = 20))
#' sc$seurat_clusters <- sample(1:3, ncol(sc), replace = TRUE)
#' sc$orig.ident <- sample(c("SampleA", "SampleB"), ncol(sc), replace = TRUE)
#'
#' # Run Fisher enrichment
#' enrichment <- CellsFisher(sc)
#' head(enrichment)
CellsFisher = function(sc) {
  cell_samples = sc[[]] %>% dplyr::pull(orig.ident) %>% unique() %>% sort()
  cell_clusters = sc[[]] %>% dplyr::pull(seurat_clusters) %>% unique() %>% sort()
  out = matrix(0+NA, nrow=0, ncol=length(cell_clusters)) %>% as.data.frame()
  for(s in cell_samples) {
    ft.list = lapply(cell_clusters, function(cl) { 
      a = sc[[]] %>% dplyr::filter(orig.ident==s, seurat_clusters==cl) %>% dplyr::count() %>% as.numeric()
      b = sc[[]] %>% dplyr::filter(orig.ident!=s, seurat_clusters==cl) %>% dplyr::count() %>% as.numeric()
      c = sc[[]] %>% dplyr::filter(orig.ident==s, seurat_clusters!=cl) %>% dplyr::count() %>% as.numeric()
      d = sc[[]] %>% dplyr::filter(orig.ident!=s, seurat_clusters!=cl) %>% dplyr::count() %>% as.numeric()
      tbl.2by2 = matrix(c(a, b, c, d), ncol=2, nrow=2, byrow=TRUE)
      ft = fisher.test(tbl.2by2, alternative="greater")
      return(c(oddsRatio=round(as.numeric(ft$estimate), 2),
               p=formatC(as.numeric(ft$p.value), format="e", digits=1)))
    })
    ft.matrix = purrr::reduce(ft.list, .f=cbind)
    rownames(ft.matrix) = paste0(s, ".", rownames(ft.matrix))
    out = rbind(out, ft.matrix)
  }
  colnames(out) = paste0("Cl_", cell_clusters)
  return(out)
}
