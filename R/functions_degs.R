#' Given the DEG contrasts table, prepares a list with DEG contrasts to do.
#' 
#' @param sc A Seurat single cell object.
#' @param degs_contrasts_list A list of contrasts. Must at least contain 'name', condition_column', 'condition_group1' and 'condition_group2'.
#' @return A list with contrasts to be analysed.
DegsSetupContrastsList = function(sc, degs_contrasts_list) {
    barcode_metadata = sc[[]]
    
    # If empty, return empty list
    if (length(degs_contrasts_list) == 0) return(list())
    
    # Convert into list, do checks and set up defaults
    degs_contrasts_list = purrr::map(seq(degs_contrasts_list), function(i) {
        contrast = degs_contrasts_list[[i]]
        
        # Trim whitespace
        contrast = purrr::map(contrast, trimws)
        
        #
        # name
        #
        assertthat::assert_that("name" %in% names(contrast), 
                                msg=FormatString("The name ('name') is missing (for comparison {i})."))
        name = contrast[["name"]]
        
        #
        # condition_column
        #
        assertthat::assert_that("condition_column" %in% names(contrast),
                                msg=FormatString("The condition column ('condition_column') is missing (for comparison {i}/{name})."))
        condition_column = contrast[["condition_column"]]
        
        #
        # condition_group1 and condition_group2
        #
        assertthat::assert_that("condition_group1" %in% names(contrast),
                                msg=FormatString("The condition group 1 ('condition_group1') is missing or empty (for comparison {i}/{name})."))
        
        assertthat::assert_that("condition_group1" %in% names(contrast),
                                msg=FormatString("The condition group 1 ('condition_group1') is missing or empty (for comparison {i}/{name})."))
        # If at least one of the condition groups is not a file, check if condition_column is part of the barcode metadata
        if (!file.exists(contrast[["condition_group1"]]) | !file.exists(contrast[["condition_group2"]])) {
            assertthat::assert_that(condition_column %in% colnames(barcode_metadata),
                                    msg=FormatString("The condition column ('condition_column') must be part of the barcode metadata if at least one of the condition groups is not a file (for comparison {i}/{name})."))
            
            assertthat::assert_that(is.factor(barcode_metadata[, condition_column]),
                                    msg=FormatString("The condition column ('condition_column') of the barcode metadata must be a factor (for comparison {i}/{name})."))
        }
        
        if (!file.exists(contrast[["condition_group1"]])) {
            condition_group1 = contrast[["condition_group1"]]
            
            # Parse condition_group1 string
            negate = grepl("^!", condition_group1)
            condition_group1 = gsub("^!", "", condition_group1) %>%
                strsplit(split="\\+") %>%
                unlist() %>%
                trimws() %>%
                unique()
            
            # If negate, get complement of condition_group1
            if (negate) {
                condition_group1 = setdiff(levels(barcode_metadata[, condition_column]), condition_group1)
            }
            
            # Make sure all levels are valid
            purrr::walk(condition_group1, function(c) {
                assertthat::assert_that(c %in% levels(barcode_metadata[, condition_column]),
                                        msg=FormatString("The condition group 1 value {c} is not level of the condition column {condition_column} of the barcode metadata (for comparison {i}/{name})."))
            })
            
            contrast[["condition_group1"]] = condition_group1
            
            # Get indices of condition_group1
            contrast[["condition_group1_idx"]] = which(barcode_metadata[, condition_column] %in% condition_group1)
        } else {
            condition_group1 = contrast[["condition_group1"]]
            
            # Sheet number appended?
            sheet = 1
            if (grepl(":\\d+$", condition_group1)) {
                sheet = gsub(pattern=".+:(\\d+)$", replacement="\\1", x=condition_group1)
                condition_group1 = gsub(pattern=":\\d+$", replacement="", x=condition_group1)
            }
            
            # Make sure file exists
            assertthat::assert_that(file.exists(condition_group1),
                                    msg=FormatString("The condition group 1 barcode file {condition_group1} does not exist (for comparison {i}/{name})."))
            
            # Decide whether it is a valid file type
            extension = tools::file_ext(gsub(pattern="\\.gz$", replacement="", x=condition_group1))
            valid_extensions = c("csv", "tsv", "xls", "xlsx")
            assertthat::assert_that(extension %in% valid_extensions,
                                    msg=FormatString("The condition group 1 barcode file must be one of: {valid_extensions*} (file can be gzipped) (for comparison {i}/{name})."))
            
            # Read file
            if (extension %in% c("csv", "tsv")) {
                barcodes = readr::read_tsv(condition_group1)
            } else if (extension %in% c("xls", "xlsx")) {
                barcodes = readxl::read_excel(condition_group1, sheet=sheet)
            }
            
            # Make sure file contains data
            assertthat::assert_that(is.data.frame(barcodes) && nrow(barcodes)>0,
                                    msg=FormatString("The condition group 1 barcode file does not contain barcodes (for comparison {i}/{name})."))
            barcodes = barcodes[, 1, drop=TRUE] %>% trimws() %>% unique()
            barcodes = barcodes[!is.na(barcodes)]
            
            # Parse barcodes with sample information
            idx = grepl("^[^:]+:.+", barcodes) %>% which()
            if (length(idx) > 0) {
                samples = gsub("^([^:]+):.+", "\\1", barcodes[idx])
                idx = idx[samples %in% levels(barcode_metadata[, "orig.ident"])]
                
                if (length(idx) > 0) {
                    orig_ident_orig_barcode = paste(barcode_metadata$orig.ident, barcode_metadata$orig_barcode, sep=":")
                    jdx = match(barcodes[idx], orig_ident_orig_barcode)
                    assertthat::assert_that(!any(is.na(jdx)),
                                            msg=FormatString("The condition group 1 barcode file contains barcodes that cannot be found in the barcode metadata (for comparison {i}/{name})."))
                    contrast[["condition_group1_idx"]] = jdx
                }
            }
        }
        
        if (!file.exists(contrast[["condition_group2"]])) {
            condition_group2 = contrast[["condition_group2"]]
            
            # Parse condition_group2 string
            negate = grepl("^!", condition_group2)
            condition_group2 = gsub("^!", "", condition_group2) %>%
                strsplit(split="\\+") %>%
                unlist() %>%
                trimws() %>%
                unique()
            
            # If negate, get complement of condition_group2
            if (negate) {
                condition_group2 = setdiff(levels(barcode_metadata[, condition_column]), condition_group2)
            }
            
            # Make sure all levels are valid
            purrr::walk(condition_group2, function(c) {
                assertthat::assert_that(c %in% levels(barcode_metadata[, condition_column]),
                                        msg=FormatString("The condition group 2 value {c} is not level of the condition column {condition_column} of the barcode metadata (for comparison {i}/{name})."))
            })
            
            contrast[["condition_group2"]] = condition_group2
            
            # Get indices of condition_group2
            contrast[["condition_group2_idx"]] = which(barcode_metadata[, condition_column] %in% condition_group2)
        } else {
            condition_group2 = contrast[["condition_group2"]]
            
            # Sheet number appended?
            sheet = 1
            if (grepl(":\\d+$", condition_group2)) {
                sheet = gsub(pattern=".+:(\\d+)$", replacement="\\1", x=condition_group2)
                condition_group2 = gsub(pattern=":\\d+$", replacement="", x=condition_group2)
            }
            
            # Make sure file exists
            assertthat::assert_that(file.exists(condition_group2),
                                    msg=FormatString("The condition group 2 barcode file {condition_group2} does not exist (for comparison {i}/{name})."))
            
            # Decide whether it is a valid file type
            extension = tools::file_ext(gsub(pattern="\\.gz$", replacement="", x=condition_group2))
            valid_extensions = c("csv", "tsv", "xls", "xlsx")
            assertthat::assert_that(extension %in% valid_extensions,
                                    msg=FormatString("The condition group 2 barcode file must be one of: {valid_extensions*} (file can be gzipped) (for comparison {i}/{name})."))
            
            # Read file
            if (extension %in% c("csv", "tsv")) {
                barcodes = readr::read_tsv(condition_group2)
            } else if (extension %in% c("xls", "xlsx")) {
                barcodes = readxl::read_excel(condition_group2, sheet=sheet)
            }
            
            # Make sure file contains data
            assertthat::assert_that(is.data.frame(barcodes) && nrow(barcodes)>0,
                                    msg=FormatString("The condition group 2 barcode file does not contain barcodes (for comparison {i}/{name})."))
            barcodes = barcodes[, 1, drop=TRUE] %>% trimws() %>% unique()
            barcodes = barcodes[!is.na(barcodes)]
            
            # Parse barcodes with sample information
            idx = grepl("^[^:]+:.+", barcodes) %>% which()
            if (length(idx) > 0) {
                samples = gsub("^([^:]+):.+", "\\1", barcodes[idx])
                idx = idx[samples %in% levels(barcode_metadata[, "orig.ident"])]
                
                if (length(idx) > 0) {
                    orig_ident_orig_barcode = paste(barcode_metadata$orig.ident, barcode_metadata$orig_barcode, sep=":")
                    jdx = match(barcodes[idx], orig_ident_orig_barcode)
                    assertthat::assert_that(!any(is.na(jdx)),
                                            msg=FormatString("The condition group 2 barcode file contains barcodes that cannot be found in the barcode metadata (for comparison {i}/{name})."))
                    contrast[["condition_group2_idx"]] = jdx
                }
            }
        }
        
        # subset_column
        if ("subset_column" %in% names(contrast)) {
            assertthat::assert_that("subset_group" %in% names(contrast),
                                    msg=FormatString("The 'subset_group' column must be used together with the subset_column column (for comparison {i}/{name})."))
            subset_column = contrast[["subset_column"]] %>% trimws()
            subset_group = contrast[["subset_group"]]
            
            if (!file.exists(subset_group)) {
                # Parse subset_group string
                subset_group = subset_group %>% strsplit(split=",") %>% unlist() %>% trimws() %>% unique()
                subset_group = purrr::map(subset_group, function(s) {
                    s = s %>%
                        strsplit(split="\\+") %>%
                        unlist() %>%
                        trimws() %>%
                        unique()
                    return(s)
                })
                purrr::walk(unlist(subset_group), function(c) {
                    assertthat::assert_that(c %in% levels(barcode_metadata[, subset_column]),
                                            msg=FormatString("The subset group value {c} is not level of the subset column {subset_column} of the barcode metadata (for comparison {i}/{name})."))
                })
                
                contrast[["subset_group"]] = subset_group
                
                # Get indices of subset_group
                contrast[["subset_group_idx"]] = purrr::map(subset_group, function(s) {
                    which(barcode_metadata[, subset_column] %in% s)
                })
            } else {
                # Sheet number appended?
                sheet = 1
                if (grepl(":\\d+$", subset_group)) {
                    sheet = gsub(pattern=".+:(\\d+)$", replacement="\\1", x=subset_group)
                    subset_group = gsub(pattern=":\\d+$", replacement="", x=subset_group)
                }
                
                # Make sure file exists
                assertthat::assert_that(file.exists(subset_group),
                                        msg=FormatString("The subset group barcode file {subset_group} does not exist (for comparison {i}/{name})."))
                
                # Decide whether it is a valid file type
                extension = tools::file_ext(gsub(pattern="\\.gz$", replacement="", x=subset_group))
                valid_extensions = c("csv", "tsv", "xls", "xlsx")
                assertthat::assert_that(extension %in% valid_extensions,
                                        msg=FormatString("The subset group barcode file must be one of: {valid_extensions*} (file can be gzipped) (for comparison {i}/{name})."))
                
                # Read file
                if (extension %in% c("csv", "tsv")) {
                    barcodes = readr::read_tsv(subset_group)
                } else if (extension %in% c("xls", "xlsx")) {
                    barcodes = readxl::read_excel(subset_group, sheet=sheet)
                }
                
                # Make sure file contains data
                assertthat::assert_that(is.data.frame(barcodes) && nrow(barcodes)>0,
                                        msg=FormatString("The subset group barcode file does not contain barcodes (for comparison {i}/{name})."))
                barcodes = barcodes[, 1, drop=TRUE] %>% trimws() %>% unique()
                barcodes = barcodes[!is.na(barcodes)]
                
                # Parse barcodes with sample information
                idx = grepl("^[^:]+:.+", barcodes) %>% which()
                if (length(idx) > 0) {
                    samples = gsub("^([^:]+):.+", "\\1", barcodes[idx])
                    idx = idx[samples %in% levels(barcode_metadata[, "orig.ident"])]
                    
                    if (length(idx) > 0) {
                        orig_ident_orig_barcode = paste(barcode_metadata$orig.ident, barcode_metadata$orig_barcode, sep=":")
                        jdx = match(barcodes[idx], orig_ident_orig_barcode)
                        assertthat::assert_that(!any(is.na(jdx)),
                                                msg=FormatString("The subset group barcode file contains barcodes that cannot be found in the barcode metadata (for comparison {i}/{name})."))
                        contrast[["subset_group_idx"]] = list(jdx)
                    }
                }
            }
        }
        
        # bulk_by
        if ("bulk_by" %in% names(contrast)) {
            # Parse bulk_by string and make sure that the columns are factors
            bulk_by = contrast[["bulk_by"]] %>% 
                strsplit(split="\\+") %>%
                unlist() %>%
                trimws() %>%
                unique()
            assertthat::assert_that(all(bulk_by %in% colnames(barcode_metadata)),
                                    msg=FormatString("The bulk by column(s) must be part of the barcode metadata (for comparison {i}/{name})."))
            purrr::walk(bulk_by, function(c) {
                assertthat::assert_that(is.factor(barcode_metadata[, c]),
                                        msg=FormatString("The bulk by column(s) must be factors (for comparison {i}/{name})."))
            })
            
            contrast[["bulk_by"]] = bulk_by
            
            # bulk_barcodes
            if ("bulk_barcodes" %in% names(contrast)) {
                contrast[["bulk_barcodes"]] = as.numeric(contrast[["bulk_barcodes"]])
            }
        }
        
        # pseudobulk_samples
        if ("pseudobulk_samples" %in% names(contrast)) {
            assertthat::assert_that(contrast[["pseudobulk_samples"]] > 1,
                                    msg=FormatString("The number of pseudobulk samples ('pseudobulk_samples') must be greater than 1 (for comparison {i}/{name})."))
            # pseudobulk_barcodes
            if ("pseudobulk_barcodes" %in% names(contrast)) {
                contrast[["pseudobulk_barcodes"]] = as.numeric(contrast[["pseudobulk_barcodes"]])
            }
        }
        
        assertthat::assert_that(!("bulk_by" %in% names(contrast) & "pseudobulk_samples" %in% names(contrast)),
                                msg=FormatString("Either 'bulk_by' or 'pseudobulk_samples' can be specified (for comparison {i}/{name})."))
        
        # test
        valid_tests = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
        if (!"test" %in% names(contrast)) contrast[["test"]] = "wilcox"
        assertthat::assert_that(contrast[["test"]] %in% valid_tests,
                                msg=FormatString("The test ('test') must be one of: {valid_tests*} (for comparison {i}/{name})."))
        
        # padj
        if (!"padj" %in% names(contrast)) contrast[["padj"]] = 0.05
        contrast[["padj"]] = as.numeric(contrast[["padj"]])
        
        # log2FC
        if (!"log2FC" %in% names(contrast)) contrast[["log2FC"]] = 0
        contrast[["log2FC"]] = as.numeric(contrast[["log2FC"]])
        
        # min_pct
        if (!"min_pct" %in% names(contrast)) contrast[["min_pct"]] = 0.05
        contrast[["min_pct"]] = as.numeric(contrast[["min_pct"]])
        
        # assay
        if (!"assay" %in% names(contrast)) contrast[["assay"]] = Seurat::DefaultAssay(sc)
        assay = contrast[["assay"]] %>% 
            strsplit(split="\\,") %>% 
            unlist() %>% 
            trimws()
        valid_assays = Seurat::Assays(sc)
        valid_reductions = Seurat::Reductions(sc)
        assertthat::assert_that(contrast[["assay"]] %in% valid_assays | contrast[["assay"]] %in% valid_reductions | all(assay %in% names(barcode_metadata)),
                                msg=FormatString("The assay ('assay') must be one of the assays {valid_assays*}, one of the reductions {valid_reductions*} or barcode metadata columns (for comparison {i}/{name})."))
        
        # slot
        if (contrast[["assay"]] %in% valid_assays) {
            if (!"slot" %in% names(contrast)) contrast[["slot"]] = "data"
            assertthat::assert_that(contrast[["slot"]] %in% c("counts", "data", "scale.data"),
                                    msg=FormatString("The slot ('slot') must be 'counts', 'data', 'scale.data' (for comparison {i}/{name})."))
        }
        
        # downsample_barcodes
        if ("downsample_barcodes" %in% names(contrast)) {
            contrast[["downsample_barcodes"]] = as.numeric(contrast[["downsample_barcodes"]])
        }
        
        # batch
        if ("batch" %in% names(contrast)) {
            batch = contrast[["batch"]] %>% 
                strsplit(split=",") %>% 
                unlist() %>% 
                trimws() %>% 
                unique()
            assertthat::assert_that(all(batch %in% colnames(barcode_metadata)),
                                    msg=FormatString("The batch column(s) must be part of the barcode metadata (for comparison {i}/{name})."))
            purrr::walk(batch, function(c) {
                assertthat::assert_that(is.factor(barcode_metadata[, c]),
                                        msg=FormatString("The batch column(s) must be factors (for comparison {i}/{name})."))
            })
            contrast[["batch"]] = batch 
        }

        # Add contrast row number
        contrast[["contrast_row"]] = i
        
        return(contrast)
    })
    
    
    # Expand subsets list so that there is now one entry per subset
    # Also get barcode indices per subset
    degs_contrasts_list = purrr::map(degs_contrasts_list, function(contrast) {
        if ("subset_group" %in% names(contrast)) {
            contrasts_expanded = purrr::map(seq(contrast[["subset_group"]]), function(i) {
                # Set up new contrast with just the subset group and barcodes indices
                con = contrast
                con[["subset_group"]] = contrast[["subset_group"]][[i]]
                con[["subset_group_idx"]] = contrast[["subset_group_idx"]][[i]]
                
                # Filter barcodes (condition_group1_idx and condition_group2_idx) so that they are in the subset group
                subset_group_barcodes = rownames(barcode_metadata[con[["subset_group_idx"]], ]) %>% unique()
                
                k = rownames(barcode_metadata[contrast[["condition_group1_idx"]], ]) %in% subset_group_barcodes
                con[["condition_group1_idx"]] = contrast[["condition_group1_idx"]][k]
                
                k = rownames(barcode_metadata[contrast[["condition_group2_idx"]], ]) %in% subset_group_barcodes
                con[["condition_group2_idx"]] = contrast[["condition_group2_idx"]][k]
                return(con)
            })
        } else {
            contrasts_expanded = list(contrast)
        }
        return(contrasts_expanded)
    }) %>% purrr::flatten()
    
    # Process further
    degs_contrasts_list = purrr::map(seq(degs_contrasts_list), function(i) {
        contrast = degs_contrasts_list[[i]]
    
        # Downsample barcodes if requested (downsample_barcodes, bulk_barcodes, pseudobulk_barcodes)
        if ("bulk_by" %in% names(contrast) & "bulk_barcodes" %in% names(contrast)) {
            # Group by bulk column(s) and sample barcodes per group, sample (at most) x rows per group, then subset
            set.seed(getOption("random_seed"))
            sampled_condition_group1_idx = barcode_metadata[contrast[["condition_group1_idx"]], ] %>% 
                dplyr::mutate(condition_group1_idx=contrast[["condition_group1_idx"]]) %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(contrast[["bulk_y"]]))) %>%
                dplyr::slice_sample(n=contrast[["bulk_barcodes"]]) %>%
                dplyr::pull(condition_group1_idx)
            k = contrast[["condition_group1_idx"]] %in% sampled_condition_group1_idx
            contrast[["condition_group1_idx"]] = contrast[["condition_group1_idx"]][k]
            
            set.seed(getOption("random_seed"))
            sampled_condition_group2_idx = barcode_metadata[contrast[["condition_group2_idx"]], ] %>% 
                dplyr::mutate(condition_group2_idx=contrast[["condition_group2_idx"]]) %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(contrast[["bulk_y"]]))) %>%
                dplyr::slice_sample(n=contrast[["bulk_barcodes"]]) %>%
                dplyr::pull(condition_group2_idx)
            k = contrast[["condition_group2_idx"]] %in% sampled_condition_group2_idx
            contrast[["condition_group2_idx"]] = contrast[["condition_group2_idx"]][k]
            
        } else if ("pseudobulk_samples" %in% names(contrast) & "pseudobulk_barcodes" %in% names(contrast)) {
            # Just sample barcodes: number of pseudosamples * barcodes per pseudosample
            downsample_n = contrast[["pseudobulk_samples"]]*contrast[["pseudobulk_barcodes"]]
            set.seed(getOption("random_seed"))
            contrast[["condition_group1_idx"]] = sample(x=contrast[["condition_group1_idx"]], 
                                                    size=min(downsample_n, length(contrast[["condition_group1_idx"]])))
            
            set.seed(getOption("random_seed"))
            contrast[["condition_group2_idx"]] = sample(x=contrast[["condition_group2_idx"]], 
                                                    size=min(downsample_n, length(contrast[["condition_group2_idx"]])))
        } else if ("downsample_n" %in% names(contrast)) {
            set.seed(getOption("random_seed"))
            contrast[["condition_group1_idx"]] = sample(x=contrast[["condition_group1_idx"]], 
                                                    size=min(contrast[["downsample_n"]], length(contrast[["condition_group2_idx"]])))
            
            set.seed(getOption("random_seed"))
            contrast[["condition_group2_idx"]] = sample(x=contrast[["cells_group2_idx"]], 
                                                    size=min(contrast[["downsample_n"]], length(contrast[["condition_group2_idx"]])))
        }
        
        return(contrast)
    })
    
    return(degs_contrasts_list)
}

#' Sorts table of differentially expressed genes per performed test. Introduces a signed p-value score calculated as follows:
#' p_val_adj_score = -log10(p_val_adj) * sign(avg_log2FC).
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function or the "Seurat::FindMarkers" function.
#' @param group Group results first by column(s) before sorting.
#' @return Sorted table with differentially expressed genes. 
DegsSort = function(degs, group=NULL) { 
  # Introduce signed p-value score and group (if requested)
  degs$p_val_adj_score = -log10(degs$p_val_adj) * sign(degs$avg_log2FC)
  if (!is.null(group)) degs = degs %>% dplyr::group_by(dplyr::across(dplyr::all_of(group)))
  
  # Now sort table
  degs = degs %>% dplyr::arrange(-p_val_adj_score, -avg_log2FC, .by_group=TRUE) %>% as.data.frame()
  
  return(degs)
}

#' Filter table of differentially expressed genes.
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" or the "Seurat::FindMarkers" functions.
#' @param cut_log2FC Log2 fold change threshold.
#' @param cut_padj Adjusted p-value threshold; not advised for filtering markers
#' @param split_by_dir Split filtered table into a table with all degs, a table with up-regulated degs and a table down-regulated degs.
#' @return If split_by_dir is set to FALSE filtered table else list of filtered tables with all, up-regulated and down-regulated degs.
DegsFilter = function(degs, cut_log2FC, cut_padj=NULL, split_by_dir=TRUE) { 
  
  # Filter differentially expressed genes based on fold change 
  filt = degs %>% 
    dplyr::filter(abs(avg_log2FC) >= cut_log2FC) %>% 
    as.data.frame()
  
  # Filter based on p-values 
  if (!is.null(cut_padj)) {
    filt = filt %>%
      dplyr::filter(p_val_adj <= cut_padj)  
  }
  
  # Separate up- and down-regulated genes if requested
  if (split_by_dir) {
    down = filt %>% 
      dplyr::filter(avg_log2FC <= -cut_log2FC)
    up = filt %>% 
      dplyr::filter(avg_log2FC >= cut_log2FC)
    filt = list(all=filt, up=up, down=down)
  }
  
  return(filt)
}

#' Display top marker genes (=up-regulated genes).
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function (requires column "cluster")
#' @param n Number of top genes to show
#' @param column_1 First column to sort genes on
#' @param column_2 Second column to sort genes on
#' @return Data.frame of top genes 
DegsUpDisplayTop = function(degs, n=5, column_1="p_val_adj_score", column_2="pct.diff") { 
  
  # Calculate difference in percentage of cells that express the gene
  degs = degs %>% dplyr::mutate(pct.diff=abs(pct.1-pct.2))
  
  # Get top 5 up-regulated markers
  top = degs %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::top_n(n=n, wt=get(column_1)) %>% 
    dplyr::top_n(n=n, wt=get(column_2)) %>% 
    dplyr::ungroup() %>% 
    dplyr::transmute(cluster=cluster,
                     gene=gene,
                     avg_log2FC=round(avg_log2FC, digits=3),
                     p_val=formatC(as.numeric(p_val), format="e", digits=1),
                     p_val_adj=formatC(as.numeric(p_val_adj), format="e", digits=1),
                     pct.1=pct.1,
                     pct.2=pct.2) %>% 
    as.data.frame()
  return(top)
}

#' Compute average gene expression data per identity class for a set of genes.
#' 
#' @param sc Seurat object.
#' @param genes Gene list for which average data are to be extracted.
#' @return A table with average RNA counts and data per identity class for each gene.
DegsAvgDataPerIdentity = function(sc, genes, assay="RNA") { 
  # The standard average log FC is derived from assay and layer="data"
  # Add average scaled data per cluster for default assay
  avg_set = list()
  avg_set[[assay]] = "counts"
  avg_set[[DefaultAssay(sc)]] = c(avg_set[[DefaultAssay(sc)]], "data")
  avg_data = matrix(NA+0, nrow=length(genes), ncol=0)
  genes_unique = unique(genes)
  identities = levels(Idents(sc))
  for (as in names(avg_set)) { 
    for (sl in avg_set[[as]]) {
      if (length(genes) > 0) {
        avg_per_id = mapply(function(id) { 
          id_cells = WhichCells(sc, idents=id)
          tmp_subset = LayerData(sc, assay=as, layer=sl, cells=id_cells , features=genes_unique)
          tmp_subset_matrix = as(tmp_subset, "dgCMatrix")
          if (sl=="data") {
            # This calculation is in accordance with what Seurat is doing
            id_avg = log2(Matrix::rowMeans(exp(tmp_subset_matrix)) + 1)[genes]
          } else if (sl=="counts") {
            id_avg = Matrix::rowMeans(tmp_subset_matrix)[genes]
          }
          return(id_avg)
        }, identities)
      } else {
        avg_per_id = matrix(NA, nrow=0, ncol=length(identities)) %>% as.data.frame()
      }
      colnames(avg_per_id) = paste0("avg_", as, "_", sl, "_id", identities)
      avg_data = cbind(avg_data, avg_per_id)
    }
  }
  return(avg_data)
}

#' Compute average gene expression data for a set of cells and a set of genes.
#' 
#' @param object Seurat assay object.
#' @param cells Cells to be used. NULL if all cells should be used.
#' @param genes Gene list for which average data are to be extracted. Can be NULL where all genes will be calculated.
#' @param slot Slot to be used (data). Can be 'counts' or 'data' or both (vector).
#' @return A table with average data for each gene.
DegsAvgData = function(object, cells=NULL, genes=NULL, slot="data") {
  if ("data" %in% slot) avg_data = as.numeric() else avg_data = NULL
  if ("counts" %in% slot) avg_counts = as.numeric() else avg_counts = NULL
  
  # If object is NULL or cells==0 or genes==0 return empty data.frame
  if (is.null(object) | (!is.null(cells) && length(cells)==0) | (!is.null(genes) && length(genes)==0)) {
    return(cbind(avg_counts, avg_data))
  }
  
  # If genes is NULL set defaults
  if (is.null(genes)) genes = rownames(object)
  
  # Subset by cells first (if requested)
  if (!is.null(cells)) object = subset(object, cells=cells)
  
  if ("data" %in% slot) {
    if (length(genes)>0) {
      # This calculation is in accordance with what Seurat is doing
      avg_data = log2(Matrix::rowMeans(exp(Seurat::GetAssayData(object, slot="data")[genes, , drop=FALSE])) + 1)
    }
  }
  
  if ("counts" %in% slot) {
    if (length(genes)>0) {
      avg_counts = Matrix::rowMeans(Seurat::GetAssayData(object, slot="data")[genes, , drop=FALSE])
    }
  }
  
  avg = as.data.frame(cbind(avg_counts, avg_data))
  if (length(genes)>0) rownames(avg) = genes
  
  return(avg)
}


#' Write differentially expressed genes or markers to an Excel file.
#' 
#' @param degs Result table of the "Seurat::FindMarkers" or "Seurat::FindAllMarkers" functions. Can also be a list of tables so that each table is written into an extra Excel tab. 
#' @param annot_ensembl Ensembl annotation for all genes with Ensembl IDs as rownames.
#' @param gene_to_ensembl Named vector for translating the Seurat gene names of the result table(s) to Ensembl IDs.
#' @param file Output file name.
#' @param additional_readme A data.frame to describe additional columns. Should contain columns 'Column' and 'Description'. Can be NULL.
#' @return Output file name.
DegsWriteToFile = function(degs, annot_ensembl, file, additional_readme=NULL) {
  # Always use list and add annotation
  if (is.data.frame(degs)) degs_lst = list(All=degs) else degs_lst = degs
  
  # Add Ensembl annotation
  for (i in seq(degs_lst)) {
    degs_lst[[i]] = cbind(degs_lst[[i]], annot_ensembl[as.character(degs_lst[[i]]$gene), ])
  }
  
  # Add README
  readme_table = data.frame(Column=c("p_val"), Description=c("Uncorrected p-value"))
  readme_table = rbind(readme_table, 
                       c("avg_log2FC", "Mean log2 fold change group 1 vs group 2"),
                       c("pct.1", "Fraction cells expressing gene in group 1"),
                       c("pct.2", "Fraction cells expressing gene in group 2"),
                       c("p_val_adj", "Adjusted p-value"),
                       c("gene", "Gene"),
                       c("cluster", "Cluster"),
                       c("p_val_adj_score", "Score calculated as follows: -log10(p_val_adj)*sign(avg_log2FC)"),
                       c("avg_RNA_counts_id1", "Average counts in group 1"),
                       c("avg_RNA_data_id1", "Average normalized expression in group 1"))
  
  if (!is.null(additional_readme)) readme_table = dplyr::bind_rows(readme_table, additional_readme)
  
  degs_lst = c(list("README"=readme_table), degs_lst)
  
  # Fix names that are more than 31bp (which is too long for Excel)
  names(degs_lst) = strtrim(names(degs_lst), 31)
  
  # Output in Excel sheet
  openxlsx::write.xlsx(degs_lst, file=file)
  return(file)
}

#' Plot the number of DEGs per test.
#' 
#' @param markers Result table of the "Seurat::FindAllMarkers" function.
#' @param group Group results by column for plotting.
#' @param title Plot title.
#' @return A ggplot object.
DegsPlotNumbers = function(degs, group=NULL, title=NULL) {
  
  degs_up = degs %>% dplyr::filter(avg_log2FC > 0)
  degs_down = degs %>% dplyr::filter(avg_log2FC < 0)
  
  if ((nrow(degs_up) > 0) | (nrow(degs_down) > 0)) {
    degs_n = rbind(data.frame(degs_up, Direction=rep("Up", nrow(degs_up))), data.frame(degs_down, Direction=rep("Down", nrow(degs_down))))
    
    if (!is.null(group)) {
      degs_group_levels = unique(levels(degs_up[, group, drop=TRUE]), levels(degs_down[, group, drop=TRUE]))
      degs_n$Identity = factor(degs_n[, group, drop=TRUE], levels=degs_group_levels)
    } else {
      degs_n$Identity = as.factor("Genes")
    }
    
    degs_n = degs_n %>% dplyr::group_by(Identity, Direction) %>% dplyr::summarise(n=length(Direction))
    p = ggplot(degs_n, aes(x=Identity, y=n, fill=Direction)) + 
      geom_bar(stat="identity") +
      xlab(group) +
      AddPlotStyle(title=title,
               fill=setNames(c("steelblue", "darkgoldenrod1"), c("Down", "Up")))
    return(p)
  }
}

#' Returns an empty deg test table with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj' and 'gene'.
#' 
#' @param col_def Additional columns.
#' @return An R data.frame with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj' and 'gene'.
DegsEmptyResultsTable = function() {
  empty_table = data.frame(p_val=as.numeric(), avg_log2FC=as.numeric(), pct.1=as.numeric(), pct.2=as.numeric(), p_val_adj=as.numeric(), gene=as.character())
  return(empty_table)
}

#' Returns an empty deg marker test table with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster' and 'gene'.
#' 
#' @param clusters: Cluster names to set factor levels of empty 'cluster' column.
#' @return An R data.frame with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster' and 'gene'.
DegsEmptyMarkerResultsTable = function(clusters) {
  empty_table = DegsEmptyResultsTable()
  empty_table$cluster = factor(as.character(), levels=clusters)
  return(empty_table[ c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene')])
}

#' Tests two sets of cells for differential expression using Seurat::FindMarkers.
#' 
#' @param object A Seurat assay object or a Seurat DimReduc object.
#' @param slot If object is a Seurat assay object, which slot to use.
#' @param cells_1 The cell names in set1.
#' @param cells_2 The cell names in set2.
#' @param is_reduction Object is a Seurat DimReduc object. 
#' @param ... Additional parameters passed on to FindMarkers.
#' @return A table with the columns p_val, avg_log2FC, pct.1, pct.2, p_val_adj and gene.
DegsTestCellSets = function(object, slot="data", cells_1=NULL, cells_2=NULL, is_reduction=FALSE, ...){
  # Additional arguments for FindMarkers in the three-dots construct
  additional_arguments = list(...)
  
  # Make sure that object, assay, slot, reduction, cells.1 and cells.2 are not in the additional_arguments list
  additional_arguments[["object"]] = NULL
  additional_arguments[["slot"]] = NULL
  additional_arguments[["cells.1"]] = NULL
  additional_arguments[["cells.2"]] = NULL

  # Create empty table to return when there are no results
  no_degs_results = DegsEmptyResultsTable()
  
  # Check that there are at least 3 cell names and that all cell names are part of the Seurat object
  if (is.null(cells_1) || length(cells_1) < 3 || any(!cells_1 %in% colnames(object))) return(no_degs_results)
  if (is.null(cells_2) || length(cells_2) < 3 || any(!cells_2 %in% colnames(object))) return(no_degs_results)
  
  # Run Seurat::FindMarkers
  if (!is_reduction) {
    arguments = c(list(object=object, slot=slot, cells.1=cells_1, cells.2=cells_2), additional_arguments)
  } else {
    arguments = c(list(object=object, cells.1=cells_1, cells.2=cells_2), additional_arguments)
  }
  
  deg_results = suppressMessages(do.call(Seurat::FindMarkers, arguments))
  if (nrow(deg_results)==0) return(no_degs_results)
  
  # Fix base 2 for log fold change
  if (!"avg_log2FC" %in% colnames(deg_results)) {
    lfc_idx = grep("avg_log\\S*FC", colnames(deg_results))
    deg_results[,lfc_idx] = deg_results[,lfc_idx] / log(2)
    col_nms = colnames(deg_results)
    col_nms[2] = "avg_log2FC"
    colnames(deg_results) = col_nms
  }
  
  # Add column gene
  deg_results$gene = rownames(deg_results)
  return(deg_results)
}


#' Returns an empty Enrichr results table.
# '
# '
#' @param overlap_split: If TRUE, then table will contain the two columns 'In.List' and 'In.Annotation' (which result from splitting 'Overlap') instead of the column 'Overlap'.
#' @return An empty Enrichr results dataframe.
EmptyEnrichrDf = function(overlap_split=FALSE) {
  if (overlap_split) {
    return(data.frame(Term=as.character(), In.List=as.numeric(), In.Annotation=as.numeric(), P.value=as.numeric(), Adjusted.P.value=as.numeric(), Odds.Ratio=as.numeric(), Combined.Score=as.numeric(), Genes=as.character())) 
  } else {
    return(data.frame(Term=as.character(), Overlap=as.character(), P.value=as.numeric(), Adjusted.P.value=as.numeric(), Odds.Ratio=as.numeric(), Combined.Score=as.numeric(), Genes=as.character())) 
  }
}

#' Tests a list of entrez gene symbols for functional enrichment via Enrichr.
# '
#' @param genes: A vector of entrez gene symbols.
#' @param databases: A vector of Enrichr databases with functional annotation.
#' @param padj: Maximum adjusted p-value (0.05).
#' @return The path to the data directory.
EnrichrTest = function(genes, databases, padj=0.05) {
  empty_enrichr_df = EmptyEnrichrDf()
  
  # Run enrichr
  if (length(genes) >= 3 & length(databases) > 0) {
    enrichr_results = suppressMessages(enrichR::enrichr(unique(unname(genes)), databases=databases))
  } else {
    enrichr_results = purrr::map(databases, function(d) return(empty_enrichr_df)) # no good gene lists
    names(enrichr_results) = databases
  }
  
  databases = setNames(databases, databases)
  enrichr_results = purrr::map(databases, function(d) {
    if (nrow(enrichr_results[[d]]) == 0) { # table is empty
      return(empty_enrichr_df)
    } else if (ncol(enrichr_results[[d]]) == 1) { # there was an actual error
      warning(paste("Enrichr returns with error: ", enrichr_results[[d]][1, 1]))
      return(empty_enrichr_df)
    } else {
      return(enrichr_results[[d]])
    }
  })
  
  # Drop Old.P.value and Old.Adjusted.P.value columns if they exist
  enrichr_results = purrr::map(enrichr_results, function(df) {
    df[["Old.P.value"]] = NULL
    df[["Old.Adjusted.P.value"]] = NULL
    return(df)
  })
  
  # Some enrichr server do not return all columns - add missing columns and set them to NA
  enrichr_results = purrr::map(enrichr_results, function(df) {
    character_cols = c("Term", "Overlap", "Genes")
    missing = !character_cols %in% colnames(df)
    if (any(missing)) {
      for (c in character_cols[missing]) df[[c]] = as.character(NA)
    }
    
    numeric_cols = c("P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score")
    missing = !numeric_cols %in% colnames(df)
    if (any(missing)) {
      for (c in numeric_cols[missing]) df[[c]] = as.numeric(NA)
    }
    
    df = df %>% dplyr::select(Term, Overlap, P.value, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes)
    return(df)
  })
  
  # Filter by adjusted p-value
  enrichr_results = purrr::map(enrichr_results, dplyr::filter, Adjusted.P.value < padj)
  
  # Split column "Overlap" into numbers
  enrichr_results = purrr::map(enrichr_results, tidyr::separate, col=Overlap, into=c("In.List", "In.Annotation"), sep="/", convert=TRUE)
  
  # Remap entrez genes to seurat rownames (if seurat rownames are the names of the genes vector)
  if (!is.null(names(genes))) {
    value2names = split(names(genes), genes)
    
    enrichr_results = purrr::map(enrichr_results, function(df) {
      values_mapped_to_names = purrr::map(strsplit(df$Genes, split=";", fixed=TRUE), function(g) {
        return(paste(unlist(value2names[g]), collapse=";")) 
      })
      df$Genes.Species = values_mapped_to_names
      return(df)
    })
  }
  
  return(enrichr_results)
}

#' Writes the enrichr results to an Excel file(s). Includes a README.
# '
#' @param enrichr_results: A list with enrichr results.
#' @param file: A file path.
#' @return The path to the file.
EnrichrWriteResults = function(enrichr_results, file) {
  
  # README
  readme_table = data.frame(Column=c("Term"), Description=c("Functional annotation in database"))
  readme_table = rbind(readme_table, c("In.List", "Number of genes in list of interest with this functional annotation"))
  readme_table = rbind(readme_table, c("In.Annotation", "Total number of genes with this functional annotation"))
  readme_table = rbind(readme_table, c("P.value", "P-value (uncorrected)"))
  readme_table = rbind(readme_table, c("Adjusted.P.value", "P-value (adjusted for multiple testing), use this one"))
  readme_table = rbind(readme_table, c("Odds.Ratio", "How to interpret whether or not the gene list is just random (<1 less than by chance, =1 equals chance, >1 more than by chance)"))
  readme_table = rbind(readme_table, c("Combined.Score", "Combined enrichr score"))
  readme_table = rbind(readme_table, c("Genes", "Enrichr genes in functional annotation"))
  readme_table = rbind(readme_table, c("Genes.Species", "Species genes in functional annotation (which are translated to Enrichr genes)"))
  enrichr_results = c(list("README"=readme_table), enrichr_results)
  
  # Fix names that are more 31bp (which is too long for Excel)
  names(enrichr_results) = strtrim(names(enrichr_results), 31)
  
  # Output in Excel sheet
  openxlsx::write.xlsx(enrichr_results, file=file)
  
  return(file)
}


#' Flattens the list of Enrichr results into one data.frame.
# '
#' @param enrichr_results: A list with enrichr results.
#' @return A data.frame with the columns reported by Enrichr as well as the db
FlattenEnrichr = function(enrichr_results) {
  enrichr_results_flat = purrr::map(names(enrichr_results), function(n) {
    return(data.frame(Database=rep(n, nrow(enrichr_results[[n]])), enrichr_results[[n]]))
  })
  enrichr_results_flat = purrr::invoke(dplyr::bind_rows, enrichr_results_flat)
  return(enrichr_results_flat)
}
