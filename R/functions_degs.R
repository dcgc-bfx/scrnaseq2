#' Sets up a list with contrasts for analysis. Makes sure that all required information is available. Sets default and does basic checks. 
#' Each contrast entry is a list either of length 1 (if there is only one comparison) or of length >1 (if the comparison is done for multiple subsets).
#' 
#' @param sc A Seurat single cell object.
#' @param contrasts_list A list of contrasts. Must at least contain 'name', condition_column', 'condition_group1' and 'condition_group2'.
#' @return A list with contrasts. Each list entry is a list either of length 1 (if there is only one comparison per contrast) or of length >1 
#' (if comparisons are done for multiple subsets for a contrast)
NewContrastsList = function(sc, contrasts_list) {
    # If empty, return empty list
    if (length(contrasts_list) == 0) return(list())
  
    # Get barcode metadata
    barcode_metadata = sc[[]]
    
    # Convert into list, do checks and set up defaults
    contrasts_list = purrr::map(seq(contrasts_list), function(i) {
        contrast = contrasts_list[[i]]
        
        # Trim whitespace
        contrast = purrr::map(contrast, trimws)
        
        # name: should only contain alphanumeric characters, underscores and hyphens
        assertthat::assert_that("name" %in% names(contrast), 
                                msg=FormatString("The name ('name') is missing (for comparison {i})."))
        name = contrast[["name"]]
        assertthat::assert_that(grepl("^[[:alnum:]_\\-]+$", name),
                                msg=FormatString("The name ('name') must only contain alphanumeric characters, underscores and hyphens (for comparison {i})."))
        
        # assay: data that is tested
        if (!"assay" %in% names(contrast)) contrast[["assay"]] = Seurat::DefaultAssay(sc)
        contrast[["assay"]] = contrast[["assay"]] %>% 
            trimws()
        valid_assays = Seurat::Assays(sc)
        valid_reductions = Seurat::Reductions(sc)
        assertthat::assert_that(length(contrast[["assay"]]) == 1 && contrast[["assay"]] %in% c(valid_assays, valid_reductions) | 
                                    all(contrast[["assay"]] %in% names(barcode_metadata)),
                                msg=FormatString("The assay ('assay') must be one of the assays {valid_assays*}, one of the reductions {valid_reductions*} or barcode metadata columns (for comparison {i}/{name})."))
        assay = contrast[["assay"]]
        
        # data_type: which type of data is tested
        if (length(assay) == 1 && assay %in% valid_assays) {
          contrast[["data_type"]] = "feature_data"
        } else if (length(assay) == 1 && assay %in% valid_reductions) {
          contrast[["data_type"]] = "reduction"
        } else {
          contrast[["data_type"]] = "barcode_metadata"
        }
        
        # Identify barcodes that are in assay/reduction. Result is a boolean vector that will be used later to get the correct barcodes.
        if (contrast[["data_type"]] %in% c("feature_data", "reduction")) {
            barcodes_in_assay = rownames(barcode_metadata) %in% SeuratObject::Cells(sc[[assay]])
        } else {
            barcodes_in_assay = rep(TRUE, nrow(barcode_metadata))
        }
        
        # condition_column
        assertthat::assert_that("condition_column" %in% names(contrast),
                                msg=FormatString("The condition column ('condition_column') is missing (for comparison {i}/{name})."))
        condition_column = contrast[["condition_column"]]
        
        # condition_group1 and condition_group2
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
            contrast[["condition_group1_idx"]] = which(barcodes_in_assay & barcode_metadata[, condition_column] %in% condition_group1)
        } else {
            condition_group1 = contrast[["condition_group1"]]
            
            # Sheet number/name appended?
            sheet = 1
            if (grepl(":[^:]+$", condition_group1)) {
                sheet = gsub(pattern=".+:([^:]+)$", replacement="\\1", x=condition_group1) %>% as.integer()
                condition_group1 = gsub(pattern=":[^:]+$", replacement="", x=condition_group1)
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
            
            # Add base file name:sheet number or sheet name as 'condition_group1'
            contrast[["condition_group1"]] = paste0(basename(contrast[["condition_group1"]]), ", sheet ", sheet)
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
            contrast[["condition_group2_idx"]] = which(barcodes_in_assay & barcode_metadata[, condition_column] %in% condition_group2)
        } else {
            condition_group2 = contrast[["condition_group2"]]
            
            # Sheet number appended?
            sheet = 1
            if (grepl(":\\d+$", condition_group2)) {
                sheet = gsub(pattern=".+:(\\d+)$", replacement="\\1", x=condition_group2) %>% as.integer()
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
            
            # Add base file name:sheet number or sheet name as 'condition_group1'
            contrast[["condition_group2"]] = paste0(basename(contrast[["condition_group2"]]), ", sheet ", sheet)
        }
        
        # subset_column and subset_group
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
                    return(which(barcodes_in_assay & barcode_metadata[, subset_column] %in% s))
                })
            } else {
                # Sheet number appended?
                sheet = 1
                if (grepl(":\\d+$", subset_group)) {
                    sheet = gsub(pattern=".+:(\\d+)$", replacement="\\1", x=subset_group) %>% as.integer()
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
                
                # Add base file name:sheet number or sheet name as 'condition_group1'
                contrast[["subset_group"]] = paste0(basename(contrast[["subset_group"]]), ", sheet ", sheet)
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
        }
        
        # pseudobulk_samples
        if ("pseudobulk_samples" %in% names(contrast)) {
            contrast[["pseudobulk_samples"]] = as.integer(contrast[["pseudobulk_samples"]])
            assertthat::assert_that(contrast[["pseudobulk_samples"]] > 1,
                                    msg=FormatString("The number of pseudobulk samples ('pseudobulk_samples') must be greater than 1 (for comparison {i}/{name})."))
        }
        
        # test
        valid_tests = c("wilcox", "bimod", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
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
        if (!"min_pct" %in% names(contrast)) contrast[["min_pct"]] = 0.01
        contrast[["min_pct"]] = as.numeric(contrast[["min_pct"]])
        
        # bulk if bulk_by or pseudobulk_samples is set
        if ("bulk_by" %in% names(contrast) | "pseudobulk_samples" %in% names(contrast)) {
            contrast[["bulk_method"]] = dplyr::case_when(
                contrast[["data_type"]] == "feature_data" ~ "aggregate",
                contrast[["data_type"]] == "reduction" ~ "average",
                contrast[["data_type"]] == "barcode_metadata" ~ "average",
                TRUE ~ "average")
        }
            
        # layer
        if (!"layer" %in% names(contrast)) {
            if (contrast[["data_type"]] == "feature_data") {
                if (contrast[["test"]] %in% c("negbinom", "poisson", "DESeq2")) {
                    contrast[["layer"]] = "counts"
                } else{
                    contrast[["layer"]] = "data"
                }
            } else if(contrast[["data_type"]] == "reduction") {
                contrast[["layer"]] = "counts"
            } else if(contrast[["data_type"]] == "barcode_metadata") {
                contrast[["layer"]] = "counts"
            }
        }
        assertthat::assert_that(contrast[["layer"]] %in% c("counts", "data", "scale.data"),
                                msg=FormatString("The layer ('layer') must be 'counts', 'data', 'scale.data' (for comparison {i}/{name})."))

        # downsample_barcodes
        if ("downsample_barcodes" %in% names(contrast)) {
            contrast[["downsample_barcodes"]] = as.integer(contrast[["downsample_barcodes"]])
        }
        
        # covariate
        if ("covariate" %in% names(contrast)) {
            covariate = contrast[["covariate"]] %>% 
                unlist() %>% 
                trimws() %>% 
                unique()
            
            # Check that covariate columns are in the barcode metadata and factors or numeric
            assertthat::assert_that(all(covariate %in% colnames(barcode_metadata)),
                                    msg=FormatString("The covariate column(s) must be part of the barcode metadata (for comparison {i}/{name})."))
            purrr::walk(covariate, function(c) {
                assertthat::assert_that(is.factor(barcode_metadata[, c]) | is.numeric(barcode_metadata[, c]),
                                        msg=FormatString("The covariate column(s) must be factors or numeric (for comparison {i}/{name})."))
            })
            
            # Check whether the test works with:
            # - covariates at all
            # - character/factor covariates
            # - numeric covariates
            test = contrast[["test"]]
            valid_tests = c("LR", "negbinom", "poisson", "MAST", "DESeq2")
            assertthat::assert_that(test %in% valid_tests,
                                    msg=FormatString("Only the following tests allow covariates: {valid_tests*}."))
            
            categorial_covariates = covariate %>% 
              purrr::keep(function(c) is.factor(barcode_metadata[, c, drop=TRUE]))
            valid_tests = c("LR", "negbinom", "poisson", "MAST", "DESeq2")
            assertthat::assert_that(test %in% valid_tests | length(categorial_covariates) == 0,
                                    msg=FormatString("Only the following tests allow character/categorial covariates: {valid_tests*}."))
            
            numeric_covariates = covariate %>% 
              purrr::keep(function(c) is.numeric(barcode_metadata[, c, drop=TRUE]))
            valid_tests = c("LR", "negbinom", "poisson", "MAST")
            assertthat::assert_that(test %in% valid_tests | length(numeric_covariates) == 0,
                                    msg=FormatString("Only the following tests allow numeric covariates: {valid_tests*}."))
            
            contrast[["covariate"]] = covariate 
        }

        # Add contrast row number
        contrast[["contrast_row"]] = i
        
        return(contrast)
    })
    
    
    # Expand subsets list so that there is now one entry per subset
    # Also get barcode indices per subset
    contrasts_list = purrr::map(contrasts_list, function(contrast) {
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
          names(contrasts_expanded) = purrr::map(contrast[["subset_group"]], paste, collapse="+") %>% unlist()
      } else {
          contrasts_expanded = list(contrast)
      }
      return(contrasts_expanded)
    })
    
    # Add names for contrasts
    contrasts_names = purrr::map_chr(contrasts_list, function(contrast) {
        return(contrast[[1]][["name"]])
    })
    names(contrasts_list) = contrasts_names
    
    return(contrasts_list)
}

#' Given a contrast configuration, prepares a Seurat object. The function extracts all 
#' relevant data and if requested bulk-aggregates and downsamples barcodes. Note that if
#' a contrast has multiple comparisons (for each subset), this function needs to be run on each of them separately.
#' 
#' @param sc A Seurat single cell object.
#' @param contrast A contrast configuration. Must have been set up with NewContrastsList.
#' @return A contrast configuration with an 'object' entry for the Seurat object.
PrepareContrast = function(sc, contrast) {
    barcode_metadata = sc[[]]
    name = contrast[["name"]]
        
    # If condition_group1_idx or condition_group2_idx is empty, return
    if (length(contrast[["condition_group1_idx"]]) == 0 | length(contrast[["condition_group2_idx"]]) == 0) return(contrast)
        
    # Get barcodes indices and barcode names
    condition_group1_idx = contrast[["condition_group1_idx"]]
    condition_group1 = barcode_metadata[condition_group1_idx, , drop=FALSE] %>% rownames()
    condition_group2_idx = contrast[["condition_group2_idx"]]
    condition_group2 = barcode_metadata[condition_group2_idx, , drop=FALSE] %>% rownames()
    
    # When creating a new Seurat, do not calculate nCounts and nFeatures (not needed); restore default on exit
    op = options(Seurat.object.assay.calcn = FALSE)
    on.exit(expr = options(op), add = TRUE)
    
    # Get barcode metadata columns
    barcode_metadata_columns = c(contrast[["condition_column"]], 
                                 contrast[["subset_column"]], 
                                 contrast[["bulk_by"]], 
                                 contrast[["covariate"]])
    barcode_metadata_columns = barcode_metadata_columns[barcode_metadata_columns %in% colnames(barcode_metadata)]
    
    # Type of data to test: feature data, reduction or barcode metadata
    data_type = contrast[["data_type"]]
    
    # Extract relevant data and save it as Seurat object
    if (data_type == "feature_data") {
        assay = contrast[["assay"]]
        barcodes_idx = unique(c(condition_group1_idx, condition_group2_idx))
        
        # Create new assay object (cannot have Seurat without one)
        assay_obj = suppressWarnings({subset(sc[[assay]],
                           cells=barcodes_idx)})
        
        # Update feature metadata to include only feature_id, feature_name, feature_type
        feature_metadata = assay_obj[[]]
        feature_metadata = feature_metadata[, c("feature_id", "feature_name", "feature_type")]
        assay_obj@meta.data = data.frame()
        assay_obj = SeuratObject::AddMetaData(assay_obj, feature_metadata)
        
        # Remove scale.data layer (never needed and saves a lot of memory)
        SeuratObject::DefaultLayer(assay_obj) = contrast[["layer"]]
        SeuratObject::LayerData(assay_obj, layer="scale.data") = NULL
        
        # Set up Seurat object for this analysis
        sc_subset = SeuratObject::CreateSeuratObject(assay_obj, 
                                                     assay=assay)
        
        # Update barcode metadata
        sc_subset = SeuratObject::AddMetaData(sc_subset,
                                              barcode_metadata[barcodes_idx, barcode_metadata_columns, drop=FALSE])
    } else if (data_type == "reduction") {
        reduction_name = contrast[["assay"]]
        assay = "reduction"

        # Get reduction and subset
        reduction = sc[[reduction_name]]
        reduction = subset(reduction, cells=unique(c(condition_group1_idx, condition_group2_idx)))

        # Create new assay object with reduction data
        assay_obj = SeuratObject::CreateAssay5Object(counts=SeuratObject::Embeddings(reduction) %>% 
                                                       as("dgCMatrix") %>% 
                                                       t())
        
        # Add dummy feature metadata
        feature_metadata = data.frame(feature_id=rownames(assay_obj), 
                                      feature_name=rownames(assay_obj), 
                                      feature_type="reduction",
                                      row.names=rownames(assay_obj))
        assay_obj = SeuratObject::AddMetaData(assay_obj, feature_metadata)

        # Set up Seurat object for this analysis
        sc_subset = SeuratObject::CreateSeuratObject(assay_obj, assay=assay, meta.data=barcode_metadata[unique(c(condition_group1_idx, condition_group2_idx)), barcode_metadata_columns, drop=FALSE])
    } else if (data_type == "barcode_metadata") {
        # Get barcode metadata columns to test
        metadata_cols = contrast[["assay"]]
        assay = "barcode_metadata"
        metadata = sc[[metadata_cols]][unique(c(condition_group1_idx, condition_group2_idx)), , drop=FALSE] %>%
            as.matrix() %>%
            Matrix::t()
        
        # CreateAssay5Object cannot deal with matrices that contain only one row, so we need to add a dummy row
        if (nrow(metadata) == 1) {
            dummy_feature = Matrix::Matrix(0, nrow=1, ncol=ncol(metadata), dimnames=list(c("-dummy-"), colnames(metadata)), sparse=TRUE)
            metadata = rbind(metadata, dummy_feature)
        }
        
        # Create new assay object (cannot have Seurat without one)
        assay_obj = SeuratObject::CreateAssay5Object(counts=metadata %>% as("dgCMatrix"))
        
        # Add dummy feature metadata
        feature_metadata = data.frame(feature_id=rownames(assay_obj), 
                                      feature_name=rownames(assay_obj), 
                                      feature_type="barcode_metadata",
                                      row.names=rownames(assay_obj))
        assay_obj = SeuratObject::AddMetaData(assay_obj, feature_metadata)

        # Set up Seurat object for this analysis
        sc_subset = SeuratObject::CreateSeuratObject(assay_obj, assay=assay, meta.data=barcode_metadata[unique(c(condition_group1_idx, condition_group2_idx)), barcode_metadata_columns, drop=FALSE])
    }
    
    # Add a condition column and update idents
    conditions = dplyr::case_when(
        SeuratObject::Cells(sc_subset) %in% condition_group1 ~ "condition1",
        SeuratObject::Cells(sc_subset) %in% condition_group2 ~ "condition2",
        TRUE ~ NA
    ) %>% factor(levels=c("condition1", "condition2"))
    
    assertthat::assert_that(!"condition_groups" %in% colnames(sc_subset[[]]),
                            msg=FormatString("The column 'condition_groups' is reserved, please use another column name (for comparison {i}/{name})."))
    sc_subset[["condition_groups"]] = conditions
    Seurat::Idents(sc_subset) = "condition_groups"
    barcode_metadata = sc_subset[[]]
    
    # Aggregate cells into bulk/pseudo-bulk samples if requested
    if ("bulk_by" %in% names(contrast) | "pseudobulk_samples" %in% names(contrast)) {
        # By default at least aggregate by conditions
        bulk_by = "condition_groups"
        
        # Additionally aggregate by these user-specified columns
        if ("bulk_by" %in% names(contrast)) {
            bulk_by = unique(c(bulk_by, contrast[["bulk_by"]]))
        }

        # Get lists of categorial and numeric covariates
        # Also aggregate by categorial covariates
        categorial_covariates = contrast[["covariate"]] %>% 
          purrr::keep(function(c) is.factor(barcode_metadata[, c, drop=TRUE]))
        numeric_covariates = contrast[["covariate"]] %>% 
          purrr::keep(function(c) is.numeric(barcode_metadata[, c, drop=TRUE]))
        if (length(categorial_covariates) > 0) bulk_by = unique(c(bulk_by, categorial_covariates))
        
        # Finally, if requested, generate x pseudo-bulk samples per bulk_by combination (to aggregate)
        if ("pseudobulk_samples" %in% names(contrast)) {
            num_pseudobulk_samples = contrast[["pseudobulk_samples"]]

            # Divide each group into pseudobulk_samples subgroups
            # Done with modulo (%%): ((row_index - 1) %% num_pseudobulk_samples) + 1
            pseudobulk_sample = barcode_metadata %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(bulk_by))) %>% 
              dplyr::mutate(pseudobulk_sample=0:(dplyr::n() - 1)) %>%
              dplyr::mutate(pseudobulk_sample=(pseudobulk_sample %% num_pseudobulk_samples) + 1) %>%
              dplyr::mutate(pseudobulk_sample=paste0("s", pseudobulk_sample)) %>%
              dplyr::pull(pseudobulk_sample)
            
            # Add pseudobulk_sample column to seurat object
            barcode_metadata$pseudobulk_sample = factor(pseudobulk_sample, levels=paste0("s", 1:num_pseudobulk_samples))
            sc_subset = SeuratObject::AddMetaData(sc_subset, metadata=barcode_metadata[,"pseudobulk_sample", drop=FALSE])
            barcode_metadata = sc_subset[[]]
            
            # Add pseudobulk_sample column to bulk_by
            bulk_by = unique(c(bulk_by, "pseudobulk_sample"))
        }
        
        # Downsample per bulk_by combination used to aggregate (if requested by 'downsample_n')
        if ("downsample_barcodes" %in% names(contrast)) {
            set.seed(getOption("random_seed"))
            sampled_barcodes = barcode_metadata %>% 
                tibble::rownames_to_column() %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(bulk_by))) %>% 
                dplyr::slice_sample(n=contrast[["downsample_barcodes"]]) %>% 
                dplyr::pull(rowname)
            sc_subset = subset(sc_subset, cells=sampled_barcodes)
            barcode_metadata = sc_subset[[]]
        }
        
        # Seurat::PseudobulkExpression will replace '_' with '-' in the columns used for aggregating
        # Keep track of original and modified values to restore original values afterwards
        bulk_by_levels = purrr::map(barcode_metadata[,bulk_by], function(column) {
          levels = levels(column) %||% unique(column)
          levels = setNames(levels, gsub("_", "-", levels))
          return(levels)
        })
        
        # When aggregating by "orig.ident" (dataset name), the column needs to be renamed to "orig_ident" since Seurat::PseudobulkExpression will overwrite it with the names of the aggregated combinations
        # When Seurat::PseudobulkExpression is done, column "orig_ident" is renamed to "orig.ident"
        if ("orig.ident" %in% bulk_by) {
          bulk_by[bulk_by == "orig.ident"] = "orig_ident"
          sc_subset$orig_ident = sc_subset$orig.ident
        }
        
        # Then aggregate counts over bulk_by columns
        sc_subset_agg = Seurat::PseudobulkExpression(object=sc_subset, 
                                                     return.seurat=TRUE, 
                                                     group.by=bulk_by,
                                                     layer="counts",
                                                     method=contrast[["bulk_method"]],
                                                     normalization.method=ifelse(data_type == "feature_data", "LogNormalize", ""),
                                                     verbose=FALSE)
        
        # Rename column orig_ident to orig.ident (see above for reason)
        if ("orig_ident" %in% bulk_by) {
          sc_subset_agg$orig.ident = sc_subset_agg$orig_ident
          sc_subset_agg$orig_ident = NULL
          bulk_by[bulk_by == "orig_ident"] = "orig.ident"
        } else{
          # Else delete
          sc_subset_agg$orig.ident = NULL
        }

        # Restore original values and factor levels for barcode metadata columns since they were modified by Seurat::PseudobulkExpression (see above for reason)
        for(column in bulk_by) {
            modified_to_original = bulk_by_levels[[column]]
            
            vals = sc_subset_agg[[]][, column, drop=TRUE] %>% as.character()
            vals = modified_to_original[vals]
            lvls = unname(modified_to_original)
            
            sc_subset_agg[[]][, column] = factor(vals, levels=lvls)
        }
        
        # Update idents to conditions
        Seurat::Idents(sc_subset_agg) = "condition_groups"
        
        # Numeric covariates need to be aggregated for each group by averaging
        barcode_metadata = sc_subset[[]]
        if (length(numeric_covariates) > 0) {
            numeric_covariate_data = barcode_metadata %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(bulk_by))) %>%
                dplyr::summarise(dplyr::across(dplyr::all_of(numeric_covariates), 
                                               ~ mean(., na.rm=TRUE)))
            
            numeric_covariate_data = dplyr::left_join(sc_subset_agg[[]], numeric_covariate_data, by=bulk_by)
            rownames(numeric_covariate_data) = numeric_covariate_data$orig.ident
            sc_subset_agg = SeuratObject::AddMetaData(sc_subset_agg, metadata=numeric_covariate_data[, numeric_covariates, drop=FALSE])
        }
        
        # Remove scale.data layer (never needed and saves memory)
        SeuratObject::LayerData(sc_subset_agg, assay=SeuratObject::DefaultAssay(sc_subset_agg), layer="scale.data") = NULL
        
        # Add feature metadata columns feature_id, feature_name, feature_type
        feature_metadata = sc_subset[[assay]][[]]
        sc_subset_agg[[assay]] = SeuratObject::AddMetaData(sc_subset_agg[[assay]], feature_metadata)
        
        sc_subset = sc_subset_agg
    } else {
      # Downsample by groups (if requested by 'downsample_barcodes')
      if ("downsample_barcodes" %in% names(contrast)) {
        barcode_metadata = sc_subset[[]]
        sampled_barcodes = barcode_metadata %>%
          tibble::rownames_to_column() %>%
          dplyr::group_by(condition_groups) %>%
          dplyr::slice_sample(n=contrast[["downsample_barcodes"]]) %>%
          dplyr::pull(rowname)
        sc_subset = subset(sc_subset, cells=sampled_barcodes)
      }
    }
        
    contrast[["sc_subset"]] = sc_subset
    
    return(contrast)
}

#' Given a contrast configuration, runs a DEG test.
#' Note that if a contrast has multiple comparisons (for each subset), this function needs to be run on each of them separately.
#' 
#' @param contrast A contrast configuration. Must have been set up with NewContrastsList followed by PrepareContrast.
#' @return A contrast with DEG results.
DegsRunTest = function(contrast) {
    # When running FindMarkers, do not calculate nCounts and nFeatures (not needed); restore default on exit
    op = options(Seurat.object.assay.calcn=FALSE)
    on.exit(expr = options(op), add = TRUE)
    
    # Name of contrast
    name = contrast[["name"]]
    
    # Set up basic arguments for FindMarkers
    arguments = list(object=contrast[["sc_subset"]],
                     test.use=contrast[["test"]], 
                     ident.1="condition1",
                     ident.2="condition2",
                     random.seed=getOption("random_seed"),
                     min.cells.group=2,
                     logfc.threshold=contrast[["log2FC"]])
    
    # Specific to what to test: assay with feature data, assay with barcode metadata or reduction
    if (contrast[["data_type"]] == "feature_data") {
      arguments[["assay"]] = assay = contrast[["assay"]]
      arguments[["slot"]] = layer = contrast[["layer"]]
      
      # If the data is stored as IterableMatrix (on-disk via BPCells), convert now to in-memory dgCMatrix unless test is wilcox
      # This is because all other tests in FindMarkers do not support IterableMatrix
      test = contrast[["test"]]
      sc_subset = contrast[["sc_subset"]]
      if (test != "wilcox" & 
          is(SeuratObject::LayerData(sc_subset, assay=assay, layer=layer), "IterableMatrix")) {
        SeuratObject::LayerData(sc_subset, assay=assay, layer=layer) = as(SeuratObject::LayerData(sc_subset, assay=assay, layer=layer), "dgCMatrix")
      }
      contrast[["sc_subset"]] = sc_subset
      arguments[["object"]] = contrast[["sc_subset"]]
      
      # How does a empty table look like
      empty_deg_results = data.frame(avg_log2FC=as.numeric(), pct.1=as.numeric(), pct.2=as.numeric())
      
      # From which layer/slot to compute the mean per group: 'data' (even though test is done using 'counts')
      arguments[["fc.slot"]] = "data"
      
      # How to compute the mean per group: Use the default function of Seurat::FindMarkers suggested for 'data'
      # Whichs is: exponentiate log-tramsformed data ('data'), calculate mean, then log-transform
      # arguments[["mean.fxn"]] = "data"
      
      # Name of the column
      arguments[["fc.name"]] = "avg_log2FC"
      
      # Minimum percentage of cells expressing
      arguments[["min.pct"]] = contrast[["min_pct"]]
      
    } else if (contrast[["data_type"]] == "reduction") {
      # Test reduction assay
      arguments[["assay"]] = assay = "reduction"
      arguments[["slot"]] = layer = "counts"
      
      # How does a empty table look like
      empty_deg_results = data.frame(avg_diff=as.numeric(), pct.1=as.numeric(), pct.2=as.numeric())
      
      # From which layer/slot to compute the mean per group: 'data' (even though test is done on 'counts')
      arguments[["fc.slot"]] = "counts"
      
      # How to compute the mean per group: apply standard mean function to unnormalized data ('counts')
      arguments[["mean.fxn"]] = rowMeans
      
      # Name of the column
      arguments[["fc.name"]] = "avg_diff"
    } else if (contrast[["data_type"]] == "barcode_metadata") {
      # Test assay with barcode metadata
      arguments[["assay"]] = assay = "barcode_metadata"
      arguments[["slot"]] = layer = "counts"
      
      # How does a empty table look like
      empty_deg_results = data.frame(avg_diff=as.numeric(), pct.1=as.numeric(), pct.2=as.numeric())
      
      # From which layer/slot to compute the mean per group: 'data' (even though test is done on 'counts')
      arguments[["fc.slot"]] = "counts"
      
      # How to compute the mean per group: apply standard mean function to unnormalized data ('counts')
      arguments[["mean.fxn"]] = rowMeans
      
      # Name of the column
      arguments[["fc.name"]] = "avg_diff"
    }
    
    # Bulk and covariate arguments
    if ("bulk_by" %in% names(contrast) | "pseudobulk_samples" %in% names(contrast)) arguments[["densify"]] = TRUE
    if ("covariate" %in% names(contrast)) {
      if (contrast[["test"]] == "DESeq2") {
        # Seurat's standard function for running DESeq2 does not allow covariates/batches
        # Replace Seurats standard function for DESeq2 with a modified version that can include a covariate as batch.
        assignInNamespace("DESeq2DETest",DESeq2DETest_covariate, ns="Seurat")
        
        # Pass batch information
        barcode_metadata = contrast[["sc_subset"]][[]]
        arguments[["batch"]] = barcode_metadata[,contrast[["covariate"]]]
        names(arguments[["batch"]]) = rownames(barcode_metadata)
      } else {
        # Other tests use the argument latent.vars (if possible)
        arguments[["latent.vars"]] = contrast[["covariate"]]
      }
    }
    
    # Run FindMarkers
    # Check that there are enough samples in each group (>=2), if not, add an error message and skip
    # Seurat object in contrast[["sc_subset"]]
    condition_counts = SeuratObject::Idents(contrast[["sc_subset"]]) %>% table()
    if (condition_counts["condition1"] >= 2 & condition_counts["condition2"] >= 2) {
        deg_results = do.call(Seurat::FindMarkers, arguments)
    } else {
        deg_results = empty_deg_results
        contrast[["message"]] = FormatString("There are fewer than two samples in at least one group for comparison {name}.")
    }
    
    # Gene rownames to column
    deg_results$gene = rownames(deg_results)
    
    # If the table is empty, add a pval and pval_adj column
    if (!"p_val" %in% colnames(deg_results)) deg_results$p_val = as.numeric(rep(NA, nrow(deg_results)))
    if (!"p_val_adj" %in% colnames(deg_results)) deg_results$p_val_adj = as.numeric(rep(NA, nrow(deg_results)))

    # Sort results
    deg_results = deg_results %>% 
        DegsSort()

    # Add means (per condition)
    # - For feature data (gene expression): use layer 'data'  (normalized and log-transformed) and calculate means of the exponentiated data followed by log.
    #   To (re-)calculate avg_log2FC from the natural log-transformed means, switch to base 2 and subtract: (mean1/log(2)) - (mean2/log(2)). See log laws.
    # - For reduction and barcode metadata: use layer 'counts (raw)' and calculate just the mean.
    avg_df = AverageCounts(contrast[["sc_subset"]], 
                           assay=assay, 
                           layer=ifelse(contrast[["data_type"]] %in% c("reduction", "barcode_metadata"), "counts", "data"))
    avg_df = as.data.frame(avg_df) %>%
      tibble::rownames_to_column(var="gene")
    deg_results = dplyr::inner_join(deg_results, avg_df, by="gene")
    
    # If bulk_by or pseudobulk_samples is set, add the bulked expression values
    if ("bulk_by" %in% names(contrast) | "pseudobulk_samples" %in% names(contrast)) {
      bulk_df = SeuratObject::GetAssayData(contrast[["sc_subset"]], 
                                           assay=assay,
                                           layer=ifelse(contrast[["data_type"]] %in% c("reduction", "barcode_metadata"), "counts", "data")) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var="gene")
      deg_results = dplyr::inner_join(deg_results, bulk_df, by="gene")
    }
    
    # Remove Seurat object
    contrast[["sc_subset"]] = NULL
    
    # Add rownames
    rownames(deg_results) = deg_results$gene
    
    # Reorder the columns
    # Note: depending what was tested there can be avg_log2FC or avg_diff. The next lines show how to use variable in tidyverse.
    fc_col = dplyr::case_match(contrast[["data_type"]], "feature_data" ~ "avg_log2FC", "reduction" ~ "avg_diff", "barcode_metadata" ~ "avg_diff")
    deg_results = deg_results %>% 
        dplyr::relocate(gene, p_val, p_val_adj, !!sym(fc_col), p_val_adj_score, pct.1, pct.2)
    
    # Add results to contrast
    contrast[["results"]] = deg_results
    
    return(contrast)
}


#' Runs an over-representation analysis (ORA) for an DEG result produced by DegsRunTest.
#' Note that if there are multiple results (e.g. for each subset), this function needs to be run on each of them separately.
#' 
#' @param deg_result A DEG result produced by DegsRunTest.
#' @param term2gene_db A data frame with genesets to be used for ORA. Must contain the columns 'gs_cat' (category of the geneset), 'gs_subcat' (sub-category of the geneset), 'gs_name' (name of the geneset), 'gene_symbol' (symbol of the gene in the geneset). See clusterProfiler documentation for more information.
#' @param genesets A character vector specifying the genesets to be used for the ORA. Each entry specifies one analysis. Entries should have the form gs_cat:gs_subcat:gs_name and can contain wildcards (e.g. C5:GO:BP:*).
#' @return A list with one or more ORA results.
DegsRunOraTest = function(deg_result, term2gene_db, genesets) {
    # Get DEGs
    degs = deg_result$results %>% 
      dplyr::filter(p_val_adj < deg_result$padj & abs(avg_log2FC) >= deg_result$log2FC) %>%
      dplyr::pull(gene) %>%
      unique()
    
    # Get universe
    universe = deg_result$results %>%
      dplyr::pull(gene) %>%
      unique()
    
    # If there are no DEGs or universe, return NULL
    if (length(degs) == 0 | length(universe) == 0) return(NULL)
    
    # Search the genesets specified by the genesets argument
    search_vector = paste(term2gene_db$gs_cat, term2gene_db$gs_subcat, term2gene_db$gs_name, sep=":")
    term2gene = purrr::map(genesets, function(term) {
      i = which(grepl(term, search_vector))
      tg = term2gene_db[i,] %>% dplyr::select(gs_name, gene_symbol)
      return(tg)
    })
    names(term2gene) = genesets
    
    # Run ORA (for one or more genesets)
    ora_result = purrr::map(term2gene, function(t2g) {
      # Set random seed so that results stay reproducible
      set.seed(getOption("random_seed"))
      
      # Run ORA
      ora = clusterProfiler::enricher(gene=degs, universe=universe, TERM2GENE=t2g, minGSSize=10, maxGSSize=500)
      
      # Sort by descreasing FoldEnrichment
      ora@result = ora@result[order(ora@result$FoldEnrichment, decreasing=TRUE),]
      
      # Fix factor levels accordingly for plots
      ora@result$ID = factor(ora@result$ID, levels=unique(ora@result$ID))
      ora@result$Description = factor(ora@result$Description, levels=unique(ora@result$Description))
      
      return(ora)
    })
    
    return(ora_result)
}

#' Runs a geneset enrichment analysis (GSEA) for a constrast produced by PrepareContrast.
#' Note that if there are multiple contrasts (e.g. for each subset), this function needs to be run on each of them separately.
#' 
#' @param contrast A contrast produced by PrepareContrast.
#' @param term2gene_db A data frame with genesets to be used for GSEA. Must contain the columns 'gs_cat' (category of the geneset), 'gs_subcat' (sub-category of the geneset), 'gs_name' (name of the geneset), 'gene_symbol' (symbol of the gene in the geneset). See clusterProfiler documentation for more information.
#' @param genesets A character vector specifying the genesets to be used for the GSEA. Each entry specifies one analysis. Entries should have the form gs_cat:gs_subcat:gs_name and can contain wildcards (e.g. C5:GO:BP:*).
#' @return A list with one or more GSEA results.
DegsRunGseaTest = function(contrast, term2gene_db, genesets) {
  # Calculate fold changes
  fold_changes = Seurat::FoldChange(object=contrast$sc_subset,
                                    ident.1="condition1",
                                    ident.2="condition2",
                                    assay=contrast$assay,
                                    slot=contrast$layer,
                                    pseudocount.use=1,
                                    base=2)
  fold_changes = setNames(fold_changes$avg_log2FC, rownames(fold_changes))
  fold_changes = fold_changes[order(fold_changes, decreasing=TRUE)]
  
  # Search the genesets specified by the genesets argument
  search_vector = paste(term2gene_db$gs_cat, term2gene_db$gs_subcat, term2gene_db$gs_name, sep=":")
  term2gene = purrr::map(genesets, function(term) {
    i = which(grepl(term, search_vector))
    tg = term2gene_db[i,] %>% dplyr::select(gs_name, gene_symbol)
    return(tg)
  })
  names(term2gene) = genesets
  
  # Run GSEA
  gsea_result = purrr::map(term2gene, function(t2g) {
    # Set random seed so that results stay reproducible
    set.seed(getOption("random_seed"))
    
    # Run GSEA
    gsea = clusterProfiler::GSEA(geneList=fold_changes, TERM2GENE=t2g, minGSSize=10, maxGSSize=500, seed=TRUE)
    return(gsea)
  })
  
  return(gsea_result)
}

#' Obtains genesets from MSigDB or user-defined input files.
#' 
#' @param msigdb_species Species name for MSigDB.
#' @param is_msigdb_species_name Is the species name a name used by MSigDB? Default is FALSE. If not, try to make it compatible with MSigDB. This will do the following: homo_sapiens => Homo sapiens.
#' @param geneset_files User-defined genesets in Excel or CSV file(s) to use instead of MSigDB. The files should have these columns: gs_cat (category), gs_subcat (subcategory), gs_name (geneset name), gene_id, gene_symbol.
#' @return A table with genesets and the following columns: gs_cat, gs_subcat, gs_name, gene_id, gene_symbol.
DegsGetGenesets = function(msigdb_species, is_msigdb_species_name=FALSE, geneset_files=NULL) {
  if (is.null(geneset_files)) {
    # No user input - get MSigDB
    
    # Format species name for MSigDB
    # If not provided explicitly, format species name for MSigDB
    if (!is_msigdb_species_name) {
      msigdb_species = strsplit(msigdb_species, "_") %>% unlist()
      if(length(msigdb_species)>1) {
        msigdb_species = msigdb_species[1:2]
      } else {
        msigdb_species = msigdb_species[1]
      }
      msigdb_species[1] = stringr::str_to_title(msigdb_species[1])
      msigdb_species = paste(msigdb_species, collapse=" ")
    }
    # Download MSigDB
    term2gene_db = msigdbr::msigdbr(species=msigdb_species)
  } else {
    # User input
    term2gene_db = purrr::map_dfr(geneset_files, function(geneset_file) {
      # Checkf file extension
      extension = tools::file_ext(gsub(pattern="\\.gz$", replacement="", x=geneset_file))
      valid_extensions = c("csv", "tsv", "xls", "xlsx", "rds")
      assertthat::assert_that(extension %in% valid_extensions,
                              msg=FormatString("Geneset file type must be: {valid_extensions*} (file can be gzipped)."))
      
      # Read table
      if (extension %in% c("xls", "xlsx")) {
        term2gene = readxl::read_excel(geneset_file, col_names=TRUE)
      } else {
        readr::read_delim(geneset_file, col_names=TRUE, comment="#", progress=FALSE, show_col_types=FALSE, col_types=readr::cols())
      }
      
      # Check that the table contains the required
      required_columns = c("gs_cat", "gs_subcat", "gs_name", "gene_id", "gene_symbol")
      assertthat::assert_that(all(required_columns %in% colnames(term2gene)),
                              msg=FormatString("Geneset file {file} is missing at least one of these columns: 'gs_cat', 'gs_subcat', 'gs_name', 'gene_id', 'gene_symbol'."))
      term2gene = dplyr::select(term2gene, gs_cat, gs_subcat, gs_name, gene_id, gene_symbol)
      
      return(term2gene)
    })
  }
  
  return(term2gene_db)
}

#' Copy of the Seurat::DESeq2DETest function with the addition of a covariate argument.
#' 
#' Needs to replaced in the Seurat namespace in the main code.
#' 
#' @param data.use Counts matrix.
#' @param cells.1 Cell names for group 1.
#' @param cells.2 Cell names for group 2.
#' @param batch Vector with batch information where names are cell names.
#' @param verbose Print progress.
#' @return A data frame with p-values.
DESeq2DETest_covariate = function (data.use, cells.1, cells.2, batch=NULL, verbose=TRUE, ...) {
  SeuratObject::CheckDots(..., fxns = "DESeq2::results")
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  
  # Improvements:
  # - Add batch variable to DESeq2 coldata
  # - Set design variable
  # - Add batch to design
  design = "~group"
  if (!is.null(batch)) {
    group.info[, "batch"] <- factor(as.character(batch[rownames(group.info)]))
    design = paste0(design, "+batch")
  }
  group.info$wellKey <- rownames(x = group.info)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(countData = data.use, 
                                         colData = group.info, design = as.formula(design))
  dds1 <- DESeq2::estimateSizeFactors(object = dds1)
  dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
  dds1 <- DESeq2::nbinomWaldTest(object = dds1)
  res <- DESeq2::results(object = dds1, contrast = c("group", 
                                                     "Group1", "Group2"), alpha = 0.05, ...)
  to.return <- data.frame(p_val = res$pvalue, row.names = rownames(res))
  return(to.return)
}
    
#' Sorts table of differentially expressed genes per performed test. Introduces a signed p-value score calculated as follows:
#' p_val_adj_score = -log10(p_val_adj) * sign(avg_log2FC).
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function or the "Seurat::FindMarkers" function.
#' @param group Group results first by column(s) before sorting.
#' @return Sorted table with differentially expressed genes. 
DegsSort = function(degs, group=NULL) { 
  # Group first (if requested)
  if (!is.null(group)) degs = degs %>% dplyr::group_by(dplyr::across(dplyr::all_of(group)))
  
  # Introduce signed p-value score and sort
  # Decide which column to use for foldchange
  if ("avg_log2FC" %in% colnames(degs)) {
    degs$p_val_adj_score = -log10(degs$p_val_adj) * sign(degs$avg_log2FC)
    degs = degs %>% dplyr::arrange(-p_val_adj_score, -avg_log2FC, .by_group=TRUE)
  } else if ("avg_diff" %in% colnames(degs)) {
    degs$p_val_adj_score = -log10(degs$p_val_adj) * sign(degs$avg_diff)
    degs = degs %>% dplyr::arrange(-p_val_adj_score, -avg_diff, .by_group=TRUE)
  }
  
  degs = as.data.frame(degs)
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

#' Creates a DEG scatterplot.
#' 
#' @param deg_result A list entry with DEG results obtained with RunDEGTests
#' @param font_size The base font size. Default is 11.
#' @return A ggplot scatterplot
DegsScatterPlot = function(deg_result, font_size=11) {
    # Get condition names. If multiple, join by '+'.
    group1 = paste(deg_result[["condition_group1"]], collapse="+")
    group2 = paste(deg_result[["condition_group2"]], collapse="+")
    
    # Get DEG table and add DEG status for significant up- and down-regulated (include log2 foldchange threshold when testing counts data)
    deg_table = deg_result[["results"]]
    padj = deg_result[["padj"]]
    
    if (deg_result[["data_type"]] == "feature_data") {
      log2FC = deg_result[["log2FC"]]
      deg_table$deg_status = dplyr::case_when(
        deg_table$p_val_adj < padj & abs(deg_table$avg_log2FC) >= log2FC & deg_table$avg_log2FC > 0 ~ "up",
        deg_table$p_val_adj < padj & abs(deg_table$avg_log2FC) >= log2FC & deg_table$avg_log2FC < 0 ~ "down",
        TRUE ~ "none"
      )
    } else {
      log2FC = -Inf
      deg_table$deg_status = dplyr::case_when(
        deg_table$p_val_adj < padj & deg_table$avg_diff > 0 ~ "up",
        deg_table$p_val_adj < padj & deg_table$avg_diff < 0 ~ "down",
        TRUE ~ "none"
      )
    }    
    deg_table$deg_status = factor(deg_table$deg_status, levels=c("none", "up", "down"))
    lims = c(min(c(deg_table$condition1, deg_table$condition2)), max(deg_table$condition1, deg_table$condition2))
    
    # Get the top 5 up- and down-regulated DEGs
    top10_deg_table = deg_table %>% 
        dplyr::filter(deg_status %in% c("up", "down")) %>%
        dplyr::group_by(deg_status) %>%
        dplyr::slice_min(order_by=abs(p_val ), n=5)
    
    # Make plot
    # Note: Font size in theme is measured in pt but font size in geom_text is measured in mm.
    # Ggplot2 provided '.pt' as conversion factor: mm = pt / .pt OR pt = mm * .pt
    # https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#text
    p = ggplot(deg_table, aes(x=condition1, y=condition2, col=deg_status)) + 
        geom_abline(slope=1, intercept=0, col="lightgrey") +
        geom_point() +
        ggrepel::geom_text_repel(data=top10_deg_table, aes(x=condition1, y=condition2, col=deg_status, label=gene), size=font_size / .pt) +
        scale_color_manual("Gene status", values=c(none="grey", up="darkgoldenrod1", down="steelblue"), 
                           labels=c(none='none', up='up', down='down')) +
        xlim(lims) + 
        ylim(lims) +
        AddPlotStyle(ylab=group1, xlab=group2, legend_position="none", font_size=font_size)
    
    # If there is a log2 threshold > 0, add lines
    if (log2FC > 0) {
        p = p + geom_abline(slope=1, intercept=c(-log2FC, log2FC), col="lightgrey", lty=2)
    }
    
    
    return(p)
}

#' Creates a DEG volcano plot
#' 
#' @param deg_result A list entry with DEG results obtained with RunDEGTests
#' @param font_size The base font size. Default is 11.
#' @return A ggplot volcano plot
DegsVolcanoPlot = function(deg_result, font_size=11) {
    # Get condition names. If multiple, join by '+'.
    group1 = paste(deg_result[["condition_group1"]], collapse="+")
    group2 = paste(deg_result[["condition_group2"]], collapse="+")
    
    # Decide whether to use column avg_log2FC or avg_diff
    fc_col = dplyr::case_match(deg_result[["data_type"]], "feature_data" ~ "avg_log2FC", "reduction" ~ "avg_diff", "barcode_metadata" ~ "avg_diff")
    
    # Get DEG table and add DEG status for significant up- and down-regulated (include log2 foldchange threshold when testing counts data)
    deg_table = deg_result[["results"]]
    padj = deg_result[["padj"]]
    
    if (deg_result[["data_type"]] == "feature_data") {
      log2FC = deg_result[["log2FC"]]
      deg_table$deg_status = dplyr::case_when(
        deg_table$p_val_adj < padj & abs(deg_table$avg_log2FC) >= log2FC & deg_table$avg_log2FC > 0 ~ "up",
        deg_table$p_val_adj < padj & abs(deg_table$avg_log2FC) >= log2FC & deg_table$avg_log2FC < 0 ~ "down",
        TRUE ~ "none"
      )
    } else {
      log2FC = -Inf
      deg_table$deg_status = dplyr::case_when(
        deg_table$p_val_adj < padj & deg_table$avg_diff > 0 ~ "up",
        deg_table$p_val_adj < padj & deg_table$avg_diff < 0 ~ "down",
        TRUE ~ "none"
      )
    }    
    deg_table$deg_status = factor(deg_table$deg_status, levels=c("none", "up", "down"))
    
    # Add -10log10(p_val) for plotting. Make sure that infinite values are replaced with 300 and negative infinite with 0
    deg_table$p_val_log10_n = -log10(deg_table$p_val)
    i = which(is.infinite(deg_table$p_val_log10_n) & deg_table$p_val_log10_n > 0)
    deg_table$p_val_log10_n[i] = 300
    i = which(is.infinite(deg_table$p_val_log10_n) & deg_table$p_val_log10_n < 0)
    deg_table$p_val_log10_n[i] = 0
    
    # Get the top 5 up- and down-regulated DEGs
    top10_deg_table = deg_table %>% 
        dplyr::filter(deg_status %in% c("up", "down")) %>%
        dplyr::group_by(deg_status) %>%
        dplyr::slice_min(order_by=abs(p_val ), n=5)
    
    # Make plot
    # Note: Font size in theme is measured in pt but font size in geom_text is measured in mm.
    # Ggplot2 provided '.pt' as conversion factor: mm = pt / .pt OR pt = mm * .pt
    # https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#text
    p = ggplot(deg_table, aes(x=!!sym(fc_col), y=p_val_log10_n, col=deg_status)) + 
        geom_point() +
        ggrepel::geom_text_repel(data=top10_deg_table, aes(x=!!sym(fc_col), y=p_val_log10_n, col=deg_status, label=gene), size=font_size / .pt) +
        scale_color_manual("Gene status", values=c(none="grey", up="darkgoldenrod1", down="steelblue"), 
                           labels=c(none='none', up='up', down='down')) +
        AddPlotStyle(xlab=ifelse(fc_col == "avg_log2FC", expression(log[2]~"foldchange"), "condition1 - condition2"),
                     ylab=expression("-"~log[10]~"p-value"), 
                     legend_position="none", font_size=font_size)
    
    # If there is a log2 threshold > 0, add vertical lines
    if (log2FC > 0) {
        p = p + geom_vline(xintercept=c(-log2FC, log2FC), linetype="dashed", color="grey")
    }
    
    return(p)
}

#' Calculate the average gene expression data for groups of cell barcodes.
#' 
#' @param sc Seurat object.
#' @param group_by How to group the barcodes. Can be a column of the barcode metadata with character or factor data or a list with groups. If NULL, the active identity will be used.
#' @param assay Assay to be used. If NULL, the default assay will be used.
#' @param layer Layer to be used. If NULL, the default layer of the assay will be used.
#' @return A table with average gene expression data per group.
AverageCounts = function(sc, group_by=NULL, assay=NULL, layer=NULL) {
  # If assay is NULL, use default assay
  if (is.null(assay)) assay = SeuratObject::DefaultAssay(sc)
  
  # If layer is NULL, use default layer
  if (is.null(layer)) layer = SeuratObject::DefaultLayer(sc[[assay]])
  
  # Get the data
  ldat = SeuratObject::LayerData(sc, assay=assay, layer=layer)
  
  # Define groups to calculate means
  if (is.null(group_by)) {
    # If group_by is NULL, use active identities
    groups = SeuratObject::Idents(sc)
    if (!is.factor(groups)) groups = factor(groups)
    groups = split(names(groups), groups)
  } else if (is.character(group_by)) {
    # Assert that the group_by column is part of the Seurat object
    assertthat::assert_that(group_by %in% colnames(sc[[]]),
                            msg="The group_by column must be part of the barcode metadata.")
    groups = sc[[]][, group_by, drop=FALSE]
    
    # Assert that groups is factor or character
    assertthat::assert_that(is.factor(groups) | is.character(groups),
                            msg="The group_by column must be a factor or character.")
    if (!is.factor(groups)) groups = factor(groups)
    groups = split(names(groups), groups)
  } else if (is.list(group_by)) {
    # If group_by is a list with barcodes, use the list
    groups = group_by
    
    # Assert that it has names
    assertthat::assert_that(!is.null(names(groups)),
                            msg="The group_by list must have names.")
  } else {
    stop("group_by must be NULL, a character or a list.")
  }
  
  # Now iterate over groups and calculate means
  group_means = purrr::map(groups, function(group) {
    if (layer %in% c("counts", "scale.data")) {
      # Counts: just calculate the mean
      means = rowMeans(ldat[, group, drop=FALSE])
    } else if (layer == "data") {
      # Data: calculate the mean of the exponentiated data, then log again
      means = (rowSums(expm1(ldat[, group, drop=FALSE])) + 1) / NCOL(ldat[, group, drop=FALSE])
      means = log(means)
    }
  })
  group_means = as.data.frame(group_means)
  return(group_means)
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
#' @param degs Table with DEG analysis results. Can also be a list of tables so that each table is written into an extra Excel tab.
#' @param file Output file name.
#' @param annotation Gene annotation to include in the tables. Will be merged using the rownames. Can be NULL.
#' @param parameter A data.frame for describing test parameter. Can be NULL.
#' @return Output file name.
DegsWriteToFile = function(degs_lst, file, annotation=NULL, parameter=NULL) {
    # Convert to list if not already
    if (is.data.frame(degs_lst)) degs_lst = list(degs_lst)
    
    # Add annotation if available (merge via rownames)
    if (!is.null(annotation)) {
        degs_lst = purrr::map(degs_lst, function(degs) {
            degs = degs %>% 
                tibble::rownames_to_column() %>%
                dplyr::left_join(annotation %>% tibble::rownames_to_column(), by=c("rowname")) %>%
                dplyr::select(-rowname)
            return(degs)
        })
    }
    
    # Add parameter
    if (!is.null(parameter)) {
      degs_lst = c(list("Parameter"=parameter), degs_lst)
    }
    
    # Add README
    readme_table = data.frame(Column=c("gene"), Description=c("Feature (e.g. gene, barcode metadata column or reduction)"))
    readme_table = rbind(readme_table, 
                         c("p_val", "Individual per-test p-value uncorrected"),
                         c("p_val_adj", "P-value adjusted for multiple testing using Bonferoni"),
                         c("avg_log2FC", "Mean log2 fold change for condition 1 vs condition 2 when testing counts data"),
                         c("avg_diff", "Difference between the means for condition 1 vs condition 2 when testing barcode metadata or reductions"),
                         c("p_val_adj_score", "Combined score defined as -log10(p_val_adj)*sign(avg_log2FC or avg_diff)"),
                         c("pct.1", "Fraction cells expressing feature in condition 1"),
                         c("pct.2", "Fraction cells expressing feature in condition 2"),
                         c("condition1", "Mean condition 1. For counts data, values are normalized and natural-log-transformed."),
                         c("condition2", "Mean condition 2. For counts data, values are normalized and natural-log-transformed."),
                         c("...", "If bulk-aggregating was done, values for aggregated samples. For counts data, values are normalized and natural-log-transformed."),
                         c("...", "Additional annotation (if provided)"))
    degs_lst = c(list("README"=readme_table), degs_lst)
    
    # Fix names that are more than 31bp (which is too long for Excel)
    names(degs_lst) = strtrim(names(degs_lst), 31)
    
    # Output in Excel sheet
    openxlsx::write.xlsx(degs_lst, file=file)
    return(file)
}

#' Write ORA (over-representation analysis) results to an Excel file.
#' 
#' @param degs Table with ORA results. Can also be a list of tables so that each table is written into an extra Excel tab.
#' @param file Output file name.
#' @param parameter A data.frame for describing test parameter. Can be NULL.
#' @return Output file name.
DegsWriteOraToFile = function(ora_lst, file, parameter=NULL) {
  # Convert to list if not already
  if (is.data.frame(ora_lst)) ora_lst = list(ora_lst)
  
  # Add parameter
  if (!is.null(parameter)) {
    ora_lst = c(list("Parameter"=parameter), ora_lst)
  }
  
  # Add README
  readme_table = data.frame(Column=c("GeneSet"), Description=c("Geneset name"))
  readme_table = rbind(readme_table, 
                       c("ID", "ID of the term tested"),
                       c("Description", "Description of the geneset term tested"),
                       c("GeneRatio", "Number of genes in the input list that are annotated to the term / Number of genes in the input list"),
                       c("BgRatio", "Number of genes in the background gene list that are annotated to the term / Number of genes in the background gene list"),
                       c("RichFactor", "Number of genes in the input list that are annotated to the term / Number of genes that are annotated to the term"),
                       c("FoldEnrichment", "Fold enrichment defined as GeneRatio / BgRatio"),
                       c("zScore", "Z-Score"),
                       c("pvalue", "Over-representation assessed using hypogeometric distribution (one sided Fisher's exact test)"),
                       c("p.adjust", "Hypergeometric p-value after correction for multiple testing (BH)"),
                       c("qvalue", "FDR adjusted p-value"),
                       c("geneID", "Gene IDs in the geneset that are annotated to the term"),
                       c("Count", "Number of genes in the geneset that are annotated to the term"),
                       c("...", "Additional annotation (if provided)"))
                       ora_lst = c(list("README"=readme_table), ora_lst)
  
  # Fix names that are more than 31bp (which is too long for Excel)
  names(ora_lst) = strtrim(names(ora_lst), 31)
  
  # Output in Excel sheet
  openxlsx::write.xlsx(ora_lst, file=file)
  return(file)
}

#' Write geneset enrichment analysis (GSEA) results to an Excel file.
#' 
#' @param gsea_lst Table with GSEA results. Can also be a list of tables so that each table is written into an extra Excel tab.
#' @param file Output file name.
#' @param parameter A data.frame for describing test parameter. Can be NULL.
#' @return Output file name.
DegsWriteGseaToFile = function(gsea_lst, file, parameter=NULL) {
  # Convert to list if not already
  if (is.data.frame(gsea_lst)) ora_lst = list(gsea_lst)
  
  # Add parameter
  if (!is.null(parameter)) {
    gsea_lst = c(list("Parameter"=parameter), gsea_lst)
  }
  
  # Add README
  readme_table = data.frame(Column=c("GeneSet"), Description=c("Geneset name"))
  readme_table = rbind(readme_table, 
                       c("ID", "ID of the term tested"),
                       c("Description", "Description of the geneset term tested"),
                       c("setSize", "Number of genes in geneset"),
                       c("enrichmentScore", "Represents the degree to which a set S is over-represented at the top (positive fold change) or bottom (negative fold change)"),
                       c("pvalue", "Significance of the enrichment score determined by a permutation test"),
                       c("p.adjust", " Permutation p-value after correction for multiple testing (BH)"),
                       c("qvalue", "FDR adjusted p-value"),
                       c("rank", "Mean rank"),
                       c("leading_edge", "Leading edge (see GSEA doc)"),
                       c("core_enrichment", "Core enrichment (see GSEA doc)"))
  gsea_lst = c(list("README"=readme_table), gsea_lst)
  
  # Fix names that are more than 31bp (which is too long for Excel)
  names(gsea_lst) = strtrim(names(gsea_lst), 31)
  
  # Output in Excel sheet
  openxlsx::write.xlsx(gsea_lst, file=file)
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
