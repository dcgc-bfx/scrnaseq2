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
        
        # assay 
        if (!"assay" %in% names(contrast)) contrast[["assay"]] = Seurat::DefaultAssay(sc)
        contrast[["assay"]] = contrast[["assay"]] %>% 
            trimws()
        valid_assays = Seurat::Assays(sc)
        valid_reductions = Seurat::Reductions(sc)
        assertthat::assert_that(contrast[["assay"]] %in% valid_assays | 
                                    contrast[["assay"]] %in% valid_reductions | 
                                    contrast[["assay"]] %in% names(barcode_metadata),
                                msg=FormatString("The assay ('assay') must be one of the assays {valid_assays*}, one of the reductions {valid_reductions*} or barcode metadata columns (for comparison {i}/{name})."))
        assay = contrast[["assay"]]
        
        # Identify barcodes that are in assay/reduction. Result is a boolean vector that will be used later to get the correct barcodes.
        if (assay %in% valid_assays | assay %in% valid_reductions) {
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
        
        # bulk if bulk_by or pseudobulk_samples is set
        if ("bulk_by" %in% names(contrast) | "pseudobulk_samples" %in% names(contrast)) {
            contrast[["bulk_method"]] = dplyr::case_when(
                all(assay %in% valid_assays) ~ "aggregate",
                all(assay %in% valid_reductions) ~ "average",
                all(assay %in% names(barcode_metadata)) ~ "average",
                TRUE ~ "average")
        }
            
        # layer
        if (!"layer" %in% names(contrast)) {
            if (contrast[["assay"]] %in% valid_assays) {
                if (contrast[["test"]] %in% c("negbinom", "poisson", "DESeq2")) {
                    contrast[["layer"]] = "counts"
                } else{
                    contrast[["layer"]] = "data"
                }
            } else if(contrast[["assay"]] %in% valid_reductions) {
                contrast[["layer"]] = "counts"
            } else {
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
                strsplit(split=",") %>% 
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
    condition_group1 = rownames(barcode_metadata[condition_group1_idx, ])
    condition_group2_idx = contrast[["condition_group2_idx"]]
    condition_group2 = rownames(barcode_metadata[condition_group2_idx, ])
    
    # When creating a new Seurat, do not calculate nCounts and nFeatures (not needed); restore default on exit
    op = options(Seurat.object.assay.calcn = FALSE)
    on.exit(expr = options(op), add = TRUE)
        
    # Extract relevant data and save it as
    if (all(contrast[["assay"]] %in% Seurat::Assays(sc))) {
        # Get assay
        assay_obj = suppressWarnings({subset(sc[[contrast[["assay"]]]],
                           cells=unique(c(condition_group1_idx, condition_group2_idx)))})
        
        # Remove scale.data layer (never needed and saves a lot of memory)
        SeuratObject::LayerData(assay_obj, layer="scale.data") = NULL
    } else if (all(contrast[["assay"]] %in% Seurat::Reductions(sc))) {
        # Convert reduction into assay
        reduction = Seurat::Embeddings(sc, reduction=contrast[["assay"]])
        reduction = reduction[unique(c(condition_group1_idx, condition_group2_idx)), , drop=FALSE] %>% 
            as.matrix() %>%
            as("dgCMatrix") %>%
            Matrix::t()

        assay_obj = SeuratObject::CreateAssay5Object(data=reduction)
        assay_obj = SeuratObject::AddMetaData(assay_obj, metadata=data.frame(feature_name=rownames(reduction), row.names=rownames(assay_obj)))
    } else {
        # Use barcode metadata
        metadata_columns = contrast[["assay"]]
        metadata = sc[[data_columns]][unique(c(condition_group1_idx, condition_group2_idx)), , drop=FALSE] %>%
            as.matrix() %>%
            as("dgCMatrix") %>%
            Matrix::t()
        
        assay_obj = SeuratObject::CreateAssay5Object(data=metadata)
        assay_obj = SeuratObject::AddMetaData(assay_obj, metadata=data.frame(feature_name=rownames(metadata), row.names=rownames(assay_obj)))
    }
        
    # Get metadata columns
    metadata_columns = c(contrast[["condition_column"]], contrast[["subset_column"]], contrast[["bulk_by"]], contrast[["covariate"]])
    metadata_columns = metadata_columns[metadata_columns %in% colnames(barcode_metadata)]

    # Create subsetted Seurat object
    sc_subset = SeuratObject::CreateSeuratObject(assay_obj, meta.data=barcode_metadata[unique(c(condition_group1_idx, condition_group2_idx)), metadata_columns, drop=FALSE])
    
    # Add a condition column and update idents
    conditions = dplyr::case_when(
        SeuratObject::Cells(sc_subset) %in% condition_group1 ~ "condition1",
        SeuratObject::Cells(sc_subset) %in% condition_group2 ~ "condition2",
        TRUE ~ NA
    ) %>% factor(levels=c("condition1", "condition2"))
    
    assertthat::assert_that(!"condition_groups" %in% names(sc_subset[[]]),
                            msg=FormatString("The column 'condition_groups' is reserved, please use another column name (for comparison {i}/{name})."))
    sc_subset[["condition_groups"]] = conditions
    Seurat::Idents(sc_subset) = "condition_groups"
        
    # Bulk/pseudo-bulk if requested
    if ("bulk_by" %in% names(contrast) | "pseudobulk_samples" %in% names(contrast)) {
        categorial_covariates = contrast[["covariate"]] %>% purrr::keep(function(c) is.factor(barcode_metadata[, c]))
        numeric_covariates = contrast[["covariate"]] %>% purrr::keep(function(c) is.numeric(barcode_metadata[, c]))
        barcode_metadata = sc_subset[[]]
        
        if ("bulk_by" %in% names(contrast)) {
            # These columns will be aggregated
            bulk_by = unique(c("condition_groups", contrast[["bulk_by"]], categorial_covariates))
            
            # Downsample by groups (if requested by 'downsample_n')
            if ("downsample_n" %in% names(contrast)) {
                sampled_barcodes = barcode_metadata %>% 
                    dplyr::group_by(dplyr::across(dplyr::all_of(bulk_by))) %>% 
                    dplyr::slice_sample(n=contrast[["downsample_n"]]) %>% 
                    rownames()
                sc_subset = subset(sc_subset, cells=sampled_barcodes)
                barcode_metadata = sc_subset[[]]
            }
            
        } else if ("pseudobulk_samples" %in% names(contrast)) {
            # Group by condition column ('condition_groups') and categorial covariate columns
            # Then divide each group into pseudobulk_samples subgroups
            num_pseudobulk_samples = contrast[["pseudobulk_samples"]]
            bulk_by = unique(c("condition_groups", categorial_covariates))
            
            # Add a column for the pseudo-bulk sample
            barcode_metadata = barcode_metadata %>%
                dplyr::mutate(pseudobulk_sample=((1:dplyr::n()) %% num_pseudobulk_samples) + 1) %>%
                dplyr::mutate(pseudobulk_sample=paste0("s", pseudobulk_sample))
            barcode_metadata$pseudobulk_sample = factor(barcode_metadata$pseudobulk_sample, levels=paste0("s", 1:num_pseudobulk_samples))
            sc_subset = SeuratObject::AddMetaData(sc_subset, metadata=barcode_metadata$pseudobulk_sample, col.name="pseudobulk_sample")
            barcode_metadata = sc_subset[[]]
            
            # These columns will be aggregated
            bulk_by = unique(c("condition_groups", categorial_covariates, "pseudobulk_sample"))
            
            # Downsample by groups (if requested by 'downsample_n')
            if ("downsample_n" %in% names(contrast)) {
                bulk_by = unique(c("condition_groups", categorial_covariates), "pseudobulk_sample")
                sampled_barcodes = barcode_metadata %>% 
                    dplyr::group_by(dplyr::across(dplyr::all_of(bulk_by))) %>%
                    dplyr::slice_sample(n=contrast[["downsample_n"]]) %>%
                    rownames()
                    
                sc_subset = subset(sc_subset, cells=rownames(sampled_barcodes))
                barcode_metadata = sc_subset[[]]
            }
        }
        
        # Then aggregate expression
        sc_subset_agg = Seurat::PseudobulkExpression(object=sc_subset, 
                                                     return.seurat=TRUE, 
                                                     group.by=bulk_by,
                                                     layer="counts",
                                                     method=contrast[["bulk_method"]],
                                                     normalization.method="LogNormalize",
                                                     verbose=FALSE)
        for(c in bulk_by) {
            sc_subset_agg[[]][, c] = factor(sc_subset_agg[[]][, c], levels=levels(barcode_metadata[, c]))
        }
        Seurat::Idents(sc_subset_agg) = "condition_groups"
        
        # Numeric covariates need to be aggregated for each group by averaging, then added
        if (length(numeric_covariates) > 0) {
            numeric_covariate_data = barcode_metadata %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(bulk_by))) %>%
                dplyr::summarise(dplyr::across(dplyr::all_of(numeric_covariates), ~ mean(., na.rm=TRUE)))
            
            numeric_covariate_data = dplyr::left_join(sc_subset_agg[[]], numeric_covariate_data, by=bulk_by)
            rownames(numeric_covariate_data) = numeric_covariate_data$orig.ident
            sc_subset_agg = SeuratObject::AddMetaData(sc_subset_agg, metadata=numeric_covariate_data[, numeric_covariates, drop=FALSE])
        }
        
        # Remove scale.data layer (never needed and saves a lot of memory)
        SeuratObject::LayerData(sc_subset_agg, assay=SeuratObject::DefaultAssay(sc_subset_agg), layer="scale.data") = NULL
        
        sc_subset = sc_subset_agg
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
    name = contrast[["name"]]
    
    # Set up arguments list for FindMarkers
    arguments = list(object=contrast[["sc_subset"]],
                     test.use=contrast[["test"]], 
                     logfc.threshold=contrast[["log2FC"]], 
                     min.pct=contrast[["min_pct"]], 
                     ident.1="condition1",
                     ident.2="condition2",
                     slot=contrast[["layer"]],
                     random.seed=getOption("random_seed"),
                     min.cells.group=2)
    if ("bulk_by" %in% names(contrast)) arguments[["densify"]] = TRUE
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
        deg_results$gene = rownames(deg_results)
    } else {
        deg_results = data.frame(p_val=as.numeric(), avg_log2FC=as.numeric(), pct.1=as.numeric(), pct.2=as.numeric(), p_val_adj=as.numeric(), gene=as.character())
        contrast[["message"]] = FormatString("There are fewer than two samples in at least one group for comparison {name}.")
    }
    
    # Sort results
    deg_results = deg_results %>% 
        DegsSort()
    
    # Add means of normalised expression values
    avg_df = Seurat::AggregateExpression(object=contrast[["sc_subset"]], verbose=FALSE, return.seurat=TRUE) %>%
        SeuratObject::LayerData(layer="data") %>% 
        as.data.frame() %>%
        tibble::rownames_to_column(var="gene")
    #colnames(avg_df) = c("gene", contrast[["condition_group1"]], contrast[["condition_group2"]])
    colnames(avg_df) = c("gene", "condition1", "condition2")
    deg_results = dplyr::inner_join(deg_results, avg_df, by="gene")
    
    # If bulk_by or pseudobulk_samples is set, add the bulked expression values
    if ("bulk_by" %in% names(contrast) | "pseudobulk_samples" %in% names(contrast)) {
      bulk_df = SeuratObject::GetAssayData(contrast[["sc_subset"]], layer="data") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var="gene")
      deg_results = dplyr::inner_join(deg_results, bulk_df, by="gene")
    }
    
    # Remove Seurat object
    contrast[["sc_subset"]] = NULL
    
    # Add rownames
    rownames(deg_results) = deg_results$gene
    
    # Reorder the columns
    deg_results = deg_results %>% 
        dplyr::relocate(gene, p_val, p_val_adj, avg_log2FC, p_val_adj_score, pct.1, pct.2, condition1, condition2)
    
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
  gsea = purrr::map(term2gene, function(t2g) return(clusterProfiler::GSEA(geneList=fold_changes, TERM2GENE=t2g, minGSSize=10, maxGSSize=500)))
  
  return(gsea)
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

#' Creates a DEG scatterplot.
#' 
#' @param result A list entry with DEG results obtained with RunDEGTests
#' @param font_size The base font size. Default is 11.
#' @return A ggplot scatterplot
DegsScatterPlot = function(result, font_size=11) {
    # Get condition names. If multiple, join by '+'.
    group1 = paste(result[["condition_group1"]], collapse="+")
    group2 = paste(result[["condition_group2"]], collapse="+")
    
    # Get DEG table and add DEG status for significant up- and down-regulated (include log2 foldchange threshold)
    deg_table = result[["results"]]
    padj = result[["padj"]]
    log2FC = result[["log2FC"]]
    deg_table$deg_status = dplyr::case_when(
        deg_table$p_val_adj < padj & abs(deg_table$avg_log2FC) >= log2FC & deg_table$avg_log2FC > 0 ~ "up",
        deg_table$p_val_adj < padj & abs(deg_table$avg_log2FC) >= log2FC & deg_table$avg_log2FC < 0 ~ "down",
        TRUE ~ "none"
    )
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
#' @param result A list entry with DEG results obtained with RunDEGTests
#' @param font_size The base font size. Default is 11.
#' @return A ggplot volcano plot
DegsVolcanoPlot = function(result, font_size=11) {
    # Get condition names. If multiple, join by '+'.
    group1 = paste(result[["condition_group1"]], collapse="+")
    group2 = paste(result[["condition_group2"]], collapse="+")
    
    # Get DEG table and add DEG status for significant up- and down-regulated (include log2 foldchange threshold)
    deg_table = result[["results"]]
    padj = result[["padj"]]
    log2FC = result[["log2FC"]]
    deg_table$deg_status = dplyr::case_when(
        deg_table$p_val_adj < padj & abs(deg_table$avg_log2FC) >= log2FC & deg_table$avg_log2FC > 0 ~ "up",
        deg_table$p_val_adj < padj & abs(deg_table$avg_log2FC) >= log2FC & deg_table$avg_log2FC < 0 ~ "down",
        TRUE ~ "none"
    )
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
    p = ggplot(deg_table, aes(x=avg_log2FC, y=p_val_log10_n, col=deg_status)) + 
        geom_point() +
        ggrepel::geom_text_repel(data=top10_deg_table, aes(x=avg_log2FC, y=p_val_log10_n, col=deg_status, label=gene), size=font_size / .pt) +
        scale_color_manual("Gene status", values=c(none="grey", up="darkgoldenrod1", down="steelblue"), 
                           labels=c(none='none', up='up', down='down')) +
        AddPlotStyle(xlab="log2FoldChange", ylab="-log10(pvalue)", legend_position="none", font_size=font_size)
    
    # If there is a log2 threshold > 0, add vertical lines
    if (log2FC > 0) {
        p = p + geom_vline(xintercept=c(-log2FC, log2FC), linetype="dashed", color="grey")
    }
    
    return(p)
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
    readme_table = data.frame(Column=c("gene"), Description=c("Gene"))
    readme_table = rbind(readme_table, 
                         c("p_val", "Individual test p-value (uncorrect)"),
                         c("p_val_adj", "Adjusted p-value"),
                         c("avg_log2FC", "Mean log2 fold change condition 1 vs condition 2"),
                         c("pct.1", "Fraction cells expressing gene in condition 1"),
                         c("pct.2", "Fraction cells expressing gene in condition 2"),
                         c("condition1", "Average normalized expression in condition 1"),
                         c("condition2", "Average normalized expression in condition 2"),
                         c("cluster", "Cluster (not always applicable)"),
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
