#' Sets up a list with contrasts for analysis. Makes sure that all required information is available. Sets default and does basic checks. 
#' Each contrast entry is a list either of length 1 (if there is only one comparison) or of length >1 (if the comparison is done for multiple subsets).
#' 
#' @param sc A Seurat single cell object.
#' @param contrasts_list A list of contrasts. Must at least contain 'name', condition_column', 'condition1' and 'condition2'.
#' @param type Type of contrast. Can be 'deg' (for changes in gene expression) or 'compositional' (for changes in cell composition). Default is 'deg'.
#' @return A list with contrasts. Each list entry is a list either of length 1 (if there is only one comparison per contrast) or of length >1 
#' (if comparisons are done for multiple subsets for a contrast)
NewContrastsList = function(sc, contrasts_list, type='deg') {
    # If empty, return empty list
    if (length(contrasts_list) == 0) return(list())
  
    # Get barcode metadata
    barcode_metadata = sc[[]]
    
    # Convert into list, do checks and set up defaults
    contrasts_list = purrr::map(seq(contrasts_list), function(i) {
        contrast = contrasts_list[[i]]
        
        # Trim whitespace
        contrast = purrr::map_depth(contrast, .depth=-1, trimws)
        
        # name: should only contain alphanumeric characters, underscores and hyphens
        assertthat::assert_that("name" %in% names(contrast), 
                                msg=FormatString("The name ('name') is missing (for comparison {i})."))
        contrast[["name"]] = name = unlist(contrast[["name"]])
        assertthat::assert_that(grepl("^[[:alnum:]_\\-]+$", name),
                                msg=FormatString("The name ('name') must only contain alphanumeric characters, underscores and hyphens (for comparison {i})."))
        
        # assay: data that is tested (not relevant for compositional contrasts)
        if (!"assay" %in% names(contrast)) contrast[["assay"]] = Seurat::DefaultAssay(sc)
        contrast[["assay"]] = unlist(contrast[["assay"]])
        valid_assays = Seurat::Assays(sc)
        valid_reductions = Seurat::Reductions(sc)
        assertthat::assert_that(length(contrast[["assay"]]) == 1 && contrast[["assay"]] %in% c(valid_assays, valid_reductions) | 
                                    all(contrast[["assay"]] %in% names(barcode_metadata)),
                                msg=FormatString("The assay ('assay') must be one of the assays {valid_assays*}, one of the reductions {valid_reductions*} or barcode metadata columns (for comparison {i}/{name})."))
        assay = contrast[["assay"]]
        
        # data_type: which type of data is tested
        if (type == "compositional") {
          # cell composition
          contrast[["data_type"]] = "compositional"
        } else if (length(assay) == 1 && assay %in% valid_assays) {
          # feature data
          contrast[["data_type"]] = "feature_data"
        } else if (length(assay) == 1 && assay %in% valid_reductions) {
          # reduction
          contrast[["data_type"]] = "reduction"
        } else {
          # barcode metadata
          contrast[["data_type"]] = "barcode_metadata"
        }
        
        # Identify barcodes that are in assay/reduction. Result is a boolean vector that will be used later to get the correct barcodes.
        if (contrast[["data_type"]] %in% c("feature_data", "reduction")) {
            barcodes_in_assay = rownames(barcode_metadata) %in% SeuratObject::Cells(sc[[assay]])
        } else {
            barcodes_in_assay = rep(TRUE, nrow(barcode_metadata))
        }
        
        # condition_column: barcode metadata column that contains the condition information
        assertthat::assert_that("condition_column" %in% names(contrast),
                                msg=FormatString("The condition column ('condition_column') is missing (for comparison {i}/{name})."))
        contrast[["condition_column"]] = condition_column = unlist(contrast[["condition_column"]])
        
        # condition1 and condition2: condition 1 vs condition 2
        assertthat::assert_that("condition1" %in% names(contrast),
                                msg=FormatString("The condition group 1 ('condition1') is missing or empty (for comparison {i}/{name})."))
        contrast[["condition1"]] = unlist(contrast[["condition1"]])
        
        assertthat::assert_that("condition2" %in% names(contrast),
                                msg=FormatString("The condition group 2 ('condition1') is missing or empty (for comparison {i}/{name})."))
        contrast[["condition2"]] = unlist(contrast[["condition2"]])
        
        # If at least one of the condition groups is not a file, check if condition_column is part of the barcode metadata
        if (!file.exists(contrast[["condition1"]]) | !file.exists(contrast[["condition2"]])) {
            assertthat::assert_that(condition_column %in% colnames(barcode_metadata),
                                    msg=FormatString("The condition column ('condition_column') must be part of the barcode metadata if at least one of the condition groups is not a file (for comparison {i}/{name})."))
            
            assertthat::assert_that(is.factor(barcode_metadata[, condition_column]),
                                    msg=FormatString("The condition column ('condition_column') of the barcode metadata must be a factor (for comparison {i}/{name})."))
        }
        
        if (!file.exists(contrast[["condition1"]])) {
            condition1 = contrast[["condition1"]]
            
            # Parse condition1 string
            negate = grepl("^!", condition1)
            condition1 = gsub("^!", "", condition1) %>%
                strsplit(split="\\+") %>%
                unlist() %>%
                trimws() %>%
                unique()
            
            # If negate (!), get complement of condition1
            if (negate) {
                condition1 = setdiff(levels(barcode_metadata[, condition_column]), condition1)
            }
            
            # Make sure all levels are valid
            purrr::walk(condition1, function(c) {
                assertthat::assert_that(c %in% levels(barcode_metadata[, condition_column]),
                                        msg=FormatString("The condition group 1 value {c} is not level of the condition column {condition_column} of the barcode metadata (for comparison {i}/{name})."))
            })
            
            contrast[["condition1"]] = condition1
            
            # Get indices of condition1
            contrast[["condition1_idx"]] = which(barcodes_in_assay & barcode_metadata[, condition_column] %in% condition1)
        } else {
            condition1 = contrast[["condition1"]]
            
            # Sheet number/name appended?
            sheet = 1
            if (grepl(":[^:]+$", condition1)) {
                sheet = gsub(pattern=".+:([^:]+)$", replacement="\\1", x=condition1) %>% as.integer()
                condition1 = gsub(pattern=":[^:]+$", replacement="", x=condition1)
            }
            
            # Make sure file exists
            assertthat::assert_that(file.exists(condition1),
                                    msg=FormatString("The condition group 1 barcode file {condition1} does not exist (for comparison {i}/{name})."))
            
            # Decide whether it is a valid file type
            extension = tools::file_ext(gsub(pattern="\\.gz$", replacement="", x=condition1))
            valid_extensions = c("csv", "tsv", "xls", "xlsx")
            assertthat::assert_that(extension %in% valid_extensions,
                                    msg=FormatString("The condition group 1 barcode file must be one of: {valid_extensions*} (file can be gzipped) (for comparison {i}/{name})."))
            
            # Read file
            if (extension %in% c("csv", "tsv")) {
                barcodes = readr::read_tsv(condition1)
            } else if (extension %in% c("xls", "xlsx")) {
                barcodes = readxl::read_excel(condition1, sheet=sheet)
            }
            
            # Make sure file contains data
            assertthat::assert_that(is.data.frame(barcodes) && nrow(barcodes)>0,
                                    msg=FormatString("The condition group 1 barcode file does not contain barcodes (for comparison {i}/{name})."))
            
            # File can contain one column (barcodes) or two columns (sample name and original barcode)
            if (ncol(barcodes) == 2) {
                # Get sample names and associated original barcodes
                sample_names = barcodes[, 1, drop=TRUE] %>% trimws() %>% unique()
                original_barcodes = barcodes[, 2, drop=TRUE] %>% trimws() %>% unique()
                i = which(!is.na(original_barcodes) & !is.na(sample_names))
                sample_names = sample_names[i]
                original_barcodes = original_barcodes[i]
                
                # Assert that all sample names are valid
                assertthat::assert_that(all(sample_names %in% levels(barcode_metadata$orig.ident)),
                                        msg=FormatString("The sample names in the condition group 1 barcode file are not levels of the 'orig.ident' column of the barcode metadata (for comparison {i}/{name})."))
                
                # Find barcodes indices for which sample name and original barcode match
                contrast[["condition1_idx"]] = which(barcode_metadata$orig.ident %in% sample_names & barcode_metadata$orig_barcode %in% original_barcodes)
            } else {
                barcodes = barcodes[, 1, drop=TRUE] %>% trimws() %>% unique()
                barcodes = barcodes[!is.na(barcodes)]
                
                # Find barcodes indices
                contrast[["condition1_idx"]] = which(rownames(barcode_metadata) %in% barcodes)
            }
            
            # Add base file name:sheet number or sheet name as 'condition1'
            contrast[["condition1"]] = paste0(basename(contrast[["condition1"]]), ", sheet ", sheet)
        }
        
        # Make sure condition1_idx is not empty
        assertthat::assert_that(length(contrast[["condition1_idx"]]) > 0,
                                msg=FormatString("The condition group 1 does not contain any barcodes (for comparison {i}/{name})."))
        
        if (!file.exists(contrast[["condition2"]])) {
            condition2 = contrast[["condition2"]]
            
            # Parse condition2 string
            negate = grepl("^!", condition2)
            condition2 = gsub("^!", "", condition2) %>%
                strsplit(split="\\+") %>%
                unlist() %>%
                trimws() %>%
                unique()
            
            # If negate (!), get complement of condition2
            if (negate) {
                condition2 = setdiff(levels(barcode_metadata[, condition_column]), condition2)
            }
            
            # Make sure all levels are valid
            purrr::walk(condition2, function(c) {
                assertthat::assert_that(c %in% levels(barcode_metadata[, condition_column]),
                                        msg=FormatString("The condition group 2 value {c} is not level of the condition column {condition_column} of the barcode metadata (for comparison {i}/{name})."))
            })
            
            contrast[["condition2"]] = condition2
            
            # Get indices of condition2
            contrast[["condition2_idx"]] = which(barcodes_in_assay & barcode_metadata[, condition_column] %in% condition2)
        } else {
            condition2 = contrast[["condition2"]]
            
            # Sheet number appended?
            sheet = 1
            if (grepl(":\\d+$", condition2)) {
                sheet = gsub(pattern=".+:(\\d+)$", replacement="\\1", x=condition2) %>% as.integer()
                condition2 = gsub(pattern=":\\d+$", replacement="", x=condition2)
            }
            
            # Make sure file exists
            assertthat::assert_that(file.exists(condition2),
                                    msg=FormatString("The condition group 2 barcode file {condition2} does not exist (for comparison {i}/{name})."))
            
            # Decide whether it is a valid file type
            extension = tools::file_ext(gsub(pattern="\\.gz$", replacement="", x=condition2))
            valid_extensions = c("csv", "tsv", "xls", "xlsx")
            assertthat::assert_that(extension %in% valid_extensions,
                                    msg=FormatString("The condition group 2 barcode file must be one of: {valid_extensions*} (file can be gzipped) (for comparison {i}/{name})."))
            
            # Read file
            if (extension %in% c("csv", "tsv")) {
                barcodes = readr::read_tsv(condition2)
            } else if (extension %in% c("xls", "xlsx")) {
                barcodes = readxl::read_excel(condition2, sheet=sheet)
            }
            
            # Make sure file contains data
            assertthat::assert_that(is.data.frame(barcodes) && nrow(barcodes)>0,
                                    msg=FormatString("The condition group 2 barcode file does not contain barcodes (for comparison {i}/{name})."))
            
            # File can contain one column (barcodes) or two columns (sample name and original barcode)
            if (ncol(barcodes) == 2) {
              # Get sample names and associated original barcodes
              sample_names = barcodes[, 1, drop=TRUE] %>% trimws() %>% unique()
              original_barcodes = barcodes[, 2, drop=TRUE] %>% trimws() %>% unique()
              i = which(!is.na(original_barcodes) & !is.na(sample_names))
              sample_names = sample_names[i]
              original_barcodes = original_barcodes[i]
              
              # Assert that all sample names are valid
              assertthat::assert_that(all(sample_names %in% levels(barcode_metadata$orig.ident)),
                                      msg=FormatString("The sample names in the condition group 2 barcode file are not levels of the 'orig.ident' column of the barcode metadata (for comparison {i}/{name})."))
              
              # Find barcodes indices for which sample name and original barcode match
              contrast[["condition2_idx"]] = which(barcode_metadata$orig.ident %in% sample_names & barcode_metadata$orig_barcode %in% original_barcodes)
            } else {
              barcodes = barcodes[, 1, drop=TRUE] %>% trimws() %>% unique()
              barcodes = barcodes[!is.na(barcodes)]
              
              # Find barcodes indices
              contrast[["condition2_idx"]] = which(rownames(barcode_metadata) %in% barcodes)
            }
            
            # Add base file name:sheet number or sheet name as 'condition1'
            contrast[["condition1"]] = paste0(basename(contrast[["condition1"]]), ", sheet ", sheet)
        }
        
        # Make sure condition2_idx is not empty
        assertthat::assert_that(length(contrast[["condition2_idx"]]) > 0,
                                msg=FormatString("The condition group 2 does not contain any barcodes (for comparison {i}/{name})."))
        
        # subset_column and subset_group: subset_column is a barcode metadata column that contains the subset information, subset_group is a string that contains the subset group(s)
        if ("subset_column" %in% names(contrast)) {
            assertthat::assert_that("subset_group" %in% names(contrast),
                                    msg=FormatString("The 'subset_group' column must be used together with the subset_column column (for comparison {i}/{name})."))
            subset_column = unlist(contrast[["subset_column"]])
            subset_group = unlist(contrast[["subset_group"]])
            
            if (!file.exists(subset_group)) {
                # Parse subset_group string
                subset_group = subset_group %>% 
                    strsplit(split=",") %>% 
                    unlist() %>% 
                    trimws() %>% 
                    unique()
                
                # Allow wilcard
                if (subset_group[1] == "*") {
                    subset_group = levels(barcode_metadata[, subset_column])
                }
                
                # Make sure all subset groups are valid
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
                
                # File can contain one column (barcodes) or two columns (sample name and original barcode)
                if (ncol(barcodes) == 2) {
                  # Get sample names and associated original barcodes
                  sample_names = barcodes[, 1, drop=TRUE] %>% trimws() %>% unique()
                  original_barcodes = barcodes[, 2, drop=TRUE] %>% trimws() %>% unique()
                  i = which(!is.na(original_barcodes) & !is.na(sample_names))
                  sample_names = sample_names[i]
                  original_barcodes = original_barcodes[i]
                  
                  # Assert that all sample names are valid
                  assertthat::assert_that(all(sample_names %in% levels(barcode_metadata$orig.ident)),
                                          msg=FormatString("The sample names in the subset group barcode file are not levels of the 'orig.ident' column of the barcode metadata (for comparison {i}/{name})."))
                  
                  # Find barcodes indices for which sample name and original barcode match
                  contrast[["subset_group_idx"]] = which(barcode_metadata$orig.ident %in% sample_names & barcode_metadata$orig_barcode %in% original_barcodes)
                } else {
                  barcodes = barcodes[, 1, drop=TRUE] %>% trimws() %>% unique()
                  barcodes = barcodes[!is.na(barcodes)]
                  
                  # Find barcodes indices
                  contrast[["subset_group_idx"]] = which(rownames(barcode_metadata) %in% barcodes)
                }
                
                # Add base file name:sheet number or sheet name as 'condition1'
                contrast[["subset_group"]] = paste0(basename(contrast[["subset_group"]]), ", sheet ", sheet)
            }
            
            # Make sure subset_group_idx is not empty
            assertthat::assert_that(length(contrast[["subset_group_idx"]]) > 0,
                                    msg=FormatString("The subset group does not contain any barcodes (for comparison {i}/{name})."))
        }
        
        # compositional_column: Column that contains the cell compositional information. By default: seurat_clusters.
        # (not relevant for deg analyses)
        if (type == "compositional") {
          if (!"compositional_column" %in% names(contrast)) contrast[["compositional_column"]] = "seurat_clusters"
          contrast[["compositional_column"]] = compositional_column = unlist(contrast[["compositional_column"]])
          assertthat::assert_that(compositional_column %in% colnames(barcode_metadata),
                                  msg=FormatString("The compositional column {compositional_column} was not found in the barcode metadata (for comparison {i}/{name})."))
        }
        
        # samples_column: Column that contains the samples (or datasets) information. By default: sample.
        # (not relevant for deg analyses)
        if (type == "compositional") {
          if (!"samples_column" %in% names(contrast)) contrast[["samples_column"]] = "sample"
          contrast[["samples_column"]] = samples_column = unlist(contrast[["samples_column"]])
          assertthat::assert_that(samples_column %in% colnames(barcode_metadata),
                                  msg=FormatString("The samples column {samples_column} was not found in the barcode metadata (for comparison {i}/{name})."))
        }
        
        # bulk_by: bulk data by a barcode metadata column
        # (ignored for compositional contrasts)
        if ("bulk_by" %in% names(contrast) & type != "compositional") {
            contrast[["bulk_by"]] = unlist(contrast[["bulk_by"]])
          
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
        
        # pseudobulk_samples: create pseudo-bulk samples
        if ("pseudobulk_samples" %in% names(contrast)) {
            contrast[["pseudobulk_samples"]] = unlist(contrast[["pseudobulk_samples"]])
            contrast[["pseudobulk_samples"]] = as.integer(contrast[["pseudobulk_samples"]])
            assertthat::assert_that(contrast[["pseudobulk_samples"]] > 1,
                                    msg=FormatString("The number of pseudobulk samples ('pseudobulk_samples') must be greater than 1 (for comparison {i}/{name})."))
        }
        
        # test: test to use
        if (type == "deg") {
          # valid deg tests
          valid_tests = c("wilcox", "bimod", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2", "DESeq2LRT")
          if (!"test" %in% names(contrast)) contrast[["test"]] = "wilcox"
        } else if (type == "compositional") {
          # valid compositional tests
          valid_tests = c("none", "fisher", "sccoda")
          if (!"test" %in% names(contrast)) contrast[["test"]] = "fisher"
        }
        contrast[["test"]] = unlist(contrast[["test"]])
        assertthat::assert_that(contrast[["test"]] %in% valid_tests,
                                msg=FormatString("The test ('test') must be one of: {valid_tests*} (for comparison {i}/{name})."))
        if (contrast[["test"]] %in% c("DESeq2", "DESeq2LRT")) {
          assertthat::assert_that("pseudobulk_samples" %in% names(contrast) | "bulk_by" %in% names(contrast),
                                  msg=FormatString("The DESeq2 tests can only be used together with 'bulk_by' and/or 'pseudobulk_samples'. (for comparison {i}/{name})."))
        }
        
        # padj
        if (!"padj" %in% names(contrast)) contrast[["padj"]] = 0.05
        contrast[["padj"]] = unlist(contrast[["padj"]])
        contrast[["padj"]] = as.numeric(contrast[["padj"]])
        
        # log2FC 
        # (not relevant for compositional contrasts)
        if (!"log2FC" %in% names(contrast)) contrast[["log2FC"]] = 0
        contrast[["log2FC"]] = unlist(contrast[["log2FC"]])
        contrast[["log2FC"]] = as.numeric(contrast[["log2FC"]])
        
        # min_pct 
        # (not relevant for compositional contrasts)
        if (!"min_pct" %in% names(contrast)) contrast[["min_pct"]] = 0.01
        contrast[["min_pct"]] = unlist(contrast[["min_pct"]])
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
        # (not relevant for compositional contrasts)
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
        contrast[["layer"]] = unlist(contrast[["layer"]])
        assertthat::assert_that(contrast[["layer"]] %in% c("counts", "data", "scale.data"),
                                msg=FormatString("The layer ('layer') must be 'counts', 'data', 'scale.data' (for comparison {i}/{name})."))

        # downsample_barcodes (ignored for compositional contrasts)
        if ("downsample_barcodes" %in% names(contrast) & type != "compositional") {
            contrast[["downsample_barcodes"]] = unlist(contrast[["downsample_barcodes"]])
            contrast[["downsample_barcodes"]] = as.integer(contrast[["downsample_barcodes"]])
        }
        
        # design
        if ("design" %in% names(contrast)) {
          # Parse design formula string
          contrast[["design"]] = design = contrast[["design"]] %>%
            unlist() %>% 
            trimws() %>% 
            unique()
          
          # Can be one or two design formulas
          assertthat::assert_that(length(design) <= 2,
                                  msg=FormatString("The parameter 'design' can specify at most two design formula: one for the full model and - when doing a LRT test - a second one needs for the reduced model (for comparison {i}/{name})."))
          # Convert to formula
          design = purrr::map(design, as.formula)
          
          # Check that at least one design formula contains the variable 'condition_groups'
          variables = purrr::map(design, all.vars) %>% unlist()
          assertthat::assert_that("condition_groups" %in% variables,
                                  msg=FormatString("One design formula must contain the variable 'condition_groups' (for comparison {i}/{name})."))
          
          # Overwrite the covariate setting
          covariate = setdiff(variables, "condition_groups")
          contrast[["covariate"]] = covariate

          contrast[["design"]] = design
        }
        
        # covariate
        if ("covariate" %in% names(contrast)) {
          contrast[["covariate"]] = covariate = contrast[["covariate"]] %>% 
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
            valid_tests = c("LR", "negbinom", "poisson", "MAST", "DESeq2", "DESeq2LRT", "sccoda")
            assertthat::assert_that(test %in% valid_tests,
                                    msg=FormatString("Only the following tests allow covariates: {valid_tests*}."))
            
            categorial_covariates = covariate %>% 
              purrr::keep(function(c) is.factor(barcode_metadata[, c, drop=TRUE]))
            valid_tests = c("LR", "negbinom", "poisson", "MAST", "DESeq2", "DESeq2LRT", "sccoda")
            assertthat::assert_that(test %in% valid_tests | length(categorial_covariates) == 0,
                                    msg=FormatString("Only the following tests allow character/categorial covariates: {valid_tests*}."))
            
            numeric_covariates = covariate %>% 
              purrr::keep(function(c) is.numeric(barcode_metadata[, c, drop=TRUE]))
            valid_tests = c("LR", "negbinom", "poisson", "MAST")
            assertthat::assert_that(test %in% valid_tests | length(numeric_covariates) == 0,
                                    msg=FormatString("Only the following tests allow numeric covariates: {valid_tests*}."))
            
            contrast[["covariate"]] = covariate 
        }
        
        # DESeq2 additional arguments to results function
        if ("deseq2_results_args" %in% names(contrast) & type != "compositional") {
          deseq2_args = purrr::map(contrast[["deseq2_results_args"]], function(p) {
            p = as.character(p)
            return(eval(parse(text=p)))
          })
          
          contrast[["deseq2_results_args"]] = deseq2_args
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
              
              # Filter barcodes (condition1_idx and condition2_idx) so that they are in the subset group
              subset_group_barcodes = rownames(barcode_metadata[con[["subset_group_idx"]], ]) %>% unique()
              
              k = rownames(barcode_metadata[contrast[["condition1_idx"]], ]) %in% subset_group_barcodes
              con[["condition1_idx"]] = contrast[["condition1_idx"]][k]
              
              k = rownames(barcode_metadata[contrast[["condition2_idx"]], ]) %in% subset_group_barcodes
              con[["condition2_idx"]] = contrast[["condition2_idx"]][k]
              return(con)
          })
          names(contrasts_expanded) = purrr::map(contrast[["subset_group"]], paste, collapse="+") %>% unlist()
          
          # Drop subset contrasts where condition1_idx or condition2_idx is empty
          contrasts_expanded = purrr::keep(contrasts_expanded, function(c) {
              return(length(c[["condition1_idx"]]) > 0 & length(c[["condition2_idx"]]) > 0)
          })
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

#' Given a contrast configuration, prepares a Seurat object for DEG analysis. The function extracts all relevant 
#' data and if requested bulk-aggregates and downsamples barcodes. Note that if a contrast has multiple comparisons (for each subset), this function needs to be run on each of them separately.
#' 
#' @param sc A Seurat single cell object.
#' @param contrast A contrast configuration. Must have been set up with NewContrastsList.
#' @return A contrast configuration with an 'object' entry for the Seurat object.
PrepareDegContrast = function(sc, contrast) {
    barcode_metadata = sc[[]]
    name = contrast[["name"]]
        
    # Get barcodes indices and names for conditions
    # Note: If the comparison is for a subset, the barcodes are already filtered for this subset. Only this data will be extracted.
    condition1_idx = contrast[["condition1_idx"]]
    condition1_barcodes = SeuratObject::Cells(sc)[condition1_idx]
    condition2_idx = contrast[["condition2_idx"]]
    condition2_barcodes = SeuratObject::Cells(sc)[condition2_idx]

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
        
        # Get relevant barcodes
        barcodes = SeuratObject::Cells(sc)
        barcodes_idx = unique(c(condition1_idx, condition2_idx))
        barcodes = barcodes[barcodes_idx]
        
        # Create new assay object (cannot have Seurat without one)
        assay_obj = suppressWarnings({subset(sc[[assay]],
                           cells=barcodes)})
        
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
                                              barcode_metadata[barcodes, barcode_metadata_columns, drop=FALSE])
    } else if (data_type == "reduction") {
        reduction_name = contrast[["assay"]]
        assay = "reduction"
        
        # Get relevant barcodes
        barcodes = SeuratObject::Cells(sc)
        barcodes_idx = unique(c(condition1_idx, condition2_idx))
        barcodes = barcodes[barcodes_idx]

        # Get reduction and subset
        reduction = sc[[reduction_name]]
        reduction = subset(reduction, cells=barcodes)

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
        sc_subset = SeuratObject::CreateSeuratObject(assay_obj, assay=assay, meta.data=barcode_metadata[barcodes, barcode_metadata_columns, drop=FALSE])
    } else if (data_type == "barcode_metadata") {
        # Get relevant barcodes
        barcodes = SeuratObject::Cells(sc)
        barcodes_idx = unique(c(condition1_idx, condition2_idx))
        barcodes = barcodes[barcodes_idx]
      
        # Get barcode metadata columns to test
        metadata_cols = contrast[["assay"]]
        assay = "barcode_metadata"
        metadata = sc[[metadata_cols]][barcodes, , drop=FALSE] %>%
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
        sc_subset = SeuratObject::CreateSeuratObject(assay_obj, assay=assay, meta.data=barcode_metadata[barcodes, barcode_metadata_columns, drop=FALSE])
    }
    
    # Add the condition_groups column with levels
    assertthat::assert_that(!"condition_groups" %in% colnames(sc_subset[[]]),
                            msg=FormatString("The column 'condition_groups' is reserved, please use another column name (for comparison {i}/{name})."))
    conditions = dplyr::case_when(
        SeuratObject::Cells(sc_subset) %in% condition1_barcodes ~ "condition1",
        SeuratObject::Cells(sc_subset) %in% condition2_barcodes ~ "condition2",
        TRUE ~ NA
    )
    conditions = factor(conditions, levels=c("condition1", "condition2"))
    sc_subset[["condition_groups"]] = conditions
    
    # Update idents
    sc_subset@active.ident = sc_subset$condition_groups
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
        sc_subset_agg@active.ident = sc_subset_agg$condition_groups
        
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

#' Given a contrast configuration, prepares a Seurat object for compositional analysis.
#' 
#' @param sc A Seurat single cell object.
#' @param contrast A contrast configuration. Must have been set up with NewContrastsList.
#' @return A contrast configuration with an 'object' entry for the Seurat object.
PrepareCompositionalContrast = function(sc, contrast) {
  barcode_metadata = sc[[]]
  name = contrast[["name"]]
  
  # If condition1_idx or condition2_idx is empty, return
  if (length(contrast[["condition1_idx"]]) == 0 | length(contrast[["condition2_idx"]]) == 0) return(contrast)
  
  # Get barcodes indices and names
  # If the comparison is for a subset, the barcodes are already filtered for this subset.
  condition1_idx = contrast[["condition1_idx"]]
  condition1_barcodes = SeuratObject::Cells(sc)[condition1_idx]
  condition2_idx = contrast[["condition2_idx"]]
  condition2_barcodes = SeuratObject::Cells(sc)[condition2_idx]
  
  # When creating a new Seurat, do not calculate nCounts and nFeatures (not needed); restore default on exit
  op = options(Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  
  # Get barcode metadata columns
  barcode_metadata_columns = c(contrast[["condition_column"]], 
                               contrast[["subset_column"]], 
                               contrast[["compositional_column"]], 
                               contrast[["covariate"]])
  barcode_metadata_columns = barcode_metadata_columns[barcode_metadata_columns %in% colnames(barcode_metadata)]
  
  # Extract the data and create a Seurat object. This depends on the test:
  # - none, fisher, sccoda: create minimal Seurat object to store the barcode metadata
  # - miloR: Seurat object with counts data, barcode metadata, dimensionality reduction and graph
  # - if the test is done on a subset, restrict to barcodes of this subset
  # - note: changes of abundance in one condition vs the other also depend on changes of abundances in all other conditions
  #   therefore include other conditions that are not part of the data as condition3
  
  # Get relevant barcodes
  barcodes = SeuratObject::Cells(sc)
  if ("subset_group_idx" %in% names(contrast)) {
    barcodes = barcodes[contrast[["subset_group_idx"]]]
  }
  
  # Create new assay object (cannot have Seurat without one)
  assay = contrast[["assay"]]
  
  if (contrast[["test"]] %in% c("none", "fisher", "sccoda")) {
    # Create minimal counts table without data
    counts = Matrix::Matrix(0, nrow=2, ncol=length(barcodes), sparse=TRUE, 
                            dimnames=list(c("-dummy1-", "-dummy2-"), barcodes)) %>% 
      as("dgCMatrix")
    assay_obj = SeuratObject::CreateAssay5Object(counts)
  } else if (contrast[["test"]] %in% c("miloR")) {
    # Get counts from assay
    assay_obj = suppressWarnings({subset(sc[[assay]],
                                         cells=barcodes)})
    
    # Update feature metadata to include only feature_id, feature_name, feature_type
    feature_metadata = assay_obj[[]]
    feature_metadata = feature_metadata[, c("feature_id", "feature_name", "feature_type")]
    assay_obj@meta.data = data.frame()
    assay_obj = SeuratObject::AddMetaData(assay_obj, feature_metadata)
    
    # Remove scale.data layer (never needed and saves a lot of memory)
    SeuratObject::DefaultLayer(assay_obj) = contrast[["layer"]]
    SeuratObject::LayerData(assay_obj, layer="scale.data") = NULL
  }

  # Store Seurat object in sc_subset
  sc_subset = SeuratObject::CreateSeuratObject(assay_obj, 
                                               assay=assay)
  # Update barcode metadata
  sc_subset = SeuratObject::AddMetaData(sc_subset,
                                        barcode_metadata[barcodes, barcode_metadata_columns, drop=FALSE])
  
  # For miloR test: add dimensionality reduction and graph TODODODO
  
  # Add the condition_groups column with levels
  # Note: condition3 is used for all barcodes not included in the analysis
  assertthat::assert_that(!"condition_groups" %in% colnames(sc_subset[[]]),
                          msg=FormatString("The column 'condition_groups' is reserved, please use another column name (for comparison {i}/{name})."))
  conditions = dplyr::case_when(
    SeuratObject::Cells(sc_subset) %in% condition1_barcodes ~ "condition1",
    SeuratObject::Cells(sc_subset) %in% condition2_barcodes ~ "condition2",
    .default = "condition3"
  )
  conditions = factor(conditions, levels=c("condition1", "condition2", "condition3"))
  sc_subset$condition_groups = conditions
  
  # Update idents
  SeuratObject::Idents(sc_subset) = "condition_groups"

  barcode_metadata = sc_subset[[]]
  
  # Group cells into pseudo-bulk samples if requested
  if ("pseudobulk_samples" %in% names(contrast)) {
    # By default at least group by conditions
    bulk_by = "condition_groups"
    
    # Get lists of categorial and numeric covariates
    # Also aggregate by categorial covariates
    categorial_covariates = contrast[["covariate"]] %>% 
      purrr::keep(function(c) is.factor(barcode_metadata[, c, drop=TRUE]))
    numeric_covariates = contrast[["covariate"]] %>% 
      purrr::keep(function(c) is.numeric(barcode_metadata[, c, drop=TRUE]))
    if (length(categorial_covariates) > 0) bulk_by = unique(c(bulk_by, categorial_covariates))
    
    # Generate x pseudo-bulk samples per bulk_by combination (to group)
    num_pseudobulk_samples = contrast[["pseudobulk_samples"]]
      
    # Divide each group into pseudobulk_samples subgroups
    # Done with modulo (%%): ((row_index - 1) %% num_pseudobulk_samples) + 1
    pseudobulk_sample = barcode_metadata %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(bulk_by))) %>% 
      dplyr::mutate(pseudobulk_sample=0:(dplyr::n() - 1)) %>%
      dplyr::mutate(pseudobulk_sample=(pseudobulk_sample %% num_pseudobulk_samples) + 1) %>%
      dplyr::mutate(pseudobulk_sample=paste(!!!syms(bulk_by), paste0("s", pseudobulk_sample), sep="-")) %>%
      dplyr::pull(pseudobulk_sample)
      
    # Add pseudobulk_sample column to seurat object
    barcode_metadata$pseudobulk_sample = pseudobulk_sample
    sc_subset = SeuratObject::AddMetaData(sc_subset, metadata=barcode_metadata[,"pseudobulk_sample", drop=FALSE])
    barcode_metadata = sc_subset[[]]
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
    
    # If data is stored as IterableMatrix (on-disk via BPCells), convert to in-memory dgCMatrix unless test is wilcox
    # This is because all other tests in FindMarkers do not support IterableMatrix
    if (contrast[["test"]] != "wilcox") {
      if ("counts" %in% SeuratObject::Layers(contrast[["sc_subset"]])) {
          SeuratObject::LayerData(contrast[["sc_subset"]], layer="counts") = SeuratObject::LayerData(contrast[["sc_subset"]], layer="counts") %>% 
        as("dgCMatrix")
      }
      if ("data" %in% SeuratObject::Layers(contrast[["sc_subset"]])) {
        SeuratObject::LayerData(contrast[["sc_subset"]], layer="data") = SeuratObject::LayerData(contrast[["sc_subset"]], layer="data") %>% 
          as("dgCMatrix")
      }
    }
    
    # Check that there are enough samples in each group (>=2), 
    # if not, add an error message and skip testing, return empty results
    condition_counts = SeuratObject::Idents(contrast[["sc_subset"]]) %>% table()
    if (min(condition_counts) >= 2) {
        # Run test
        if (contrast[["test"]] %in% c("DESeq2", "DESeq2LRT")) {
            # design and reduced argument
            design = reduced = NULL
            if ("design" %in% names(contrast)) {
                design = contrast[["design"]][1] %>% as.character()
              
              if (contrast[["test"]] == "DESeq2LRT" & length(contrast[["design"]]) == 2) {
                reduced = contrast[["design"]][2] %>% as.character()
              }
            }        
  
            # Run DegsRunDESeq2
            deg_results = DegsRunDESeq2(contrast[["sc_subset"]],
                                        assay=SeuratObject::DefaultAssay(contrast[["sc_subset"]]), 
                                        ident_1="condition1",
                                        ident_2="condition2",
                                        test=dplyr::case_match(contrast[["test"]],
                                                               "DESeq2" ~ "Wald",
                                                               "DESeq2LRT" ~ "LRT"), 
                                        design=design, 
                                        reduced=reduced,
                                        random_seed=getOption("random_seed"),
                                        results_args=contrast[["deseq2_results_args"]])
        } else {
            # How should FindMarkers calculate the means:
            # - for feature/counts data (gene expression) use default (difference in the log of the average exponentiated data, with pseudocount)
            # - for barcode metadata and reductions use the mean
            if (contrast[["data_type"]] == "feature_data") {
                mean_fxn = NULL
            } else {
                mean_fxn = rowMeans
            }
            
            # Run FindMarkers
            deg_results = Seurat::FindMarkers(contrast[["sc_subset"]],
                                              ident.1="condition1",
                                              ident.2="condition2",
                                              test.use=contrast[["test"]],
                                              assay=SeuratObject::DefaultAssay(contrast[["sc_subset"]]),
                                              slot=contrast[["layer"]],
                                              fc.slot=dplyr::case_match(contrast[["data_type"]],
                                                                        "feature_data" ~ "data",
                                                                        "reduction" ~ "counts",
                                                                        "barcode_metadata" ~ "counts"),
                                              fc.name=dplyr::case_match(contrast[["data_type"]],
                                                                         "feature_data" ~ "avg_log2FC",
                                                                         "reduction" ~ "avg_diff",
                                                                         "barcode_metadata" ~ "avg_diff"),
                                              mean.fxn=mean_fxn,
                                              random.seed=getOption("random_seed"),
                                              min.cells.group=2,
                                              logfc.threshold=contrast[["log2FC"]],
                                              min.pct=contrast[["min_pct"]],
                                              densify="bulk_by" %in% names(contrast) | 
                                                "pseudobulk_samples" %in% names(contrast),
                                              latent.vars=contrast[["latent.vars"]])
        }
    } else {
        deg_results = data.frame()
        contrast[["message"]] = FormatString("There are fewer than two samples in at least one group for comparison {name}.")
    }
    
    # If there are no DEGs, create an empty data frame
    if (nrow(deg_results) == 0) {
      empty_deg_results = data.frame(pct.1=as.numeric(), pct.2=as.numeric(), p_val=as.numeric(), p_val_adj=as.numeric())
      if (contrast[["data_type"]] == "feature_data") {
        empty_deg_results$avg_log2FC = as.numeric() 
      } else if (contrast[["data_type"]] %in% c("reduction", "barcode_metadata")) {
        empty_deg_results$avg_diff = as.numeric()
      }
    }
    
    # Gene rownames to column
    deg_results$gene = rownames(deg_results)

    # Sort results
    deg_results = deg_results %>% 
        DegsSort()

    # Add means (per condition)
    # - For feature data (gene expression): use layer 'data'  (normalized and log-transformed) and calculate means of the exponentiated data followed by log.
    #   To (re-)calculate avg_log2FC from the natural log-transformed means, switch to base 2 and subtract: (mean1/log(2)) - (mean2/log(2)). See log laws.
    # - For reduction and barcode metadata: use layer 'counts (raw)' and calculate just the mean.
    avg_df = AverageCounts(contrast[["sc_subset"]], 
                           assay=SeuratObject::DefaultAssay(contrast[["sc_subset"]]), 
                           layer=dplyr::case_match(contrast[["data_type"]],
                                                   "feature_data" ~ "data",
                                                   "reduction" ~ "counts",
                                                   "barcode_metadata" ~ "counts"))
    avg_df = as.data.frame(avg_df) %>%
      tibble::rownames_to_column(var="gene")
    deg_results = dplyr::inner_join(deg_results, avg_df, by="gene")
    
    # If bulk_by or pseudobulk_samples is set, add the bulked expression values
    if ("bulk_by" %in% names(contrast) | "pseudobulk_samples" %in% names(contrast)) {
      bulk_df = SeuratObject::GetAssayData(contrast[["sc_subset"]], 
                                           assay=SeuratObject::DefaultAssay(contrast[["sc_subset"]]),
                                           layer=dplyr::case_match(contrast[["data_type"]],
                                                                   "feature_data" ~ "data",
                                                                   "reduction" ~ "counts",
                                                                   "barcode_metadata" ~ "counts")) %>%
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
    library(clusterProfiler)
  
    # Get DEGs
    degs = deg_result$results %>% 
      dplyr::filter(p_val_adj < deg_result$padj & abs(avg_log2FC) >= deg_result$log2FC) %>%
      dplyr::pull(gene) %>%
      unique()
    
    # Get universe
    universe = deg_result$results %>%
      dplyr::pull(gene) %>%
      unique()
    
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
        
      # Define an empty ORA result to return in case we cannot do the ORA analysis
      empty_ora = new("enrichResult",
                    result = data.frame(ID=as.character(NULL), Description=as.character(NULL), 
                                        GeneRatio=as.character(NULL), BgRatio=as.character(NULL), 
                                        RichFactor=as.numeric(NULL), FoldEnrichment=as.numeric(NULL), zScore=as.numeric(NULL), 
                                        pvalue=as.numeric(NULL), p.adjust=as.numeric(NULL), qvalue=as.numeric(NULL), 
                                        geneID=as.character(NULL), Count=as.numeric(NULL)),
                    readable = FALSE, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                    qvalueCutoff = 0.2, organism = "UNKNOWN", ontology = "UNKNOWN",
                    gene = degs, keytype = "UNKNOWN", universe = universe,
                    gene2Symbol = as.character(NULL), geneSets = split(t2g$gene_symbol, t2g$gs_name))
        
      # If there are no DEGs or universe, return NULL
      if (length(degs) == 0 | length(universe) == 0) return(empty_ora)
      
      # Run ORA
      ora = clusterProfiler::enricher(gene=degs, universe=universe, TERM2GENE=t2g, minGSSize=10, maxGSSize=500)
      if (is.null(ora)) return(empty_ora)
      
      # Older clusterProfiler versions do not have RichFactor, FoldEnrichment and zScore. Set it to NA.
      if (!"RichFactor" %in% colnames(ora@result)) {
        ora@result$RichFactor = as.numeric(NA)
      }
      if (!"FoldEnrichment" %in% colnames(ora@result)) {
        ora@result$FoldEnrichment = as.numeric(NA)
      }
      if (!"zScore" %in% colnames(ora@result)) {
        ora@result$zScore = as.numeric(NA)
      }
      
      # Sort by descreasing FoldEnrichment
      ora@result = ora@result[order(ora@result$FoldEnrichment, decreasing=TRUE),]
      
      # Fix factor levels accordingly for plots
      ora@result$ID = factor(ora@result$ID, levels=unique(ora@result$ID))
      ora@result$Description = factor(ora@result$Description, levels=unique(ora@result$Description))
      
      # Order columns
      ora@result = ora@result %>%
        dplyr::relocate(ID, Description, GeneRatio, BgRatio, RichFactor, FoldEnrichment, zScore, pvalue, p.adjust, qvalue, geneID, Count)
      
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
    
    # In recent msigdbr package versions gs_cat is now gs_collection, gs_subcat gs_subcollection and gene_id gene_symbol - add columns
    term2gene_db$gs_cat = term2gene_db$gs_collection
    term2gene_db$gs_subcat = term2gene_db$gs_subcollection
    term2gene_db$gene_id = term2gene_db$gene_symbol
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
      
      # Check that the table contains the required columns
      required_columns = c("gs_cat", "gs_subcat", "gs_name", "gene_id", "gene_symbol")
      assertthat::assert_that(all(required_columns %in% colnames(term2gene)),
                              msg=FormatString("Geneset file {file} is missing at least one of these columns: 'gs_cat', 'gs_subcat', 'gs_name', 'gene_id', 'gene_symbol'."))
      term2gene = dplyr::select(term2gene, gs_cat, gs_subcat, gs_name, gene_id, gene_symbol)
      
      return(term2gene)
    })
  }
  
  return(term2gene_db)
}

#' Run DESeq2 on bulk-aggregated single-cell dataset. Operates on raw counts and does not do pre-filtering of genes.
#' 
#' By default it compares two conditions (Wald test). Alternatively, it can do a log ratio test where it compares a full model to a reduced model to see which fits better. This is done for all levels of the condition_groups column together.
#'
#' @param object Seurat object. Single-cell counts must have been aggregated into bulk samples. Object must have a 'condition_groups' column. 
#' @param ident_1 Condition 1. Must be part of the 'condition_groups' column. Ignored for LRT tests.
#' @param ident_2 Condition 2. Must be part of the 'condition_groups' column. Ignored for LRT tests.
#' @param assay Assay to use for the analysis (default: default assay of Seurat object)
#' @param test Use Wald or LRT test (default: Wald)
#' @param design Design formula of the full model to test (default: ~condition_groups).
#' @param reduced When doing an LRT test, design of the reduced model to compare (default: null).
#' @param random_seed Random seed for reproducibility (default: 1)
#' @param results_args Additional arguments that are passed to the results function of DESeq2.
#' @return A data frame with DESeq2 results. Contains columns gene, p_val, p_val_adj, avg_log2FC, pct.1, pct.2.
DegsRunDESeq2 = function(object, ident_1, ident_2, assay=NULL, test="Wald", design="~condition_groups", reduced=NULL, random_seed=1, results_args=NULL) {
  
  # Need these packages preloaded
  library(DESeq2)
  library(IHW)
  
  # Get assay
  if (is.null(assay)) assay = SeuratObject::DefaultAssay(object)
  
  # Get counts and convert them to standard matrix
  counts_data = SeuratObject::GetAssayData(object, assay=assay, layer="counts")
  counts_data = as.matrix(counts_data)
  
  # Get condition table
  col_data = object[[]]
  assertthat::assert_that("condition_groups" %in% colnames(col_data), 
                          msg="Column 'condition_groups' not found in the barcode metadata of the Seurat object.")
  
  # Set up design formula
  if (is.null(design)) {
    design = "~condition_groups"
  }
  
  # Run DESeq2
  dds = DESeq2::DESeqDataSetFromMatrix(countData=counts_data, colData=col_data, design=as.formula(design))
  if (test == "Wald") {
    # Run Wald test
    dds = DESeq2::DESeq(dds, test="Wald")
  } else if(test == "LRT") {
    # Set up reduced formula 
    if (is.null(reduced)) reduced = "~1"
    
    # Run LRT test
    dds = DESeq2::DESeq(dds, test="LRT", reduced=as.formula(reduced))
  }
  
  # Get results
  if (is.null(results_args)) {
    results_args = list(
      contrast=c("condition_groups", ident_1, ident_2), 
      alpha=0.05, 
      filterFun=ihw
    )
  }
  dds_results = do.call(DESeq2::results, 
                        c(list(dds), results_args))
  dds_results = as.data.frame(dds_results)
  
  # Fix NAs for pvalue and padj
  dds_results$pvalue[is.na(dds_results$pvalue)] = 1
  dds_results$padj[is.na(dds_results$padj)] = 1
  
  # Add percentage of cells that express the gene in condition 1 and condition 2
  cts = counts_data[, col_data$condition_groups == ident_1]
  dds_results$pct.1 = rowSums(cts > 0) / ncol(cts)
    
  cts = counts_data[, col_data$condition_groups == ident_2]
  dds_results$pct.2 = rowSums(cts > 0) / ncol(cts)
  
  # Convert to data.frame and rename columns
  dds_results = dds_results %>% 
    as.data.frame() %>% 
    dplyr::select(avg_log2FC=log2FoldChange, p_val=pvalue, p_val_adj=padj, pct.1, pct.2)
  
  return(dds_results)
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
    group1 = paste(deg_result[["condition1"]], collapse="+")
    group2 = paste(deg_result[["condition2"]], collapse="+")
    
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
        dplyr::slice_min(order_by=abs(p_val ), n=5, with_ties=FALSE)
    
    # Make plot
    # Note: Font size in theme is measured in pt but font size in geom_text is measured in mm.
    # Ggplot2 provided '.pt' as conversion factor: mm = pt / .pt OR pt = mm * .pt
    # https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#text
    p = ggplot(deg_table, aes(x=condition1, y=condition2, col=deg_status)) + 
        geom_abline(slope=1, intercept=0, col="lightgrey") +
        geom_point() +
        ggrepel::geom_text_repel(data=top10_deg_table, aes(x=condition1, y=condition2, col=deg_status, label=gene), size=font_size / .pt, colour="black") +
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
    group1 = paste(deg_result[["condition1"]], collapse="+")
    group2 = paste(deg_result[["condition2"]], collapse="+")
    
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
        dplyr::slice_min(order_by=abs(p_val ), n=5, with_ties=FALSE)
    
    # Make plot
    # Note: Font size in theme is measured in pt but font size in geom_text is measured in mm.
    # Ggplot2 provided '.pt' as conversion factor: mm = pt / .pt OR pt = mm * .pt
    # https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#text
    p = ggplot(deg_table, aes(x=!!sym(fc_col), y=p_val_log10_n, col=deg_status)) + 
        geom_point() +
        ggrepel::geom_text_repel(data=top10_deg_table, aes(x=!!sym(fc_col), y=p_val_log10_n, col=deg_status, label=gene), size=font_size / .pt, colour="black") +
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
