
#' Glue Transformer for Quote and Collapse Operations
#'
#' A transformer function for the glue package that optionally quotes variables
#' and handles multiple values with automatic collapsing.
#'
#' @param sep Character. Separator to use when collapsing multiple values.
#'   Default is \code{", "}.
#' @param quote Logical. Whether to quote variables. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{glue_collapse}.
#'
#' @return A transformer function to be used with \code{glue::glue}.
#'
#' @details
#' Use \code{*} suffix in glue strings to indicate that a variable may contain
#' multiple values that should be collapsed. For example, \code{{variable*}}
#' will collapse all values of \code{variable} with the specified separator.
#'
#' @seealso \code{\link{FormatString}}
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#' # Used internally by FormatString
#' glue::glue("Values: {x*}", x = c("a", "b", "c"),
#'   .transformer = GlueTransformer_quote_collapse())
#' # Returns: "Values: 'a', 'b', 'c'"
#' }
GlueTransformer_quote_collapse = function(sep=", ", quote=TRUE, ...) {
  function(text, envir) {
    collapse = grepl("[*]$", text)
    if (collapse) {
      text = sub("[*]$", "", text)
    }
    res = glue::identity_transformer(text, envir)
    if (quote) {
      res = glue::single_quote(res)
    }
    if (collapse) {
      glue::glue_collapse(res, sep, ...)  
    } else {
      res
    }
  }
}

#' Format Strings with Variable Interpolation
#'
#' A wrapper around \code{glue::glue} that provides convenient string formatting
#' with optional quoting and automatic handling of multiple values.
#'
#' @param x Character. The string template with variables in curly braces.
#' @param quote Logical. Whether to quote interpolated variables.
#'   Default is \code{TRUE}.
#' @param sep Character. Separator for collapsing multiple values when using
#'   the \code{*} suffix. Default is \code{", "}.
#'
#' @return A glue object (character string) with interpolated values.
#'
#' @details
#' Variables can be inserted using \code{{variable}} syntax. Add \code{*} suffix
#' to collapse multiple values: \code{{variable*}}.
#'
#' @importFrom glue glue single_quote glue_collapse identity_transformer
#' @export
#'
#' @examples
#' \dontrun{
#' # Simple interpolation
#' name <- "World"
#' FormatString("Hello, {name}!")
#' # Returns: "Hello, 'World'!"
#'
#' # Multiple values with collapse
#' items <- c("apple", "banana", "cherry")
#' FormatString("Fruits: {items*}")
#' # Returns: "Fruits: 'apple', 'banana', 'cherry'"
#'
#' # Without quoting
#' FormatString("Value: {x}", x = 42, quote = FALSE)
#' # Returns: "Value: 42"
#' }
FormatString = function(x, quote=TRUE, sep=", ") {
  return(glue::glue(x, .transformer=GlueTransformer_quote_collapse(quote=quote, sep=sep), .envir=parent.frame()))
}

#' Generate Callout Box for Documentation
#'
#' Creates a Quarto-style callout box for displaying notes, tips, warnings,
#' or other important information in rendered documents.
#'
#' @param x Character. The message content. Will be formatted using
#'   \code{\link{FormatString}}.
#' @param type Character. Type of callout box. Options: \code{"note"},
#'   \code{"tip"}, \code{"important"}, \code{"caution"}, \code{"warning"}.
#' @param print Logical. If \code{TRUE}, prints the callout box using \code{cat}.
#'   If \code{FALSE}, returns the string. Default is \code{TRUE}.
#' @param quote Logical. Whether to quote expanded variables in the message.
#'   Default is \code{TRUE}.
#'
#' @return If \code{print = TRUE}, invisibly returns \code{NULL} and prints
#'   the callout box. If \code{print = FALSE}, returns the callout box as a
#'   character string.
#'
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' # Print a note
#' CalloutBox("Remember to save your work!", type = "note")
#'
#' # Get warning as string
#' warning_text <- CalloutBox("This may take a while.", type = "warning", print = FALSE
#' }
CalloutBox = function(x, type, print=TRUE, quote=TRUE) {
  valid_types = c("note", "tip", "important", "caution", "warning")
  assertthat::assert_that(type %in% valid_types,
                          msg=FormatString("Callout box typ {type} but must be one of: {valid_types*}."))
  
  
  
  x = paste0("\n\n::: callout-", type, "\n", FormatString(x, quote=quote), "\n:::\n\n")
  if (print) {
    cat(x)
  } else {
    return(x)
  }
}

#' Get Profile YAML Configuration
#'
#' Reads and returns the content of the current Quarto profile YAML file.
#' The profile is determined by the \code{QUARTO_PROFILE} environment variable.
#'
#' @return A nested list containing the parsed YAML configuration.
#'
#' @details
#' The function expects a file named \code{_quarto-<profile>.yml} in the
#' current working directory, where \code{<profile>} is the value of the
#' \code{QUARTO_PROFILE} environment variable.
#'
#' @importFrom yaml read_yaml
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' Sys.setenv(QUARTO_PROFILE = "development")
#' config <- GetProfileYaml()
#' }
GetProfileYaml = function() {
  profile = Sys.getenv("QUARTO_PROFILE")
  assertthat::assert_that(nchar(profile) > 0,
                          msg="Environment variable 'QUARTO_PROFILE' must be set to the current profile.")
  
  profile_yml = yaml::read_yaml(paste0("_quarto-", profile, ".yml"), eval.expr=TRUE)
  return(profile_yml)
}

#' Get Previous Module Directory
#'
#' Given the current module directory, returns the module that was run before
#' this module according to the workflow defined in the current profile YAML.
#'
#' @param current_module_dir Character. The current module directory.
#'
#' @return Character string with the previous module directory, or \code{NULL}
#'   if there is no previous module (i.e., this is the first module).
#'
#' @details
#' The workflow order is determined by the \code{chapters} parameter in the
#' \code{book} section of the profile YAML file.
#'
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' prev_module <- PreviousModuleDir("02_filtering")
#' }
PreviousModuleDir = function(current_module_dir) {
  current_module_dir = module_dir
  
  # Get profile yaml parameter 'chapters' in section 'book'
  profile_yml = GetProfileYaml()
  assertthat::assert_that("book" %in% names(profile_yml),
                          msg="Profile yaml does not contain the parameter 'book'.")
  assertthat::assert_that("chapters" %in% names(profile_yml$book),
                          msg="Profile yaml does not contain the parameter 'chapters' in the section 'book'.")
  chapters = dirname(profile_yml$book$chapters)
  
  # Find position
  idx = match(current_module_dir, chapters)
  idx = idx[1]
  assertthat::assert_that(!is.na(idx),
                          msg="Module is not part of parameter 'chapters' in the section 'book'.")
  
  if (idx == 1) {
    return(NULL)
  } else {
    return(chapters[(idx-1)])
  }
}

#' Access Workflow Parameters
#'
#' Retrieves workflow parameters by merging module YAML parameters with
#' general and module-specific parameters from the profile YAML.
#'
#' @param p Character or \code{NULL}. Parameter name to retrieve. If \code{NULL},
#'   returns all parameters as a list. Default is \code{NULL}.
#'
#' @return If \code{p} is specified, returns the value of that parameter.
#'   If \code{p} is \code{NULL}, returns a list of all parameters.
#'
#' @details
#' Parameter precedence (later overrides earlier):
#' \enumerate{
#'   \item Module YAML parameters (document \code{params})
#'   \item Profile YAML general parameters
#'   \item Profile YAML module-specific parameters
#' }
#'
#' This function works both interactively in RStudio and during Quarto rendering.
#'
#' @importFrom purrr list_assign
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' # Get all parameters
#' all_params <- param()
#'
#' # Get a specific parameter
#' project_id <- param("project_id")
#' }
param = function(p=NULL) {

  # Read module parameter (document params yaml) and get module name
  assertthat::assert_that("module" %in% names(params),
                          msg="Module does not contain the document yaml parameter 'module' with the module name.")
  module_name = params[["module"]]
  param_set = params
  
  # Read profile params from the profile yaml file
  profile_yml = GetProfileYaml()
  if ("params" %in% names(profile_yml)) {
    profile_params = profile_yml[["params"]]
    
    # Get general parameter if available
    if ("general" %in% names(profile_params)) {
      param_set = purrr::list_assign(param_set, !!!profile_params[["general"]])
    }
    
    # Get module-specific parameter if available
    if ("modules" %in% names(profile_params)) {
      profile_module_params = profile_params[["modules"]]
      
      if (module_name %in% names(profile_module_params)) {
        param_set = purrr::list_assign(param_set, !!!profile_module_params[[module_name]])
      }
    }
  }
  
  if(is.null(p)) {
    return(param_set)
  } else {
    return(param_set[[p]])
  }
}

#' Get Ensembl BioMart
#'
#' Creates a connection to an Ensembl BioMart database for a specific species
#' and Ensembl version.
#'
#' @param species Character. Latin species name in \code{genus_species} format
#'   (e.g., \code{"homo_sapiens"}, \code{"mus_musculus"}).
#' @param ensembl_version Character or numeric. Ensembl version number
#'   (e.g., \code{"98"} or \code{98}).
#'
#' @return A biomaRt \code{Mart} object connected to the specified database.
#'
#' @importFrom biomaRt listEnsemblArchives useMart listDatasets useDataset
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' # Connect to human Ensembl version 98
#' mart <- GetBiomaRt("homo_sapiens", 98)
#'
#' # Connect to mouse Ensembl
#' mart <- GetBiomaRt("mus_musculus", "104")
#' }
GetBiomaRt = function(species, ensembl_version) {
  # Check if we can find find an Ensembl database for this species and annotation version
  ensembl_archives = biomaRt::listEnsemblArchives()
  assertthat::assert_that(
    ensembl_version %in% ensembl_archives$version,
    msg = FormatString("Could not find or access Ensembl version {ensembl_version}.")
  )
  
  # Get mart and check if the species is part of ensembl
  ensembl_mart = biomaRt::useMart(biomart = "ensembl", host = ensembl_archives[match(ensembl_version, ensembl_archives$version), "url"])
  ensembl_datasets = biomaRt::listDatasets(ensembl_mart)
  
  if (species == "heterocephalus_glaber") {
    # Hack for heterocephalus_glaber which does not fit in the usual Ensembl way of naming datasets
    species_dataset_name = paste0("hgfemale", "_gene_ensembl")
  } else {
    species_dataset_name = paste0(gsub("^(.)[a-z]+_", "\\1", species), "_gene_ensembl")
  }
  
  idx = which(ensembl_datasets$dataset == species_dataset_name)
  assertthat::assert_that(
    length(idx) == 1,
    msg = FormatString(
      "Could not find species {species} dataset (name: {species_dataset_name}) for Ensembl version {ensembl_annotation_version}."
    )
  )
  
  # Get dataset
  ensembl_mart = biomaRt::useDataset(ensembl_datasets[idx[1], "dataset"], ensembl_mart)
  
  return(ensembl_mart)
}

#' Fetch Gene Information from Ensembl
#'
#' Retrieves gene annotations from Ensembl BioMart for a list of gene IDs
#' or symbols.
#'
#' @param ids Character vector. Ensembl gene IDs or gene symbols to query.
#' @param symbols Logical. If \code{TRUE}, \code{ids} are interpreted as gene
#'   symbols instead of Ensembl IDs. Default is \code{FALSE}.
#' @param species Character. Species name in \code{genus_species} format.
#' @param ensembl_version Character or numeric. Ensembl version to use.
#' @param mart_attributes Named character vector. BioMart attributes to fetch.
#'   Names become column names in the output. Default includes gene ID, symbol,
#'   biotype, description, chromosome, and strand information.
#' @param useCache Logical. Use local cache for faster querying.
#'   Default is \code{TRUE}.
#'
#' @return A data frame with requested gene information. IDs not found in
#'   Ensembl are included but with \code{NA} values.
#'
#' @importFrom biomaRt getBM
#' @importFrom dplyr left_join
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' # Fetch information for Ensembl IDs
#' gene_info <- EnsemblFetchGeneInfo(
#'   ids = c("ENSG00000141510", "ENSG00000157764"),
#'   species = "homo_sapiens",
#'   ensembl_version = 98
#' )
#'
#' # Fetch by gene symbols
#' gene_info <- EnsemblFetchGeneInfo(
#'   ids = c("TP53", "BRAF"),
#'   symbols = TRUE,
#'   species = "homo_sapiens",
#'   ensembl_version = 98
#' )
#' }
EnsemblFetchGeneInfo = function(ids, symbols=FALSE, species, ensembl_version, mart_attributes=c(ensembl_id="ensembl_gene_id", ensembl_symbol="external_gene_name", ensembl_biotype="gene_biotype", ensembl_description="description", ensembl_chr="chromosome_name", ensembl_start_position="start_position", ensembl_end_position="end_position", ensembl_strand="strand"), useCache=TRUE) {
  # Get species mart
  species_mart = GetBiomaRt(species, ensembl_version)
  
  # Check that we have Ensembl ids
  if (!symbols) {
    assertthat::assert_that(any(grepl("^ENS", ids)),
                          msg="None of the ids in this dataset is Ensembl. Cannot fetch gene information with this method.")
    id_column = "ensembl_gene_id"
  } else {
    id_column = "external_gene_name"
  }
  
  if (!id_column %in% mart_attributes) {
    mart_attributes = c(setNames(id_column, id_column), mart_attributes)
  }
  
  # Fetch attributes from Ensembl (use cache to allow multiple fetches)
  species_annotation = biomaRt::getBM(mart=species_mart, 
                                      filters=id_column, 
                                      values=ids, 
                                      attributes=mart_attributes, 
                                      useCache=useCache)
  if (!is.null(names(mart_attributes))) {
    colnames(species_annotation) = names(mart_attributes)
  }
  assertthat::assert_that(nrow(species_annotation)>0,
    msg = FormatString(
      "Could not find fetch any gene information for this dataset. Is the species {species} correct?"))
  
  # Add rows for ids that were not found
  idx = which(mart_attributes == id_column)
  first_y_col = names(mart_attributes)[idx]
  x_df = data.frame(ids)
  colnames(x_df) = first_y_col
  annotation_ensembl = dplyr::left_join(x_df, species_annotation, by=first_y_col)
  
  return(annotation_ensembl)
}

#' Fetch Orthologues from Ensembl
#'
#' Retrieves orthologue information between two species from Ensembl using
#' BioMart's cross-species linking.
#'
#' @param ids Character vector. Ensembl gene IDs or gene symbols from species1.
#' @param symbols Logical. If \code{TRUE}, \code{ids} are gene symbols.
#'   Default is \code{FALSE}.
#' @param species1 Character. First species in \code{genus_species} format.
#' @param species2 Character. Second species in \code{genus_species} format.
#' @param ensembl_version Character or numeric. Ensembl version to use.
#' @param mart_attributes1 Named character vector. Attributes to fetch for
#'   species1. Default includes gene ID and symbol.
#' @param mart_attributes2 Named character vector. Attributes to fetch for
#'   species2. Default includes gene ID and symbol.
#' @param useCache Logical. Use local cache. Default is \code{TRUE}.
#'
#' @return A data frame with orthologue mappings including a \code{one_to_one}
#'   column indicating whether the orthologue relationship is one-to-one.
#'   IDs without orthologues are included with \code{NA} values.
#'
#' @importFrom biomaRt getLDS
#' @importFrom dplyr left_join group_by summarise filter pull mutate
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' # Get human-mouse orthologues
#' orthologs <- EnsemblFetchOrthologues(
#'   ids = c("ENSG00000141510", "ENSG00000157764"),
#'   species1 = "homo_sapiens",
#'   species2 = "mus_musculus",
#'   ensembl_version = 98
#' )
#' }
EnsemblFetchOrthologues = function(ids, symbols=FALSE, species1, species2, ensembl_version, mart_attributes1=c(ensembl_id1="ensembl_gene_id", ensembl_symbol1="external_gene_name"), mart_attributes2=c(ensembl_id2="ensembl_gene_id", ensembl_symbol2="external_gene_name"), useCache=TRUE) {
  # Get species marts
  species1_mart = GetBiomaRt(species1, ensembl_version)
  species2_mart = GetBiomaRt(species2, ensembl_version)
  
  # Check that we have Ensembl ids
  if (!symbols) {
    assertthat::assert_that(any(grepl("^ENS", ids)),
                            msg = "None of the ids in this dataset is Ensembl. Cannot fetch gene information with this method.")
    id_column = "ensembl_gene_id"
  } else {
    id_column = "external_gene_name"
  }
  
  if (!id_column %in% mart_attributes1) {
    mart_attributes1 = c(setNames(id_column, id_column), mart_attributes1)
  }
  
  ortholog_annotation = biomaRt::getLDS(mart=species1_mart,
                            martL=species2_mart,
                            filters=id_column,
                            values=ids,
                            attributes=mart_attributes1,
                            attributesL=mart_attributes2, 
                            uniqueRows=TRUE)
  if (!is.null(names(mart_attributes1)) & !is.null(names(mart_attributes2))) {
    colnames(ortholog_annotation) = c(names(mart_attributes1), names(mart_attributes2))
  }
  assertthat::assert_that(nrow(ortholog_annotation)>0,
                          msg = FormatString(
                            "Could not find fetch any ortholog information for this dataset. Is the species {species1} correct?"))

  
  # Identify one-to-one orthologues
  idx = which(mart_attributes1 == id_column)
  first_y_col = names(mart_attributes1)[idx]
  one_to_one = ortholog_annotation %>% 
    dplyr::group_by(!!sym(first_y_col)) %>%
    dplyr::summarise(n=dplyr::n()) %>%
    dplyr::filter(n==1) %>%
    dplyr::pull(!!sym(first_y_col))
  
  ortholog_annotation = ortholog_annotation %>%
    dplyr::mutate(one_to_one = !!sym(first_y_col) %in% one_to_one)
  
  # Add rows for ids that were not found
  idx = which(mart_attributes1 == id_column)
  first_y_col = names(mart_attributes1)[idx]
  x_df = data.frame(ids)
  colnames(x_df) = first_y_col
  ortholog_annotation = dplyr::left_join(x_df, ortholog_annotation, by=first_y_col)
  
  return(ortholog_annotation)
  
}  

#' Make Names Syntactically Valid
#'
#' Converts names to syntactically valid R identifiers by replacing invalid
#' characters with underscores.
#'
#' @param x Character vector. Names to make valid.
#'
#' @return A character vector with syntactically valid names.
#'
#' @details
#' Uses \code{make.names} internally and then cleans up multiple consecutive
#' dots, leading/trailing dots, and replaces remaining dots with underscores.
#'
#' @export
#'
#' @examples
#' MakeNamesValid(c("my-var", "123abc", "a.b.c"))
#' # Returns: c("my_var", "X123abc", "a_b_c")
MakeNamesValid = function(x) {
   x = make.names(x) %>%  
        gsub("\\.\\.+", ".", .) %>% 
        gsub("\\.$", "", .) %>% 
        gsub("^\\.", "", .) %>% 
        gsub("\\.", "_", .)
   return(x)
}

#' Add Feature Metadata to Seurat Object
#'
#' Adds or updates feature (gene) metadata in a Seurat object or Assay.
#' This function works around limitations in \code{Seurat::AddMetaData} for
#' feature metadata in Seurat v5.
#'
#' @param obj A Seurat object or Assay5/Assay object.
#' @param assay Character or \code{NULL}. Assay to modify. Ignored if \code{obj}
#'   is an Assay. Default is \code{NULL}.
#' @param metadata Data frame. Feature metadata to add, with feature names as
#'   row names.
#'
#' @return The updated Seurat or Assay object.
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join select
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' # Add gene annotations to a Seurat object
#' gene_info <- data.frame(
#'   gene_type = c("protein_coding", "lncRNA"),
#'   row.names = c("Gene1", "Gene2")
#' )
#' seurat_obj <- AddFeatureMetadata(seurat_obj, assay = "RNA", metadata = gene_info)
#' }
AddFeatureMetadata = function(obj, assay=NULL, metadata) {
  # Checks
  valid_objs = c("Seurat", "Assay5", "Assay")
  assertthat::assert_that(class(obj) %in% valid_objs,
                          msg = "AddFeatureMetadata works only for 'Seurat v5' or 'Assay5' objects")
  
  # Prepare metadata tables for join
  metadata = metadata %>% tibble::rownames_to_column()
  
  if (is(obj, "Seurat")) {
    feature_metadata = obj[[assay]]@meta.data
  } else if (is(obj, "Assay5") | is(obj, "Assay")) {
    feature_metadata = obj@meta.data
  }
  feature_metadata = feature_metadata %>% tibble::rownames_to_column()
 
  # Join
  feature_metadata = dplyr::left_join(feature_metadata, metadata, by="rowname") %>% as.data.frame()
  rownames(feature_metadata) = feature_metadata$rowname
  feature_metadata = feature_metadata %>% dplyr::select(-rowname)
  
  # Add again
  if (is(obj, "Seurat")) {
    feature_metadata = feature_metadata[rownames(obj[[assay]]), , drop=FALSE]
    obj[[assay]]@meta.data = feature_metadata
  } else if (is(obj, "Assay5") | is(obj, "Assay")) {
    feature_metadata = feature_metadata[rownames(obj), , drop=FALSE]
    obj@meta.data = feature_metadata
  }
  
  return(obj)
}

#' Evaluate Knitr Chunk Code
#'
#' Extracts and evaluates R code from knitr-style code chunks.
#'
#' @param x Character. String containing knitr code chunks (text between
#'   \code{```{r}} and \code{```}).
#'
#' @return The result of evaluating the last expression in the extracted code.
#'
#' @importFrom stringr str_extract_all regex
#' @export
#'
#' @examples
#' \dontrun{
#' code <- "```{r}\nx <- 1 + 1\nx\n```"
#' result <- EvalKnitrChunk(code)  # Returns 2
#' }
EvalKnitrChunk = function(x) {
  chunks = unlist(stringr::str_extract_all(string=x, 
                           pattern=stringr::regex(pattern="```\\s*\\{r\\}.*?\\n```", dotall=TRUE)
                           ))
  chunks = gsub("```", "#", chunks)
  r_code = paste(chunks, collapse="\n")
  r_code = parse(text=r_code)
  return(eval(r_code))
}

#' Prepare Barcode Filter
#'
#' Processes filter configurations from YAML, applying general and sample-specific
#' filters to barcode metadata.
#'
#' @param filter List. Filter configuration from YAML. Can contain general
#'   filters and sample-specific filters (named by sample).
#' @param orig_idents Character vector. Sample identifiers in the analysis.
#' @param metadata Data frame. Barcode metadata to validate filter columns against.
#'
#' @return A named list with filter entries for each sample, or \code{NULL}
#'   if no filters are defined.
#'
#' @details
#' Filters can be specified as:
#' \itemize{
#'   \item Numeric vectors of length 2 for \code{c(min, max)} thresholds
#'   \item Character vectors for categorical values to keep
#' }
#'
#' @importFrom purrr map_depth list_modify
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' filter <- list(
#'   nCount_RNA = c(500, 50000),  # Applied to all samples
#'   Sample1 = list(pMito_RNA = c(NA, 20))  # Sample-specific
#' )
#' prepared <- PrepareBarcodeFilter(filter, c("Sample1", "Sample2"), metadata)
#' }
PrepareBarcodeFilter = function(filter, orig_idents, metadata) {
  if (is.null(filter) | length(filter) == 0) {
    return(NULL)
  }
  
  sample_specific_filter = filter[names(filter) %in% orig_idents]
  general_filter = filter[!names(filter) %in% orig_idents]
  filter = list()
  
  # Apply general filter to all samples
  if (length(general_filter) > 0) {
    filter = rep(list(general_filter), length(orig_idents))
    names(filter) = orig_idents
  }
  
  # Apply sample_specific filter (overwrite or add)
  filter = purrr::list_modify(filter, !!!sample_specific_filter)
  
  # Deal with filter that are all NA - assume that they are numeric
  filter = purrr::map_depth(filter, -1, function(f) {
    if (all(is.na(f))) {
      f = as.numeric(f)
    }
    return(f)
  })
  
  # Make sure filter columns exist in metadata
  filter_cols = purrr::map(filter, names) %>% 
    unlist() %>% 
    unique()
  assertthat::assert_that(all(filter_cols %in% colnames(metadata)),
                          msg=FormatString("Filter columns {filter_cols*} not found in barcode metadata."))
  
  return(filter)
}

#' Prepare Feature Filter
#'
#' Processes feature filter configurations from YAML, applying general and
#' sample-specific filters.
#'
#' @param filter List. Filter configuration from YAML.
#' @param orig_idents Character vector. Sample identifiers in the analysis.
#'
#' @return A named list with filter entries for each sample, or \code{NULL}
#'   if no filters are defined.
#'
#' @importFrom purrr map_depth list_modify
#' @export
#'
#' @examples
#' \dontrun{
#' filter <- list(
#'   min_counts = 10,
#'   Sample1 = list(min_cells = 5)
#' )
#' prepared <- PrepareFeatureFilter(filter, c("Sample1", "Sample2"))
#' }
PrepareFeatureFilter = function(filter, orig_idents) {
  if (is.null(filter) | length(filter) == 0) {
    return(NULL)
  }
  
  sample_specific_filter = filter[names(filter) %in% orig_idents]
  general_filter = filter[!names(filter) %in% orig_idents]
  filter = list()
  
  # Apply general filter to all samples
  if (length(general_filter) > 0) {
    filter = rep(list(general_filter), length(orig_idents))
    names(filter) = orig_idents
  }
  
  # Apply sample_specific filter (overwrite or add)
  filter = purrr::list_modify(filter, !!!sample_specific_filter)
  
  # Deal with filter that are all NA - assume that they are numeric
  filter = purrr::map_depth(filter, -1, function(f) {
    if (all(is.na(f))) {
      f = as.numeric(f)
    }
    return(f)
  })
  
  return(filter)
}

#' Add Lists to Seurat Misc Slot
#'
#' Stores named lists (e.g., gene lists, signatures) in the misc slot of a
#' Seurat object for later retrieval.
#'
#' @param sc A Seurat object.
#' @param lists Named list. Lists to store. Each element should be a character
#'   vector.
#' @param lists_slot Character. Name of the misc slot to use.
#'   Default is \code{"gene_lists"}.
#' @param add_to_list Logical. If \code{TRUE} and a list with the same name
#'   exists, append to it instead of overwriting. Default is \code{FALSE}.
#' @param make_unique Logical. If \code{TRUE}, remove duplicates from lists
#'   after storing. Default is \code{FALSE}.
#'
#' @return The Seurat object with updated lists in the misc slot.
#'
#' @importFrom Seurat Misc
#' @importFrom purrr map_chr
#' @export
#'
#' @examples
#' \dontrun{
#' # Add gene lists
#' gene_lists <- list(
#'   markers = c("Gene1", "Gene2"),
#'   housekeeping = c("GAPDH", "ACTB")
#' )
#' seurat_obj <- ScAddLists(seurat_obj, gene_lists)
#' }
ScAddLists = function(sc, lists, lists_slot='gene_lists', add_to_list=FALSE, make_unique=FALSE) {
  stored_lists = Seurat::Misc(sc, slot=lists_slot)
  if (is.null(stored_lists)) stored_lists = list()
  
  # Check that the lists have names, otherwise set to 'list_1', 'list_2' etc
  lists_names = purrr::map_chr(seq(lists), function(i) {
    n = names(lists[i])
    if (is.null(n) || nchar(n) == 0) return(paste0("list_",as.character(i))) else return(names(lists[i]))
  })
  names(lists) = lists_names
  
  for(n in names(lists)) {
    if (n %in% names(stored_lists) & add_to_list == TRUE) {
      stored_lists[[n]] = c(stored_lists[[n]], unlist(lists[[n]]))
    } else {
      stored_lists[[n]] = unlist(lists[[n]])
    }
    
    if (make_unique) stored_lists[[n]] = unique(stored_lists[[n]])
  }
  
  # Add to Seurat object
  suppressWarnings({Seurat::Misc(sc, slot=lists_slot) = stored_lists})
  return(sc)
}

#' Get Lists from Seurat Misc Slot
#'
#' Retrieves one or more named lists from the misc slot of a Seurat object.
#'
#' @param sc A Seurat object.
#' @param lists Character vector or \code{NULL}. Names of lists to retrieve.
#'   If \code{NULL}, returns all lists. Default is \code{NULL}.
#' @param lists_slot Character or \code{NULL}. Name of the misc slot containing
#'   the lists. If \code{NULL}, retrieves from the top level of misc.
#'
#' @return If a single list is requested, returns that list directly.
#'   If multiple lists are requested, returns a named list of lists.
#'
#' @importFrom Seurat Misc
#' @importFrom assertthat assert_that
#' @export
#'
#' @examples
#' \dontrun{
#' # Get a specific list
#' markers <- ScLists(seurat_obj, "markers", lists_slot = "gene_lists")
#'
#' # Get all lists
#' all_lists <- ScLists(seurat_obj, lists_slot = "gene_lists")
#' }
ScLists = function(sc, lists=NULL, lists_slot=NULL) {
  stored_lists = Seurat::Misc(sc, slot=lists_slot)
  assertthat::assert_that(!is.null(stored_lists), 
                          msg=FormatString("No lists found in misc slot of Seurat object (list slot: {lists_slot})."))
  
  if (is.null(lists)) lists = names(stored_lists)
  assertthat::assert_that(all(lists %in% names(stored_lists)), 
                          msg=FormatString("List(s) {lists} not found in misc slot of Seurat object (list slot: {lists_slot})."))
  
  if (length(lists) > 1) {
    return(stored_lists[lists])
  } else {
    return(stored_lists[[lists]])
  }
}

#' Add Colours to Seurat Object
#'
#' Stores colour palettes for categorical variables in the misc slot of a
#' Seurat object. This is a wrapper around \code{\link{ScAddLists}}.
#'
#' @param sc A Seurat object.
#' @param colours Named list. Colour palettes where names are category values
#'   and values are colour codes.
#' @param colours_slot Character. Name of the misc slot to use.
#'   Default is \code{"colour_lists"}.
#' @param add_to_list Logical. If \code{TRUE}, append to existing palettes.
#'   Default is \code{FALSE}.
#' @param make_unique Logical. If \code{TRUE}, remove duplicate colour entries.
#'   Default is \code{FALSE}.
#'
#' @return The Seurat object with updated colour palettes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add cluster colours
#' cluster_colours <- list(
#'   seurat_clusters = c("0" = "red", "1" = "blue", "2" = "green")
#' )
#' seurat_obj <- ScAddColours(seurat_obj, cluster_colours)
#' }
ScAddColours = function(sc, colours, colours_slot='colour_lists', add_to_list=FALSE, make_unique=FALSE) {
  return(ScAddLists(sc, colours, lists_slot=colours_slot, add_to_list=add_to_list, make_unique=make_unique))
}

#' Get Colours from Seurat Object
#'
#' Retrieves colour palettes for one or more categories from the misc slot
#' of a Seurat object.
#'
#' @param sc A Seurat object.
#' @param categories Character vector or \code{NULL}. Category names to retrieve
#'   colours for. If \code{NULL}, returns all colour palettes.
#' @param colours_slot Character. Name of the misc slot containing colours.
#'   Default is \code{"colour_lists"}.
#'
#' @return Named character vector(s) of colours.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get colours for clusters
#' colours <- ScColours(seurat_obj, "seurat_clusters")
#' }
ScColours = function(sc, categories=NULL, colours_slot="colour_lists") {
  return(ScLists(sc, lists=categories, lists_slot=colours_slot))
}

#' Get Default Dimensionality Reduction
#'
#' Retrieves the name of the default dimensionality reduction for an assay.
#'
#' @param sc A Seurat object.
#' @param assay Character or \code{NULL}. Assay to query. If \code{NULL},
#'   uses the default assay.
#'
#' @return Character string with the name of the default dimensionality reduction.
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom SeuratObject Misc
#' @export
#'
#' @examples
#' \dontrun{
#' default_red <- DefaultReduct(seurat_obj, assay = "RNA")
#' }
DefaultReduct = function(sc, assay=NULL) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  return(SeuratObject::Misc(sc[[assay]], slot="default.dimred"))
}

#' Set Default Dimensionality Reduction
#'
#' Sets the name of the default dimensionality reduction for an assay.
#'
#' @param sc A Seurat object.
#' @param assay Character or \code{NULL}. Assay to modify. If \code{NULL},
#'   uses the default assay.
#' @param value Character. Name of the dimensionality reduction to set as default.
#'
#' @return The Seurat object with updated default dimensionality reduction.
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom SeuratObject Misc
#' @export
#'
#' @examples
#' \dontrun{
#' DefaultReduct(seurat_obj, assay = "RNA") <- "pca"
#' }
"DefaultReduct<-" <- function(sc, assay=NULL, value) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  suppressWarnings({SeuratObject::Misc(sc[[assay]], slot="default.dimred") = value})
  return(sc)
}

#' Get Default Visualization Method
#'
#' Retrieves the name of the default visualization method (e.g., "umap", "tsne")
#' for an assay.
#'
#' @param sc A Seurat object.
#' @param assay Character or \code{NULL}. Assay to query. If \code{NULL},
#'   uses the default assay.
#'
#' @return Character string with the name of the default visualization method.
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom SeuratObject Misc
#' @export
#'
#' @examples
#' \dontrun{
#' default_viz <- DefaultVisualization(seurat_obj)
#' }
DefaultVisualization = function(sc, assay=NULL) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  return(SeuratObject::Misc(sc[[assay]], slot="default.visualization"))
}

#' Set Default Visualization Method
#'
#' Sets the name of the default visualization method for an assay.
#'
#' @param sc A Seurat object.
#' @param assay Character or \code{NULL}. Assay to modify. If \code{NULL},
#'   uses the default assay.
#' @param value Character. Name of the visualization method to set as default.
#'
#' @return The Seurat object with updated default visualization method.
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom SeuratObject Misc
#' @export
#'
#' @examples
#' \dontrun{
#' DefaultVisualization(seurat_obj) <- "umap"
#' }
"DefaultVisualization<-" <- function(sc, assay=NULL, value) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  suppressWarnings({SeuratObject::Misc(sc[[assay]], slot="default.visualization") = value})
  return(sc)
}

#' Format Chunk Option
#'
#' Formats a single knitr/quarto chunk option as a YAML comment.
#'
#' @param name Character. Option name.
#' @param value Any. Option value (logical, character, numeric, or list).
#'
#' @return Character string with the formatted chunk option.
#'
#' @keywords internal
#' @export
FormatChunkOption = function(name, value) {
    # Convert value to correct format
    if (is.logical(value)) {
        value = ifelse(value, "true", "false")
    } else if (is.character(value)) {
        value = paste0("'", value, "'")
    } else if (is.null(value)) {
        value = ifelse(is.null(value), "null", value)
    }
    
    # Expand lists
    if (length(value) > 1) {
        value = paste0(c("", value), collapse = "\n#|\t- ")
    }
    
    return(paste0("#| ", name, ": ", value))
}

#' Generate Code Chunks
#'
#' Recursively generates knitr/quarto code chunks from a nested list of
#' objects (plots or tables).
#'
#' @param olist A list of objects to generate chunks for.
#' @param olist_name Character or \code{NULL}. Name of the list variable.
#'   Usually auto-detected.
#' @param indices Integer vector. Internal parameter for tracking recursion.
#' @param chunk_label Character or \code{NULL}. Label for the chunk.
#' @param chunk_caption Character or \code{NULL}. Caption for the chunk.
#' @param chunk_opts List or \code{NULL}. Additional chunk options.
#'
#' @return Character vector of chunk code.
#'
#' @keywords internal
#' @export
GenerateChunks = function(olist, olist_name=NULL, indices=c(), chunk_label=NULL, chunk_caption=NULL, chunk_opts=NULL) {
    # Note: This function will recurse through a nested list. All arguments except for 'indices' 
    # only reflect the current recursion. The 'indices' argument is used to track the recursion process so
    # that in the chunk code the correct entry of the nested list is used.
    
    # Use this to debug
    # browser()
    
    # Check if input is a single object; not trivial since some objects are lists
    class_lst = class(olist)
    is_single_object = !(length(class_lst) == 1 && class_lst[1] == "list")
    
    # This is the name of the nested list that is provided as input
    # Do not set it since it will be set automatically
    if (is.null(olist_name)) {
        olist_name = deparse(substitute(olist))
    }
    
    if (is_single_object) {
        # If input is not a list but a single object, end of recursion - generate chunk code
        chunk = c("```{r}")
        
        # Decide if object is a table or a plot
        type = ifelse(colnames(olist) %>% length() > 0, "table", "plot")
        
        # Add chunk label
        if (!is.null(chunk_label)) {
            if (type == "plot") {
                chunk_label = paste0("fig-", chunk_label)
            } else if (type == "table") {
                chunk_label = paste0("tbl-", chunk_label)
            }
            chunk = c(chunk, GenerateChunkOption("label", chunk_label))
        }
        
        # Add chunk caption
        if (!is.null(chunk_caption)) {
            if (type == "plot") {
                chunk = c(chunk, GenerateChunkOption("fig-cap", chunk_caption))
            } else if (type == "table") {
                chunk = c(chunk, GenerateChunkOption("tbl-cap", chunk_caption))
            }
        }
        
        # Parse and add additional chunk options
        chunk = c(chunk,
                  purrr::map_chr(names(chunk_opts), function(n) {
                      v = chunk_opts[[n]]
                      return(GenerateChunkOption(n, v))
                  })
        )
        
        # Print plot or table
        chunk = c(chunk, "")
        chunk = c(chunk, paste0("print(", olist_name, "[[c(", paste(indices, collapse=", "), ")]])"))
        
        # Finish chunk and collapse
        chunk = c(chunk, "```")
        chunks = chunk
    } else {
        # If input is list, this means that there is still a list of chunks (or a list of lists of chunks) to process
        # However, these chunks or list of chunks can have a header or belong to a panel tabset.
        chunks = c()
        
        # Get attributes for this list
        list_attributes = list()
        if ("attributes" %in% names(olist)) {
            list_attributes = olist[["attributes"]]
            olist = olist[setdiff(names(olist), "attributes")]   
        }
        
        # Get names for this list
        list_names = names(olist)
        
        # Do the elements of this list are part of a panel tabset? (incompatible with headers)
        panel_tabset = "panel-tabset" %in% names(list_attributes) && list_attributes[["panel-tabset"]]
        if (panel_tabset) chunks = c(chunks, "::: panel-tabset")
        
        # Do the elements of this list have headers? (incompatible with panel tabset)
        headers = "headers" %in% names(list_attributes) && list_attributes[["headers"]]
        
        # Now run next recursion - run GenerateChunks on each element
        chunks = c(chunks, purrr::map(seq_along(olist), function(i) {
            chunk = c()
            
            # If element of a panel tabset, add a panel tabset header
            #if (panel_tabset) chunk = c(chunk, "######")
            
            # Keep track of indices
            new_indices = c(indices, i)
            
            # Start another recursion
            chunk = c(chunk, GenerateChunks(olist[[i]],
                                            olist_name=olist_name,
                                            indices=new_indices,
                                            chunk_label=chunk_label[[new_indices]],
                                            chunk_caption=chunk_label[[new_indices]],
                                            chunk_opts=chunk_opts[[new_indices]]))
            return(chunk)
        }))
        if (panel_tabset) chunks = c(chunks, ":::")
    }
    
    return(chunks)
}

######################
######################

#' Test if Values Convert to Numbers
#'
#' Checks whether values can be converted to numeric type.
#'
#' @param x A vector of values.
#'
#' @return Logical vector indicating which values can be converted.
#'
#' @export
#'
#' @examples
#' converts_to_number(c("1", "2.5", "abc"))
#' # Returns: c(TRUE, TRUE, FALSE)
converts_to_number = function(x) {
  return(suppressWarnings(!is.na(as.numeric(na.omit(x)))))
}

#' Test if Values Convert to Logical
#'
#' Checks whether values can be converted to logical type.
#'
#' @param x A vector of values.
#'
#' @return Logical vector indicating which values can be converted.
#'
#' @export
#'
#' @examples
#' converts_to_logical(c("TRUE", "FALSE", "abc"))
#' # Returns: c(TRUE, TRUE, FALSE)
converts_to_logical = function(x) {
  return(suppressWarnings(!is.na(as.logical(na.omit(x)))))
}


#' Format First N Elements as String
#'
#' Given a vector, reports at most n elements as a concatenated string,
#' adding "..." if there are more elements.
#'
#' @param x A vector.
#' @param n Integer. Maximum number of elements to include. Default is \code{5}.
#' @param sep Character. Separator for concatenation. Default is \code{","}.
#'
#' @return A character string with at most n elements concatenated.
#'
#' @export
#'
#' @examples
#' first_n_elements_to_string(letters, n = 3)
#' # Returns: "a,b,c,..."
first_n_elements_to_string = function(x, n=5, sep=",") {
  s = paste(x[1:min(n,length(x))], collapse=sep)
  if (length(x) > n) s = paste(s, "...", sep=sep)
  return(s)
}

#' Get Git Repository Version
#'
#' Retrieves the current commit SHA of a git repository.
#'
#' @param path_to_git Character. Path to the git repository.
#'
#' @return Character string with the commit SHA, or \code{"NA"} if not available.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' version <- GitRepositoryVersion(".")
#' }
GitRepositoryVersion = function(path_to_git) {
  repo = tryCatch({system(paste0("git --git-dir=", path_to_git, "/.git rev-parse HEAD"), intern=TRUE)},
                  warning = function(war) {return("NA")})
  return(repo)
}

#' Get Container Version Information
#'
#' Retrieves container information from environment variables if the code
#' is running inside a container.
#'
#' @return Character string with container git name, commit ID, and build date,
#'   or \code{"NA"} if not running in a container.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' container_info <- ContainerVersion()
#' }
ContainerVersion = function() {
  container_info = c(Sys.getenv("CONTAINER_GIT_NAME"), Sys.getenv("CONTAINER_VERSION"), Sys.getenv("CONTAINER_GIT_COMMIT_ID"), Sys.getenv("CONTAINER_BUILD_DATE"))
  container_info = container_info[nchar(container_info)>0]
  if (length(container_info) > 0) {
    container_info = paste(container_info, collapse=", ")
  } else {
    container_info = c("NA")
  }
  return(container_info)
}

#' Get scrnaseq Session Information
#'
#' Collects comprehensive session information for reproducibility, including
#' git version, container info, R version, and loaded packages.
#'
#' @param path_to_git Character. Path to the scrnaseq2 git repository.
#'   Default is \code{"."}.
#'
#' @return A matrix with two columns (Name, Value) containing session information.
#'
#' @importFrom sessioninfo package_info
#' @export
#'
#' @examples
#' \dontrun{
#' session_info <- ScrnaseqSessionInfo()
#' knitr::kable(session_info)
#' }
ScrnaseqSessionInfo = function(path_to_git=".") {
  out = matrix(NA, nrow=0, ncol=2)
  colnames(out) = c("Name", "Value")
  
  # Run on
  out = rbind(out, c("Run on:", format(Sys.time(), "%a %b %d %X %Y")))
  
  # Git
  repo = GitRepositoryVersion(path_to_git)
  out = rbind(out, c("dcgc-bfx/scrnaseq2", repo))
  
  # Container (if available)
  out = rbind(out, c("Container", ContainerVersion()))
  
  # System
  info_session = sessionInfo()
  out = rbind(out, c("R", info_session$R.version$version.string))
  out = rbind(out, c("Platform", info_session$platform))
  out = rbind(out, c("Operating system", info_session$running))
  out = rbind(out, c("Host name", unname(Sys.info()["nodename"])))
  out = rbind(out, c("Host OS", paste0(unname(Sys.info()["version"]), " (", unname(Sys.info()["release"]), ") ")))
  
  info_pkgs = sessioninfo::package_info()
  out = rbind(out, c("Packages", paste(paste(info_pkgs$package, info_pkgs$loadedversion, sep=""), collapse=", ")))
  
  return(out)
}

#' Subsample Cells from Seurat Object
#'
#' Randomly samples a specified number of cells from a Seurat object,
#' optionally stratified by a grouping variable.
#'
#' @param sc A Seurat v5 object.
#' @param n Integer. Total number of barcodes to sample. Default is \code{500}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{1}.
#' @param group Character or \code{NULL}. If provided, samples equally from each
#'   group defined by this barcode metadata column. Default is \code{NULL}.
#'
#' @return Character vector of sampled barcode names.
#'
#' @importFrom purrr map flatten
#' @export
#'
#' @examples
#' \dontrun{
#' # Sample 500 cells total
#' sampled_cells <- SubsampleSC(seurat_obj, n = 500)
#'
#' # Sample equally from each cluster
#' sampled_cells <- SubsampleSC(seurat_obj, n = 500, group = "seurat_clusters")
#' }
SubsampleSC = function(sc, n=500, seed=1, group=NULL) {
  barcode_metadata = sc[[]]
  barcodes = rownames(barcode_metadata)
  if (!is.null(group)) {
    barcodes = split(barcodes, barcode_metadata[, group, drop=TRUE])
  } else {
    barcodes = list(barcodes)
  }
  
  n = round(n/length(barcodes))
  
  barcodes = purrr::map(barcodes, function(x) {
    set.seed(seed)
    return(sample(x, min(n, length(x))))
  })
  
  return(unlist(purrr::flatten(barcodes)))
}



#' Get Named Vector of Names
#'
#' Returns a named vector where both names and values are the original names.
#'
#' @param x A list or vector with names.
#'
#' @return A named vector with names as both names and values.
#'
#' @export
#'
#' @examples
#' list_names(list(a = 1, b = 2))
#' # Returns: c(a = "a", b = "b")
list_names = function(x) {
  return(setNames(names(x), names(x)))
}

#' Set Values as Names
#'
#' Returns a vector with its values used as names.
#'
#' @param x A vector.
#'
#' @return A vector with its values as names.
#'
#' @export
#'
#' @examples
#' values_to_names(c("a", "b", "c"))
#' # Returns: c(a = "a", b = "b", c = "c")
values_to_names = function(x) {
  return(setNames(x,x))
}

#' Get Named Vector of Indices
#'
#' Returns a named vector of indices, where names are from the input and
#' values are sequential integers.
#'
#' @param x A list or vector with names.
#'
#' @return A named vector with names as names and indices as values.
#'
#' @export
#'
#' @examples
#' list_indices(list(a = 1, b = 2, c = 3))
#' # Returns: c(a = 1, b = 2, c = 3)
list_indices = function(x) {
  return(setNames(seq(x), names(x)))
}

#' Generate Colours from a Palette
#'
#' Generates a specified number of colours from a palette function.
#' If more colours are needed than the palette provides, uses different
#' alpha values to extend the palette.
#'
#' @param num_colours Integer. Number of colours to generate.
#' @param names Character vector or \code{NULL}. Names to assign to colours.
#' @param palette Character. Name of the palette function (e.g.,
#'   \code{"ggsci::pal_igv"}). Default is \code{"ggsci::pal_igv"}.
#' @param palette_options List. Additional arguments (besides alpha) to pass
#'   to the palette function.
#' @param alphas Numeric vector. Alpha values to use. If more colours are
#'   needed than the palette provides, subsequent alpha values are used.
#'   Default is \code{c(1, 0.7, 0.3)}.
#'
#' @return Character vector of colour codes.
#'
#' @importFrom purrr flatten_chr map
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate 10 colours
#' colours <- GenerateColours(10)
#'
#' # Generate named colours
#' colours <- GenerateColours(3, names = c("A", "B", "C"))
#' }
GenerateColours = function(num_colours, names=NULL, palette="ggsci::pal_igv", alphas=c(1,0.7,0.3), palette_options=list()) {
  palette = tryCatch({eval(parse(text=palette))}, error=function(cond) return(NULL))
  if (is.null(palette)) stop("GenerateColours: Could not find specified palette!")
  
  colours = purrr::flatten_chr(purrr::map(alphas, function(a) {
    palette_options[["alpha"]] = a
    cols = suppressWarnings(do.call(do.call(palette, palette_options), list(100)))
    cols[!is.na(cols)]
  }))
  
  if (num_colours > length(colours)) {
    stop("GenerateColours: Cannot generate the requested number of colours. Please change palette or add alpha values.")
  }
  
  colours = colours[1:num_colours]
  if (!is.null(names)) colours = setNames(colours, names)
  return(colours)
}

#' Format Parameters as Table
#'
#' Converts a parameter list to a two-column table for display.
#'
#' @param params Named list. Parameters to format.
#'
#' @return A matrix with columns "Name" and "Value".
#'
#' @export
#'
#' @examples
#' \dontrun{
#' params <- list(a = 1, b = c(1, 2, 3), c = list(x = 1, y = 2))
#' knitr::kable(ScrnaseqParamsInfo(params))
#' }
ScrnaseqParamsInfo = function(params) { 
  
  # Initialize output table
  out = matrix(NA, nrow=0, ncol=2)
  colnames(out) = c("Name", "Value")
  
  # Convert a (named) vector to a string
  VectorToString = function(x) { 
    if (is.null(names(x))) return(toString(x))
    if (sum(names(x) != "") != length(x)) return(toString(x))
    return(paste(names(x), x, sep="=", collapse=", "))   
  }
  
  # Convert list to table
  for (i in seq(params)) {
    x = params[[i]]
    # List, function or simple vector?
    if(is.list(x)) {
      y = paste(names(x), sapply(names(x), function(j) VectorToString(x[[j]])), sep=":", collapse="; ")
    } else if (is.function(x)) {
      y = "function"
    } else {
      y = VectorToString(x)
    }
    out = rbind(out, c(names(params)[i], y))
  }
  
  return(out)
}

#' Check scrnaseq Workflow Parameters
#'
#' Validates and converts parameters for the scrnaseq workflow.
#'
#' @param param Named list. Parameters to validate.
#'
#' @return The validated parameter list with an \code{error_messages} element
#'   if any validation errors occurred.
#'
#' @export
check_parameters_scrnaseq = function(param) {
  error_messages = c()
  param[["error_messages"]] = NULL
  
  # Check project_id ###
  if (!"project_id" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'project_id' is missing!")
  } else {
    param$project_id = as.character(param$project_id)
  }
  
  # Check path_data ###
  if (!"path_data" %in% names(param) | 
      !is.data.frame(param$path_data) | 
      !ncol(param$path_data) >= 4 | 
      !nrow(param$path_data) > 0 | 
      any(!c("name", "type", "path", "stats") %in% colnames(param$path_data))) {
    
    error_messages = c(error_messages, "The parameter 'path_data' needs to be a non-empty data.frame with at least four columns - 'name' (dataset name), 'type' (10x or smartseq2), 'path' (path to counts directory or file), 'stats' (path to 10x metrics summary file; can be NA) - and an optional fifth column 'suffix' (suffix to cell names)")
  } else {
    # name
    param$path_data$name = as.character(param$path_data$name)
    if (any(duplicated(param$path_data$name))) {
      error_messages = c(error_messages, "The column 'name' of the parameter 'path_data' must contain unique values!")
    }
    
    # type
    param$path_data$type = as.character(param$path_data$type)
    if (any(!param$path_data$type %in% c("10x","smartseq2"))) {
      error_messages = c(error_messages, "The column 'type' of the parameter 'path_data' should be either '10x' or 'smartseq2'!")
    }
    
    # path
    param$path_data$path = file.path(param$path_data$path)
    datasets_10x = param$path_data %>% dplyr::filter(type=="10x")
    is_valid = purrr::map_lgl(datasets_10x$path, function(p){
      return(dir.exists(p) & file.exists(file.path(p,"barcodes.tsv.gz")) & file.exists(file.path(p,"features.tsv.gz")) & file.exists(file.path(p,"matrix.mtx.gz")))
    })
    if (length(is_valid) > 0 & any(!is_valid)) {
      error_messages = c(error_messages, "At least one 10x dataset 'path' of the parameter 'path_data' is not a directory or misses at least one of the following files: 'barcodes.tsv', 'features.tsv.gz', 'matrix.mtx.gz'!")
    }
  
    datasets_smartseq2 = param$path_data %>% dplyr::filter(type=="smartseq2")
    is_valid = purrr::map_lgl(datasets_smartseq2$path, file.exists)
    if (length(is_valid) > 0 & any(!is_valid)) {
      error_messages = c(error_messages, "For at least one smartseq2 dataset, the 'path' of the parameter 'path_data' could not be found!")
    }
  
    # stats
    param$path_data$stats = ifelse(is.na(param$path_data$stats), NA, file.path(param$path_data$stats))
    is_valid = purrr::map_lgl(param$path_data$stats, function(p){
      if (is.na(p)) return(TRUE) else return(file.exists(p))
    })
    if (length(is_valid) > 0 && any(!is_valid)) {
      error_messages = c(error_messages, "At least one 'stats' file of the parameter 'path_data' could not be found. If not available, please set to NA!")
    }
    
    # suffix
    if (!"suffix" %in% colnames(param$path_data)) {
      param$path_data$suffix = paste0("-", 1:nrow(param$path_data))
    }
  }
  
  # Check downsample_cells_n ###
  if (("downsample_cells_n" %in% names(param)) & (!is.null(param$downsample_cells_n))) {
    if (!converts_to_number(param$downsample_cells_n)) {
      error_messages = c(error_messages, "The parameter 'path_out' (path for output) is missing!")
    } else {
      param$downsample_cells_n = as.numeric(param$downsample_cells_n)
    }
  } else {
    param["downsample_cells_n"] = list(NULL)
  }
  
  # Check path_out ###
  if (!"path_out" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'path_out' (path for output) is missing!")
  } else {
    param$path_out = file.path(param$path_out)
  }
  
  # Check file_known_markers ###
  if (!is.null(param$file_known_markers)){
    param$file_known_markers = file.path(param$file_known_markers)
    if (!file.exists(param$file_known_markers)) error_messages = c(error_messages, "The parameter 'file_known_markers' (Excel file with markers) is set but the file cannot be found. If not available, please set to NULL!")
  }
  
  # Check mart_dataset ###
  if (!"mart_dataset" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'mart_dataset' (Biomart dataset name) is missing!")
  } else {
    param$mart_dataset = as.character(param$mart_dataset)
  }
  
  # Check annot_version ###
  if (!"annot_version" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'annot_version' (Ensembl version) is missing!")
  } else {
    param$annot_version = as.character(param$annot_version)
  }
  
  # Check annot_main ###
  if (!"annot_main" %in% names(param) |
      any(!c("ensembl", "symbol", "entrez") %in% names(param$annot_main))) {
    error_messages = c(error_messages, "The parameter 'annot_main' is missing or is not a named vector with names 'ensembl', 'symbol' and 'entrez' as well as corresponding values!")
  } else {
    param$annot_main = setNames(as.character(param$annot_main), names(param$annot_main))
  }
  
  # Check file_annot ###
  if (!is.null(param$file_annot)) {
    param$file_annot = file.path(param$file_annot)
  }
  
  # Check mart_attributes ###
  if (!"mart_attributes" %in% names(param) || !is.vector(param$mart_attributes) || length(param$mart_attributes)==0) {
    error_messages = c(error_messages, "The parameter 'mart_attributes' (Biomart attributes) is missing or is not a non-empty vector!")
  } else {
    param$mart_attributes = as.character(param$mart_attributes)
  }
  
  # Check biomart_mirror ###
  if (!is.null(param$biomart_mirror)) {
    if (!param$biomart_mirror %in% c("www", "useast", "uswest", "asia")) error_messages = c(error_messages, "The parameter 'biomart_mirror' (Biomart mirror) is set but does not contain one of the following values: 'www', 'uswest', 'useast' ,'asia'!")
    param$biomart_mirror = as.character(param$biomart_mirror)
  } else {
    param$biomart_mirror = "www"
  }
  
  # Check mt ###
  if (!"mt" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'mt' (prefix of mitochondrial genes) is missing!")
  } else {
    param$mt = as.character(param$mt)
  }
  
  # Check cell_filter: can be numeric with length 2 or character or factor with any length; can also contain sublists per sample with the same criteria ###
  if ("cell_filter" %in% names(param) && length(param$cell_filter) > 0) {
    is_valid = TRUE
    for (i in seq(param$cell_filter)) {
      f = param$cell_filter[[i]]
      
      if (is.list(f)) {
        # Sample-specific values
        for (j in seq(f)) {
          f_s = f[[j]]
          if (length(f_s) == 2 & all(converts_to_number(f_s))) {
            f_s = as.numeric(f_s)
          } else if ( (is.character(f_s) | is.factor(f_s)) & length(f_s) > 0) {
            f_s = as.character(f_s)
          } else {
            is_valid = FALSE
          }
          f[[j]] = f_s
        }
      } else {
        # Global values
        if (length(f) == 2 & all(converts_to_number(f))) {
          f = as.numeric(f)
        } else if ( (is.character(f) | is.factor(f)) & length(f) > 0) {
          f = as.character(f)
        } else {
          is_valid = FALSE
        }
      }
      param$cell_filter[[i]] = f
    }
    
    if (!is_valid) {
      error_messages = c(error_messages, "The parameter 'cell_filter' should contain filters for cell properties with the following structures: a) lists of length 2 with minimum and maxium for numeric properties (set NA if not applicable/no min/no max), b) character/factor vectors for categorial properties. This can also specified by sample using sublists with the sample names as used by the script!")
    }
  } else {
    param$cell_filter = list()
  }
  
  # Check feature_filter: so far only contains: min_counts and min_cells; can also contain sublists per sample with the same criteria  ###
  if ("feature_filter" %in% names(param) && length(param$feature_filter) > 0) {
    valid = TRUE
    
    # Always set a global default for the minimum number of counts; but sample-specific values will overrule it
    if (!"min_counts" %in% names(param$feature_filter)) {
      param$feature_filter[["min_counts"]] = 1
    } else {
      if (converts_to_number(param$feature_filter[["min_counts"]])) {
        param$feature_filter[["min_counts"]] = as.numeric(param$feature_filter[["min_counts"]])
      } else {
        valid = FALSE
      }
    }
    
    # Always set a global default for the minimum number of cells; but sample-specific values will overrule it
    if (!"min_cells" %in% names(param$feature_filter)) {
      param$feature_filter[["min_cells"]] = 1
    } else {
      if (converts_to_number(param$feature_filter[["min_cells"]])) {
        param$feature_filter[["min_cells"]] = as.numeric(param$feature_filter[["min_cells"]])
      } else {
        valid = FALSE
      }
    }
    
    # Check sample-specific values
    valid = TRUE
    for (n in setdiff(names(param$feature_filter), c("min_counts", "min_cells"))) {
      f = param$feature_filter[[n]]
      
      if ("min_counts" %in% names(f)) {
        if (converts_to_number(f[["min_counts"]])) {
          f[["min_counts"]] = as.numeric(f[["min_counts"]])
        } else {
          valid = FALSE
        }
      }
      
      if ("min_cells" %in% names(f)) {
        if (converts_to_number(f[["min_cells"]])) {
          f[["min_cells"]] = as.numeric(f[["min_cells"]])
        } else {
          valid = FALSE
        }
      }
      
      param$feature_filter[[n]] = f
    }
    
    if (!valid) error_messages =  c(error_messages, paste("The parameter 'feature_filter' can contain: a) 'min_counts' for the minimum counts for a gene to be considered expressed,",
                   "b) 'min_cells' for the minimum number of cells in which a gene must be expressed. This can also specified by sample using sublists with the sample names as used by the script!"))
    
  } else {
    param$feature_filter = list(min_counts = 1, min_cells = 1)
  }
  
  # Check samples_to_drop  ###
  if ("samples_to_drop" %in% names(param)) {
    param$samples_to_drop = as.character(param$samples_to_drop)
  } else {
    param$samples_to_drop = as.character(NULL)
  }
  
  # Check samples_min_cells  ###
  if ("samples_min_cells" %in% names(param)) {
    if (converts_to_number(param$samples_min_cells)) {
      param$samples_min_cells = as.numeric(param$samples_min_cells)
    } else{
      error_messages = c(error_messages, "The parameter 'samples_min_cells' must be a number specifying the minimum number of cells a sample must have. Please set to NULL if there is no minimum!")
    }
  } else {
    param$samples_min_cells = 10
  }
  
  # Check cc_remove  ###
  if ("cc_remove" %in% names(param)) {
    if (converts_to_logical(param$cc_remove)) {
      param$cc_remove = as.logical(param$cc_remove)
    } else {
      error_messages = c(error_messages, "The parameter 'cc_remove' (correct for cell cycle) is missing or is not a logical value!")
    }
  } else {
    param$cc_remove = FALSE
  }

  # Check cc_remove_all  ###
  if ("cc_remove_all" %in% names(param)) {
    if (converts_to_logical(param$cc_remove_all)) {
      param$cc_remove_all = as.logical(param$cc_remove_all)
    } else {
      error_messages = c(error_messages, "The parameter 'cc_remove_all' (remove all cell cycle) is missing or is not a logical value!")
    }
    
    if (param$cc_remove_all && !param$cc_remove) {
      error_messages = c(error_messages, "The parameter 'cc_remove_all' (remove all cell cycle)  cannot be set to TRUE while the parameter 'cc_remove' is set to FALSE!")
    }
    
  } else {
    param$cc_remove_all = FALSE
  }
  
  # Check cc_rescore_after_merge  ###
  if ("cc_rescore_after_merge" %in% names(param)) {
    if (converts_to_logical(param$cc_rescore_after_merge)) {
      param$cc_rescore_after_merge = as.logical(param$cc_rescore_after_merge)
    } else {
      error_messages = c(error_messages, "The parameter 'cc_rescore_after_merge' (rescore after merging/integrating multiple samples) is missing or is not a logical value!")
    }
    
    if (param$cc_remove_all && !param$cc_remove) {
      error_messages = c(error_messages, "The parameter 'cc_rescore_after_merge' (rescore after merging/integrating multiple samples)  cannot be set to TRUE while the parameter 'cc_remove' is set to FALSE!")
    }
  } else {
    param$cc_rescore_after_merge = FALSE
  }
  
  # Check vars_to_regress  ###
  if ("vars_to_regress" %in% names(param)) {
    param$vars_to_regress = as.character(param$vars_to_regress)
  } else {
    param["vars_to_regress"] = NULL
  }
  
  # Check latent_vars  ###
  if ("latent_vars" %in% names(param)) {
    param$latent_vars = as.character(param$latent_vars)
  } else {
    param["latent_vars"] = NULL
  }
  
  # Check integrate_samples
  if ("integrate_samples" %in% names(param) && is.list(param$integrate_samples)) {
    if (!"method" %in% names(param$integrate_samples) || !param$integrate_samples$method %in% c("single", "merge", "integrate")) {
      error_messages = c(error_messages, "The parameter 'integrate_samples' misses a 'method' entry that is one of: 'single' (only one dataset), 'merge' (just merge) or 'integrate' (integrate)!")
    }
    
    if ("method" %in% names(param$integrate_samples) && param$integrate_samples$method=="integrate") {
      if (!"dimensions" %in% names(param$integrate_samples)) {
        error_messages = c(error_messages, "The parameter 'integrate_samples' misses a 'dimensions' entry. Please specify the number of dimensions to include for integration!")
      } else if (!converts_to_number(param$integrate_samples$dimensions)) {
        error_messages = c(error_messages, "The 'dimensions' entry of the parameter 'integrate_samples' must be numeric!")
      } else {
        param$integrate_samples$dimensions = as.numeric(param$integrate_samples$dimensions)
      }
      
      if (!"reference" %in% names(param$integrate_samples)) {
        param$integrate_samples["reference"] = NULL
      }
      
      if (!"use_reciprocal_pca" %in% names(param$integrate_samples)) {
        param$integrate_samples[["use_reciprocal_pca"]] = FALSE
      }
    }
  } else {
    param$integrate_samples = list(method="merge")
  }

  # Check norm (normalisation)
  if ("norm" %in% names(param)) {
    if (!param$norm %in% c("RNA", "SCT")) {
      error_messages = c(error_messages, "The parameter 'norm' (normalisation method to use) must be one of: 'RNA', 'SCT'!")
    }
  } else {
    param$norm = "RNA"
  }

  # Check pc_n
  if ("pc_n" %in% names(param)) {
    if (converts_to_number(param$pc_n)) {
      param$pc_n = as.numeric(param$pc_n)
    } else {
      error_messages = c(error_messages, "The parameter 'pc_n' is not a numeric value!")
    }
  } else {
    param$pc_n = 10
  }
  
  # Check cluster_k
  if ("cluster_k" %in% names(param)) {
    if (converts_to_number(param$cluster_k)) {
      param$cluster_k = as.numeric(param$cluster_k)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_k' is not a numeric value!")
    }
  } else {
    param$cluster_k = 20
  }
  
  # Check umap_k
  if ("umap_k" %in% names(param)) {
    if (converts_to_number(param$umap_k)) {
      param$umap_k = as.numeric(param$umap_k)
    } else {
      error_messages = c(error_messages, "The parameter 'umap_k' is not a numeric value!")
    }
  } else {
    param$umap_k = 30
  }
  
  # Check cluster_resolution_test
  if ("cluster_resolution_test" %in% names(param)) {
    if (all(purrr::map_lgl(param$cluster_resolution_test, converts_to_number))) {
      param$cluster_resolution_test = as.numeric(param$cluster_resolution_test)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_resolution_test' does not contain numeric values!")
    }
  } else {
    param$cluster_resolution_test = c()
  }
  
  # Check cluster_resolution
  if ("cluster_resolution_test" %in% names(param)) {
    if (all(purrr::map_lgl(param$cluster_resolution_test, converts_to_number))) {
      param$cluster_resolution_test = as.numeric(param$cluster_resolution_test)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_resolution_test' does not contain numeric values!")
    }
  } else {
    param$cluster_resolution_test = c()
  }
  
  if ("cluster_resolution" %in% names(param)) {
    if (converts_to_number(param$cluster_resolution)) {
      param$cluster_resolution = as.numeric(param$cluster_resolution)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_resolution' does not contain a numeric value!")
    }
  } else {
    param$cluster_resolution = 0.5
  }

  # Check marker_padj
  if ("marker_padj" %in% names(param)) {
    if (converts_to_number(param$marker_padj)) {
      param$marker_padj = as.numeric(param$marker_padj)
    } else {
      error_messages = c(error_messages, "The parameter 'marker_padj' is not a numeric value!")
    }
  } else {
    param$marker_padj = 0.05
  }

  # Check marker_log2FC
  if ("marker_log2FC" %in% names(param)) {
    if (converts_to_number(param$marker_log2FC)) {
      param$marker_log2FC = as.numeric(param$marker_log2FC)
    } else {
      error_messages = c(error_messages, "The parameter 'log2fc' is not a numeric value!")
    }
  }
  
  # Check deg_contrasts
  if ("deg_contrasts" %in% names(param) & !is.null(param$deg_contrasts)) {
    if (is.character(param$deg_contrasts) && file.exists(param$deg_contrasts)) {
      deg_contrasts_table = openxlsx::read.xlsx(param$deg_contrasts)
    } else {
      deg_contrasts_table = param$deg_contrasts
    }
    
    if (!is.data.frame(deg_contrasts_table)) {
      error_messages = c(error_messages, "The parameter 'deg_contrasts' must be either a data.frame or an existing Excel file with a valid table in the first sheet!")
    } else {
      if (!all(c("condition_column","condition_group1", "condition_group2") %in% colnames(deg_contrasts_table))){
        error_messages = c(error_messages, "The table specified by parameter 'deg_contrasts' must contain the following columns: 'condition_column', 'condition_group1' and 'condition_group2'!")
      }
      
    }
  }

  # Check enrichr_padj
  if (!"enrichr_padj" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'p_enrichr' is missing!")
  }

  # Check col
  if (!"col" %in% names(param) || !param$col %in% colors()) {
    error_messages = c(error_messages, "The parameter 'col' is missing or not a valid colour!")
  }

  # Check col_palette_samples
  if (!"col_palette_samples" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'col_palette_samples' is missing!")
  }
  
  # Check col_palette_clusters
  if (!"col_palette_clusters" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'col_palette_clusters' is missing!")
  }
  
  # Add to param object and return
  if (length(error_messages) > 0) param[["error_messages"]] = error_messages
  return(param)
  
  # Check
}

#' Check Python Availability
#'
#' Verifies that Python is available and required modules are installed.
#'
#' @return Character vector of error messages, or empty vector if all checks pass.
#'
#' @importFrom reticulate py_available py_config py_module_available
#' @export
check_python = function() {
  error_messages = c()
  
  if (!reticulate::py_available(initialize = TRUE) || is.null(reticulate::py_config())) {
    return("Python is not installed on this system or not found in the specified path!")
  }
  
  python_modules = c("leidenalg", "anndata", "scipy")
  is_available = purrr::map_lgl(python_modules, reticulate::py_module_available)
  if (any(!is_available)) {
    error_messages = c(error_messages, paste0("The following python packages are missing: ", paste(python_modules[!is_available], sep=", "),"!"))
  }

  return(error_messages)
}

#' Check Pandoc Availability
#'
#' Verifies that Pandoc is available for document rendering.
#'
#' @return Character vector of error messages, or empty vector if Pandoc is found.
#'
#' @importFrom rmarkdown find_pandoc
#' @export
check_pandoc = function() {
  error_messages = c()
  
  if (length(rmarkdown::find_pandoc()) == 0) {
    return("Pandoc is not installed on this system or not found in the specified path!")
  }

  return(error_messages)
}

#' Check EnrichR Availability
#'
#' Verifies that EnrichR is available and specified databases exist.
#'
#' @param databases Character vector. EnrichR databases to check.
#' @param site Character. EnrichR site to use. Default is \code{"Enrichr"}.
#'
#' @return Character vector of error messages, or empty vector if all checks pass.
#'
#' @importFrom enrichR setEnrichrSite listEnrichrDbs
#' @export
check_enrichr = function(databases, site="Enrichr") {
  if(is.null(databases) || length(databases)==0) return(c())
  
  # Is enrichR live
  if (is.null(options("enrichR.live")) || !options("enrichR.live")[[1]]) {
    return("EnrichR is not available or cannot connect to the databases")
  }
  
  # Set Enrichr site
  suppressMessages(enrichR::setEnrichrSite(param$enrichr_site))

  # Are databases available at all
  available_databases = tryCatch({ enrichR::listEnrichrDbs()[,"libraryName"] }, error=function(e) {return(NULL) })
  if (is.null(available_databases)) return("Could not list databases available at enrichR! Please check the enrichR vignette!")
  
  
  # Are the requested databases available
  are_valid = databases %in% available_databases
  if (any(!are_valid)) {
    return(paste0("The following enrichR databases are not available: ",paste(databases[!are_valid],collapse=", "),"!"))
  }
  
  return(c())
}

#' Check if Packages are Installed
#'
#' Checks whether specified R packages are installed.
#'
#' @param packages Character vector. Package names to check.
#'
#' @return Logical vector indicating which packages are installed.
#'
#' @export
#'
#' @examples
#' packages_installed(c("ggplot2", "nonexistent_package"))
packages_installed = function(packages) {
  return(packages %in% installed.packages()[ , "Package"])
}


#' Check Required Packages for scrnaseq Workflow
#'
#' Checks that all required packages for the scrnaseq workflow are installed.
#'
#' @return Character vector of error messages, or empty vector if all packages
#'   are installed.
#'
#' @export
check_installed_packages_scrnaseq = function() {
  required_packages = c("Seurat", "ggplot2", "patchwork", "magrittr",
                        "reticulate", "enrichR", "future", "knitr",
                        "dplyr", "tidyr", "purrr", "stringr", "sctransform", 
                        "Matrix", "kableExtra", "DT", "ggsci",
                        "openxlsx", "readr", "R.utils", "biomaRt",
                        "MAST", "enrichR", "sessioninfo", "cerebroApp",
                        "knitcitations", "sceasy")
  
  is_installed = packages_installed(packages=required_packages)
  if(any(!is_installed)) {
    return(paste0("The R packages '", required_packages[!is_installed],"' are not installed!"))
  } else {
    return(c())
  }
}

#' Check Ensembl Availability
#'
#' Checks whether Ensembl is available for downloading annotations and
#' cell cycle markers.
#'
#' @param biomart Character. BioMart database name.
#' @param dataset Character. Dataset name.
#' @param mirror Character. Ensembl mirror to use.
#' @param version Character. Ensembl version.
#' @param attributes Character vector. Attributes to verify.
#' @param file_annot Character or \code{NULL}. Path to existing annotation file.
#' @param file_cc_markers Character or \code{NULL}. Path to existing cell cycle
#'   markers file.
#'
#' @return Character vector of error messages, or empty vector if all checks pass.
#'
#' @importFrom biomaRt listAttributes
#' @export
check_ensembl = function(biomart, dataset, mirror, version, attributes, file_annot=NULL, file_cc_markers=NULL) {
  error_messages = c()
  
  if (is.null(file_annot) || !file.exists(file_annot) || is.null(file_cc_markers) || !file.exists(file_cc_markers)) {
    # See if mart is available
    annot_mart = suppressWarnings(GetBiomaRt(biomart, dataset, mirror, version))
    if (is.null(annot_mart)) {
      if (is.null(mirror)) {
        return(paste0("Cannot download Ensembl annotation for dataset '",dataset,"', version '",version ,"' using biomaRt!"))
      } else {
        return(paste0("Cannot download Ensembl annotation for dataset '",dataset,"', version '",version ,"' at mirror '",mirror,"' using biomaRt!"))
      }
    }
  
    # See if attributes are valid
    attributes = unique(attributes)
    available_attributes = biomaRt::listAttributes(annot_mart)[,1]
    is_available = attributes %in% available_attributes
    if (any(!is_available)) {
      error_messages = c(error_messages, paste0("The following Ensembl attributes could not be found using biomaRt: ", paste(attributes[!is_available], sep=", "),"!"))
    }
  } else {
    file_annot_tbl = read.delim(file_annot, nrows=50)
    is_available = attributes %in% colnames(file_annot_tbl)
    if (any(!is_available)) {
      error_messages = c(error_messages, paste0("The existing annotation file '", file_annot,"' misses the following Ensembl attributes: ", paste(attributes[!is_available], sep=", "),"!"))
    }
  }
  
  
  return(error_messages)
}

#' Set Error Handler to Print Traceback
#'
#' Configures R to print a traceback on error instead of starting a debugger.
#' Useful for non-interactive use.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @importFrom rlang caller_env trace_back
#' @export
on_error_just_print_traceback = function(x) {
  options(rlang_trace_top_env = rlang::caller_env())
  options(error = function() {
    sink()
    print(rlang::trace_back(bottom = sys.frame(-1)), simplify = "none")
  })
  return(invisible(NULL))
}

#' Set Error Handler to Start Terminal Debugger
#'
#' Configures R to start a debugger on the terminal on error.
#' Useful for interactive use without X11.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @export
on_error_start_terminal_debugger = function(x) {
  options(error = function() {
    sink()
    recover()
  })
  return(invisible(NULL))
}

#' Set Default Error Handler
#'
#' Keeps R's default error handling behavior.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @export
on_error_default_debugging = function(x) {
  return(invisible(NULL))
}

#' Cite a Reference
#'
#' Wrapper around \code{knitcitations::citet} and \code{knitcitations::citep}
#' with error handling for connection problems.
#'
#' @param reference Reference to cite (DOI, URL, or BibTeX key).
#' @param type Character. Citation type: \code{"citet"} for textual citation
#'   or \code{"citep"} for parenthetical citation. Default is \code{"citet"}.
#'
#' @return Formatted citation string, or the original reference if citation
#'   fails.
#'
#' @importFrom knitcitations citet citep
#' @export
#'
#' @examples
#' \dontrun{
#' Cite("10.1038/s41592-019-0654-x")
#' }
Cite = function(reference, type="citet") {
  formatted = tryCatch({
    if (type=="citet") knitcitations::citet(reference) else knitcitations::citep(reference)
  },
  error=function(cond) {
    return(NULL)
  })
  
  if (is.null(formatted)) formatted = reference
  return(formatted)
}

