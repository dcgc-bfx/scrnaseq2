#!/usr/bin/env Rscript

#QUARTO_PROJECT_RENDER_ALL	Set to “1” if this is a render of all files in the project (as opposed to an incremental render or a render for preview). This unset if Quarto is not rendering all files.
#QUARTO_PROJECT_OUTPUT_DIR	Output directory
#QUARTO_PROJECT_INPUT_FILES	Newline separated list of all input files being rendered (passed only to pre-render)
#QUARTO_PROJECT_OUTPUT_FILES	Newline separated list of all output files rendered (passed only to post-render).

# Libraries
library(magrittr)

# Splits a path
split_path = function(x) {
  if (dirname(x)==x) return(x) else return(c(split_path(dirname(x)), basename(x)))
}

# The environment variable QUARTO_PROJECT_OUTPUT_FILES should contains a newline-separated list of all quarto output files rendered 
quarto_output_files = Sys.getenv("QUARTO_PROJECT_OUTPUT_FILES")
quarto_project_output_dir = Sys.getenv("QUARTO_PROJECT_OUTPUT_DIR")

if (nchar(quarto_output_files) > 0 & nchar(quarto_project_output_dir) > 0) {
  quarto_output_files = stringi::stri_split_lines(quarto_output_files, omit_empty=TRUE) %>% unlist()

  # Find all modules that were rendered
  modules_rendered = purrr::map(quarto_output_files, function(f) {
    m = NULL
    
    # Split path and check whether it was produced as part of a module
    f = split_path(f)
    i = which(f == "modules")
    if (length(i) == 0) return(NULL)
    i = i[1]
    j = length(f) - 1
    
    # Return the path to the module directory
    m = do.call(file.path, as.list(f[i:j]))
    return(m)
  }) %>% purrr::flatten_chr() %>% unique()
  
  # Clear and make 'results' directory
  results_dir = file.path(quarto_project_output_dir, "results")
  unlink(results_dir, recursive=TRUE)
  dir.create(results_dir, showWarnings=FALSE)
  
  # Make 'figures' directory
  figures_dir = file.path(results_dir, "figures")
  dir.create(figures_dir, showWarnings=FALSE)
  
  # Now copy files located in the 'results'  directories of the rendered modules
  # Also copy figures
  for (m in modules_rendered) {
    module_results_files = list.files(file.path(m, "results"), full.names=TRUE)
    if (length(module_results_files) > 0) {
      for(f in module_results_files) file.copy(f, results_dir, recursive=TRUE)
    }
    
    module_figures_files = list.files(file.path(quarto_project_output_dir, m), full.names=TRUE, pattern="(png|pdf)$", recursive=TRUE)
    if (length(module_figures_files) > 0) {
      for(f in module_figures_files) file.copy(f, figures_dir)
    }
  }
}
