---
# Module-specific parameters (that are not already in the profile yaml)
params:
  # Name of the module. Must be the same as name of qmd file.
  module: "index"

# Module execution
execute:
  # Should this module be frozen and never re-rendered?
  # - auto: if the code has not changed (default)
  # - true/false: freeze/re-render
  # Does not apply in interactive mode or when explicitly rendering this document via in rstudio
  # Set to false for this page since it can always be rendered
  freeze: false
---



## Introduction {.unnumbered}


::: {.cell}

```{.r .cell-code}
# If running code interactively in rstudio, set profile here
# When rendering with quarto, profile is already set and should not be overwritten
if (nchar(Sys.getenv("QUARTO_PROFILE")) == 0) {
  Sys.setenv("QUARTO_PROFILE" = "default")
}

# Source general configurations (always)
source("R/general_configuration.R")

# Source required R functions
source("R/functions_util.R")
source("R/functions_io.R")

# Load libraries
library(magrittr)
suppressMessages(library(Seurat))

# Get module name and directory (needed to access files within the module directory)
module_name <- param("module")
module_dir <- file.path("modules", module_name)
```
:::


This main page ...

**Important**: When starting a project, please fill out the details in the next section. Furthermore, use the section to track all tasks.

## Scientific project {.unnumbered}

#### Description {.unnumbered}

For this tutorial, we will be analyzing two datasets of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics.

-   two datasets, 1k and 5k cells

-   10x 3' assay

-   blood cells

#### Questions {.unnumbered}

-   How many clusters do we find?

-   What are the most abundant cell types?

-   How many T-cells do we find?

#### TODOs {.unnumbered}

-   ~~Rerun clustering with relaxed parameter~~

-   Map clusters to common cell types

-   Find DEGs

## Content {.unnumbered}

## Common questions {.unnumbered}

-   Common question1
-   Common question2
-   Common question3

## About the workflow {.unnumbered}

[**scrnaseq**](https://github.com/ktrns/scrnaseq) contains bioinformatics workflows for the analysis for single-cell data. The workflows are written in R and largely based on the R package Seurat. They follow best practices recommendations [@Luecken_2019; @Heumos_2023] and contain additional visualisations, tables and documentation to better understand the analysis. They support single-cell data from one or more samples processed with 10X Genomics, SmartSeq, Parse Biosciences and Scale Bio.

The workflow generates an extensive HTML report. Are you curious about what that looks like and if it would be useful for your own data? If so, you can download the GitHub repository and open scrnaseq.html. This report has been generated for test data as mentioned below.

If you are a researcher, and you would like to start analysing your own data, the workflow can be your starting point. If you work in a bioinformatics core facility and frequently support other researchers with bioinformatics analyses, the workflow can be run in a standardised fashion both interactively in Rstudio and on command line. We typically first run the workflow with default parameters, and communicate with our collaborators. We then optimise the parameters in further rounds to improve the results.

## How to run {.unnumbered}

Do this to run it by yourself.

## How to cite {.unnumbered}

If scrnaseq was used to analyse single-cell data that is part of your publication, please cite it by mentioning the Dresden-concept Genome Center URL "https://genomecenter.tu-dresden.de". ZENODO!