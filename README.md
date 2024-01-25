# Introduction

**scrnaseq2** is a bioinformatics workflow to analyze single-cell data, and it is written by the bioinformatics team of the DRESDEN-concept Genome Center. 

The workflow is written in R  and is largely based on the [Seurat v5](https://satijalab.org/seurat/) R package. It follows best practices recommendations ([Luecken and Theis, 2019](https://www.embopress.org/doi/full/10.15252/msb.20188746), and [Heumos et al., 2023](https://www.nature.com/articles/s41576-023-00586-w)) and contains additional visualisations, tables and documentation for a better understanding. The workflow supports RNA sequencing data from one or more samples processed with 10X Genomics, SmartSeq, Parse Biosciences and Scale Bio. 

The workflow generates an extensive HTML report. Are you curious about what that looks like and if it would be useful for your own data? If so, you can  download the GitHub repository and open `docs/example_scrnaseq2.html`. This report has been generated for test data as mentioned [below](#quick_start). 

If you are a researcher, and you would like to start analysing your own data, the workflow can be your starting point. If you work in a *bioinformatics core facility* and frequently support other researchers with bioinformatics analyses, the workflow can be run in a standardised fashion both interactively in `Rstudio` and on command line. We typically first run the workflow with default parameters, and communicate with our collaborators. We then optimise the parameters in further rounds to improve the results. 

# News
<a name="news"></a>

(2024-01-15) We started working intensively on this new workflow. Stay tuned, updates keep coming in the next days. 

# Workflow summary
<a name="workflow_summary"/>

## Workflow: Single-cell RNA-seq analysis 
UPDATE

# Quick start
UPDATE

The workflow is inialised for test data in `test_datasets`. First, navigate to the respective test dataset folder(s), and download the test dataset(s) by running the `download.R` script(s). Once all test data is downloaded, you can knit the workflow to HTML. 

The repository provides several other useful test data that you can use to get to know the functionality of the workflow. To run the workflow for another than the initial dataset, you need to adapt the `project_parameters` code chunk and provide all relevant paths and parameters. 


# Documentation 
UPDATE 

The scrnaseq workflow comes with a good amount of documentation, found in the `docs/` directory:
 
[Installation](docs/installation.md)   
[Running the workflow](docs/usage_workflow.md)   
[Running the pre-workflow](docs/usage_preworkflow.md)   

# Credits

The [Seurat Vignettes](https://satijalab.org/seurat/vignettes.html) were initially used as templates for this workflow. 

The workflow was developed by [Katrin Sameith](https://github.com/ktrns), [Andreas Petzold](https://github.com/andpet0101) and [Ulrike Friedich](https://github.com/ulrikefriedrich) at the [DRESDEN-concept Genome Center (Dresden, Germany)](https://genomecenter.tu-dresden.de/about-us). **scrnaseq2** is based on [**scrnaseq**](https://github.com/ktrns/scrnaseq). 

# Contributions and Support

If you would like to contribute to this workflow, please first create your own fork of the GitHub repository. You can then work on your master branch, or create feature branches for developement that you later merge into your master branch. Once your code is finalised and working, you can create a pull request. 

# Citation

If you used **scrnaseq2** to analyze your data, please cite it by mentioning the DRESDEN-concept Genome Center URL "https://genomecenter.tu-dresden.de". (TODO: Zenodo!)

