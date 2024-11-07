# Dataset provided by 10x: Xenium In Situ preview data generated using the Mouse Tissue Atlassing panel on a one day old mouse pup
# https://www.10xgenomics.com/datasets/mouse-pup-preview-data-xenium-mouse-tissue-atlassing-panel-1-standard


# Important: Run this script in its directory
files = c("analysis_summary.html", "analysis.zarr.zip", "cell_boundaries.csv.gz", "cell_boundaries.parquet", "cell_feature_matrix.h5",
          "cell_feature_matrix.zarr.zip", "cells.csv.gz", "cells.parquet", "cells.zarr.zip", "experiment.xenium", "gene_panel.json",
          "metrics_summary.csv", "morphology_focus.ome.tif", "morphology_mip.ome.tif", "morphology.ome.tif", "nucleus_boundaries.csv.gz",
          "nucleus_boundaries.parquet", "transcripts.csv.gz", "transcripts.parquet", "transcripts.zarr.zip")
unlink(files)
dirs = c("analysis", "cell_feature_matrix") 
unlink(dirs, recursive=TRUE)

# download and unzip
url = 'https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/1.6.0/Xenium_V1_mouse_pup/Xenium_V1_mouse_pup_outs.zip'
curl::curl_download(url=url, destfile=basename(path=url))
utils::unzip(basename(path=url))
unlink(basename(path=url))

untar(tarfile = "analysis.tar.gz")
unlink(basename(path=url))