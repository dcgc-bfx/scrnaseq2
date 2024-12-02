# Dataset provided by 10x: 1.3 Million Brain Cells from E18 Mice (v3 chemistry)
# https://www.10xgenomics.com/datasets/1-3-million-brain-cells-from-e-18-mice-2-standard-1-2-0

# Important: Run this script in its directory

unlink("filtered_feature_bc_matrix.h5", recursive=T)

# download and untar
url = 'https://cf.10xgenomics.com/samples/cell-exp/1.2.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5'
curl::curl_download(url=url, destfile="filtered_feature_bc_matrix.h5")

# get metrics summary
url = 'https://cf.10xgenomics.com/samples/cell-exp/1.2.0/1M_neurons/1M_neurons_reanalyze.csv'
curl::curl_download(url=url, destfile='metrics_summary.csv')

# get some analysis results and use "analysis/clustering/graphclust/clusters.csv" as barcode metadata
url = 'https://cf.10xgenomics.com/samples/cell-exp/1.2.0/1M_neurons/1M_neurons_analysis.tar.gz'
curl::curl_download(url=url, destfile=basename(path='analysis.tar.gz'))
untar(tarfile = basename(path=url))
unlink(basename(path=url))

