# Dataset provided by 10x: 1.3 million cells from cortex, hippocampus and subventricular zone of two E18 mice (v3 chemistry)
# https://www.10xgenomics.com/datasets/1-3-million-brain-cells-from-e-18-mice-2-standard-1-2-0

# Important: Run this script in its directory

unlink("filtered_feature_bc_matrix.h5", recursive=T)

# download
url = 'https://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5'
curl::curl_download(url=url, destfile="filtered_feature_bc_matrix.h5")
