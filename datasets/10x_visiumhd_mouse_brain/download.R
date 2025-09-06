# Dataset provided by 10x: Visium HD Spatial Gene Expression Library, Mouse Brain (FFPE)
# https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he-v4

# Important: Run this script in its directory


# download and untar
url = 'https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_binned_outputs.tar.gz'
curl::curl_download(url=url, destfile='Visium_HD_Mouse_Brain_binned_outputs.tar.gz')
untar('Visium_HD_Mouse_Brain_binned_outputs.tar.gz')
unlink('Visium_HD_Mouse_Brain_binned_outputs.tar.gz')

url = 'https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_segmented_outputs.tar.gz'
curl::curl_download(url=url, destfile='Visium_HD_Mouse_Brain_segmented_outputs.tar.gz')
untar('Visium_HD_Mouse_Brain_segmented_outputs.tar.gz')
unlink('Visium_HD_Mouse_Brain_segmented_outputs.tar.gz')

url = 'https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_metrics_summary.csv'
curl::curl_download(url=url, destfile='metrics_summary.csv')
