# Readme

## Frequency of co-expression networks: 
You will need to install OutDeCo and download the networks. Note, these are large (GBs)! As this is currently in development, there may be some changes in the future. 
```
install.packages("devtools")
devtools::install_github("sarbal/OutDeCo") 

download_network_hdf5(network_type="generic") 
download_network_hdf5(network_type="blood") 
download_network_hdf5(network_type="brain") 
```

## GIANT networks
Run scripts: get_GIANT_networks.sh to download ~143 tissue specific networks
Note: These files are large (GBs!) and this will take some time.

## RNA-seq count data from TAF1 pedigrees
Download the counts_data.Rdata file for the counts mapped and annotated with GENCODE v22 (processed_count_data.txt). 
For version v24, the counts text file is given:  processed_count_data.GENCODEV24.txt
Or, you can access the raw data from: GSE84891 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84891).

## Gene annotation file 
Contains the gene names and IDs of the genes used in the original analysis (GENCODE v22).
gene_annotations_v22.Rdata

## GTEX
Data obtained from version 8 of GTEX. 
These files were downloaded and parsed:  
Phenotypes: 
https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx
https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx
https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

Expression data: 
https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz




