# RWRF and RWRNF

## “Multi-dimensional Data Integration Algorithm Based on Random Walk with Restart”


Abstract：The accumulation of various multi-omics data and computational approaches for data integration can accelerate the development of precision medicine. 
However, the algorithm devel-opment for multi-omics data integration remains a pressing challenge. We propose a multi-omics data integration algorithm based on random walk with restart (RWR) on heterogeneous network. We call the resulting methodology RWRF (Random Walk with Restart for multi-dimensional data Fusion). RWRF uses similarity network of samples as the basis for integration. It constructs the similarity network for each data type and then connects corresponding samples of multiple similarity networks to create a heterogeneous sample network. By applying RWR on the heterogeneous network, RWRF uses stationary probability distribution to fuse similarity networks. We applied RWRF to TCGA data to identify subtypes in different cancer data set. Three types of data (mRNA expression, DNA methylation, and microRNA expression data) are integrated and network clustering is conducted. Experiment results show that RWRF performs better than sin-gle data type analysis and previous integrative methods.

# Get Started

## Example Datasets

To get started, you need to download example datasets (3 types of adrenocortical carcinoma data: mRNA expression, DNA methylation, and microRNA (miRNA) expression) from TCGA:  

[TCGA](https://portal.gdc.cancer.gov)

## Run Example

example.py: Three types of data (mRNA expression, DNA methylation, and microRNA expression data) are integrated and network clustering is conducted. Dunn and P value for the log-rank test of survival analysis is calculated.





