# demo_R
Re-analysis of RNASeq dataset combining UMAP, k-means and co-expression network analysis

Whole-transcriptome sequencing technologies, such as microarray and RNA-Seq, provide an effective and comprehensive description of gene expression profiles over time in each tissue.
Analysis of differentially expressed genes is typically centered on identifying genes with differential expression, but no study has explored how the genes are interconnected.
Identification of genes with correlated expression can shed more light on their possible functions, as genes with similar expression patterns might be related.
Moreover, DGE analysis only within given groups without knowledge of sample-to-sample heterogeneity within those groups can often lead to erroneous or biased results.
Therefore, it is imperative to analyze the sample-to-sample heterogeneity within groups in order to identify outliers or subgroups. Only when such information is available can analytical methods be used to correct batch effects, remove outliers, and distinguish subgroups.
Here, a new analysis was conducted on the publicly available GSE98212 dataset, using a combination of UMAP and k-mean clustering to unable to shed light on hidden structures inside the dataset. Then, the result was integrated into a weighted gene co-expression network analysis approach, in order to show the network of genes related to samples groups. 
