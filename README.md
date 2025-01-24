This directory contains Python scripts and Jupyter notebooks I wrote for simulations and data analysis.
## random_walker.ipynb
This notebook simulates and visualizes random walk processes for multiple entities, analyzing the resulting position distributions through Probability Density Functions (PDFs) and Cumulative Probability Functions (CPFs).
## vaccine_simulation.ipynb
This notebook simulates the effectiveness of vaccination using permutation testing to assess the likelihood of observing a statistically significant reduction in infection rates in a treated group compared to a control group.
## motif_finding.ipynb
This notebook implements various bioinformatics algorithms, including functions for analyzing genetic sequences. It covers tasks like calculating Hamming distance, translating RNA into proteins, finding DNA motifs, RNA splicing, identifying common substrings or subsequences, and constructing shortest common supersequences.
## human_muscle_stem_cell_sc_rna_seq.ipynb
This notebook implements a computational pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data of muscle stem cells (MuSCs) using ```scanpy``` and ```scvi-tools```.  
The pipeline starts by loading datasets, filtering the data to focus on MuSCs, and storing raw count data for reference. It then normalizes and log-transforms the gene expression data to ensure comparability across cells.  
Next, the code identifies the top 10,000 highly variable genes using the Seurat v3 method and prepares the data for latent representation analysis. A Variational Autoencoder (VAE) is trained using ```scVI``` to learn a 30-dimensional latent space representation of the data, enabling downstream clustering and visualization.  
Clustering is performed using the Leiden algorithm, and dimensionality reduction is achieved with UMAP to create a low-dimensional representation of the data. The processed data is saved for future use.  
Finally, the pipeline generates visualizations, including a UMAP plot to display cell clusters and a dot plot to examine the expression of key marker genes across specified cell groups. This workflow provides a comprehensive approach to studying single-cell transcriptomics data of MuSCs.
