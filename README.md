# Single cell RNAseq data analysis bundle

This repository contains a series of pipelines that can be used for processing, analysing and exploratory analysis of single cell RNA sequencing transcriptomics data using a server running on a Grid Engine. The scripts can also be used on personal machines and/or in interactive mode or by adapting the shel scripts to be used on other cluster softwares. 

These pipelines were developed for the study __Decoding the development of the blood and immune systems during human fetal liver haematopoiesis__ and due to their wider application they are being used and further extended for other projects.

The bundle includes:
* tools for building data sets from multiple CellRanger count tables
* data annotation aids
* training machine learning cell type classifiers
* doublet removal
* computing data reduction coordinates like tSNE, UMAP, force directed graph and diffusion map
* trajectory analysis
* interactive plots
* batch correction
* signature genes
* cell type comparison metrics and plots
* creating animated force directed graphs
* handling big data sets as Seurat objects
 
To create mode advanced data exploration apps and web portal visit https://github.com/DoruMP/Fast-data-portals-for-scRNAseq-data

To run the pipelines you must download the entire bundle and transfer to a server/personal computer. The folder structure must be kept as it is.

## A proposed standard work flow (the order of steps is mandatory):

* create the seurat object (seurat_from_count_tables.sh)
* train a doublet detector (train_doublets_SVM.sh)
* update doublet assignment in data (can use Apply_Classifier_On_Seurat_Object function in a customised script or, the easier approach is to run again the seurat object creation process with seurat_from_count_tables.sh and this time setting the argument identify.doublets = TRUE)
* subset the seurat object into singlets and doublets (split_seurat_by_category.sh). It is advised to keep the doublets data (even if you have not further use for them) for future reference.
* compute dimensionality reduction coordinates (add_dr.sh)
* make an annotation template and annotate (make_cell_annotation_template.sh)
* optional: annotation can be made a lot easier by plotting annotation assisting word clouds (wordclouds.sh) and making an interactive heat map (interactive heat map tool portal_tools/interactive_heatmap_dotplot).
* important notice: it is highly advisable to use only alphanumeric characters in naming cell population. Some characters (e.g. "\", "/") can create problems and raise errors with some pipelines. While many of these issues are solved for, it is still advisable as good practice to avoid fancy characters in naming. This is because it is imposible to predict all possible issues created by non-alphanumeric characters and even if when they do trigger errors, the error messages are particularly vague in such situations. Here alpha-numeric is defined as the collection of Latin letters and Arabic digits.
* update the annotation (update_annotation.sh)
* make all the apps that allow easy data exploration:
* plot_dr.sh for interactive UMAP, FDG, tSNE and AGA; 
   * portal_tools/interactive_heatmap_dotplot interactive heatmap but this time with the labels, not clusters; 
   * portal_tools/web_portal web portal tool for gene expression (you must have access to web server or alternatively set an Apache server on your local machine)
   * portal_tools/gene_grouping for gene expression patterns
   * super_markers.sh get cell types signatures. You will understand the power of these signature when you input them in the interactive heat map (if you make one). This will be very useful for annotating new data, for supporting data annotation or for checking expression patterns in new data sets
   * optional: you could run again the word clouds because this also show DEGs as word clouds
   * optional: you could train a cell type classifier for fast integration of new data sets. However it is highly recommended that any published conclusions should be made on whole data annotation, not on classifiers results from this bundle (or any other type of label propagations/projections made with tools that are not part of this bundle)
* it is recommended that you treat portals, doublet detectors, cell type classifiers and gene signatures as resources not as results. You should only share such resources with relevant people otherwise you might risk leaking results to others before publication.
* next steps are project specific

## Prerequisites

Python version 3.6
R version 3.4.2

