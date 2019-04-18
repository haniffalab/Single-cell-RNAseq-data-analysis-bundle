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
* interactive plots as html pages
* batch correction
* signature genes
* cell type comparison metrics and plots
* creating animated force directed graphs
* handling big data sets as Seurat objects
 
To create mode advanced data exploration apps and web portal visit https://github.com/DoruMP/Fast-data-portals-for-scRNAseq-data. This repository will be referred in this documentation as __Fast Portals__.

To run the pipelines you must download the entire bundle and transfer to a server/personal computer. The folder structure must be kept as it is.

## A proposed standard work flow (the order of steps is mandatory):

* create the seurat object (_seurat\_from\_count\_tables.sh_)
* train a doublet detector (_train\_doublets\_SVM.sh_)
* update doublet assignment in data (can use _Apply\_Classifier\_On\_Seurat\_Object(..)_ function in a customised script or, the easier approach is to run again the seurat object creation process with _seurat\_from\_count\_tables.sh_ and this time setting the argument _identify.doublets = TRUE_)
* subset the seurat object into singlets and doublets (_split\_seurat\_by\_category.sh_). It is advised to keep the doublets data (even if you have not further use for them) for future reference.
* compute dimensionality reduction coordinates (_add\_dr.sh_)
* make an annotation template and annotate (_make\_cell\_annotation\_template.sh_)
* optional: annotation can be made a lot easier by plotting annotation assisting word clouds (_wordclouds.sh_) and making an interactive heat map (interactive heat map tool at  __Fast Portals__ _interactive\_heatmap\_dotplot_).
* important notice: it is highly advisable to use only alphanumeric characters in naming cell population. Some characters (e.g. "\\", "/") can create problems and raise errors with some pipelines. While many of these issues are solved for, it is still advisable as good practice to avoid fancy characters in naming. This is because it is imposible to predict all possible issues created by non-alphanumeric characters and even when they do trigger errors, the error messages are particularly vague in such situations. Here alpha-numeric is defined as the collection of Latin letters and Arabic digits.
* update the annotation (_update\_annotation.sh_)
* make all the apps that allow easy data exploration:
   * _plot\_dr.sh_ for interactive UMAP, FDG, tSNE and AGA; 
   * __Fast Portals__ _interactive\_heatmap\_dotplot_ interactive heatmap but this time with the labels, not clusters; 
   * __Fast Portals__ _web\_portal_ web portal tool for gene expression (you must have access to web server or alternatively set an Apache server on your local machine)
   * __Fast Portals__ _gene\_grouping_ for gene expression patterns
   * _super\_markers.sh_ gets cell types signatures. You will understand the power of these signature when you input them in the interactive heat map (if you make one). This will be very useful for annotating new data, for supporting data annotation, for exploring expression patterns in new data sets or for designing new flow panels
   * optional: you could run again the word clouds because this also show DEGs as word clouds
   * optional: you could train a cell type classifier for fast integration of new data sets. However it is highly recommended that any published conclusions should be made on whole data annotation, not on classifiers results from this bundle (or any other type of machine learning classification/label propagations/projections made with tools that are not part of this bundle)
* it is recommended that you treat portals, doublet detectors, cell type classifiers and gene signatures as resources not as results. You should only share such resources with relevant people otherwise you might risk leaking results to others before publication.
* next steps are project specific

## Prerequisites

Python version 3.6
>pandas 0.22.0
>pptx 0.6.9
>patsy 0.5.0
>scanpy 1.2.2
>sklearn 0.19.1
>numpy 1.14.2
>scipy 1.0.0
>umap 0.2.3
 
R version 3.4.2
>Seurat 2.3.4
>dplyr 0.7.6
>reshape 0.8.8
>plyr 1.8.4
>ggplot2 3.0.0
>RColorBrewer 1.1.2
>BiocParallel 1.12.0
>gridExtra 2.3
>grid 3.4.2
>sva 3.26.0
>destiny 2.6.2
 ggplot2 3.0.0
 monocle 2.6.4
 harmony 0.1.0
 methods 3.4.2
 utils 3.4.2
 wordcloud 2.6
 
 ## Structure
 
 The main folder is called _single\_cell\_data\_analysis\_bundle_ and must contain:
 * data folder where all seurat object as kept as RDS file and scanpy objects as h5ad files
 * output folder where jobs save their output. This is where the user can get the results of running a pipeline
 * pipelines folder contain a folder for each pipeline
 * resources folder containing sample key, colour keys, cell type classifier, doublet detectors, options file etc.
 * tools folder - there is no reason why the user should be concern with this folder. It contains programs for running force-directed graph, AGA, UMAP, classifier prediction, doublet detection, pseudotime 3D viewer app builder and a file with lots of utilities. These are never required to be called directly by the user. 
 * the tools folder includes the file bunddle_utils.R. 
 * if you want to use the portal tools and fast gene expression explorer you must have an additional folder named portal_tool where these tools are stored and can be used as pipelines

## Colour keys

 Color keys compatible with the single cell analysis bundle can be generated using the interactive tool _color\_management.html_. Instructions can be found inside the interactive tool if opened in a browser (recommended: Chrome, Firefox; to avoid if possible: all versions of Internet Explorer and all versions of Microsoft Edge).

## Utilities

The functions in the _bunddle\_utils.R_ can be used in new pipelines or in customized scripts. If this is required check the parameter description for the required function and ensure access to required scripts that it need to call.

The _bunddle\_utils.R_:
* declares the tool.addr and python.addr variables. Change the python.addr if you want to use a different python version for your work but make sure first that the version you are trying to use has installed all the required packages.
* the functions _runFDG_, _RunUMAP_ are used to compute force-directed graph and UMAP coordinates for a seurat object. Although currently Seurat package has a function to compute UMAP which goes by the same name, the function in this bundle was created before Seurat published its umap computing function. Both the in-house and the Seurat RunUMAP functions do the same thing but because the bundle was build before Seurat had the ability to compute UMAP it is recomended to use the RunUMAP from the _bunddle\_utils.R_ script with the current bundle for compatibility reasons.
* __runFDG(pca.df, snn, iterations = 600, tool\_addr, python.addr)__
    * computes force directed coordinates on a seurat object. This function requires time and computation resources for big data sets.
    * @param `pca.df` the input data frame with variables as columns. In the pipeline this is used on the pca coordinates of a seurat object which can be retrieved at `seurat.obj@dr$pca@cell.embeddings`. The function is very flexible due to this input and can used on many types data formats not limited to a seurat object. Further flexibility of this function comes from the fact that other embeddings can be used besides pca (e.g. batch corrected pca stored at `seurat.obj@dr$harmony@cell.embeddings`).
    * @param `snn` shared nearest neighbor graph. In a seurat object this is available at `seurat.obj@snn`. If you want to use this function outside a pipeline you must make sure that the snn has been computed on the seurat object or if you are not using a seurat object you must computed using other available tools. You must pay attention to seurat object subsetting which as of writing this documentation does not recompute the snn.
    * @param `tool_addr` folder name were the bundle tools are stored. this function uses _force\_abstract\_graph/make\_fdg.py_ for the actual computations. If you need to use this function outside the pipelines you must make sure you have this script and set the tool.addr argument properly. 
    * @param `python.addr` python address. This is pre-set in all pipelines but having this as an argument allows the user to re-use the function in other scripts and choose the python version
* __RunUMAP(pca.df, tool\_addr, python.addr)__
    * computes umap coordinates. Currently Seurat packages has published a function with same name that computs UMAP coordinates on a seurat object. However the RunUMAP in this bundle is more flexible and can be used not just on a seurat object. This function is the fastest dimension reduction method (compared to PCA, tSNE and FDG).
    * @param `pca.df` the input data frame with variables as columns. In the pipeline this is used on the pca coordinates of a seurat object which can be retrieved at `seurat.obj@dr$pca@cell.embeddings`. The function is very flexible due to this input and can used on many types data formats not limited to a seurat object. Further flexibility of this function comes from the fact that other embeddings can be used besides pca (e.g. batch corrected pca stored at `seurat.obj@dr$harmony@cell.embeddings`).
    * @param `tool_addr` folder name were the bundle tools are stored. this function uses _umap/umap\_compute.py_ stored in the tools folder to compute the umap coordinates.
    * @param `python.addr` python address. This is pre-set in all pipelines but having this as an argument allows the user to re-use the function in other scripts and choose the python version
* __Apply\_Classifier\_On\_Seurat\_Object(seurat.obj, classifier.fname, tool\_addr, python.addr)__
    * this function applies an SVM classifier to the seurat objects and outputs predictions as a character vector;
    * the predictions can be added to the meta-data of the seurat objects
    * can be used to predict cell types or doublets. Make sure that the cell type classifier or doublet detector is relevant for the data you want to predict (e.g. do not use a thymus cell type classifier on spleen)
    * @param `seurat.obj` name of seurat object - data must be normalized before applying the function
    * @param `classifier.fname` folder name were the classifier is stored. cell type classifier are trained with the _train\_classifier.sh_ pipeline and doublet detectors are obtained with _train\_doublets\_classifier.sh pipeline_;
    * @param `tool_addr` folder name were the bundle tools are stored. This function calls the _predict\_by_classifier/predict.py_ script to load the classifier and must be present in the tool.addr folder. Inside the pipelines this argument is already set but if you want to use this function in other scripts you must set the proper tool.addr
    * @param `python.addr` python address. This is pre-set in all pipelines but having this as an argument allows the user to re-use the function in other scripts and choose the python version
* __make\_3D\_interactive\_page(data\_frame\_3D, tool\_addr, python.addr, save.to)__
    * outdated since the web portal was created (see __Fast Portals__)
    * this function was used to create a 3D visualising html page. 
    * this function is the ancestor of the pseudotime web portal
    * this function is used by the _plot\_dr.sh_ pipeline to make a interactive tool for visualising diffusion map coordinates.
    * It can also used outside the pipeline to visualize other types of three-dimensional data sets (3D UMAP, 3D FDG, 3D tSNE) but this will be to be pre-computed
    * @param `data_frame_3D` data frame with the 3 dimensions in the first 3 columns, colours as hexdecimals in the forth columns and cell labels in 5th columns which must be named "Labels". The names of the other columns are not important. Make sue you do not have non-alphanumeric character in the labels (e.g. /, \, @ etc.) which can cause issues with the output.
    * @param `tool_addr` folder name were the bundle tools are stored. This function uses _interactive\_3D_viewer/html_WebGL\_3D\_viewer.py_ script which must be found in the tools folder. The function can be used outside the pipelines by setting the tool\_addr argument.
    * @param `python.addr` python address. This is pre-set in all pipelines but having this as an argument allows the user to re-use the function in other scripts and choose the python version
    * @param `save.to` address to save the resulted interactive page
* __make\_2D\_interactive\_page(data\_frame\_2D, tool\_addr, python.addr, save.to="./")__
    * this functions is similar to _make\_3D\_interactive\_page_ and it produces and interactive html page for vizualizing 2D data (UMAP, tSNE, FDG)
    * @param `data_frame_2D` has similar formating with the parameter data\_frame\_3D in function make\_3D\_interactive\_page the difference being that it features only two dimensions
    * @param `tool_addr` folder name were the bundle tools are stored. This function uses interactive\_2D\_viewer/html\_WebGL\_2D\_viewer.py script which be found in the tools folder.
    * @param `python.addr` python address. This is pre-set in all pipelines but having this as an argument allows the user to re-use the function in other scripts and choose the python version
    * @param `save.to` address to save the resulted interactive page
* __create\_gene\_expression\_viewer\_apps(seurat.obj, dim.type = 'umap', save.to, tool\_addr, python.addr, categories.colours=NA)__
    * deprecated
    * this has been replaced by the web portal creating tool (see __Fast Portals__)
    * this can still be used by the pipeline _seurat\_to\_interactive\_gene\_expression.R_
    * this creates html interactive pages that allow data exploration for dimensional reduction plots and gene expression. However the data is embedded in the page so the results is heavy and will require time to load in a browser. It is not recommended to used the output on a web server. The output is useful for internal distribution of data. For other situation I recommend you use the web portal building tool.
* __plot.indexed.legend(label.vector, color.vector, ncols = 2, left.limit = 3.4, symbol.size = 8, text.size = 10, padH = 1, padV = 1, padRight = 0)__
    * function that creates a legend pdf file to be used with dimension reduction plots
    * this function is called in _plot\_dr.sh_
    * @param `label.vector` a character vector with the cell labels
    * @param `color.vector` a character vector with the colors for the labels written in hexdecimal format. hexdecimal colour keys can be created using the interactive tool _color\_management.html_
    * @param `ncols` number of columns to arrange the labels in the legend
    * @param `left.limit` left padding
    * @param `symbol.size` symbol size
    * @param `text.size` text size
    * @param `padH` horizontal padding
    * @param `padV` vertical padding
    * @param `padRight` right padding
* __dr.plot(point.labels, dr1, dr2, dr1.name, dr2.name, no.legend = F, plt.lb.sz = 5, txt.lb.size = 3, pt.size = .2, random\_state = 2, use.cols = NULL, use.labels = NULL, limits = NULL, annotate.plot = T, index.map = NA)__
    * function used to plot dimensionality reduced coordinates (e.g. tSNE, UMAP)
    * @param `point.labels` character vector of all labels in the data set
    * @param `dr1` numeric vector of first dimension coordinates
    * @param `dr2` numeric vector of second dimension coordinates
    * @param `no.legend` boolean indicating if to plot the legend inside the final plot. Sometimes it is desirable to plot the legend separatly to allow for more flexibility in arranging figure panels. In this case set this argument to FALSE and use plot.indexed.legend to create a separate legend figure
    * @param `plt.lb.sz` label size
    * @param `txt.lb.size` text label size
    * @param `pt.size` point size
    * @param `random_state` used for generating repeatable random colours when colours are not available
    * @param `use.cols` if null, then generate colours randomly
    * @param `use.labels` vector of labels
    * @param `limits` plot limits
    * @param `annotate.plot` boolean to indicate if plot should be annotated by with indices
    * @param `index.map` either a NA type or a numeric vector to set the indices of each label
 