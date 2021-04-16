################
## exploreDDS ##
################
exploreDDS <- function(countMatrix, targets, cmp=cmp[[1]], preFilter=NULL, transformationMethod="raw", blind=TRUE) {
  ## A few validations ##
  if (!transformationMethod %in% c("raw", "rlog", "vst"))
    stop("Supported methods include 'raw', 'rlog' and 'vst'")
  if(is.data.frame(countMatrix)){
    countMatrix <- as.matrix(countMatrix)
  } else if(is.matrix(countMatrix)){
    countMatrix <- countMatrix
  } else {
    stop("countMatrix needs to be assigned an object of class 'data.frame' OR 'matrix'")
  }
  if (!is.data.frame(targets)) stop("targets needs to be assignes an object of class 'data.frame'")
  if (all(class(cmp) != "matrix" & length(cmp)==2)) cmp <- t(as.matrix(cmp))
  ## Samples
  samples <- as.character(targets$Factor); names(samples) <- paste(as.character(targets$SampleName), "", sep="")
  ## Create full DESeqDataSet object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=countMatrix,
                                        colData=data.frame(condition=samples), design = ~ condition)
  ## Pre-filtering
  if(!is.null(preFilter)){
    if (!is.numeric(preFilter)) stop ("'preFilter' needs to be numeric value.")
    keep <- rowSums(DESeq2::counts(dds)) >= preFilter
    dds <- dds[keep,]
  }
  ## Estimate of (i) size factors, (ii) dispersion, (iii) negative binomial GLM fitting and (iv) Wald statistics
  dds_deseq2 <- DESeq2::DESeq(dds, quiet=TRUE)
  ## Count data transformations
  if (transformationMethod == "rlog") {
    normdata <- DESeq2::rlog(dds_deseq2, blind=TRUE)
  } else if (transformationMethod == "vst") {
    normdata <- DESeq2::varianceStabilizingTransformation(dds_deseq2, blind = TRUE)
  } else if (transformationMethod == "raw") {
    normdata <- dds
  }
  return(normdata)
}

####################
## exploreDDSplot ##
####################
exploreDDSplot <- function(countMatrix, targets, cmp=cmp[[1]], preFilter=NULL, samples, blind=TRUE,
                           scattermatrix = FALSE, plotly = FALSE, savePlot=FALSE, filePlot=NULL) {
  ## Validations
  SampleName <- targets$SampleName; names(SampleName) <- targets$SampleName
  if(is.numeric(samples)){
    samples <- SampleName[samples]
    if(!all(samples %in% SampleName)) stop(paste("samples position can be assigned from the following options", paste0(1:length(SampleName), collapse=", "), sep = " "))
  } else if(is.character(samples)){
    if(all(samples=="ALL")){
      samples <- SampleName
      if(!scattermatrix=="TRUE") stop("'scattermatrix' argument needs to set as TRUE in the case of ALL the samples selected.")
    } else {
      samples <- SampleName[samples]
      if(!all(samples %in% SampleName)) stop(paste("samples names can be assigned from the following options", paste0((SampleName), collapse=", "), sep = " "))
    }
  }
  transformation <- NULL
  ## Calculate the data transformations
  suppressWarnings({
    vst <- exploreDDS(countMatrix, targets, cmp=cmp, preFilter=preFilter, transformationMethod="vst", blind=blind)
    rlog <- exploreDDS(countMatrix, targets, cmp=cmp, preFilter=preFilter, transformationMethod="rlog", blind=blind)
    dss <- exploreDDS(countMatrix, targets, cmp=cmp, preFilter=preFilter, transformationMethod="raw")
    dss <- DESeq2::estimateSizeFactors(dss)})
  ## create dataframe with transformed values
  transform_df <- dplyr::bind_rows(
    dplyr::as_tibble(log2(DESeq2::counts(dss, normalized = TRUE)[, samples] + 1)) %>%
      dplyr::mutate(transformation = "log2(x + 1)"),
    dplyr::as_tibble(SummarizedExperiment::assay(vst)[, samples]) %>% dplyr::mutate(transformation = "vst"),
    dplyr::as_tibble(SummarizedExperiment::assay(rlog)[, samples]) %>% dplyr::mutate(transformation = "rlog"))
  names <- colnames(transform_df)[1:2]
  lvl <- levels(factor(transform_df$transformation))
  ## plot
  if (scattermatrix==TRUE){
    plot <- GGally::ggpairs(transform_df, title="Scatterplot of transformed counts", ggplot2::aes_string(colour="transformation"))
  } else {
    plot <- ggplot2::ggplot(transform_df, ggplot2::aes(x = .data[[names[1]]], y = .data[[names[2]]])) +
      ggplot2::geom_hex(bins = 80) +
      ggplot2::coord_fixed() + ggplot2::facet_grid( . ~transformation) +
      ggplot2::xlab(names[1]) + ggplot2::ylab(names[2])
  }
  if (savePlot == TRUE) {
    ggplot2::ggsave(filePlot, scale = 0.8)
  }
  ## Return
  if (plotly == TRUE) {
    plot <- transform_df %>%
      dplyr::group_by(transformation) %>%
      dplyr::do(p=plotly::plot_ly(., x = .data[[names[1]]], y = .data[[names[2]]], color = ~transformation, type = "scatter",
                                  name= ~transformation, showlegend=TRUE, legendgroup = ~transformation)) %>%
      plotly::subplot(nrows = 1, shareX = TRUE, shareY = TRUE)
  }
  return(plot)
}

################
## hclustplot ##
################
hclustplot <- function(exploredds, method = "spearman", plotly = FALSE, savePlot = FALSE, filePlot = NULL) {
  ## Validations
  if (!class(exploredds) == "DESeqTransform") stop("'exploredds' needs to be assignes an object of class 'DESeqTransform'. For more information check 'help(exploreDDS)'.")
  ## cor() computes the correlation coefficient
  d <- stats::cor(SummarizedExperiment::assay(exploredds), method = method)
  ## Hierarchical cluster analysis
  hc <- stats::hclust(stats::dist(1 - d))
  ## plot phylogenetic trees
  plot <- ggtree::ggtree(ape::as.phylo(hc), color="blue") + ggtree::geom_tiplab() +
    ggplot2::coord_cartesian(clip = 'off') + ggtree::theme_tree(plot.margin=ggplot2::margin(6, 60, 6, 6))
  if (savePlot == TRUE){
    ggplot2::ggsave(filePlot, scale = 0.8)
  }
  ##Return
  if (plotly == TRUE){
    return(plotly::ggplotly(plot)) }
  return(plot)
}

################
## heatMaplot ##
################
heatMaplot <- function(exploredds, clust, DEGlist = NULL, plotly = FALSE, savePlot = FALSE, filePlot = NULL, ...) {
  ## Validations
  if (!class(exploredds) == "DESeqTransform") stop("'exploredds' needs to be assignes an object of class 'DESeqTransform'. For more information check 'help(exploreDDS)'")
  anno <- as.data.frame(exploredds$condition); colnames(anno) <- "Condition"
  ## sample-to-sample distances
  if (clust == "samples") {
    sampleDists <- stats::dist(t(SummarizedExperiment::assay(exploredds)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(anno) <- colnames(sampleDistMatrix)
    if (plotly == FALSE) {
      pheatPlot <- pheatmap::pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
                                      clustering_distance_cols = sampleDists, annotation_col = anno)
    } else if (plotly == TRUE) {
      plot <- plotly::plot_ly(x = colnames(sampleDistMatrix), y = rownames(sampleDistMatrix),
                              z = sampleDistMatrix, type = "heatmap")
    }
  } else if (clust == "ind") {
    ## Hierarchical clustering on the transformed expression matrix subsetted by the DEGs identified in differential expression analysis.
    if (any(is.null(DEGlist) | !is.character(DEGlist))) stop("Provide a character vector with the gene names identified in differential expression analysis.")
    dist <- SummarizedExperiment::assay(exploredds)[DEGlist, ]
    rownames(anno) <- colnames(dist)
    if (plotly == FALSE) {
      pheatPlot <- pheatmap::pheatmap(dist, scale = "row", clustering_distance_rows = "correlation",
                                      clustering_distance_cols = "correlation", annotation_col = anno)
    } else if (plotly == TRUE) {
      plot <- plotly::plot_ly(x = colnames(dist), y = rownames(dist), z = dist, type = "heatmap")
    }
  } else {stop("Supported clust include 'samples' and 'ind'") }
  if (savePlot == TRUE) {
    ggplot2::ggsave(plot = pheatPlot, filename = filePlot)
  }
  ##Return
  if (plotly == TRUE) {
    return(plot) }
  return(pheatPlot)
}

#############
## PCAplot ##
#############
PCAplot <- function(exploredds, plotly = FALSE, savePlot = FALSE, filePlot = NULL) {
  ## Validations
  if (!class(exploredds) == "DESeqTransform") {
    warning("'exploredds' needs to be assignes an object of class 'DESeqTransform'.
    Here we are converting the object into a 'DESeqTransform'class for
    downstream analysis. For more information check 'help(exploreDDS)'")
    exploredds <- DESeq2::DESeqTransform(exploredds) }
  ## Plot
  pcaData <- DESeq2::plotPCA(exploredds, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  Sample <- exploredds$condition
  plot <- ggplot2::ggplot(pcaData, ggplot2::aes_string("PC1", "PC2", color = Sample)) +
    ggplot2::geom_point(size=3) +
    ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggplot2::coord_fixed() + ggplot2::ggtitle("Principal Component Analysis (PCA)")
  ## Save plot
  if (savePlot == TRUE){
    ggplot2::ggsave(plot = plot, filename = filePlot)
  }
  ## Return
  if (plotly == TRUE){
    return(plotly::ggplotly(plot)) }
  return(plot)
}

#############
## GLMplot ##
#############
GLMplot <- function(exploredds, plotly = FALSE, savePlot = FALSE, filePlot = NULL, ...) {
  ## Add validation, need to be counts reads
  if (is.data.frame(exploredds)) {
    count_mat <- exploredds
  } else if (class(exploredds) == "DESeqDataSet") {
    count_mat <- DESeq2::counts(exploredds)
  } else if (!class(exploredds) == "DESeqDataSet") {
    stop("'exploredds' needs to be assignes an object of class 'DESeqDataSet'.
                For more information check 'help(exploreDDS)', and select the transformationMethod='raw'")
  }
  ##glmpca is performed on raw counts
  nozero <- count_mat[which(rowSums(count_mat) > 0),]
  gpca <- glmpca::glmpca(nozero, L=2, ...)
  gpca.dat <- gpca$factors
  gpca.dat$condition <- exploredds$condition
  Samples <- as.character(exploredds$condition)
  plot <- ggplot2::ggplot(gpca.dat, ggplot2::aes_string("dim1", "dim2")) +
    ggplot2::geom_point(size = 3, ggplot2::aes(color=Samples)) + ggplot2::coord_fixed() +
    ggplot2::ggtitle("Generalized PCA (GLM-PCA)")
  ## Save plot
  if (savePlot == TRUE){
    ggplot2::ggsave(plot = plot, filename = filePlot)
  }
  ## Return
  if (plotly == TRUE){
    return(plotly::ggplotly(plot)) }
  return(plot)
}

#############
## MDSplot ##
#############
MDSplot <- function(exploredds, method="spearman", plotly = FALSE, savePlot = FALSE, filePlot = NULL) {
  ## Add validation
  if (!class(exploredds) == "DESeqTransform") {
    warning("'exploredds' needs to be assignes an object of class 'DESeqTransform'.
    Here we are converting the object into a 'DESeqTransform'class for
    downstream analysis. For more information check 'help(exploreDDS)'")
    exploredds <- DESeq2::DESeqTransform(exploredds) }
  ## transformation to a distance matrix
  d <- stats::cor(SummarizedExperiment::assay(exploredds), method = method)
  distmat <- stats::dist(1 - d)
  ## perform MDS
  mdsData <- data.frame(stats::cmdscale(distmat))
  mds <- cbind(mdsData, as.data.frame(SummarizedExperiment::colData(exploredds)))
  Sample <- exploredds$condition
  ## plot
  plot <- ggplot2::ggplot(mds, ggplot2::aes_string("X1", "X2", color=Sample)) +
    ggplot2::geom_point(size=3) + ggplot2::scale_y_reverse() +
    ggplot2::ggtitle("Multidimensional Scaling (MDS)")
  ## Save plot
  if (savePlot == TRUE){
    ggplot2::ggsave(plot = plot, filename = filePlot)
  }
  ## Return
  if (plotly == TRUE){
    return(plotly::ggplotly(plot)) }
  return(plot)
}

###############
## tSNEplot ##
###############
tSNEplot <- function(countMatrix, targets, plotly = FALSE, savePlot = FALSE, filePlot = NULL, ...) {
  ## Validations
  if(is.data.frame(countMatrix)) {
    countMatrix <- as.matrix(countMatrix)
  } else if(is.matrix(countMatrix)) {
    countMatrix <- countMatrix
  } else {
    stop("countMatrix needs to be assigned an object of class 'data.frame' OR 'matrix'")
  }
  if (!is.data.frame(targets)) stop("targets needs to be assignes an object of class 'data.frame'")
  ## data manipulation
  countDF_uni <- t(unique(countMatrix)) # removes duplicates and transpose matrix, samples perspective
  set.seed(42)
  tsne_out <- Rtsne::Rtsne(countDF_uni, dims = 2, theta = 0.0, ...)
  targets <- data.frame(targets)
  Sample <- targets$Factor
  plotdata <- data.frame(tsne_x = tsne_out$Y[,1], tsne_y = tsne_out$Y[,2])
  ## Plot
  plot <- ggplot2::ggplot(plotdata, ggplot2::aes_string(x = "tsne_x", y = "tsne_y")) +
    ggplot2::geom_point(size = 3, ggplot2::aes(color = Sample)) + ggplot2::ggtitle("t-SNE")
  ## Save plot
  if (savePlot == TRUE) {
    ggplot2::ggsave(plot = plot, filename = filePlot)
  }
  ## Return
  if (plotly == TRUE) {
    return(plotly::ggplotly(plot)) }
  return(plot)
}

#############
## MAplot ##
#############
MAplot <- function(exploredds, lfcShrink= FALSE, padj.cutoff = 0.05, plotly = FALSE, savePlot = FALSE, filePlot = NULL) {
  ## Add validation, need to be raw, counts
  padj <- NULL
  if (!class(exploredds) == "DESeqDataSet") {
    stop("'exploredds' needs to be assignes an object of class 'DESeqDataSet'.
                For more information check 'help(exploreDDS)'")
  }
  ## lfcShrink
  if (lfcShrink == FALSE) {
    res <- DESeq2::results(DESeq2::DESeq(exploredds))
  } else if (lfcShrink == TRUE) {
    resLFC <- DESeq2::lfcShrink(DESeq2::DESeq(exploredds))
  }
  results <- as.data.frame(res)
  if (any(is.na(results$padj))) {
    print("removing NA from the results")
    results[is.na(results)] = 0.99
  }
  ## plot
  plot <- ggplot2::ggplot(results, ggplot2::aes_string(x = "baseMean", y = "log2FoldChange")) +
    ggplot2::geom_point(ggplot2::aes(colour = padj < padj.cutoff), size = 0.5) +
    ggplot2::scale_colour_manual(name = paste0('padj < ', padj.cutoff),
                                 values = stats::setNames(c('red','grey'), c(TRUE, FALSE))) +
    ggplot2::scale_x_continuous(trans = "log10", limits = c(0.1,300000))
  #ggplot2::geom_smooth(colour = "red")
  ## Save plot
  if (savePlot == TRUE) {
    ggplot2::ggsave(plot = plot, filename = filePlot)
  }
  ## Return
  if (plotly == TRUE) {
    return(suppressWarnings(suppressMessages(plotly::ggplotly(plot)))) }
  return(suppressWarnings(suppressMessages(print(plot))))
}

##################
## volcanoplot ##
##################
volcanoplot <- function(degseqDF, comparison = "M12-A12", filter = c(Fold = 2, FDR = 10),
                        genes="NULL", plotly = FALSE, savePlot = FALSE, filePlot = NULL) {
  ## Validations
  ##TODO
  ## Selecting comparison
  table <- degseqDF %>%
    dplyr::select(paste0(comparison,"_FDR"), paste0(comparison,"_logFC")) %>%
    tibble::rownames_to_column(var = "names") %>%
    dplyr::mutate(significant = .data[[paste0(comparison,"_FDR")]] <= filter["FDR"]/100 &
                    abs(as.numeric(.data[[paste0(comparison,"_logFC")]])) > filter["Fold"])
  table$significant[is.na(table$significant)] <- FALSE
  if(!is.null(genes)){
    genes <- table %>%
      dplyr::filter(names %in% genes)
  }
  ## plot
  plot <- ggplot2::ggplot(table, ggplot2::aes(x=.data[[paste0(comparison,"_logFC")]],
                                              y=-log10(as.numeric(.data[[paste0(comparison,"_FDR")]])), label=names)) +
    ggplot2::geom_point(ggplot2::aes_string(color = "significant")) +
    ggplot2::geom_vline(xintercept = c(-filter["Fold"], filter["Fold"]), linetype=2) +
    #  ggplot2::geom_hline(yintercept = -log10(filter["FDR"]/100), linetype = 2) +
    ggplot2::ggtitle(comparison) +
    ggplot2::xlab("log2 fold change") +
    ggplot2::ylab("-log10(p-value)") +
    ggrepel::geom_text_repel(data=genes) +
    #scale_y_continuous(limits = c(0,50)) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0.5),
                   axis.title = ggplot2::element_text(size = ggplot2::rel(1.25))) +
    ggplot2::scale_color_manual(values=c("#524b4b", "#e11f28")) +
    ggplot2::theme_bw()
  ## Save plot
  if (savePlot == TRUE) {
    ggplot2::ggsave(plot = plot, filename = filePlot)
  }
  ## Return
  if (plotly == TRUE) {
    return(suppressWarnings(suppressMessages(plotly::ggplotly(plot)))) }
  return(suppressWarnings(suppressMessages(print(plot))))
}
