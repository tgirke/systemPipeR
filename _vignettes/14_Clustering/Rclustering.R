## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(ggplot2)
    library(ape)
})

## ----hclust_heatmap_example, eval=TRUE, warning=FALSE, message=FALSE-----
library(gplots) 
y <- matrix(rnorm(500), 100, 5, dimnames=list(paste("g", 1:100, sep=""), paste("t", 1:5, sep=""))) 
heatmap.2(y) # Shortcut to final result

## ----hclust_heatmap_example_setpwise, eval=TRUE, warning=FALSE, message=FALSE----
## Row- and column-wise clustering 
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") 
## Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
## Plot heatmap 
mycol <- colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) 

## ----kmeans_example, eval=TRUE, warning=FALSE, message=FALSE-------------
km <- kmeans(t(scale(t(y))), 3)
km$cluster 

## ----fanny_example, eval=TRUE, warning=FALSE, message=FALSE--------------
library(cluster) # Loads the cluster library.
fannyy <- fanny(y, k=4, metric = "euclidean", memb.exp = 1.2)
round(fannyy$membership, 2)[1:4,]
fannyy$clustering 

## ----pca_example, eval=TRUE, warning=FALSE, message=FALSE----------------
pca <- prcomp(y, scale=T)
summary(pca) # Prints variance summary for all principal components
plot(pca$x, pch=20, col="blue", type="n") # To plot dots, drop type="n"
text(pca$x, rownames(pca$x), cex=0.8)

## ----mds_example, eval=TRUE, warning=FALSE, message=FALSE----------------
loc <- cmdscale(eurodist) 
plot(loc[,1], -loc[,2], type="n", xlab="", ylab="", main="cmdscale(eurodist)")
text(loc[,1], -loc[,2], rownames(loc), cex=0.8) 

## ----jaccard_index, eval=TRUE, warning=FALSE, message=FALSE--------------
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/clusterIndex.R") 
library(cluster); y <- matrix(rnorm(5000), 1000, 5, dimnames=list(paste("g", 1:1000, sep=""), paste("t", 1:5, sep=""))); clarax <- clara(y, 49); clV1 <- clarax$clustering; clarax <- clara(y, 50); clV2 <- clarax$clustering 
ci <- cindex(clV1=clV1, clV2=clV2, self=FALSE, minSZ=1, method="jaccard")
ci[2:3] # Returns Jaccard index and variables used to compute it 

## ----jaccard_index_clustering, eval=TRUE, warning=FALSE, message=FALSE----
clVlist <- lapply(3:12, function(x) clara(y[1:30, ], k=x)$clustering); names(clVlist) <- paste("k", "=", 3:12)
d <- sapply(names(clVlist), function(x) sapply(names(clVlist), function(y) cindex(clV1=clVlist[[y]], clV2=clVlist[[x]], method="jaccard")[[3]]))
hv <- hclust(as.dist(1-d))
plot(as.dendrogram(hv), edgePar=list(col=3, lwd=4), horiz=T, main="Similarities of 10 Clara Clustering Results for k: 3-12") 

## ----scaling, eval=TRUE, warning=FALSE, message=FALSE--------------------
## Sample data set
set.seed(1410)
y <- matrix(rnorm(50), 10, 5, dimnames=list(paste("g", 1:10, sep=""), 
            paste("t", 1:5, sep="")))
dim(y)
## Scaling
yscaled <- t(scale(t(y))) # Centers and scales y row-wise
apply(yscaled, 1, sd)

## ----euclidean_distance_matrices, eval=TRUE, warning=FALSE, message=FALSE----
dist(y[1:4,], method = "euclidean")

## ----correlation_similarity_matrices, eval=TRUE, warning=FALSE, message=FALSE----
c <- cor(t(y), method="pearson") 
as.matrix(c)[1:4,1:4]

## ----correlation_distance_matrices, eval=TRUE, warning=FALSE, message=FALSE----
d <- as.dist(1-c)
as.matrix(d)[1:4,1:4]

## ----hclust1, eval=TRUE, warning=FALSE, message=FALSE--------------------
hr <- hclust(d, method = "complete", members=NULL)
names(hr)
par(mfrow = c(1, 2)); plot(hr, hang = 0.1); plot(hr, hang = -1) 

## ----hclust_plot_tree1, eval=TRUE, warning=FALSE, message=FALSE----------
plot(as.dendrogram(hr), edgePar=list(col=3, lwd=4), horiz=T) 

## ----hclust_plot_tree2, eval=TRUE, warning=FALSE, message=FALSE----------
library(ape) 
plot.phylo(as.phylo(hr), type="p", edge.col=4, edge.width=2, 
           show.node.label=TRUE, no.margin=TRUE)

## ----hclust_object, eval=TRUE, warning=FALSE, message=FALSE--------------
hr
## Print row labels in the order they appear in the tree
hr$labels[hr$order] 

## ----cutree, eval=TRUE, warning=FALSE, message=FALSE---------------------
mycl <- cutree(hr, h=max(hr$height)/2)
mycl[hr$labels[hr$order]] 

## ----heatmap2a, eval=TRUE, warning=FALSE, message=FALSE------------------
library(gplots)
heatmap.2(y, col=redgreen(75))

## ----pheatmap, eval=TRUE, warning=FALSE, message=FALSE-------------------
library(pheatmap); library("RColorBrewer")
pheatmap(y, color=brewer.pal(9,"Blues"))

## ----heatmap2_custom, eval=TRUE, warning=FALSE, message=FALSE------------
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete")
mycol <- colorpanel(40, "darkblue", "yellow", "white")
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(mycl))

## ----kmeans2, eval=TRUE, warning=FALSE, message=FALSE--------------------
library(cluster)
pamy <- pam(d, 4)
(kmcol <- pamy$clustering)
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol,
          scale="row", density.info="none", trace="none", 
          RowSideColors=as.character(kmcol))

## ----kmeans_fuzzy, eval=TRUE, warning=FALSE, message=FALSE---------------
library(cluster)
fannyy <- fanny(d, k=4, memb.exp = 1.5)
round(fannyy$membership, 2)[1:4,]
fannyy$clustering 
## Returns multiple cluster memberships for coefficient above a certain 
## value (here >0.1)
fannyyMA <- round(fannyy$membership, 2) > 0.10 
apply(fannyyMA, 1, function(x) paste(which(x), collapse="_"))

## ----cmdscale2, eval=TRUE, warning=FALSE, message=FALSE------------------
loc <- cmdscale(eurodist) 
## Plots the MDS results in 2D plot. The minus is required in this example to 
## flip the plotting orientation.
plot(loc[,1], -loc[,2], type="n", xlab="", ylab="", main="cmdscale(eurodist)")
text(loc[,1], -loc[,2], rownames(loc), cex=0.8) 

## ----pca2, eval=TRUE, warning=FALSE, message=FALSE-----------------------
library(scatterplot3d)
pca <- prcomp(y, scale=TRUE)
names(pca)
summary(pca) # Prints variance summary for all principal components.
scatterplot3d(pca$x[,1:3], pch=20, color="blue") 

## ----sessionInfo---------------------------------------------------------
sessionInfo()

