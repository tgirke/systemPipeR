---
title: 3. Clustering Algorithms
last_updated: Wed May 24 16:24:26 2017
sidebar: mydoc_sidebar
permalink: mydoc_Rclustering_3.html
---

## Hierarchical Clustering

### Overview of algorithm 

1. Identify clusters (items) with closest distance  
2. Join them to new clusters
3. Compute distance between clusters (items)
4. Return to step 1

#### Hierarchical clustering: agglomerative Approach
<center><img title="hierarchical clustering" src="./pages/mydoc/Rclustering_files/hierarchical.png" ></center>

#### Hierarchical Clustering with Heatmap
<center><img title="heatmap" src="./pages/mydoc/Rclustering_files/heatmap.png" ></center>

- A heatmap is a color coded table. To visually identify patterns, the rows and columns of a heatmap are often sorted by hierarchical clustering trees.  
- In case of gene expression data, the row tree usually represents the genes, the column tree the treatments and the colors in the heat table represent the intensities or ratios of the underlying gene expression data set.

### Hierarchical Clustering Approaches

1. Agglomerative approach (bottom-up)
    - R functions: `hclust()` and `agnes()`

2. Divisive approach (top-down)
    - R function: `diana()`

### Tree Cutting to Obtain Discrete Clusters

1. Node height in tree
2. Number of clusters
3. Search tree nodes by distance cutoff


### Examples

#### Using `hclust` and `heatmap.2`


```r
library(gplots) 
y <- matrix(rnorm(500), 100, 5, dimnames=list(paste("g", 1:100, sep=""), paste("t", 1:5, sep=""))) 
heatmap.2(y) # Shortcut to final result
```

<img src="./pages/mydoc/Rclustering_files/hclust_heatmap_example-1.png" width="672" />

#### Stepwise Approach with Tree Cutting


```r
## Row- and column-wise clustering 
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") 
## Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
## Plot heatmap 
mycol <- colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) 
```

<img src="./pages/mydoc/Rclustering_files/hclust_heatmap_example_setpwise-1.png" width="672" />

## K-Means Clustering

### Overview of algorithm 

1. Choose the number of k clusters   
2. Randomly assign items to the k clusters
3. Calculate new centroid for each of the k clusters
4. Calculate the distance of all items to the k centroids
5. Assign items to closest centroid
6. Repeat until clusters assignments are stable

<center><img title="kmeans" src="./pages/mydoc/Rclustering_files/kmeans.png" ></center>
	
### Examples


```r
km <- kmeans(t(scale(t(y))), 3)
km$cluster 
```

```
##   g1   g2   g3   g4   g5   g6   g7   g8   g9  g10  g11  g12  g13  g14  g15  g16  g17  g18  g19  g20 
##    1    1    2    3    1    2    2    2    3    2    2    3    3    3    3    2    3    2    3    3 
##  g21  g22  g23  g24  g25  g26  g27  g28  g29  g30  g31  g32  g33  g34  g35  g36  g37  g38  g39  g40 
##    3    2    2    3    3    1    3    2    1    2    2    2    1    2    1    3    3    1    3    3 
##  g41  g42  g43  g44  g45  g46  g47  g48  g49  g50  g51  g52  g53  g54  g55  g56  g57  g58  g59  g60 
##    1    1    1    3    2    2    2    3    2    3    3    1    2    3    2    3    3    3    3    1 
##  g61  g62  g63  g64  g65  g66  g67  g68  g69  g70  g71  g72  g73  g74  g75  g76  g77  g78  g79  g80 
##    3    3    2    1    3    1    3    1    2    2    2    2    3    3    2    3    2    1    1    2 
##  g81  g82  g83  g84  g85  g86  g87  g88  g89  g90  g91  g92  g93  g94  g95  g96  g97  g98  g99 g100 
##    1    2    2    3    3    2    2    3    3    2    2    1    1    2    3    3    1    1    1    3
```

## Fuzzy C-Means Clustering

- In contrast to strict (hard) clustering approaches, fuzzy (soft) clustering methods allow multiple cluster memberships of the clustered items (Hathaway et al., 1996). 
- This is commonly achieved by assigning to each item a weight of belonging to each cluster. 
- Thus, items at the edge of a cluster, may be in a cluster to a lesser degree than items at the center of a cluster. 
- Typically, each item has as many coefficients (weights) as there are clusters that sum up for each item to one.

### Examples

#### Fuzzy Clustering with `fanny`


```r
library(cluster) # Loads the cluster library.
fannyy <- fanny(y, k=4, metric = "euclidean", memb.exp = 1.2)
round(fannyy$membership, 2)[1:4,]
```

```
##    [,1] [,2] [,3] [,4]
## g1 0.85 0.06 0.01 0.08
## g2 0.10 0.69 0.16 0.05
## g3 0.02 0.01 0.95 0.02
## g4 0.01 0.02 0.01 0.96
```

```r
fannyy$clustering 
```

```
##   g1   g2   g3   g4   g5   g6   g7   g8   g9  g10  g11  g12  g13  g14  g15  g16  g17  g18  g19  g20 
##    1    2    3    4    2    3    1    1    2    3    3    2    2    4    4    1    2    3    4    4 
##  g21  g22  g23  g24  g25  g26  g27  g28  g29  g30  g31  g32  g33  g34  g35  g36  g37  g38  g39  g40 
##    4    3    3    2    4    4    4    3    1    3    3    3    1    3    4    2    4    2    2    2 
##  g41  g42  g43  g44  g45  g46  g47  g48  g49  g50  g51  g52  g53  g54  g55  g56  g57  g58  g59  g60 
##    1    2    3    2    1    3    3    2    1    2    1    2    1    2    3    4    2    2    2    1 
##  g61  g62  g63  g64  g65  g66  g67  g68  g69  g70  g71  g72  g73  g74  g75  g76  g77  g78  g79  g80 
##    4    4    1    2    2    4    2    4    3    3    4    3    2    2    1    2    2    2    1    1 
##  g81  g82  g83  g84  g85  g86  g87  g88  g89  g90  g91  g92  g93  g94  g95  g96  g97  g98  g99 g100 
##    3    3    3    2    4    3    3    4    4    3    3    3    2    3    4    1    1    1    2    4
```
	
## Principal Component Analysis (PCA)

Principal components analysis (PCA) is a data reduction technique that allows to simplify multidimensional data sets to 2 or 3 dimensions for plotting purposes and visual variance analysis.

### Basic Steps

- Center (and standardize) data
- First principal component axis
    - Across centroid of data cloud
	- Distance of each point to that line is minimized, so that it crosses the maximum variation of the data cloud
- Second principal component axis 
    - Orthogonal to first principal component
	- Along maximum variation in the data
- First PCA axis becomes x-axis and second PCA axis y-axis 
- Continue process until the necessary number of principal components is obtained 

<center><img title="pca" src="./pages/mydoc/Rclustering_files/pca.png" ></center>

### Example


```r
pca <- prcomp(y, scale=T)
summary(pca) # Prints variance summary for all principal components
```

```
## Importance of components%s:
##                           PC1    PC2    PC3    PC4    PC5
## Standard deviation     1.1165 1.0418 1.0237 0.9311 0.8679
## Proportion of Variance 0.2493 0.2170 0.2096 0.1734 0.1506
## Cumulative Proportion  0.2493 0.4664 0.6760 0.8494 1.0000
```

```r
plot(pca$x, pch=20, col="blue", type="n") # To plot dots, drop type="n"
text(pca$x, rownames(pca$x), cex=0.8)
```

<img src="./pages/mydoc/Rclustering_files/pca_example-1.png" width="672" />
1st and 2nd principal components explain x% of variance in data.

## Multidimensional Scaling (MDS)

- Alternative dimensionality reduction approach
- Represents distances in 2D or 3D space
- Starts from distance matrix (PCA uses data points)

### Example

The following example performs MDS analysis with `cmdscale` on the geographic distances among European cities.


```r
loc <- cmdscale(eurodist) 
plot(loc[,1], -loc[,2], type="n", xlab="", ylab="", main="cmdscale(eurodist)")
text(loc[,1], -loc[,2], rownames(loc), cex=0.8) 
```

<img src="./pages/mydoc/Rclustering_files/mds_example-1.png" width="672" />

## Biclustering

Finds in matrix subgroups of rows and columns which are as similar as possible to each other and as different as possible to the remaining data points.

<center><img title="pca" src="./pages/mydoc/Rclustering_files/biclust.png"  ></center>
<center> Unclustered --------------------------> Clustered</center>

## Similarity Measures for Clusters

- Compare the numbers of identical and unique item pairs appearing in cluster sets
- Achieved by counting the number of item pairs found in both clustering sets _(a)_ as well as the pairs appearing only in the first _(b)_ or the second _(c)_ set. 
- With this a similarity coefficient, such as the Jaccard index, can be computed. The latter is defined as the size of the intersect divided by the size of the union of two sample sets: _a/(a+b+c)_. 
- In case of partitioning results, the Jaccard Index measures how frequently pairs of items are joined together in two clustering data sets and how often pairs are observed only in one set. 
- Related coefficient are the Rand Index and the Adjusted Rand Index. These indices also consider the number of pairs _(d)_ that are not joined together in any of the clusters in both sets. 

### Example: 

#### Jaccard index for cluster sets

The following imports the `cindex()` function and computes the Jaccard Index for two sample clusters.


```r
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/clusterIndex.R") 
library(cluster); y <- matrix(rnorm(5000), 1000, 5, dimnames=list(paste("g", 1:1000, sep=""), paste("t", 1:5, sep=""))); clarax <- clara(y, 49); clV1 <- clarax$clustering; clarax <- clara(y, 50); clV2 <- clarax$clustering 
ci <- cindex(clV1=clV1, clV2=clV2, self=FALSE, minSZ=1, method="jaccard")
ci[2:3] # Returns Jaccard index and variables used to compute it 
```

```
## $variables
##    a    b    c 
## 9245 4238 3951 
## 
## $Jaccard_Index
## [1] 0.5302856
```

#### Clustering cluster sets with Jaccard index

The following example shows how one can cluster entire cluster result sets. First, 10 sample cluster results are created with Clara using k-values from 3 to 12. The results are stored as named clustering vectors in a list object. Then a nested sapply loop is used to generate a similarity matrix of Jaccard Indices for the clustering results. After converting the result into a distance matrix, hierarchical clustering is performed with `hclust`.}


```r
clVlist <- lapply(3:12, function(x) clara(y[1:30, ], k=x)$clustering); names(clVlist) <- paste("k", "=", 3:12)
d <- sapply(names(clVlist), function(x) sapply(names(clVlist), function(y) cindex(clV1=clVlist[[y]], clV2=clVlist[[x]], method="jaccard")[[3]]))
hv <- hclust(as.dist(1-d))
plot(as.dendrogram(hv), edgePar=list(col=3, lwd=4), horiz=T, main="Similarities of 10 Clara Clustering Results for k: 3-12") 
```

<img src="./pages/mydoc/Rclustering_files/jaccard_index_clustering-1.png" width="672" />

- Remember: there are many additional clustering algorithms.
- Additional details can be found in the Clustering Section of the [R/Bioconductor Manual](http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual\#TOC-Clustering-and-Data-Mining-in-R).

<br><br><center><a href="mydoc_Rclustering_2.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_Rclustering_4.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
