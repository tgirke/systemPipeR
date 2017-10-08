## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----scatter_plot_basic, eval=TRUE---------------------------------------
set.seed(1410)
y <- matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))
plot(y[,1], y[,2]) 

## ----scatter_plot_allpairs, eval=TRUE------------------------------------
pairs(y) 

## ----scatter_plot_labels, eval=TRUE--------------------------------------
plot(y[,1], y[,2], pch=20, col="red", main="Symbols and Labels")
text(y[,1]+0.03, y[,2], rownames(y))

## ----scatter_plot_row_names, eval=FALSE----------------------------------
## plot(y[,1], y[,2], type="n", main="Plot of Labels")
## text(y[,1], y[,2], rownames(y))

## ----plot_parameters, eval=FALSE-----------------------------------------
## grid(5, 5, lwd = 2)
## op <- par(mar=c(8,8,8,8), bg="lightblue")
## plot(y[,1], y[,2], type="p", col="red", cex.lab=1.2, cex.axis=1.2,
##      cex.main=1.2, cex.sub=1, lwd=4, pch=20, xlab="x label",
##      ylab="y label", main="My Main", sub="My Sub")
## par(op)

## ----scatter_plot_regress, eval=TRUE-------------------------------------
plot(y[,1], y[,2])
myline <- lm(y[,2]~y[,1]); abline(myline, lwd=2) 
summary(myline) 

## ----scatter_plot_log_scale, eval=TRUE-----------------------------------
plot(y[,1], y[,2], log="xy") 

## ----scatter_plot_math, eval=TRUE----------------------------------------
plot(y[,1], y[,2]); text(y[1,1], y[1,2], 
     expression(sum(frac(1,sqrt(x^2*pi)))), cex=1.3) 

## ----iris_structure, eval=TRUE-------------------------------------------
class(iris)
iris[1:4,]
table(iris$Species)

## ----exercise1_solution, eval=FALSE, echo=FALSE--------------------------
## plot(iris[,1], iris[,2], col=iris$Species, lwd=2, pch=19)
## plot(iris[,1], iris[,2], col=iris$Species, lwd=2, pch=19, xlim=c(4,16), ylim=c(2,8))

## ----line_plot_single, eval=TRUE-----------------------------------------
plot(y[,1], type="l", lwd=2, col="blue") 

## ----line_plot_many, eval=TRUE-------------------------------------------
split.screen(c(1,1)) 
plot(y[,1], ylim=c(0,1), xlab="Measurement", ylab="Intensity", type="l", lwd=2, col=1)
for(i in 2:length(y[1,])) { 
	screen(1, new=FALSE)
	plot(y[,i], ylim=c(0,1), type="l", lwd=2, col=i, xaxt="n", yaxt="n", ylab="", 
             xlab="", main="", bty="n") 
}
close.screen(all=TRUE) 

## ----bar_plot_basic, eval=TRUE-------------------------------------------
barplot(y[1:4,], ylim=c(0, max(y[1:4,])+0.3), beside=TRUE, 
        legend=letters[1:4]) 
text(labels=round(as.vector(as.matrix(y[1:4,])),2), x=seq(1.5, 13, by=1)
     +sort(rep(c(0,1,2), 4)), y=as.vector(as.matrix(y[1:4,]))+0.04) 

## ----bar_plot_error_bar, eval=TRUE---------------------------------------
bar <- barplot(m <- rowMeans(y) * 10, ylim=c(0, 10))
stdev <- sd(t(y))
arrows(bar, m, bar, m + stdev, length=0.15, angle = 90)

## ----bar_plot_mirrored, eval=TRUE----------------------------------------
df <- data.frame(group = rep(c("Above", "Below"), each=10), x = rep(1:10, 2), y = c(runif(10, 0, 1), runif(10, -1, 0)))
plot(c(0,12),range(df$y),type = "n")
barplot(height = df$y[df$group == "Above"], add = TRUE,axes = FALSE)
barplot(height = df$y[df$group == "Below"], add = TRUE,axes = FALSE)

## ----hist_plot, eval=TRUE------------------------------------------------
hist(y, freq=TRUE, breaks=10)

## ----density_plot, eval=TRUE---------------------------------------------
plot(density(y), col="red")

## ----pie_chart, eval=TRUE------------------------------------------------
pie(y[,1], col=rainbow(length(y[,1]), start=0.1, end=0.8), clockwise=TRUE)
legend("topright", legend=row.names(y), cex=1.3, bty="n", pch=15, pt.cex=1.8, 
col=rainbow(length(y[,1]), start=0.1, end=0.8), ncol=1) 

## ----color_palette, eval=TRUE--------------------------------------------
palette()
palette(rainbow(5, start=0.1, end=0.2))
palette()
palette("default")

## ----gray_palette, eval=TRUE---------------------------------------------
gray(seq(0.1, 1, by= 0.2))

## ----color_panel, eval=FALSE---------------------------------------------
## library(gplots)
## colorpanel(5, "darkblue", "yellow", "white")

## ----par_mfrow, eval=TRUE------------------------------------------------
par(mfrow=c(2,3)); for(i in 1:6) { plot(1:10) } 

## ----layout_plot, eval=TRUE----------------------------------------------
nf <- layout(matrix(c(1,2,3,3), 2, 2, byrow=TRUE), c(3,7), c(5,5), 
             respect=TRUE)
# layout.show(nf)
for(i in 1:3) { barplot(1:10) } 

## ----save_pdf, eval=FALSE------------------------------------------------
## pdf("test.pdf"); plot(1:10, 1:10); dev.off()

## ----save_svg, eval=FALSE------------------------------------------------
## svg("test.svg"); plot(1:10, 1:10); dev.off()

## ----structure_iris2, eval=TRUE------------------------------------------
class(iris)
iris[1:4,]
table(iris$Species)

## ----exercise2_solution, eval=FALSE, echo=FALSE--------------------------
## mMA <- sapply(colnames(iris[,1:4]), function(x) tapply(iris[,x], iris[,5], mean))
## barplot(mMA, beside=T, legend=rownames(mMA))

## ----help_lattice, eval=FALSE--------------------------------------------
## library(help=lattice)

## ----param_lattice, eval=FALSE-------------------------------------------
## ?lattice.options
## ?trellis.device

## ----scatter_plot_lattice, eval=TRUE-------------------------------------
library(lattice)
p1 <- xyplot(1:8 ~ 1:8 | rep(LETTERS[1:4], each=2), as.table=TRUE) 
plot(p1)

## ----line_plot_lattice, eval=TRUE----------------------------------------
library(lattice)
p2 <- parallelplot(~iris[1:4] | Species, iris, horizontal.axis = FALSE, 
              layout = c(1, 3, 1))  
plot(p2)

## ----qplot_scatter_plot_sample, eval=TRUE--------------------------------
library(ggplot2)
x <- sample(1:10, 10); y <- sample(1:10, 10); cat <- rep(c("A", "B"), 5)

## ----qplot_scatter_plot, eval=TRUE---------------------------------------
qplot(x, y, geom="point")

## ----qplot_scatter_plot_dot_param, eval=TRUE-----------------------------
qplot(x, y, geom="point", size=x, color=cat, 
      main="Dot Size and Color Relative to Some Values")

## ----qplot_scatter_plot_no_legend, eval=TRUE-----------------------------
qplot(x, y, geom="point", size=x, color=cat) + 
      theme(legend.position = "none")

## ----qplot_scatter_plot_shapes, eval=TRUE--------------------------------
qplot(x, y, geom="point", size=5, shape=cat)

## ----qplot_scatter_plot_colored_groups, eval=TRUE------------------------
p <- qplot(x, y, geom="point", size=x, color=cat, 
            main="Dot Size and Color Relative to Some Values") + 
     theme(legend.position = "none")
print(p)

## ----qplot_scatter_regression_line, eval=TRUE----------------------------
set.seed(1410)
dsmall <- diamonds[sample(nrow(diamonds), 1000), ]
p <- qplot(carat, price, data = dsmall) +
           geom_smooth(method="lm")
print(p)

## ----qplot_scatter_regression_loess, eval=TRUE---------------------------
p <- qplot(carat, price, data=dsmall, geom=c("point", "smooth")) 
print(p) # Setting se=FALSE removes error shade

## ----ggplot_background_color, eval=FALSE---------------------------------
## ... + theme(panel.background=element_rect(fill = "white", colour = "black"))

## ----ggplot_store_plot, eval=FALSE---------------------------------------
## p <- ggplot(dsmall, aes(carat, price)) + geom_point()
## p # or print(p)

## ----ggplot_summary, eval=FALSE------------------------------------------
## summary(p)

## ----ggplot_dot_sizes, eval=FALSE----------------------------------------
## bestfit <- geom_smooth(methodw = "lm", se = F, color = alpha("steelblue", 0.5), size = 2)
## p + bestfit # Plot with custom regression line

## ----ggplot_pass_data, eval=FALSE----------------------------------------
## p %+% diamonds[sample(nrow(diamonds), 100),]

## ----ggplot_save_plot, eval=FALSE----------------------------------------
## ggsave(p, file="myplot.pdf")

## ----ggplot_scatter_plot1, eval=TRUE-------------------------------------
p <- ggplot(dsmall, aes(carat, price, color=color)) + 
            geom_point(size=4)
print(p) 

## ----ggplot_regression_line, eval=TRUE-----------------------------------
p <- ggplot(dsmall, aes(carat, price)) + geom_point() + 
            geom_smooth(method="lm", se=FALSE) +
    	    theme(panel.background=element_rect(fill = "white", colour = "black"))
print(p) 

## ----ggplot_many_regression_lines, eval=TRUE-----------------------------
p <- ggplot(dsmall, aes(carat, price, group=color)) + 
            geom_point(aes(color=color), size=2) + 
            geom_smooth(aes(color=color), method = "lm", se=FALSE) 
print(p) 

## ----ggplot_loess_regression, eval=TRUE----------------------------------
p <- ggplot(dsmall, aes(carat, price)) + geom_point() + geom_smooth() 
print(p) # Setting se=FALSE removes error shade

## ----ggplot_line_plot, eval=TRUE-----------------------------------------
p <- ggplot(iris, aes(Petal.Length, Petal.Width, group=Species, 
            color=Species)) + geom_line() 
print(p) 

## ----ggplot_line_plot_faceting, eval=TRUE--------------------------------
p <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) + 
    	    geom_line(aes(color=Species), size=1) + 
            facet_wrap(~Species, ncol=1)
print(p) 

## ----structure_iris_data3, eval=TRUE-------------------------------------
class(iris)
iris[1:4,]
table(iris$Species)

## ----exercise3_solution, eval=FALSE, echo=FALSE--------------------------
## ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_point(aes(color=Species), size=3) + xlim(4,15) + ylim(2,12)
## ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_line(aes(color=Species), size=3) + ggtitle("Line Plot") + facet_wrap(~Species, ncol=1)

## ----iris_mean, eval=TRUE------------------------------------------------
iris_mean <- aggregate(iris[,1:4], by=list(Species=iris$Species), FUN=mean) 

## ----iris_sd, eval=TRUE--------------------------------------------------
iris_sd <- aggregate(iris[,1:4], by=list(Species=iris$Species), FUN=sd) 

## ----iris_mean_melt, eval=TRUE-------------------------------------------
library(reshape2) # Defines melt function
df_mean <- melt(iris_mean, id.vars=c("Species"), variable.name = "Samples", value.name="Values")

## ----iris_sd_melt, eval=TRUE---------------------------------------------
df_sd <- melt(iris_sd, id.vars=c("Species"), variable.name = "Samples", value.name="Values")

## ----iris_limits, eval=TRUE----------------------------------------------
limits <- aes(ymax = df_mean[,"Values"] + df_sd[,"Values"], ymin=df_mean[,"Values"] - df_sd[,"Values"])

## ----iris_mean_bar_plot, eval=TRUE---------------------------------------
p <- ggplot(df_mean, aes(Samples, Values, fill = Species)) + 
	    geom_bar(position="dodge", stat="identity")
print(p) 

## ----iris_factor_order, eval=FALSE---------------------------------------
## df_mean$Species <- factor(df_mean$Species, levels=unique(df_mean$Species), ordered=TRUE)

## ----iris_mean_bar_plot_sideways, eval=TRUE------------------------------
p <- ggplot(df_mean, aes(Samples, Values, fill = Species)) + 
            geom_bar(position="dodge", stat="identity") + coord_flip() + 
            theme(axis.text.y=element_text(angle=0, hjust=1))
print(p) 

## ----iris_mean_bar_plot_facetting, eval=FALSE----------------------------
## p <- ggplot(df_mean, aes(Samples, Values)) + geom_bar(aes(fill = Species), stat="identity") +
##             facet_wrap(~Species, ncol=1)
## print(p)

## ----iris_mean_bar_plot_error, eval=TRUE---------------------------------
p <- ggplot(df_mean, aes(Samples, Values, fill = Species)) + 
	    geom_bar(position="dodge", stat="identity") + geom_errorbar(limits, position="dodge") 
print(p) 

## ----ggplot2_mirrored_barplot, eval=TRUE---------------------------------
df <- data.frame(group = rep(c("Above", "Below"), each=10), x = rep(1:10, 2), y = c(runif(10, 0, 1), runif(10, -1, 0)))
p <- ggplot(df, aes(x=x, y=y, fill=group)) + 
	    geom_bar(stat="identity", position="identity")
print(p) 

## ----ggplot2_color1, eval=TRUE-------------------------------------------
library(RColorBrewer)
# display.brewer.all() 
p <- ggplot(df_mean, aes(Samples, Values, fill=Species, color=Species)) +
            geom_bar(position="dodge", stat="identity") + geom_errorbar(limits, position="dodge") + 
            scale_fill_brewer(palette="Blues") + scale_color_brewer(palette = "Greys") 
print(p) 

## ----ggplot2_color2, eval=TRUE-------------------------------------------
p <- ggplot(df_mean, aes(Samples, Values, fill=Species, color=Species)) + 
            geom_bar(position="dodge", stat="identity") + geom_errorbar(limits, position="dodge") + 
            scale_fill_manual(values=c("red", "green3", "blue")) + 
            scale_color_manual(values=c("red", "green3", "blue")) 
print(p) 

## ----iris_structure4, eval=TRUE------------------------------------------
class(iris)
iris[1:4,]
table(iris$Species)

## ----exercise4_solution, eval=FALSE, echo=FALSE--------------------------
## iris_mean <- aggregate(iris[,1:4], by=list(Species=iris$Species), FUN=mean)
## library(reshape2)
## df_mean <- melt(iris_mean, id.vars=c("Species"), variable.name = "Samples", value.name="Values")
## ggplot(df_mean, aes(Samples, Values, fill=Species)) + geom_bar(position="stack", stat="identity")
## ggplot(df_mean, aes(Samples, Values, fill=Species)) + geom_bar(position="dodge", stat="identity")

## ----ggplot_melt_data, eval=TRUE-----------------------------------------
y <- matrix(rnorm(500), 100, 5, dimnames=list(paste("g", 1:100, sep=""), paste("Sample", 1:5, sep="")))
y <- data.frame(Position=1:length(y[,1]), y)
y[1:4, ] # First rows of input format expected by melt()
df <- melt(y, id.vars=c("Position"), variable.name = "Samples", value.name="Values")
p <- ggplot(df, aes(Position, Values)) + geom_line(aes(color=Samples)) + facet_wrap(~Samples, ncol=1)
print(p)

## ----ggplot_melt_box_plot, eval=FALSE------------------------------------
## ggplot(df, aes(Samples, Values, fill=Samples)) + geom_boxplot()

## ----ggplot_jitter_plot, eval=TRUE---------------------------------------
p <- ggplot(dsmall, aes(color, price/carat)) + 
            geom_jitter(alpha = I(1 / 2), aes(color=color))
print(p) 

## ----ggplot_box_plot2, eval=TRUE-----------------------------------------
p <- ggplot(dsmall, aes(color, price/carat, fill=color)) + geom_boxplot()
print(p) 

## ----ggplot_density_plot_linecol, eval=TRUE------------------------------
p <- ggplot(dsmall, aes(carat)) + geom_density(aes(color = color))
print(p) 

## ----ggplot_density_plot_areacol, eval=TRUE------------------------------
p <- ggplot(dsmall, aes(carat)) + geom_density(aes(fill = color))
print(p) 

## ----ggplot_histogram, eval=TRUE-----------------------------------------
p <- ggplot(iris, aes(x=Sepal.Width)) + geom_histogram(aes(y = ..density.., 
            fill = ..count..), binwidth=0.2) + geom_density()  
print(p) 

## ----ggplot_pie_chart, eval=TRUE-----------------------------------------
df <- data.frame(variable=rep(c("cat", "mouse", "dog", "bird", "fly")), 
                 value=c(1,3,3,4,2)) 
p <- ggplot(df, aes(x = "", y = value, fill = variable)) + 
            geom_bar(width = 1, stat="identity") + 
            coord_polar("y", start=pi / 3) + ggtitle("Pie Chart") 
print(p) 

## ----ggplot_windrose_pie_chart, eval=TRUE--------------------------------
p <- ggplot(df, aes(x = variable, y = value, fill = variable)) + 
       geom_bar(width = 1, stat="identity") + coord_polar("y", start=pi / 3) + 
       ggtitle("Pie Chart") 
print(p) 

## ----ggplot_arrange_graphics, eval=TRUE----------------------------------
library(grid)
a <- ggplot(dsmall, aes(color, price/carat)) + geom_jitter(size=4, alpha = I(1 / 1.5), aes(color=color))
b <- ggplot(dsmall, aes(color, price/carat, color=color)) + geom_boxplot()
c <- ggplot(dsmall, aes(color, price/carat, fill=color)) + geom_boxplot() + theme(legend.position = "none")
grid.newpage() # Open a new page on grid device
pushViewport(viewport(layout = grid.layout(2, 2))) # Assign to device viewport with 2 by 2 grid layout 
print(a, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(b, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(c, vp = viewport(layout.pos.row = 2, layout.pos.col = 2, width=0.3, height=0.3, x=0.8, y=0.8))

## ----ggplot_insert_graphics, eval=TRUE-----------------------------------
library(grid)
print(a)
print(b, vp=viewport(width=0.3, height=0.3, x=0.8, y=0.8))

## ----specgraph_venn, eval=TRUE-------------------------------------------
library(systemPipeR)
setlist5 <- list(A=sample(letters, 18), B=sample(letters, 16), C=sample(letters, 20), D=sample(letters, 22), E=sample(letters, 18))
OLlist5 <- overLapper(setlist=setlist5, sep="_", type="vennsets")
vennPlot(OLlist5, mymain="", mysub="", colmode=2, ccol=c("blue", "red"))

## ----specgraph_structure, eval=TRUE--------------------------------------
library(ChemmineR)
data(sdfsample)
plot(sdfsample[1], print=FALSE)

## ----ROCR_example, eval=TRUE, warning=FALSE, message=FALSE---------------
# install.packages("ROCR") # Install if necessary on your laptop
library(ROCR)
data(ROCR.simple)
ROCR.simple
pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels)
perf <- performance( pred, "tpr", "fpr" )
plot(perf)

## ----trees_ape1, eval=TRUE-----------------------------------------------
library(ape)
tree1 <- rtree(40)
tree2 <- rtree(20)
association <- cbind(tree2$tip.label, tree2$tip.label)
cophyloplot(tree1, tree2, assoc = association,
            length.line = 4, space = 28, gap = 3)

## ----ggbio_align, eval=TRUE, warning=FALSE, message=FALSE----------------
library(ggbio)
df1 <- data.frame(time = 1:100, score = sin((1:100)/20)*10)
p1 <- qplot(data = df1, x = time, y = score, geom = "line")
df2 <- data.frame(time = 30:120, score = sin((30:120)/20)*10, value = rnorm(120-30 +1))
p2 <- ggplot(data = df2, aes(x = time, y = score)) + geom_line() + geom_point(size = 2, aes(color = value))
tracks(time1 = p1, time2 = p2) + xlim(1, 40) + theme_tracks_sunset()

## ----ggbio_granges, eval=TRUE--------------------------------------------
library(GenomicRanges)
set.seed(1); N <- 100; gr <- GRanges(seqnames = sample(c("chr1", "chr2", "chr3"), size = N, replace = TRUE), IRanges(start = sample(1:300, size = N, replace = TRUE), width = sample(70:75, size = N,replace = TRUE)), strand = sample(c("+", "-"), size = N, replace = TRUE), value = rnorm(N, 10, 3), score = rnorm(N, 100, 30), sample = sample(c("Normal", "Tumor"), size = N, replace = TRUE), pair = sample(letters, size = N, replace = TRUE))
autoplot(gr, aes(color = strand, fill = strand), facets = strand ~ seqnames)

## ----ggbio_coverage, eval=TRUE-------------------------------------------
autoplot(gr, aes(color = strand, fill = strand), facets = strand ~ seqnames, stat = "coverage")

## ----ggbio_mirrored_coverage, eval=TRUE----------------------------------
pos <- sapply(coverage(gr[strand(gr)=="+"]), as.numeric)
pos <- data.frame(Chr=rep(names(pos), sapply(pos, length)), Strand=rep("+", length(unlist(pos))), Position=unlist(sapply(pos, function(x) 1:length(x))), Coverage=as.numeric(unlist(pos)))
neg <- sapply(coverage(gr[strand(gr)=="-"]), as.numeric)
neg <- data.frame(Chr=rep(names(neg), sapply(neg, length)), Strand=rep("-", length(unlist(neg))), Position=unlist(sapply(neg, function(x) 1:length(x))), Coverage=-as.numeric(unlist(neg)))
covdf <- rbind(pos, neg)
p <- ggplot(covdf, aes(Position, Coverage, fill=Strand)) + 
	    geom_bar(stat="identity", position="identity") + facet_wrap(~Chr)
p

## ----ggbio_circular1, eval=TRUE------------------------------------------
ggplot(gr) + layout_circle(aes(fill = seqnames), geom = "rect")

## ----ggbio_circular2, eval=FALSE-----------------------------------------
## seqlengths(gr) <- c(400, 500, 700)
## values(gr)$to.gr <- gr[sample(1:length(gr), size = length(gr))]
## idx <- sample(1:length(gr), size = 50)
## gr <- gr[idx]
## ggplot() + layout_circle(gr, geom = "ideo", fill = "gray70", radius = 7, trackWidth = 3) +
##   layout_circle(gr, geom = "bar", radius = 10, trackWidth = 4,
##                 aes(fill = score, y = score)) +
##   layout_circle(gr, geom = "point", color = "red", radius = 14,
##                 trackWidth = 3, grid = TRUE, aes(y = score)) +
##   layout_circle(gr, geom = "link", linked.to = "to.gr", radius = 6, trackWidth = 1)

## ----ggbio_align_variants, eval=TRUE, warning=FALSE, message=FALSE-------
library(rtracklayer); library(GenomicFeatures); library(Rsamtools); library(GenomicAlignments); library(VariantAnnotation)
ga <- readGAlignments("./data/SRR064167.fastq.bam", use.names=TRUE, param=ScanBamParam(which=GRanges("Chr5", IRanges(4000, 8000))))
p1 <- autoplot(ga, geom = "rect")
p2 <- autoplot(ga, geom = "line", stat = "coverage")
vcf <- readVcf(file="data/varianttools_gnsap.vcf", genome="ATH1")
p3 <- autoplot(vcf[seqnames(vcf)=="Chr5"], type = "fixed") + xlim(4000, 8000) + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y=element_blank())
txdb <- makeTxDbFromGFF(file="./data/TAIR10_GFF3_trunc.gff", format="gff3")
p4 <- autoplot(txdb, which=GRanges("Chr5", IRanges(4000, 8000)), names.expr = "gene_id")
tracks(Reads=p1, Coverage=p2, Variant=p3, Transcripts=p4, heights = c(0.3, 0.2, 0.1, 0.35)) + ylab("")

## ----systempiper_symbolic_links, eval=FALSE------------------------------
## library("systemPipeR")
## symLink2bam(sysargs=args, htmldir=c("~/.html/", "somedir/"),
##             urlbase="http://myserver.edu/~username/",
##             urlfile="IGVurl.txt")

## ----control_igv, eval=FALSE, warning=FALSE, message=FALSE---------------
## library(SRAdb)
## startIGV("lm")
## sock <- IGVsocket()
## session <- IGVsession(files=myurls,
##                       sessionFile="session.xml",
##                       genome="A. thaliana (TAIR10)")
## IGVload(sock, session)
## IGVgoto(sock, 'Chr1:45296-47019')

