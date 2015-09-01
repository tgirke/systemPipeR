\name{plotfeaturetypeCounts}
\alias{plotfeaturetypeCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot feature type counts
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotfeaturetypeCounts(x, graphicsfile, graphicsformat = "pdf", scales = "fixed", anyreadlength = FALSE, drop_N_total_aligned = TRUE, scale_count_val = 10^6, scale_length_val = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{graphicsfile}{
%%     ~~Describe \code{graphicsfile} here~~
}
  \item{graphicsformat}{
%%     ~~Describe \code{graphicsformat} here~~
}
  \item{scales}{
%%     ~~Describe \code{scales} here~~
}
  \item{anyreadlength}{
%%     ~~Describe \code{anyreadlength} here~~
}
  \item{drop_N_total_aligned}{
%%     ~~Describe \code{drop_N_total_aligned} here~~
}
  \item{scale_count_val}{
%%     ~~Describe \code{scale_count_val} here~~
}
  \item{scale_length_val}{
%%     ~~Describe \code{scale_length_val} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, graphicsfile, graphicsformat = "pdf", scales = "fixed", 
    anyreadlength = FALSE, drop_N_total_aligned = TRUE, scale_count_val = 10^6, 
    scale_length_val = 1000) 
{
    if (class(x) != "data.frame") 
        stop("x needs to be object of class 'data.frame'")
    if (any(colnames(x)[1:3] != c("SampleName", "Strand", "Featuretype"))) 
        stop("First three column names need to be: 'SampleName', 'Strand', 'Featuretype'")
    if (!is.numeric(as.matrix(x[, -c(1:3)]))) 
        stop("Following columns (after 3rd one) need to be numeric.")
    featuretypelength <- tapply(x$Featuretypelength, x$Featuretype, 
        unique)
    featuretypelength <- featuretypelength[!is.na(featuretypelength)]
    x <- x[, !colnames(x) \%in\% "Featuretypelength"]
    N_total_aligned_DF <- x[x$Featuretype == "N_total_aligned", 
        ]
    if (colnames(N_total_aligned_DF)[4] != "anyreadlength") 
        N_total_aligned_DF <- data.frame(N_total_aligned_DF[, 
            1:3], anyreadlength = rowSums(N_total_aligned_DF[, 
            -c(1:3)]))
    N_total_aligned <- tapply(as.numeric(as.vector(N_total_aligned_DF[, 
        4])), N_total_aligned_DF$SampleName, sum)
    if (drop_N_total_aligned == TRUE) 
        x <- x[x$Featuretype != "N_total_aligned", ]
    if (colnames(x)[4] != "anyreadlength" & anyreadlength == 
        TRUE) 
        x <- data.frame(x[, 1:3], anyreadlength = rowSums(x[, 
            -c(1:3)]))
    if (is.numeric(scale_count_val[1])) {
        if (is.numeric(scale_length_val[1])) {
            correctv <- scale_length_val/featuretypelength[as.character(x$Featuretype)]
            for (i in seq_along(colnames(x))[-c(1:3)]) x[, i] <- as.numeric(as.vector(x[, 
                i])) * (scale_count_val/N_total_aligned[x$SampleName]) * 
                correctv
        }
        else if (length(scale_length_val) == 0) {
            for (i in seq_along(colnames(x))[-c(1:3)]) x[, i] <- as.numeric(as.vector(x[, 
                i])) * (scale_count_val/N_total_aligned[x$SampleName])
        }
        else {
            stop("'scale_length_val' needs to be assinged NULL or numeric value.")
        }
    }
    else if (is.null(scale_count_val[1])) {
        if (is.numeric(scale_length_val[1])) {
            correctv <- scale_length_val/featuretypelength[as.character(x$Featuretype)]
            for (i in seq_along(colnames(x))[-c(1:3)]) x[, i] <- as.numeric(as.vector(x[, 
                i])) * correctv
        }
        else if (length(scale_length_val) == 0) {
            x <- x
        }
        else {
            stop("'scale_length_val' needs to be assinged NULL or numeric value.")
        }
    }
    else {
        stop("'scale_count_val' needs to be assinged NULL or numeric value.")
    }
    x[, 1] <- factor(x[, 1], levels = unique(x[, 1]), ordered = TRUE)
    if (colnames(x)[4] == "anyreadlength") {
        x[, 3] <- factor(x[, 3], levels = sort(as.character(unique(x[, 
            3])), decreasing = TRUE), ordered = TRUE)
        colnames(x)[3:4] <- c("Feature", "Counts")
        myplot <- ggplot(x, aes(x = Feature, y = Counts)) + geom_bar(aes(fill = Strand), 
            position = "stack", stat = "identity") + facet_wrap(~SampleName) + 
            theme(axis.text.x = element_text(angle = -90, hjust = 0, 
                vjust = 0.4)) + coord_flip()
        get(graphicsformat)(graphicsfile)
        print(myplot)
        dev.off()
        cat("Generated graphics file", graphicsfile, "\n")
        return(myplot)
    }
    else {
        .convert2DFlist <- function(x) {
            featureCountslist <- split(x, x$SampleName, drop = TRUE)
            featureCountslist <- featureCountslist[unique(as.character(x$SampleName))]
            for (j in seq_along(featureCountslist)) {
                df <- featureCountslist[[j]]
                tmpDF <- data.frame()
                for (i in seq_along(df[, 4:length(df[1, ])])) {
                  tmpDF <- rbind(tmpDF, data.frame(df[, 1:3], 
                    Counts = df[, i + 3], Length = colnames(df[, 
                      i + 3, drop = FALSE])))
                }
                featureCountslist[[j]] <- tmpDF
            }
            return(featureCountslist)
        }
        featureCountslist <- .convert2DFlist(x)
        tmp <- data.frame(paste(x[, 1], x[, 3], sep = "_"), x[, 
            -c(1:3)])
        mymax <- max(t(sapply(unique(tmp[, 1]), function(x) colSums(tmp[tmp[, 
            1] == x, -1]))))
        myplotlist <- sapply(names(featureCountslist), function(x) NULL, 
            simplify = FALSE)
        for (i in seq(along = featureCountslist)) {
            myplot <- ggplot(featureCountslist[[i]], aes(x = Length, 
                y = Counts)) + geom_bar(aes(fill = Strand), position = "stack", 
                stat = "identity") + facet_wrap(~Featuretype, 
                ncol = 1, scales = scales) + theme(legend.position = "bottom") + 
                theme(axis.text.x = element_text(angle = -90, 
                  hjust = 0, vjust = 0.4)) + ggtitle(featureCountslist[[i]][1, 
                1])
            if (scales == "fixed") 
                myplot <- myplot + ylim(0, mymax)
            myplotlist[[i]] <- myplot
        }
        if (tolower(graphicsformat) == "pdf") 
            mydim <- c(height = 20, width = 6 * length(myplotlist))
        if (tolower(graphicsformat) \%in\% c("png", "jpeg")) 
            mydim <- c(height = 20 * 96, width = 6 * length(myplotlist) * 
                96)
        get(graphicsformat)(graphicsfile, width = mydim["width"], 
            height = mydim["height"])
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(1, length(myplotlist))))
        for (i in seq(along = myplotlist)) print(myplotlist[[i]], 
            vp = viewport(layout.pos.row = 1, layout.pos.col = i))
        dev.off()
        cat("Generated graphics file", graphicsfile, "\n")
        return(myplotlist)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
