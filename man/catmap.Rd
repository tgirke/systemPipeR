\name{catmap}
\alias{catmap}
\alias{catmap-methods}
\alias{catlist}
\alias{catlist-methods}
\alias{idconv}
\alias{idconv-methods}
\title{ catDB accessor methods}

\description{
	Methods to access information from \code{catDB} object.
}
\usage{
	catmap(x)
}
\arguments{
  \item{x}{
	object of class \code{catDB}
}
}
\value{
	\code{catmap}: \code{list} of \code{data.frame}s with the gene-to-GO
	mappings \code{D_BP}, \code{D_CC} and \code{D_MF} (columns \code{GOID} and
	\code{GeneID}, plus \code{GOCAT} when built from \code{myfile} or
	\code{Count} when built from \code{lib}), and \code{GO_DF} with the GO term
	descriptions (columns \code{GOID}, \code{Term}, \code{Ont}).

	\code{catlist}: \code{list} with components \code{L_BP}, \code{L_CC} and
	\code{L_MF}, each a named \code{list} giving for every GO node the IDs of
	the genes assigned to it or to any of its offspring nodes.

	\code{idconv}: the ID conversion data provided to \code{makeCATdb} via
	its \code{idconv} argument, or \code{NULL} if none was provided.
}
\author{
Thomas Girke
}
\examples{
\dontrun{
## Obtain annotations from BioMart
m <- useMart("ENSEMBL_MART_PLANT"); listDatasets(m) 
m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
listAttributes(m) # Choose data types you want to download
go <- getBM(attributes=c("go_accession", "tair_locus", 
                            "go_namespace_1003"), mart=m)
go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
write.table(go, "GOannotationsBiomart_mod.txt", quote=FALSE, 
            row.names=FALSE, col.names=FALSE, sep="\t")

## Create catDB instance (takes a while but needs to be done only once)
catdb <- makeCATdb(myfile="GOannotationsBiomart_mod.txt", lib=NULL, 
                    org="", colno=c(1,2,3), idconv=NULL)
catdb

## Access methods for catDB
catmap(catdb)$D_MF[1:4,]
catlist(catdb)$L_MF[1:4]
idconv(catdb)

}
}
\keyword{ utilities }

