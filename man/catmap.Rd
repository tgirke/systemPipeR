\name{catmap}
\alias{catmap}
\alias{catmap-methods}
\alias{catlist}
\alias{catlist-methods}
\alias{idconv}
\alias{idconv-methods}
\docType{data}
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
	various outputs
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
go <- getBM(attributes=c("go_accession", "tair_locus", "go_namespace_1003"), mart=m)
go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
write.table(go, "GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

## Create catDB instance (takes a while but needs to be done only once)
catdb <- makeCATdb(myfile="GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
catdb

## Access methods for catDB
catmap(catdb)$D_MF[1:4,]
catlist(catdb)$L_MF[1:4]
idconv(catdb)

}
}
\keyword{ utilities }

