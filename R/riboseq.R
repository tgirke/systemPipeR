######################################################
## Generate various feature types from TxDb objects ##
######################################################
genFeatures <- function(txdb, featuretype="all", reduce_ranges, upstream=1000, downstream=0, verbose=TRUE) {
    ## Check validity of inputs
    supported_features <- c("tx_type", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic")
    if(tolower(featuretype[1])=="all") featuretype <- supported_features
    if(class(txdb)!="TxDb") stop("Argument 'txdb' needs to be assigned object of class TxDb!")
    if(any(!featuretype %in% supported_features)) stop("featuretype supports any of the following values: ", paste(supported_features, collapse=", "))
	featuresGRl <- GRangesList()
	
    ## Gene/transcript id mappings (required for some features)
    ids <- mcols(transcripts(txdb,  columns=c("tx_id", "tx_name", "gene_id")))
    ge_id <- unstrsplit(ids$gene_id)
    names(ge_id) <- ids$tx_id
	
    ## Transcript ranges: each 'tx_type' as separate GRanges object (reduced by gene)
    if("tx_type" %in% featuretype) {
        tx <- transcripts(txdb, columns=c("tx_name", "gene_id", "tx_type"))
        tx_type <- unique(mcols(tx)$tx_type)
        for(i in seq_along(tx_type)) {
            tmp <- tx[mcols(tx)$tx_type==tx_type[i]]
            if(reduce_ranges==TRUE) {
                tmp <- reduce(split(tmp, as.character(mcols(tmp)$gene_id))) 
                tmp <- unlist(tmp)
                feature_id <- paste(names(tmp), ":", tx_type[i], "_red", sep="")
		        mcols(tmp) <- DataFrame(feature_by=names(tmp), featuretype_id=feature_id, featuretype=paste0(tx_type[i], "_red"))
	            names(tmp) <- seq_along(tmp)
		        featuresGRl <- c(featuresGRl, GRangesList(tmp=tmp))
                names(featuresGRl)[length(featuresGRl)] <- paste0(tx_type[i], "_red")
            } else {
		        mcols(tmp) <- DataFrame(feature_by=mcols(tmp)$gene_id, featuretype_id=mcols(tmp)$tx_name, featuretype=tx_type[i])
                featuresGRl <- c(featuresGRl, GRangesList(tmp=tmp))
	            names(tmp) <- seq_along(tmp)
                names(featuresGRl)[length(featuresGRl)] <- tx_type[i]
            }
        }
        if(verbose==TRUE) cat("Created feature ranges:", paste(tx_type, collapse=", "), "\n")
    }

    ## Create promoter ranges reduced by gene
	if("promoter" %in% featuretype) {
        if(reduce_ranges==TRUE) {
		    mypromoters <- suppressWarnings(unlist(reduce(promoters(transcriptsBy(txdb, "gene"), upstream, downstream))))
		    mypromoters <- trim(mypromoters)
            feature_id <- paste0(names(mypromoters), ":P_red")
            mcols(mypromoters) <- DataFrame(feature_by=names(mypromoters), featuretype_id=feature_id, featuretype="promoter_red")
	        names(mypromoters) <- seq_along(mypromoters)
		    featuresGRl <- c(featuresGRl, GRangesList("promoter_red"=mypromoters))
        } else {
            mypromoters <- suppressWarnings(promoters(txdb, upstream, downstream, columns=c("tx_name", "gene_id", "tx_type")))
		    mypromoters <- trim(mypromoters)
            mcols(mypromoters) <- DataFrame(feature_by=as(mcols(mypromoters)$gene_id, "CharacterList"), featuretype_id=mcols(mypromoters)$tx_name, featuretype="promoter")
		    featuresGRl <- c(featuresGRl, GRangesList("promoter"=mypromoters))
        }
        if(verbose==TRUE) cat("Created feature ranges: promoter", "\n")
	}
    
    ## Create intron ranges reduced by gene
	if("intron" %in% featuretype) {
        ## Introns by transcript
		myintrons <- intronsByTranscript(txdb)
        ## Convert to introns by genes
        intron_count <- sapply(width(myintrons), length) 
        names(intron_count) <- ge_id[names(intron_count)]
        intron_count <- intron_count[intron_count!=0]
        ge_id_intron <- rep(names(intron_count), intron_count)
        myintrons <- unlist(myintrons)
        myintrons <- split(myintrons, ge_id_intron) # Introns by gene
        if(reduce_ranges==TRUE) {
            ## Reduce introns by gene, use gene ID for naming and append counter
            myintrons <- reduce(myintrons) 
            intron_red_count <- sapply(width(myintrons), length) 
            ge_id_red_intron <- rep(names(intron_red_count), intron_red_count)
            feature_id <- paste(ge_id_red_intron, sprintf("I%03d_red", unlist(sapply(as.integer(intron_red_count), function(x) seq(from=1, to=x)))), sep=":")
            myintrons <- unlist(myintrons) 
		    mcols(myintrons) <- DataFrame(feature_by=ge_id_red_intron, featuretype_id=feature_id, featuretype="intron_red")
	        names(myintrons) <- seq_along(myintrons)
		    featuresGRl <- c(featuresGRl, GRangesList("intron_red"=myintrons))
        } else {
            intron_red_count <- sapply(width(myintrons), length) 
            ge_id_red_intron <- rep(names(intron_red_count), intron_red_count)
            feature_id <- paste(ge_id_red_intron, sprintf("I%03d", unlist(sapply(as.integer(intron_red_count), function(x) seq(from=1, to=x)))), sep=":")
            myintrons <- unlist(myintrons) 
		    mcols(myintrons) <- DataFrame(feature_by=as(ge_id_red_intron, "CharacterList"), featuretype_id=feature_id, featuretype="intron")
	        names(myintrons) <- seq_along(myintrons)
		    featuresGRl <- c(featuresGRl, GRangesList("intron"=myintrons))
        }
        if(verbose==TRUE) cat("Created feature ranges: intron", "\n")
	}
    
    ## Create exon ranges reduced by gene 
	if("exon" %in% featuretype) {
        ## exons by gene
		myexons <- exonsBy(txdb, "gene")
        if(reduce_ranges==TRUE) {
	        myexons <- reduce(myexons)
            exon_red_count <- sapply(width(myexons), length) 
            ge_id_red_exon <- rep(names(exon_red_count), exon_red_count)
            feature_id <- paste(ge_id_red_exon, sprintf("E%03d_red", unlist(sapply(as.integer(exon_red_count), function(x) seq(from=1, to=x)))), sep=":")
            myexons <- unlist(myexons) 
		    mcols(myexons) <- DataFrame(feature_by=ge_id_red_exon, featuretype_id=feature_id, featuretype="exon_red")
	        names(myexons) <- seq_along(myexons)
		    featuresGRl <- c(featuresGRl, GRangesList("exon_red"=myexons))
        } else {
            exon_red_count <- sapply(width(myexons), length) 
            ge_id_red_exon <- rep(names(exon_red_count), exon_red_count)
            feature_id <- paste(ge_id_red_exon, sprintf("E%03d", unlist(sapply(as.integer(exon_red_count), function(x) seq(from=1, to=x)))), sep=":")
            myexons <- unlist(myexons) 
		    mcols(myexons) <- DataFrame(feature_by=as(ge_id_red_exon, "CharacterList"), featuretype_id=feature_id, featuretype="exon")
	        names(myexons) <- seq_along(myexons)
		    featuresGRl <- c(featuresGRl, GRangesList("exon"=myexons))
        
        }
        if(verbose==TRUE) cat("Created feature ranges: exon", "\n")
    }
    
    ## Create CDS ranges reduced by gene 
	if("cds" %in% featuretype) {
        ## CDS by gene
		mycds <- cdsBy(txdb, "gene")
        if(reduce_ranges==TRUE) {
	        mycds <- reduce(mycds)
            cds_red_count <- sapply(width(mycds), length) 
            ge_id_red_cds <- rep(names(cds_red_count), cds_red_count)
            feature_id <- paste(ge_id_red_cds, sprintf("CDS%03d_red", unlist(sapply(as.integer(cds_red_count), function(x) seq(from=1, to=x)))), sep=":")
            mycds <- unlist(mycds) 
		    mcols(mycds) <- DataFrame(feature_by=ge_id_red_cds, featuretype_id=feature_id, featuretype="cds_red")
	        names(mycds) <- seq_along(mycds)
		    featuresGRl <- c(featuresGRl, GRangesList("cds_red"=mycds))
        } else {
            cds_red_count <- sapply(width(mycds), length) 
            ge_id_red_cds <- rep(names(cds_red_count), cds_red_count)
            feature_id <- paste(ge_id_red_cds, sprintf("CDS%03d", unlist(sapply(as.integer(cds_red_count), function(x) seq(from=1, to=x)))), sep=":")
            mycds <- unlist(mycds) 
		    mcols(mycds) <- DataFrame(feature_by=as(ge_id_red_cds, "CharacterList"), featuretype_id=feature_id, featuretype="cds")
	        names(mycds) <- seq_along(mycds)
		    featuresGRl <- c(featuresGRl, GRangesList("cds"=mycds))
        }
        if(verbose==TRUE) cat("Created feature ranges: cds", "\n")
    }

    ## Create 5'UTR ranges reduced by gene
	if("fiveUTR" %in% featuretype) {
        ## 5'UTRs by transcript
		myfiveutr <- fiveUTRsByTranscript(txdb)
        ## Convert to 5'UTRs by genes
        utr_count <- sapply(width(myfiveutr), length) 
        names(utr_count) <- ge_id[names(utr_count)]
        utr_count <- utr_count[utr_count!=0]
        ge_id_utr <- rep(names(utr_count), utr_count)
        myfiveutr <- unlist(myfiveutr)
        myfiveutr <- split(myfiveutr, ge_id_utr) # UTRs by gene
        if(reduce_ranges==TRUE) {
            ## Reduce introns by gene, use gene ID for naming and append counter
            myfiveutr <- reduce(myfiveutr) 
            utr_red_count <- sapply(width(myfiveutr), length) 
            ge_id_red_utr <- rep(names(utr_red_count), utr_red_count)
            feature_id <- paste(ge_id_red_utr, sprintf("fiveUTR%03d_red", unlist(sapply(as.integer(utr_red_count), function(x) seq(from=1, to=x)))), sep=":")
            myfiveutr <- unlist(myfiveutr) 
		    mcols(myfiveutr) <- DataFrame(feature_by=ge_id_red_utr, featuretype_id=feature_id, featuretype="fiveUTR_red")
	        names(myfiveutr) <- seq_along(myfiveutr)
		    featuresGRl <- c(featuresGRl, GRangesList("fiveUTR_red"=myfiveutr))
        } else {
            utr_red_count <- sapply(width(myfiveutr), length) 
            ge_id_red_utr <- rep(names(utr_red_count), utr_red_count)
            feature_id <- paste(ge_id_red_utr, sprintf("fiveUTR%03d", unlist(sapply(as.integer(utr_red_count), function(x) seq(from=1, to=x)))), sep=":")
            myfiveutr <- unlist(myfiveutr) 
		    mcols(myfiveutr) <- DataFrame(feature_by=as(ge_id_red_utr, "CharacterList"), featuretype_id=feature_id, featuretype="fiveUTR")
	        names(myfiveutr) <- seq_along(myfiveutr)
		    featuresGRl <- c(featuresGRl, GRangesList("fiveUTR"=myfiveutr))
        }    
        if(verbose==TRUE) cat("Created feature ranges: fiveUTR", "\n")
	}
    
    ## Create 3'UTR ranges reduced by gene
	if("threeUTR" %in% featuretype) {
        ## 3'UTRs by transcript
		mythreeutr <- threeUTRsByTranscript(txdb)
        ## Convert to 3'UTRs by genes
        utr_count <- sapply(width(mythreeutr), length) 
        names(utr_count) <- ge_id[names(utr_count)]
        utr_count <- utr_count[utr_count!=0]
        ge_id_utr <- rep(names(utr_count), utr_count)
        mythreeutr <- unlist(mythreeutr)
        mythreeutr <- split(mythreeutr, ge_id_utr) # UTRs by gene
        if(reduce_ranges==TRUE) {
            ## Reduce introns by gene, use gene ID for naming and append counter
            mythreeutr <- reduce(mythreeutr) 
            utr_red_count <- sapply(width(mythreeutr), length) 
            ge_id_red_utr <- rep(names(utr_red_count), utr_red_count)
            feature_id <- paste(ge_id_red_utr, sprintf("threeUTR%03d_red", unlist(sapply(as.integer(utr_red_count), function(x) seq(from=1, to=x)))), sep=":")
            mythreeutr <- unlist(mythreeutr) 
		    mcols(mythreeutr) <- DataFrame(feature_by=ge_id_red_utr, featuretype_id=feature_id, featuretype="threeUTR_red")
	        names(mythreeutr) <- seq_along(mythreeutr)
		    featuresGRl <- c(featuresGRl, GRangesList("threeUTR_red"=mythreeutr))
        } else {
            utr_red_count <- sapply(width(mythreeutr), length) 
            ge_id_red_utr <- rep(names(utr_red_count), utr_red_count)
            feature_id <- paste(ge_id_red_utr, sprintf("threeUTR%03d", unlist(sapply(as.integer(utr_red_count), function(x) seq(from=1, to=x)))), sep=":")
            mythreeutr <- unlist(mythreeutr) 
		    mcols(mythreeutr) <- DataFrame(feature_by=as(ge_id_red_utr, "CharacterList"), featuretype_id=feature_id, featuretype="threeUTR")
	        names(mythreeutr) <- seq_along(mythreeutr)
		    featuresGRl <- c(featuresGRl, GRangesList("threeUTR"=mythreeutr))
        }
        if(verbose==TRUE) cat("Created feature ranges: threeUTR", "\n")
	}

    ## Create intergenic ranges
	if("intergenic" %in% featuretype) {
        if(verbose==TRUE & any(is.na(seqlengths(txdb)))) warning("seqlengths missing for: ", paste(head(names(seqlengths(txdb))[is.na(seqlengths(txdb))]), collapse=", "), ". Thus, corresponding chromosome end ranges will be missing in intergenic results.")
        ge <- genes(txdb)
        mynames <- names(ge)
        strand(ge) <- "*"
		ge <- reduce(ge, with.revmap=TRUE)
        myids <- character(length = length(ge))
        revmaplist <- as.list(mcols(ge)$revmap)
        for(i in seq(along=ge)) myids[i] <- paste(mynames[revmaplist[[i]]], collapse="_")
		myintergenics <- gaps(ge)
        myintergenics <- myintergenics[strand(myintergenics)=="*"] # Removes chromosome ranges created by gaps when seqlengths available
        index <- findOverlaps(ge, myintergenics, minoverlap=0)
        myids <- tapply(myids[as.matrix(index)[,1]], as.matrix(index)[,2], paste, collapse="__")
        if(reduce_ranges==TRUE) {
		    mcols(myintergenics) <- DataFrame(feature_by=sprintf("INTER%08d", seq_along(myids)), featuretype_id=as.character(myids), featuretype="intergenic")
	    } else {
		    mcols(myintergenics) <- DataFrame(feature_by=as(sprintf("INTER%08d", seq_along(myids)), "CharacterList"), featuretype_id=as.character(myids), featuretype="intergenic")
        }
        names(myintergenics) <- seq_along(myintergenics)
		featuresGRl <- c(featuresGRl, GRangesList("intergenic"=myintergenics))
        if(verbose==TRUE) cat("Created feature ranges: intergenic", "\n")
	}

    ## Append results
    return(featuresGRl)
}
## Usage:
# file <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
# txdb <- makeTxDbFromGFF(file=file, format="gff3", organism="Arabidopsis")
# feat <- genFeatures(txdb, featuretype="all", reduce_ranges=TRUE, upstream=1000, downstream=0, verbose=TRUE)

