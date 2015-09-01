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

#############################################################
## Compute and Plot Read Distribution Across Feature Types ##
#############################################################
## Compute distribution of reads across feature types
featuretypeCounts <- function(bfl, grl, singleEnd=TRUE, readlength=NULL, type="data.frame") {
    ## Check for valid inputs
    if(!(is.null(readlength[1]) | is.integer(readlength))) stop("readlength needs to be assigned NULL or integer vector")
    if(length(type) != 1) stop("Argument needs to be character vector of length 1.")
    if(!any(c("list", "data.frame") %in% type)) stop("Argument 'type' needs to be assigned one of 'list' or 'data.frame'")
    
    ## List container for storing results
    resultlist <- sapply(names(bfl), function(x) list(NULL))
    
    ## Vector to store counts of total aligned reads
    totalcounts <- numeric(length(bfl)); names(totalcounts) <- names(bfl)

    ## Total length (in bps) of each reduced feature type
    feattypelength <- sapply(names(grl), function(x) sum(width(reduce(grl[[x]]))))

    ## Read count for every feature given in gr (coming from gff)
    for(i in names(grl)) strand(grl[[i]][strand(grl[[i]])=="*"]) <- "+" # Sets "*" to "+" to avoid double-counting of reads in this strand-specific method
    feature <- names(grl)
    featurewithtotal <- c(names(grl), "N_total_aligned")
    for(j in seq(along=bfl)) {
        myMA <- matrix(rep(0, length(featurewithtotal)*length(readlength)), nrow=length(featurewithtotal), ncol=length(readlength), dimnames=list(featurewithtotal, readlength))
        if(is.null(readlength[1])) colnames(myMA) <- "anyreadlength"
        myMAsense <- myMA
        bf <- open(bfl[[j]])
        while(isIncomplete(bf)) {
            if(singleEnd==TRUE) {
                aligns <- readGAlignments(bf, use.names=FALSE)
            } else if(singleEnd==FALSE) {
                aligns <- readGAlignmentPairs(bf, use.names=FALSE)
            } else {
                stop("singleEnd needs to be TRUE or FALSE")
            }
            for(i in feature) {
                gr <- grl[[i]]
                counter <- which(feature==i)
                ## Non-strand-specific counts
                alignssub <- subsetByOverlaps(aligns, gr, ignore.strand=TRUE)
                if(singleEnd==TRUE) {
                    if(is.null(readlength[1])) {
                        counts <- length(alignssub)
                        if(counter==1) totalcounts <- length(aligns)
                    } else {
                        counts <- table(qwidth(alignssub))[colnames(myMA)]
                        if(counter==1) totalcounts <- table(qwidth(aligns))[colnames(myMA)]
                    } 
                } else { # The following reports for PE reads the mean read length for each pair
                    if(is.null(readlength[1])) {
                        counts <- length(alignssub)
                        if(counter==1) totalcounts <- length(aligns)
                    } else {
                        mywidth <- round((qwidth(last(alignssub)) + qwidth(first(alignssub)))/2)
                        counts <- table(mywidth)[colnames(myMA)]
                        if(counter==1) { 
                            totalmywidth <- round((qwidth(last(aligns)) + qwidth(first(aligns)))/2)
                            totalcounts <- table(totalmywidth)[colnames(myMA)]
                        }
                    } 
                }
                counts[is.na(counts)] <- 0
                myMA[i,] <- myMA[i,] + counts
                if(counter==1) {
                    totalcounts[is.na(totalcounts)] <- 0
                    myMA["N_total_aligned",] <- myMA["N_total_aligned",] + totalcounts
                }
                ## Strand-specific counts
                alignssub <- subsetByOverlaps(aligns, gr, ignore.strand=FALSE)
                if(singleEnd==TRUE) {
                    if(is.null(readlength[1])) {
                        sensecounts <- length(alignssub)
                        if(counter==1) totalsensecounts <- length(aligns[strand(aligns)=="+"])
                    } else {
                        sensecounts <- table(qwidth(alignssub))[colnames(myMAsense)]
                        if(counter==1) totalsensecounts <- table(qwidth(aligns[strand(aligns)=="+"]))[colnames(myMAsense)]
                    }
                } else { # The following reports for PE reads the mean read length for each pair
                    if(is.null(readlength[1])) {
                        sensecounts <- length(alignssub)
                        if(counter==1) totalsensecounts <- length(aligns[strand(aligns)=="+"])
                    } else {
                        mywidth <- round((qwidth(last(alignssub)) + qwidth(first(alignssub)))/2)
                        sensecounts <- table(mywidth)[colnames(myMA)]
                        if(counter==1) { 
                            totalmywidth <- round((qwidth(last(aligns)) + qwidth(first(aligns)))/2)
                            totalsensecounts <- table(totalmywidth)[colnames(myMA)]
                        }
                    } 
                }
                sensecounts[is.na(sensecounts)] <- 0
                myMAsense[i,] <- myMAsense[i,] + sensecounts
                if(counter==1) {
                    totalsensecounts[is.na(totalsensecounts)] <- 0
                    myMAsense["N_total_aligned",] <- myMAsense["N_total_aligned",] + totalsensecounts
                }
            }
        }
        close(bf)
        myMAantisense <- myMA-myMAsense
        myMAantisense[myMAantisense < 0] <- 0; myMAsense[myMAsense <- 0] <- 0
        myMAsense <- data.frame(SampleName=names(bfl)[j], Strand="Sense", Featuretype=row.names(myMAsense), Featuretypelength=as.numeric(feattypelength[row.names(myMAsense)]), myMAsense, check.names=FALSE)
        myMAantisense <- data.frame(SampleName=names(bfl)[j], Strand="Antisense", Featuretype=row.names(myMAantisense), Featuretypelength=as.numeric(feattypelength[row.names(myMAantisense)]), myMAantisense, check.names=FALSE)
        resultlist[[j]] <- list(Sense=myMAsense, Antisense=myMAantisense)
        cat(paste0("Processed sample ", j, ": ", names(bfl)[j]), "\n")
    }
    
    ## Object types to return
    if(type=="list") {
        return(resultlist)
    } else if(type=="data.frame") {
        resultdf <- do.call("rbind", unlist(resultlist, recursive=FALSE))
        row.names(resultdf) <- NULL
        return(resultdf)
    } else {
        stop("Argument 'type' needs to be assigned one of 'list' or 'data.frame'")
    }
}

## Usage: 
# featureCounts <- featuretypeCounts(bfl=BamFileList(outpaths(args), yieldSize=50000), grl=feat, singleEnd=TRUE, readlength=c(74:76,99:102), type="data.frame")
# write.table(featureCounts, "results/featureCounts.xls", quote=FALSE, row.names=FALSE, sep="\t")
# featureCounts2 <- featuretypeCounts(bfl=BamFileList(outpaths(args), yieldSize=50000), grl=feat, singleEnd=TRUE, readlength=NA, type="data.frame")
# write.table(featureCounts2, "results/featureCounts2.xls", quote=FALSE, row.names=FALSE, sep="\t")

## Plot distribution of reads across feature types
plotfeaturetypeCounts <- function(x, graphicsfile, graphicsformat="pdf", scales="fixed", anyreadlength=FALSE, drop_N_total_aligned=TRUE, scale_count_val=10^6, scale_length_val=10^3) { 
    ## Input validity checks
    if(class(x)!="data.frame") stop("x needs to be object of class 'data.frame'")
    if(any(colnames(x)[1:3]!=c("SampleName", "Strand", "Featuretype"))) stop("First three column names need to be: 'SampleName', 'Strand', 'Featuretype'")
    if(!is.numeric(as.matrix(x[,-c(1:3)]))) stop("Following columns (after 3rd one) need to be numeric.")
    
    ## Store Featuretypelength and remove it from x
    featuretypelength <- tapply(x$Featuretypelength, x$Featuretype, unique)
    featuretypelength <- featuretypelength[!is.na(featuretypelength)]
    x <- x[, !colnames(x) %in% "Featuretypelength"]
        
    ## Get numbers of total aligned reads per sample
    N_total_aligned_DF <- x[x$Featuretype=="N_total_aligned",]
    if(colnames(N_total_aligned_DF)[4]!="anyreadlength") N_total_aligned_DF <- data.frame(N_total_aligned_DF[,1:3], anyreadlength=rowSums(N_total_aligned_DF[,-c(1:3)]))
    N_total_aligned <- tapply(as.numeric(as.vector(N_total_aligned_DF[,4])), N_total_aligned_DF$SampleName, sum)

    ## Eliminate total read counts to not show them in plots
    if(drop_N_total_aligned==TRUE) x <- x[x$Featuretype!="N_total_aligned",]
    
    ## Sum up read length counts if 'anyreadlength=TRUE'  
    if(colnames(x)[4]!="anyreadlength" & anyreadlength==TRUE) x <- data.frame(x[,1:3], anyreadlength=rowSums(x[,-c(1:3)]))
    
    ## Scale per x reads (e.g. per million reads)
    if(is.numeric(scale_count_val[1])) {
        if(is.numeric(scale_length_val[1])) {
            correctv <- scale_length_val/featuretypelength[as.character(x$Featuretype)]
            for(i in seq_along(colnames(x))[-c(1:3)]) x[,i] <- as.numeric(as.vector(x[,i])) * (scale_count_val/N_total_aligned[x$SampleName]) * correctv
        } else if(length(scale_length_val) == 0) {
            for(i in seq_along(colnames(x))[-c(1:3)]) x[,i] <- as.numeric(as.vector(x[,i])) * (scale_count_val/N_total_aligned[x$SampleName])
        } else {
            stop("'scale_length_val' needs to be assinged NULL or numeric value.")
        }
    } else if(is.null(scale_count_val[1])) {
        if(is.numeric(scale_length_val[1])) {
            correctv <- scale_length_val/featuretypelength[as.character(x$Featuretype)]
            for(i in seq_along(colnames(x))[-c(1:3)]) x[,i] <- as.numeric(as.vector(x[,i])) * correctv
        } else if(length(scale_length_val) == 0) {
            x <- x
        } else {
            stop("'scale_length_val' needs to be assinged NULL or numeric value.")
        }
    } else {
        stop("'scale_count_val' needs to be assinged NULL or numeric value.")
    }

    ## Define plotting order of samples
    x[,1] <- factor(x[,1], levels=unique(x[,1]), ordered=TRUE) 
    
    ## Plot without read length resolution
    if(colnames(x)[4]=="anyreadlength") {
        x[,3] <- factor(x[,3], levels=sort(as.character(unique(x[,3])), decreasing=TRUE), ordered=TRUE) # Defines plotting order of bars!!!
        colnames(x)[3:4] <- c("Feature", "Counts") # Assign nicer names for plot
        myplot <- ggplot(x, aes(x=Feature, y=Counts)) + 
                        geom_bar(aes(fill=Strand), position="stack", stat="identity") + 
                        facet_wrap(~SampleName) + 
                        theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.4)) +
                        coord_flip()
                        
        get(graphicsformat)(graphicsfile)
            print(myplot)
        dev.off()
        cat("Generated graphics file", graphicsfile, "\n")
        return(myplot)
    ## Plot with read length resolution
    } else {
    ## Split data on samples (SampleName column) and convert to format expected by ggplot2
        .convert2DFlist <- function(x) {
            featureCountslist <- split(x, x$SampleName, drop=TRUE)
            featureCountslist <- featureCountslist[unique(as.character(x$SampleName))]
            for(j in seq_along(featureCountslist)) {
                df <- featureCountslist[[j]]
                tmpDF <- data.frame()
                for(i in seq_along(df[,4:length(df[1,])])) {
                    tmpDF <- rbind(tmpDF, data.frame(df[,1:3], Counts=df[,i+3], Length=colnames(df[,i+3,drop=FALSE])))
                }
                featureCountslist[[j]] <- tmpDF
            }
            return(featureCountslist)
        }
        featureCountslist <- .convert2DFlist(x)
    
        ## Get maximum read counts in data set after summing up strand spec counts
        tmp <- data.frame(paste(x[,1], x[,3], sep="_"), x[,-c(1:3)])
        mymax <- max(t(sapply(unique(tmp[,1]), function(x) colSums(tmp[tmp[,1]==x,-1]))))
        ## Generate plotting instructions
        myplotlist <- sapply(names(featureCountslist), function(x) NULL, simplify=FALSE)
        for(i in seq(along=featureCountslist)) {
            myplot <- ggplot(featureCountslist[[i]], aes(x=Length, y=Counts)) + 
                             geom_bar(aes(fill=Strand), position="stack", stat="identity") + 
                             facet_wrap(~Featuretype, ncol=1, scales=scales) + 
                             theme(legend.position="bottom") +
                             theme(axis.text.x=element_text(angle=-90, hjust=0, vjust=0.4)) +
                             ggtitle(featureCountslist[[i]][1,1])
            ## Assure same scale for all panels if scales="fixed"
            if(scales=="fixed") myplot <- myplot + ylim(0, mymax)
            myplotlist[[i]] <- myplot
        }
    
        ## Generate graphics and write to file in pdf, png or jpeg format.
        if(tolower(graphicsformat)=="pdf") mydim <- c(height=20, width=6*length(myplotlist))
        if(tolower(graphicsformat) %in% c("png", "jpeg")) mydim <- c(height=20*96, width=6*length(myplotlist)*96)
        get(graphicsformat)(graphicsfile, width=mydim["width"], height=mydim["height"])
            grid.newpage() # Open a new page on grid device
            pushViewport(viewport(layout = grid.layout(1, length(myplotlist)))) # Assign to device viewport with 1 by 2 grid layout 
            for(i in seq(along=myplotlist)) print(myplotlist[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i))
        dev.off()
        cat("Generated graphics file", graphicsfile, "\n")
        return(myplotlist)
    }
}

## Usage: 
# library(ggplot2); library(grid)
# featureCounts <- read.delim("results/featureCounts.xls", check.names=FALSE)
# myplots <- plotfeaturetypeCounts(x=featureCounts, graphicsfile="results/featureCounts.pdf", graphicsformat="pdf", scales="fixed", anyreadlength=FALSE, drop_N_total_aligned=TRUE, scale_count_val=10^6, scale_length_val=10^3) 
# featureCounts2 <- read.delim("results/featureCounts2.xls", check.names=FALSE)
# myplots <- plotfeaturetypeCounts(x=featureCounts2, graphicsfile="results/featureCounts2.pdf", graphicsformat="pdf", scales="fixed", anyreadlength=TRUE, drop_N_total_aligned=TRUE, scale_count_val=10^6, scale_length_val=10^3) 

