######################################################
## Generate various feature types from TxDb objects ##
######################################################
genFeatures <- function(txdb, featuretype="all", reduce_ranges, upstream=1000, downstream=0, verbose=TRUE) {
    ## Check validity of inputs
    supported_features <- c("tx_type", "promoter", "intron", "exon", "cds", "fiveUTR", "threeUTR", "intergenic")
    if(tolower(featuretype[1])=="all") featuretype <- supported_features
    if(class(txdb)!="TxDb") stop("Argument 'txdb' needs to be assigned object of class TxDb!")
    if(any(!featuretype %in% supported_features)) stop("featuretype supports any of the following values: ", paste(supported_features, collapse=", "))
	
    ## Empty GRangesList to store results
    featuresGRl <- GRangesList()
    GenomeInfoDb::seqinfo(featuresGRl) <- GenomeInfoDb::seqinfo(txdb)
    
    ## Empty GRanges object for feature types that don't exist
    gr_empty <- GRanges()
    mcols(gr_empty) <- DataFrame(feature_by=character(), featuretype_id=character(), featuretype=character())
    GenomeInfoDb::seqinfo(gr_empty) <- GenomeInfoDb::seqinfo(txdb)

    ## Gene/transcript id mappings (required for some features)
    ids <- mcols(transcripts(txdb,  columns=c("tx_id", "tx_name", "gene_id")))
    ge_id <- unstrsplit(ids$gene_id)
    names(ge_id) <- ids$tx_id
	
    ## Transcript ranges: each 'tx_type' as separate GRanges object (reduced by gene)
    if("tx_type" %in% featuretype) {
        tx <- transcripts(txdb, columns=c("tx_name", "gene_id", "tx_type"))
        if(length(tx)==0) { # If range set is empty, append 'gr_empty'.
            featuresGRl <- c(featuresGRl, GRangesList("transcript"=gr_empty))
        } else {
            mcols(tx)$tx_type <- ifelse(is.na(mcols(tx)$tx_type), "transcript", mcols(tx)$tx_type)
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
                    featuresGRl <- c(tmp=featuresGRl, GRangesList(tmp=tmp))
                    names(tmp) <- seq_along(tmp)
                    names(featuresGRl)[length(featuresGRl)] <- tx_type[i]
                }
            }
            if(verbose==TRUE) cat("Created feature ranges:", paste(tx_type, collapse=", "), "\n")
        }
    }

    ## Create promoter ranges reduced by gene
	if("promoter" %in% featuretype) {
        mycheck <- suppressWarnings(promoters(transcriptsBy(txdb, "gene"), upstream, downstream))
        if(length(mycheck)==0) { # If range set is empty, append 'gr_empty'.
            featuresGRl <- c(featuresGRl, GRangesList("promoter"=gr_empty))
        } else {
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
    }
    
    ## Create intron ranges reduced by gene
	if("intron" %in% featuretype) {
        ## Introns by transcript
		myintrons <- intronsByTranscript(txdb)
        if(length(myintrons)==0) { # If range set is empty, append 'gr_empty'.
            featuresGRl <- c(featuresGRl, GRangesList("intron"=gr_empty))
        } else {
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
    }
    
    ## Create exon ranges reduced by gene 
	if("exon" %in% featuretype) {
        ## exons by gene
		myexons <- exonsBy(txdb, "gene")
        if(length(myexons)==0) { # If range set is empty, append 'gr_empty'.
            featuresGRl <- c(featuresGRl, GRangesList("exon"=gr_empty))
        } else {
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
    }
    
    ## Create CDS ranges reduced by gene 
	if("cds" %in% featuretype) {
        ## CDS by gene
		mycds <- cdsBy(txdb, "gene")
        if(length(mycds)==0) { # If range set is empty, append 'gr_empty'.
            featuresGRl <- c(featuresGRl, GRangesList("cds"=gr_empty))
        } else {
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
    }

    ## Create 5'UTR ranges reduced by gene
	if("fiveUTR" %in% featuretype) {
        ## 5'UTRs by transcript
		myfiveutr <- fiveUTRsByTranscript(txdb)
        if(length(myfiveutr)==0) { # If range set is empty, append 'gr_empty'.
            featuresGRl <- c(featuresGRl, GRangesList("fiveUTR"=gr_empty))
        } else {
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
    }
    
    ## Create 3'UTR ranges reduced by gene
	if("threeUTR" %in% featuretype) {
        ## 3'UTRs by transcript
		mythreeutr <- threeUTRsByTranscript(txdb)
        if(length(mythreeutr)==0) { # If range set is empty, append 'gr_empty'.
            featuresGRl <- c(featuresGRl, GRangesList("threeUTR"=gr_empty))
        } else {
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
    }

    ## Create intergenic ranges
	if("intergenic" %in% featuretype) {
        # if(verbose==TRUE & any(is.na(seqlengths(txdb)))) warning("seqlengths missing for: ", paste(head(names(seqlengths(txdb))[is.na(seqlengths(txdb))]), collapse=", "), ". Thus, corresponding chromosome end ranges will be missing in intergenic results.")
        ge <- genes(txdb)
        if(length(ge)==0) { # If range set is empty, append 'gr_empty'.
            featuresGRl <- c(featuresGRl, GRangesList("intergenic"=gr_empty))
        } else {
            myseqinfo <- GenomeInfoDb::seqinfo(ge)
            GenomeInfoDb::seqlengths(ge) <- NA # Note: the following ignores chromosome end ranges on the right end
            mynames <- names(ge)
            strand(ge) <- "*"
            ge <- reduce(ge, with.revmap=TRUE)
            myids <- character(length = length(ge))
            revmaplist <- as.list(mcols(ge)$revmap)
            for(i in seq(along=ge)) myids[i] <- paste(mynames[revmaplist[[i]]], collapse="_")
            myintergenics <- gaps(ge)
            myintergenics <- myintergenics[strand(myintergenics)=="*"] # Removes chromosome ranges created by gaps when seqlengths available
            # index <- findOverlaps(ge, myintergenics, minoverlap=0L) 
            index <- findOverlaps(ge, myintergenics, maxgap=0L, minoverlap=0L) # change necessary since default in IRanges changed from maxgap=0L to maxgap=-1L
            myids <- tapply(myids[as.matrix(index)[,1]], as.matrix(index)[,2], paste, collapse="__")
            if(reduce_ranges==TRUE) {
                mcols(myintergenics) <- DataFrame(feature_by=sprintf("INTER%08d", seq_along(myids)), featuretype_id=as.character(myids), featuretype="intergenic")
            } else {
                mcols(myintergenics) <- DataFrame(feature_by=as(sprintf("INTER%08d", seq_along(myids)), "CharacterList"), featuretype_id=as.character(myids), featuretype="intergenic")
            }
            names(myintergenics) <- seq_along(myintergenics)
            GenomeInfoDb::seqinfo(myintergenics) <- myseqinfo
            featuresGRl <- c(featuresGRl, GRangesList("intergenic"=myintergenics))
            if(verbose==TRUE) cat("Created feature ranges: intergenic", "\n")
	    }
    }

    ## Return results
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
    ## global functions or variables
    readGAlignments <- readGAlignmentPairs <- qwidth <- last <- strand <- NULL
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
        if(is.null(readlength[1])) {
            myMA <- matrix(rep(0, length(featurewithtotal)*1), nrow=length(featurewithtotal), ncol=1, dimnames=list(featurewithtotal, "anyreadlength"))
        } else {
            myMA <- matrix(rep(0, length(featurewithtotal)*length(readlength)), nrow=length(featurewithtotal), ncol=length(readlength), dimnames=list(featurewithtotal, readlength))
        }
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
# outpaths <- subsetWF(args , slot="output", subset=1, index=1)
# file.exists(outpaths)
# featureCounts <- featuretypeCounts(bfl=BamFileList(outpaths, yieldSize=50000), grl=feat, singleEnd=TRUE, readlength=c(74:76,99:102), type="data.frame")
# write.table(featureCounts, "results/featureCounts.xls", quote=FALSE, row.names=FALSE, sep="\t")
# featureCounts2 <- featuretypeCounts(bfl=BamFileList(outpaths), grl=feat, singleEnd=TRUE, readlength=NULL, type="data.frame")
# write.table(featureCounts2, "results/featureCounts2.xls", quote=FALSE, row.names=FALSE, sep="\t")

######################################################
## Plot distribution of reads across feature types ##
###################################################### 
plotfeaturetypeCounts <- function(x, graphicsfile, graphicsformat="pdf", scales="fixed", anyreadlength=FALSE, drop_N_total_aligned=TRUE, scale_count_val=10^6, scale_length_val=NULL) { 
    ## Input validity checks
    if(class(x)!="data.frame") stop("x needs to be object of class 'data.frame'")
    if(any(colnames(x)[1:3]!=c("SampleName", "Strand", "Featuretype"))) stop("First three column names need to be: 'SampleName', 'Strand', 'Featuretype'")
    if(!is.numeric(as.matrix(x[,-c(1:3)]))) stop("Following columns (after 3rd one) need to be numeric.")
    ## global functions or variables
    Strand <- Feature <- Counts <- Length <- NULL
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
    
    ## Label for count axis
    if(length(scale_count_val)==0 & length(scale_length_val)==0) axis_label <- "Raw counts"
    if(length(scale_count_val)>0 & length(scale_length_val)>0) axis_label <- paste("Counts normalized per", scale_count_val, "reads", "and", scale_length_val, "bp feature length")
    if(length(scale_count_val)>0 & length(scale_length_val)==0) axis_label <- paste("Counts normalized per", scale_count_val, "reads") 
    if(length(scale_count_val)==0 & length(scale_length_val)>0) axis_label <- paste("Counts normalized per", scale_length_val, "bp feature length")

    ## Plot without read length resolution
    if(colnames(x)[4]=="anyreadlength") {
        x[,3] <- factor(x[,3], levels=sort(as.character(unique(x[,3])), decreasing=TRUE), ordered=TRUE) # Defines plotting order of bars!!!
        colnames(x)[3:4] <- c("Feature", "Counts") # Assign nicer names for plot
        myplot <- ggplot(x, aes(x=Feature, y=Counts)) + 
                        geom_bar(aes(fill=Strand), position="stack", stat="identity") + 
                        facet_wrap(~SampleName) + 
                        scale_y_continuous(axis_label) +
                        theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.4)) +
                        coord_flip() +
            scale_x_discrete(limits = unique(x$Feature))
        get(graphicsformat)(graphicsfile)
            print(myplot)
        grDevices::dev.off()
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
                             ggtitle(featureCountslist[[i]][1,1]) +
                scale_x_discrete(limits = unique(featureCountslist[[i]]$Length))
            ## Assure same scale for all panels if scales="fixed"
            if(scales=="fixed") {
                myplot <- myplot + scale_y_continuous(axis_label, limits=c(0, mymax)) 
           } else {
                myplot <- myplot + scale_y_continuous(axis_label) 
           }
            myplotlist[[i]] <- myplot
        }
    
        ## Generate graphics and write to file in pdf, png or jpeg format.
        if(tolower(graphicsformat)=="pdf") mydim <- c(height=20, width=6*length(myplotlist))
        if(tolower(graphicsformat) %in% c("png", "jpeg")) mydim <- c(height=20*96, width=6*length(myplotlist)*96)
        get(graphicsformat)(graphicsfile, width=mydim["width"], height=mydim["height"])
            grid::grid.newpage() # Open a new page on grid device
            grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(myplotlist)))) # Assign to device viewport with 1 by 2 grid layout 
            for(i in seq(along=myplotlist)) print(myplotlist[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
        grDevices::dev.off()
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

##############################################
## Compute and plot coverage along features ## 
##############################################
## Computes coverage along single and multi component features, such as exons/cds by transcripts. 
featureCoverage <- function(bfl, grl, resizereads=NULL, readlengthrange=NULL, Nbins=20, method=mean, fixedmatrix, resizefeatures, upstream, downstream, outfile, overwrite=FALSE) {
    ## global functions or variables
    readGAlignments <- qwidth <- NULL
    ## Input validity checks
    if(!is.null(outfile)) {
        if(file.exists(outfile) & overwrite==FALSE) stop(paste("File", outfile, "exists. Delete it or set 'overwrite=TRUE'"))
        if(file.exists(outfile) & overwrite==TRUE) unlink(outfile)
    }
    if(class(bfl)!="BamFileList") stop("'bfl' needs to be of class 'BamFileList'")
    if(!is(grl, "GRangesList")) stop("'grl' needs to be a GRangesList object.")
    if(!is.null(readlengthrange)) {
        if(is.numeric(readlengthrange) & length(readlengthrange)!=2) stop("'readlengthrange' needs to be assigned NULL or numeric vector of length 2")
        if(!is.numeric(readlengthrange)) stop("'readlengthrange' needs to be assigned NULL or numeric vector of length 2")
        if((is.numeric(readlengthrange) & length(readlengthrange)==2) & readlengthrange[1] > readlengthrange[2]) stop("Second value in read length range needs to be equal or larger then first one.") 
    }  
    if(fixedmatrix==TRUE & resizefeatures==FALSE) stop("If 'fixedmatrix=TRUE' then 'resizefeatures' needs to be set to 'TRUE' too.")
    
    ## If Nbins is not NULL, then remove ranges with width < Nbins
    if(!is.null(Nbins)) grl <- grl[sum(width(reduce(grl))) >= Nbins]

    ## Resize features by upstream and downstream value
    if(resizefeatures==TRUE) {
        grl <- .resizeFeature(grl=grl, upstream=upstream, downstream=upstream, component_resort=TRUE)
    }

    ## Following codes operates mostly on GRanges object
    gr <- unlist(grl)
    
    ## Define empty buckets for collcecting results generated in master loop
    binDFresult <- data.frame()
    fixedMAresult <- data.frame()
    fixedbinMAresult <- data.frame()
    rle_result <- list()

    for(mysample in seq_along(bfl)) {
        ## Compute coverage for ranges in stored in gr from BAM file. Note, coverage() works 
        ## directly on BAM files but for resizing read lengths and strand specific analysis 
        ## one needs to import the corresponding alignment sections first.
        
        ## Import alignment chunk
        aligns <- readGAlignments(bfl[[mysample]], param=ScanBamParam(what=scanBamWhat(), which=gr))
        ## Remove duplicated mappings generated by readGAlignment for nearby/overlapping ranges; reduce() doesn't fix that
        aligns <- aligns[!duplicated(paste(as.character(seqnames(aligns)), start(aligns), end(aligns), as.character(strand(aligns)), mcols(aligns)$qname, sep="_"))]   
     
        ## Subset to specific read length range specified under readlengthrange
        if(is.numeric(readlengthrange)) {
                aligns <- aligns[qwidth(aligns) >= readlengthrange[1] & qwidth(aligns) <= readlengthrange[2]]
        } 
       
        ## Compute coverage 
        if(is.numeric(resizereads[1])) {
            cov <- coverage(resize(granges(aligns), resizereads[1]))
            cov_pos <- coverage(resize(granges(subsetByOverlaps(aligns, gr, ignore.strand=FALSE)), resizereads[1]))
        } else if(is.null(resizereads)) {
            cov <- coverage(aligns)
            cov_pos <- coverage(subsetByOverlaps(aligns, gr, ignore.strand=FALSE))
        } else {
            stop("'resizereads' needs to be assigned NULL or positive integer of length 1.")   
        }
        cov_neg <- cov - cov_pos    

        ## Extract coverage for gr components. 
        ## Sense coverage
        # cov_posreg <- suppressWarnings(Views(cov_pos, as(gr, "IntegerRangesList"))) # delete 
        cov_posreg <- suppressWarnings(Views(cov_pos[names(as(gr, "IntegerRangesList"))], as(gr, "IntegerRangesList"))) # Update 22-Nov-15: cov_pos needs to be subsetted by seqnames in IntegerRangesList. This is relevant if txdb/grl was created from gff with scaffolds not containing any genes
        cov_posreg <- cov_posreg[sapply(cov_posreg, length) > 0] # Removes empty components (chr) 
        ## Antisense coverage
        # cov_negreg <- suppressWarnings(Views(cov_neg, as(gr, "IntegerRangesList"))) # delete
        cov_negreg <- suppressWarnings(Views(cov_neg[names(as(gr, "IntegerRangesList"))], as(gr, "IntegerRangesList"))) # Update 22-Nov-15: cov_neg needs to be subsetted by seqnames in IntegerRangesList. This is relevant if txdb/grl was created from gff with scaffolds not containing any genes
        cov_negreg <- cov_negreg[sapply(cov_negreg, length) > 0] # Removes empty components (chr) 
        
        ## Collapse (splice) multicomponent coverage ranges, e.g. cds exons to full cds   
        ## Sense coverage
        mystrand <- strand(grl); mystrand <- sapply(names(mystrand), function(x) as.character(mystrand[[x]])[1])
        cov_posList <- as.list(mystrand)
        ## Features on "+" strand
        for(i in seq_along(cov_posreg)) {
            for(j in names(mystrand[mystrand=="+"])) {
                myview <- cov_posreg[[i]][names(cov_posreg[[i]]) %in% j]
                if(length(myview)!=0) {
                    cov_posList[[j]] <- Rle(as.numeric(na.omit(as.vector(t(as.matrix(myview))))))
                }
            }
        }
        ## Features on "-" strand: note coverage of component features needs to be reversed here prior to collapsing to single feature!
        for(i in seq_along(cov_posreg)) {
            for(j in names(mystrand[mystrand=="-"])) {
                myview <- cov_posreg[[i]][names(cov_posreg[[i]]) %in% j]
                if(length(myview)!=0) {
                    cov_posList[[j]] <- Rle(as.numeric(na.omit(as.vector(t(as.matrix(reverse(myview)))))))
                }
            }
        }
        ## Antisense coverage
        cov_negList <- as.list(mystrand)
        ## Features on "+" strand
        for(i in seq_along(cov_negreg)) {
            for(j in names(mystrand[mystrand=="+"])) {
                myview <- cov_negreg[[i]][names(cov_negreg[[i]]) %in% j]
                if(length(myview)!=0) {
                    cov_negList[[j]] <- Rle(as.numeric(na.omit(as.vector(t(as.matrix(myview))))))
                }
            }
        }
        ## Features on "-" strand: note coverage of component features needs to be reversed here prior to collapsing to single feature!
        for(i in seq_along(cov_negreg)) {
            for(j in names(mystrand[mystrand=="-"])) {
                myview <- cov_negreg[[i]][names(cov_negreg[[i]]) %in% j]
                if(length(myview)!=0) {
                    cov_negList[[j]] <- Rle(as.numeric(na.omit(as.vector(t(as.matrix(reverse(myview)))))))
                }
            }
        }
        
        ## Get total number of aligned reads from BAM file and add to output for normalization
        param <- ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE, isUnmappedQuery=FALSE))
        aligned_reads <- countBam(bfl[mysample], param=param)
        
        ## If Nbins is integer, compute relative bin coverage for fixed number of bins (intervals) of each feature
        if(!is.null(Nbins)) {
            ## Function to bin coverage
            .interComp <- function(x, slice, method) {
                interval <- cut(1:length(x), slice, labels=FALSE)
                resultvec <- tapply(x, as.factor(interval), method)
                return(resultvec)
            }
            ## Sense coverage
            myMApos <- t(sapply(cov_posList, function(y) .interComp(x=y, slice=Nbins, method=method)))
            ## Antisense coverage
            myMAneg <- t(sapply(cov_negList, function(y) .interComp(x=y, slice=Nbins, method=method)))
            ## Combine results
            myMAneg <- myMAneg[rownames(myMApos),, drop=FALSE] # Assure proper sorting
            myDFpos <- data.frame(SampleName=names(bfl[mysample]), N_total_aligned=aligned_reads$records, IDs=rownames(myMApos), Strand="Sense", myMApos, check.names=FALSE)
            myDFneg <- data.frame(SampleName=names(bfl[mysample]), N_total_aligned=aligned_reads$records, IDs=rownames(myMAneg), Strand="Antisense", myMAneg, check.names=FALSE)
            myDFbins <- rbind(myDFpos, myDFneg, make.row.names=FALSE)
            row.names(myDFbins) <- NULL                  
            binDFresult <- rbind(binDFresult, myDFbins) 
            exportDF <- myDFbins
        } 
        ## Generate fixed matrix output
        if(fixedmatrix==TRUE) {
            mycolnames <- c(-upstream:-1, 0, 1:downstream)
            fixedMAstartpos <- t(sapply(names(cov_posList), function(x) as.numeric(cov_posList[[x]])[1:(upstream+downstream+1)]))
            colnames(fixedMAstartpos) <- mycolnames 
            fixedMAendpos <- t(sapply(names(cov_posList), function(x) as.numeric(cov_posList[[x]])[(length(as.numeric(cov_posList[[x]]))-(upstream+downstream)):length(as.numeric(cov_posList[[x]]))]))
            colnames(fixedMAendpos) <- mycolnames 
            fixedMAstartneg <- t(sapply(names(cov_negList), function(x) as.numeric(cov_negList[[x]])[1:(upstream+downstream+1)]))
            colnames(fixedMAstartneg) <- mycolnames 
            fixedMAendneg <- t(sapply(names(cov_negList), function(x) as.numeric(cov_negList[[x]])[(length(as.numeric(cov_negList[[x]]))-(upstream+downstream)):length(as.numeric(cov_negList[[x]]))]))
            colnames(fixedMAendneg) <- mycolnames 
            fixedMA <- rbind(data.frame(SampleName=names(bfl[mysample]), N_total_aligned=aligned_reads$records, IDs=rownames(fixedMAstartpos), Strand="Sense", fixedMAstartpos, check.names=FALSE),
                             data.frame(SampleName=names(bfl[mysample]), N_total_aligned=aligned_reads$records, IDs=rownames(fixedMAstartneg), Strand="Antisense", fixedMAstartneg, check.names=FALSE))
            tmpMA <- rbind(data.frame(fixedMAendpos[rownames(fixedMAstartpos),, drop=FALSE], check.names=FALSE),
                             data.frame(fixedMAendneg[row.names(fixedMAstartneg),,drop=FALSE], check.names=FALSE))
            fixedMAout <- cbind(fixedMA, "<S|E>"="|", tmpMA)
            row.names(fixedMAout) <- NULL                  
            fixedMAresult <- rbind(fixedMAresult, fixedMAout)
            exportDF <- fixedMAout
        }
        ## Combine fixed and binned coverage matrix
        if(!is.null(Nbins) & fixedmatrix==TRUE) {
            fixedbinMAout <- data.frame(fixedMA, "<S|B>"="|", myDFbins[,-c(1:4), drop=FALSE], "<B|E>"="|", tmpMA, check.names=FALSE)
            row.names(fixedbinMAout) <- NULL 
            fixedbinMAresult <- rbind(fixedbinMAresult, fixedbinMAout)
            exportDF <- fixedbinMAout
        } 
        ## Organize rle results in list
        if(fixedmatrix==FALSE & is.null(Nbins)) {
            rle_out <- list(list(N_total_aligned=aligned_reads$records, sense=as(cov_posList, "RleList"), antisense=as(cov_negList, "RleList")))
            names(rle_out)  <- names(bfl[mysample])
            rle_result <- c(rle_result, rle_out)
        } 
        
        ## Write tabular data to file in append mode
        if(!is.null(outfile) & (fixedmatrix==TRUE | !is.null(Nbins))) {
            if(!file.exists(outfile)) {    
                write.table(exportDF, outfile, append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
                cat("Processed sample:", names(bfl[mysample]), "and wrote results to file", outfile, "\n")
            } else {
                write.table(exportDF, outfile, append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
                cat("Processed sample:", names(bfl[mysample]), "and appended results to file", outfile, "\n")
            }
        } else {
            cat("Processed sample:", names(bfl[mysample]), "\n")
        }
    }
    
    ## Return proper final result 
    if(!is.null(Nbins) & fixedmatrix==FALSE) {
        return(binDFresult)
    } else if(is.null(Nbins) & fixedmatrix==TRUE) {
        return(fixedMAresult)
    } else if(!is.null(Nbins) & fixedmatrix==TRUE) {
        return(fixedbinMAresult)
    } else {
        return(rle_result)    
    }
}
## Usage:
# grl <- cdsBy(txdb, "tx", use.names=TRUE)
# fcov <- featureCoverage(bfl=BamFileList(outpaths[1]), grl=grl[1:4], resizereads=NULL, readlengthrange=NULL, Nbins=NULL, method=mean, fixedmatrix=TRUE, resizefeatures=TRUE, upstream=20, downstream=20, outfile="results/zzz.xls", overwrite=TRUE)

## Helper function to extend single and multi component ranges
## Extends featureBy GRangeList objects containing component features 
## such as exons in CDSs or transcripts so that only the first and last components 
## get extended. Single component features will be extended the same way.
.resizeFeature <- function(grl, upstream, downstream, component_resort=TRUE) { 
    ## global functions or variables
    end <- start <- mcols <- NULL
    if(!is(grl, "GRangesList")) stop("'grl' needs to be a GRangesList object.")
    gr <- unlist(grl)    
    if(!all(names(grl) %in% unique(names(gr)))) stop("None or not all components in grl are named.")
    
    ## Add sort_index
    mcols(gr) <- DataFrame(mcols(gr), sort_index=seq_along(gr))

    ## Ranges on + strand
    gr_pos <- gr[strand(gr)=="+" | strand(gr)=="*"] # Note: treats unstranded * the same as pos +
    if(component_resort==TRUE) {
        gr_pos <- gr_pos[order(as.character(seqnames(gr_pos)), names(gr_pos), start(gr_pos))] # Assures proper sorting
    }
    index_pos_first <- !duplicated(names(gr_pos), fromLast=FALSE)
    suppressWarnings(start(gr_pos[index_pos_first]) <- start(gr_pos[index_pos_first]) - upstream)
    gr_pos <- trim(gr_pos)
    index_pos_last <- !duplicated(names(gr_pos), fromLast=TRUE)
    suppressWarnings(end(gr_pos[index_pos_last]) <- end(gr_pos[index_pos_last]) + downstream) 
    gr_pos <- trim(gr_pos)

    ## Ranges on - strand
    gr_neg <- gr[strand(gr)=="-"]
    if(component_resort==TRUE) {
        gr_neg <- gr_neg[order(as.character(seqnames(gr_neg)), names(gr_neg), -start(gr_neg))]
    } 
    index_neg_first <- !duplicated(names(gr_neg), fromLast=FALSE)
    suppressWarnings(end(gr_neg[index_neg_first]) <- end(gr_neg[index_neg_first]) + upstream)
    gr_neg <- trim(gr_neg)
    index_neg_last <- !duplicated(names(gr_neg), fromLast=TRUE)
    suppressWarnings(start(gr_neg[index_neg_last]) <- start(gr_neg[index_neg_last]) - downstream)
    gr_neg <- trim(gr_neg)
   
    ## Return output in same format as input 
    gr_mod <- c(gr_pos, gr_neg)
    gr_mod <- gr_mod[order(mcols(gr_mod)$sort_index)]
    gr_mod <- gr_mod[, seq_along(colnames(mcols(grl[[1]])))]
    myfactor <- names(gr_mod)
    if(is.null(names(grl[[1]]))) names(gr_mod) <- NULL
    grl_mod <- split(gr_mod, myfactor)
    grl_mod <- grl_mod[names(grl)]
    return(grl_mod)
}
## Usage:
# grl_mod <- .resizeFeature(grl=grl, upstream=20, downstream=20, component_resort=TRUE)

###########################
## Plot feature coverage ##
###########################
## Plots tabular coverage data generated by featureCoverage()
plotfeatureCoverage <- function(covMA, method=mean, scales="fixed", extendylim=2, scale_count_val=10^6) {
    ## Some input validity checks
    if(class(covMA) != "data.frame") stop("'covMA' needs to be assigned an object of class 'data.frame'.")
    expectedcol <- c("SampleName", "N_total_aligned", "IDs", "Strand")
    if(any(!colnames(covMA)[1:4] %in% expectedcol)) stop(paste("The first 4 columns in 'covMA' need to be named:", paste(expectedcol, collapse=", ")))
    ## global functions or variables
    Coverage <- Strand <- NULL
    ## Determine split type required for provided input 
    splitloc <- which(grepl("(<S\\|B>)|(<B\\|E>)|(<S\\|E>)", colnames(covMA)))
    
    ## Get numbers of total aligned reads per sample
    N_total_aligned <- covMA$N_total_aligned; names(N_total_aligned) <- covMA$SampleName
    N_total_aligned <- N_total_aligned[!duplicated(names(N_total_aligned))]    
    
    ## Define function to convert tabular input into ggplot2 friendly format
    .convertDFlist <- function(malist) {
        for(j in seq_along(malist)[-1]) {
            df <- malist[[j]]
            tmpDF <- data.frame() 
            for(i in seq_along(colnames(df))[-c(1,2)]) {
                tmpDF <- rbind(tmpDF, data.frame(df[,1:2], Position=colnames(df)[i], Coverage=df[,i]))
            }
            malist[[j]] <- tmpDF
        }
        return(malist)
    }
   
    ## Construct informative title: print ID if only one or their count if more than 1 
    ti <- as.character(unique(covMA[,"IDs"]))
    if(length(ti) > 1) ti <- paste(length(ti), "features")
     
    ## Split and compute summary stats in format expected by ggplot2
    ## For binned matrix 
    if(length(splitloc)==0) {
        malist <- sapply(c("title", "bin"), function(x) NULL)
        malist[["title"]] <- c(titlebin=paste("CDS of", ti), x="Bins", y="Coverage")
        bin <- covMA
        malist[["bin"]] <- aggregate(bin[,-c(1:4)], by=list(SampleName=bin$SampleName, Strand=bin$Strand), FUN=method)
        malist <- .convertDFlist(malist)
    } 
    ## For fixed matrix around start and stop codons 
    if(length(splitloc)==1) {
        malist <- sapply(c("title", "start", "stop"), function(x) NULL)
        malist[["title"]] <- c(titlestart=paste("Start of", ti), titlestop=paste("Stop of", ti), x="Bins", y="Coverage")
        start <- covMA[ , 1:(splitloc[1] - 1)]
        malist[["start"]] <- aggregate(start[,-c(1:4)], by=list(SampleName=start$SampleName, Strand=start$Strand), FUN=method)
        stop <- covMA[ , c(1:4, (splitloc[1] + 1):length(colnames(covMA)))]  
        malist[["stop"]] <- aggregate(stop[,-c(1:4)], by=list(SampleName=stop$SampleName, Strand=stop$Strand), FUN=method)
        malist <- .convertDFlist(malist)
    } 
    ## For fixed matrix and binned matrix 
    if(length(splitloc)==2) {
        malist <- sapply(c("title","start","bin","stop"), function(x) NULL)
        malist[["title"]] <- c(titlestart=paste("Start of", ti), titlebin=paste("CDS of", ti), titlestop=paste("Stop of", ti), x="Bins", y="Coverage")
        start <- covMA[ , 1:(splitloc[1] - 1)]
        malist[["start"]] <- aggregate(start[,-c(1:4)], by=list(SampleName=start$SampleName, Strand=start$Strand), FUN=method)
        bin <- covMA[ , c(1:4, (splitloc[1] + 1):(splitloc[2] - 1))]
        malist[["bin"]] <- aggregate(bin[,-c(1:4)], by=list(SampleName=bin$SampleName, Strand=bin$Strand), FUN=method)
        stop <- covMA[ , c(1:4, (splitloc[2] + 1):length(colnames(covMA)))]  
        malist[["stop"]] <- aggregate(stop[,-c(1:4)], by=list(SampleName=stop$SampleName, Strand=stop$Strand), FUN=method)
        malist <- .convertDFlist(malist)
    }
    
    ## Scale per x reads (e.g. per million reads)
    if(is.numeric(scale_count_val[1])) {
        for(i in seq_along(malist)[-1]) malist[[i]][,"Coverage"] <- malist[[i]][,"Coverage"] * (scale_count_val/N_total_aligned[malist[[i]][,"SampleName"]])
    } else if(length(scale_count_val) == 0) {
            malist <- malist
    } else {
            stop("'scale_count_val' needs to be assinged NULL or numeric value.")
    }
    
    ## Label for coverage axis
    if(length(scale_count_val)==0) axis_label <- "Raw Coverage"
    if(length(scale_count_val)>0) axis_label <- paste("Coverage normalized per", scale_count_val, "reads")
    
    ## Get maximum coverage in data sets after summing up strand spec counts
    mymax <- numeric()
    for(i in seq_along(malist)[-1]) {
        myfactor <- paste(malist[[i]][,1], malist[[i]][,2], malist[[i]][,3], sep="_")
        mymax <- max(c(mymax, max(tapply(malist[[i]][,"Coverage"], myfactor, sum))))
    }
    mymax <- mymax * extendylim
    
    ## Generate plotting instructions
    myplotlist <- sapply(names(malist[-1]), function(x) NULL, simplify=FALSE)
    for(i in seq_along(myplotlist)) {
        myplot <- ggplot(malist[[names(myplotlist)[i]]], aes(x=Position, y=Coverage)) + 
                         geom_bar(aes(fill=Strand), position="stack", stat="identity") + 
                         facet_wrap(~SampleName, ncol=1, scales=scales) +
                         theme(legend.position="bottom") + 
                         ggtitle(malist[["title"]][i]) +
            scale_x_discrete(limits = unique(malist[[names(myplotlist)[i]]]$Position))
        ## Assure same scale for all panels if scales="fixed"
        if(scales=="fixed") {
            myplot <- myplot + scale_y_continuous(axis_label, limits=c(0, mymax)) 
        } else {
            myplot <- myplot + scale_y_continuous(axis_label) 
        }
        myplotlist[[i]] <- myplot
    }
    
    ## Generate graphics 
    grid::grid.newpage() # Open a new page on grid device
    grid::pushViewport(viewport(layout = grid::grid.layout(1, length(myplotlist)))) # Assign to device viewport with 1 by 2 grid layout 
    for(i in seq(along=myplotlist)) print(myplotlist[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
}

## Usage:
# grl <- cdsBy(txdb, "tx", use.names=TRUE)
# fcov <- featureCoverage(bfl=BamFileList(outpaths[1:2]), grl=grl[1:4], resizereads=NULL, readlengthrange=NULL, Nbins=20, method=mean, fixedmatrix=TRUE, resizefeatures=TRUE, upstream=20, downstream=20, outfile="results/zzz.xls", overwrite=TRUE)
# plotfeatureCoverage(covMA=fcov, method=mean, scales="fixed", extendylim=2, scale_count_val=10^6)

##################
## Predict ORFs ##
##################
## Function to predict ORFs in DNA sequences provided as DNAString/DNAStringSet objects
predORF <- function(x, n=1, type="grl", mode="orf", strand="sense", longest_disjoint=FALSE, startcodon="ATG", stopcodon=c("TAA", "TAG", "TGA")) {
	## Check input validity 
    if(any(nchar(c(startcodon, stopcodon))!=3)) stop("startcodon and stopcodons can only contain 3-letter strings.")
    if(!toupper(mode) %in% c("ORF", "CDS")) stop("'mode' can only be assigned one of: 'orf' or 'cds'")
    if(length(names(x))==0 | any(duplicated(names(x)))) stop("Sequence name slot of x need be populated with unique names.")
    
    ## Remove sequence with less than 6 nucleotides in length
    x <- x[width(x)>=6]

    ## Function for predicting ORFs/CDSs on single sequence
    .predORF <- function(x, n, mode, strand, ...) {
        ## start/stop codon assignment
        if(tolower(strand)=="sense") {
		    mystrand <- 1 # for +
            startcodon <- startcodon
            stopcodon <- stopcodon
		} else if(tolower(strand)=="antisense") {
            ## Reverse and complement start/stop plus swap their assignment
		    mystrand <- 2 # for -
		    stopcodon_temp <- as.character(reverseComplement(DNAStringSet(stopcodon)))
            startcodon_temp <- as.character(reverseComplement(DNAStringSet(startcodon)))
		    stopcodon <- startcodon_temp 
            startcodon <- stopcodon_temp
        
        } else {
            stop("strand can only be assigned 'sense', 'antisense' or 'both'")
        }
        
        ## Sequences containing N are not processed
        if(alphabetFrequency(x)["N"] > 0) {
            orfRanges <- cbind(subject_id=numeric(), start=numeric(), end=numeric(), width=numeric(), strand=numeric(), inframe2end=numeric())
            warning("Skipped sequence containing Ns.")
            return(orfRanges)
        }

        ## Tripletize x for each frame 
        c1 <- as.character(suppressWarnings(codons(x)))
		c2 <- as.character(suppressWarnings(codons(x[2:length(x)])))
		c3 <- as.character(suppressWarnings(codons(x[3:length(x)])))
        
        ## Identify position of start/stop in tripletized x
		startpos1 <- which(c1 %in% startcodon)
		stoppos1 <- which(c1 %in% stopcodon)
		startpos2 <- which(c2 %in% startcodon)
		stoppos2 <- which(c2 %in% stopcodon)
		startpos3 <- which(c3 %in% startcodon)
		stoppos3 <- which(c3 %in% stopcodon)
		
		## Make sure subsequent code also finds coding sequence frames (CDS) rather than just strict open reading frames (ORFs) 
		if(mode=="cds") {
			stoppos1 <- unique(c(0, stoppos1, length(c1)+1)); stoppos2 <- unique(c(0, stoppos2, length(c2)+1)); stoppos3 <- unique(c(0, stoppos3, length(c3)+1))
			startpos1 <- stoppos1; startpos2 <- stoppos2; startpos3 <- stoppos3
		}
        
        ## Map tripletized matches back to sequence of x 
        if(tolower(strand)=="sense") {
		    orfpos1 <- t(sapply(seq(along=startpos1), function(x) c((startpos1[x] * 3) -2, stoppos1[stoppos1 > startpos1[x]][1] * 3)))
		    if(length(orfpos1)==0) orfpos1 <- matrix(nrow = 0, ncol = 2)
		    orfpos2 <- t(sapply(seq(along=startpos2), function(x) c((startpos2[x] * 3) -1, (stoppos2[stoppos2 > startpos2[x]][1] * 3) + 1)))
		    if(length(orfpos2)==0) orfpos2 <- matrix(nrow = 0, ncol = 2)
		    orfpos3 <- t(sapply(seq(along=startpos3), function(x) c((startpos3[x] * 3) -0, (stoppos3[stoppos3 > startpos3[x]][1] * 3) + 2)))
		    if(length(orfpos3)==0) orfpos3 <- matrix(nrow = 0, ncol = 2)
		}
        if(tolower(strand)=="antisense") { # Note - for antisense predictions startpos contain stop codons and stoppos contain start codons
		    startpos1 <- sort(startpos1, decreasing = TRUE); startpos2 <- sort(startpos2, decreasing = TRUE); startpos3 <- sort(startpos3, decreasing = TRUE)
		    stoppos1 <- sort(stoppos1, decreasing = TRUE); stoppos2 <- sort(stoppos2, decreasing = TRUE); stoppos3 <- sort(stoppos3, decreasing = TRUE)
		    orfpos1 <- t(sapply(seq(along=stoppos1), function(x) c((startpos1[startpos1 < stoppos1[x]][1] * 3) -2, stoppos1[x] * 3)))
		    if(length(orfpos1)==0) orfpos1 <- matrix(nrow = 0, ncol = 2)
		    orfpos2 <- t(sapply(seq(along=stoppos2), function(x) c((startpos2[startpos2 < stoppos2[x]][1] * 3) -1, (stoppos2[x] * 3) + 1)))
		    if(length(orfpos2)==0) orfpos2 <- matrix(nrow = 0, ncol = 2)
		    orfpos3 <- t(sapply(seq(along=stoppos3), function(x) c((startpos3[startpos3 < stoppos3[x]][1] * 3) -0, (stoppos3[x] * 3) + 2)))
		    if(length(orfpos3)==0) orfpos3 <- matrix(nrow = 0, ncol = 2)
		}
        orfRanges <- rbind(orfpos1, orfpos2, orfpos3)
		orfRanges <- na.omit(orfRanges)
		if(mode=="cds") { # Modifications required for mode=="cds"
			orfRanges[,1] <- orfRanges[,1] + 3
			orfRanges[orfRanges[,2] > length(x), 2] <- length(x)
		}
		orfRanges <- IRanges::IRanges(start=orfRanges[,1], end=orfRanges[,2])
		orfRanges <- orfRanges[rev(order(width(orfRanges)))]
		
        ## Organize results in data.frame and also add info about frame of predicted ORF to downstream ORF e.g. prediction is uORF of 5'-UTR
	    orfRanges <- as.data.frame(orfRanges)
	 	inframe <- (length(x) - orfRanges$end) / 3; inframe2 <- inframe
		inframe2[abs((inframe - as.integer(inframe))) == 0] <- 1
		inframe2[abs(round((inframe - as.integer(inframe)),2)) == round(1/3, 2)] <- 2
		inframe2[abs(round((inframe - as.integer(inframe)),2)) == round(2/3, 2)] <- 3
		if(nrow(orfRanges)>0) {
            orfRanges <- cbind(subject_id=1:length(orfRanges[,1]), orfRanges, strand=mystrand, inframe2end=inframe2)
        } else {
            orfRanges <- cbind(subject_id=numeric(), start=numeric(), end=numeric(), width=numeric(), strand=numeric(), inframe2end=numeric())
        }
        ## Return all ORFs 
        if(n=="all") {
			## Subset to non-overlapping ORF set containing longest ORF
            if(nrow(orfRanges) > 0 & longest_disjoint==TRUE) {
                tmpgr <- GRanges(seqnames="dummy", IRanges::IRanges(orfRanges[,2], orfRanges[,3]), strand="+")
                orfRanges <- orfRanges[disjointBins(tmpgr)==1, , drop=FALSE]
			    orfRanges[,1] <- 1:nrow(orfRanges)
            }
            return(orfRanges)
        ## Return only as many ORFs as specified under n sorted decreasingly by length
		} else if(is.numeric(n)) {
            ## Make sure subsetting does not exceed number of records in orfRanges
            if(nrow(orfRanges) < n & nrow(orfRanges) != 0) { 
                upperlimit <- length(orfRanges[,1]) 
            } else if(nrow(orfRanges) > 0) { 
                upperlimit <- n
            } else {
                upperlimit <- NULL
            }
            if(is.null(upperlimit)) {
                return(orfRanges[upperlimit, , drop=FALSE])
		    } else {
                return(orfRanges[1:upperlimit, , drop=FALSE])
            }
        } else {
            stop("n needs to be assigned positive integer or 'all'")
        }
	}
	
    ## Run .predORF
    if(class(x)=="DNAString") {
        if(tolower(strand) == "sense" | tolower(strand) == "antisense") {
		    tmp <- .predORF(x, n, mode, strand, longest_disjoint, startcodon, stopcodon)
		    tmp[,"strand"] <- ifelse(as.numeric(tmp[,"strand"])==1, "+", "-")
        } else if(tolower(strand) == "both") {
		    tmp_pos <- .predORF(x, n, mode, strand="sense", longest_disjoint, startcodon, stopcodon)
		    tmp_neg <- .predORF(x, n, mode, strand="antisense", longest_disjoint, startcodon, stopcodon)
            tmp <- rbind(tmp_pos, tmp_neg)
		    tmp[,"strand"] <- ifelse(as.numeric(tmp[,"strand"])==1, "+", "-")
        } else {
            stop("strand can only be assigned 'sense', 'antisense' or 'both'")
        }
	    if(type=="df") return(tmp)
	    if(type=="gr") return(makeGRangesFromDataFrame(data.frame(seqnames="unknown", tmp), keep.extra.columns=TRUE))
    }
	if(class(x)=="DNAStringSet") {
        if(tolower(strand) == "sense" | tolower(strand) == "antisense") {
		    tmp <- lapply(x, function(y) .predORF(y, n=n, mode=mode, strand, longest_disjoint, startcodon, stopcodon))
		    names(tmp) <- names(x)
	        tmpdf <- do.call("rbind", tmp)
            rownames(tmpdf) <- NULL
            tmpdf <- data.frame(seqnames=rep(names(tmp), sapply(tmp, nrow)), tmpdf)
            tmpdf[,"strand"] <- ifelse(as.numeric(tmpdf[,"strand"])==1, "+", "-")
        } else if(tolower(strand) == "both") {
		    ## for sense strand
            tmp_pos <- lapply(x, function(y) .predORF(y, n=n, mode=mode, strand="sense", longest_disjoint, startcodon, stopcodon))
		    names(tmp_pos) <- names(x)
	        tmpdf_pos <- do.call("rbind", tmp_pos)
            rownames(tmpdf_pos) <- NULL
		    ## for antisense strand
            tmpdf_pos <- data.frame(seqnames=rep(names(tmp_pos), sapply(tmp_pos, nrow)), tmpdf_pos)
		    tmp_neg <- lapply(x, function(y) .predORF(y, n=n, mode=mode, strand="antisense", longest_disjoint, startcodon, stopcodon))
		    names(tmp_neg) <- names(x)
	        tmpdf_neg <- do.call("rbind", tmp_neg)
            rownames(tmpdf_neg) <- NULL
            tmpdf_neg <- data.frame(seqnames=rep(names(tmp_neg), sapply(tmp_neg, nrow)), tmpdf_neg)
            ## combine and sort results
            tmpdf <- rbind(tmpdf_pos, tmpdf_neg)
            tmpdf[,"strand"] <- ifelse(as.numeric(tmpdf[,"strand"])==1, "+", "-")
            tmpdf <- tmpdf[order(tmpdf$seqnames, tmpdf$subject_id), ]
        } else {
            stop("strand can only be assigned 'sense', 'antisense' or 'both'")
        }
	    
        ## Return results in format defined by type
        if(tolower(type)=="df") {
                rownames(tmpdf) <- NULL
                return(tmpdf)
        } else if(tolower(type)=="gr") {
                tmpdf <- makeGRangesFromDataFrame(tmpdf, keep.extra.columns=TRUE)
		        return(tmpdf)
        } else if(tolower(type)=="grl") {
                tmpdf <- makeGRangesFromDataFrame(tmpdf, keep.extra.columns=TRUE)
		        tmpdf <- split(tmpdf, as.character(seqnames(tmpdf)))
                return(tmpdf)
        } else {
            stop("type can only be assigned 'df' of 'gr'")
        }
	}
}
# library(Biostrings)
# file <- system.file("extdata", "someORF.fa", package="Biostrings")
# dna <- readDNAStringSet(file)
# orf <- predORF(dna[1:4], n=1, type="df", mode="orf", strand="antisense", startcodon="ATG", stopcodon=c("TAA", "TAG", "TGA"))

######################################################
## Scale feature level mappings to genome positions ##
######################################################
## Function to scale ranges within genomic features (query)
## to corresponding genome positions (subject). The function
## also accounts for intron insertions of compound features 
## (e.g. between exons of transcribed regions) that are absent in 
## the query ranges, but present in the corresponding subject ranges.
scaleRanges <- function(subject, query, type="custom", verbose=TRUE) {
    ## global functions or variables
    mcols <- NULL
    ## Both input objects need to be of class GRangesList
    if(!is(subject, "GRangesList") | !is(query, "GRangesList")) stop("Both subject and query need to be GRangesList objects.")
    ## All names(query) need to be present in names(subject)
    if(any(!names(query) %in% names(subject))) stop("All 'names(query)' need to be present in 'names(subject)'.")

    ## Perform scaling on single subject/query pair each containing on entry
    .scaleRanges <- function(subject, query, returntype="df") {
        ## global functions or variables
        mcols <- NULL
        ## Check for validity of query
        if(length(query)>1) warning("Only the first range in 'query' will be used.")
        query <- query[1]
        if(sum(width(subject)) < width(query)) stop("Sum of width of subject ranges cannot be smaller than width of query range.")
    
        ## Check for validity of subject
        subjectstrand <- unique(as.character(strand(subject)))
        if(length(subjectstrand) != 1) stop("More than one orientation detected. There can only be one.")
        querystrand <- unique(as.character(strand(query)))
    
        ## Check for validity of seqnames
        myseqname <- unique(as.character(seqnames(subject)))
        if(length(myseqname) != 1) stop("More than one seqname detected. There can only be one.")
    
        ## Scale query range to range of subject ranges (of genomic feature) using interval trees from IRanges 
        subject <- subject[order(start(subject))]
        rangev <- paste("c(", paste(paste(start(subject), end(subject), sep=":"), collapse=", "),")", sep="")
        rangev <- eval(parse(text=rangev)) 
        names(rangev) <- seq_along(rangev)
        if(subjectstrand=="+") rangev <- rangev[start(query):end(query)]
        if(subjectstrand=="-") rangev <- rev(rangev)[start(query):end(query)]
        ir <- reduce(IRanges::IRanges(rangev, rangev))
    
        ## Set orientation properly
        if(querystrand=="-" & subjectstrand=="-") mystrand <- "+"
        if(querystrand=="+" & subjectstrand=="+") mystrand <- "+"
        if(querystrand=="+" & subjectstrand=="-") mystrand <- "-"
        if(querystrand=="-" & subjectstrand=="+") mystrand <- "-"
    
        ## Organize result as GRanges 
        if(mystrand=="+") { # Proper exon ranking for exons on +/- strand
            gr <- GRanges(myseqname, ir[order(start(ir), decreasing=FALSE)], mystrand)
        } else {
            gr <- GRanges(myseqname, ir[order(start(ir), decreasing=TRUE)], mystrand)  
        }
        mcols(gr) <- data.frame(type=type)
        
        ## Return results as data.frame or GRanges
        if(returntype=="gr") return(gr)
        if(returntype=="df") return(as.data.frame(gr))
    }

    ## Run .scaleRanges on two GRangesLists as input
    subject <- subject[names(query)] # Subset subject to entries in query
    mygrl <- unlist(query)
    myids1 <- rep(names(query), sapply(start(query), length))
    myids2 <- paste(myids1, ":", start(mygrl), "_", end(mygrl), sep="")
    mcols(mygrl) <- data.frame(type=type)
    mygrl <- split(mygrl, seq_along(mygrl))
    names(mygrl) <- myids2
    mylist <- as.list(character(length=length(mygrl))); names(mylist) <- names(mygrl) 
    for(i in seq_along(mygrl)) {
        myentryid <- myids1[i]
        suppressWarnings(mylist[[i]] <- .scaleRanges(subject[[myentryid]], mygrl[[i]]))
        if(verbose==TRUE) {
            cat(paste0("Scaled range ", i, " of ", length(mygrl), " ranges: ", names(mygrl[i])), "\n")
        }
    }
    ## If loop gets interrupted, output only completed results!
    mylist <- mylist[1:i] 
    gr <- makeGRangesFromDataFrame(do.call("rbind", mylist), keep.extra.columns=FALSE)
    mcols(gr) <- DataFrame(feature_by=as(rep(names(mylist), sapply(mylist, nrow)), "CharacterList"), featuretype_id=names(gr), featuretype=type)
    names(gr) <- NULL
    grl <- split(gr, as.character(mcols(gr)$feature_by))
    return(grl)
}
## Usage for simple example
# subject <- GRanges(seqnames="Chr1", IRanges(c(5,15,30),c(10,25,40)), strand="+")
# query <- GRanges(seqnames="myseq", IRanges(1, 9), strand="+")
# scaleRanges(GRangesList(myid1=subject), GRangesList(myid1=query), type="test")

## Usage for more complex example
# library(GenomicFeatures)
# gff <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
# txdb <- makeTxDbFromGFF(file=gff, format="gff3", organism="Arabidopsis")
# futr <- fiveUTRsByTranscript(txdb, use.names=TRUE)
# genome <- system.file("extdata/annotation", "tair10.fasta", package="systemPipeRdata")
# dna <- extractTranscriptSeqs(FaFile(genome), futr)
# uorf <- predORF(dna, n="all", mode="orf", longest_disjoint=TRUE, strand="sense")
# grl_scaled <- scaleRanges(subject=futr, query=uorf, type="uORF", verbose=TRUE)
# rtracklayer::export.gff3(unlist(grl_scaled), "uorf.gff")

#######################################
## Translational Efficiency Analysis ##
#######################################
## Currently, this step is covered entire by the code given in the Ribo-Seq
## workflow vignette using DESeq2 out of the box. The approach ended up being
## very similar to some of time course analysis steps given in the DESeq2
## vignette as well as some community discussions:
##      https://support.bioconductor.org/p/61509/ also here
##      https://support.bioconductor.org/p/67455/
## To-dos: add convenience function similar to run_DESeq2 (or extend it)
## to handle this analysis step for many experiment sets in a more convenient manner.

