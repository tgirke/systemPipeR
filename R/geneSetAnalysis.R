#######################
## Gene Set Analysis ##
#######################

############################################
## Class and Method Definitions for catDB ##
############################################
## Define catDB class
setClass("catDB", representation(catmap="list",	catlist="list", idconv="ANY"))

## Methods to return catDB components as lists 
setGeneric(name="catmap", def=function(x) standardGeneric("catmap"))
setMethod(f="catmap", signature="catDB", definition=function(x) {return(x@catmap)})
setGeneric(name="catlist", def=function(x) standardGeneric("catlist"))
setMethod(f="catlist", signature="catDB", definition=function(x) {return(x@catlist)})
setGeneric(name="idconv", def=function(x) standardGeneric("idconv"))
setMethod(f="idconv", signature="catDB", definition=function(x) {return(x@idconv)})

## Constructor methods
## List to catDB with: as(mylist, "catDB")
setAs(from="list", to="catDB",  
        def=function(from) {
		new("catDB", catmap=from$catmap,
		              catlist=from$catlist,
		              idconv=from$idconv)
})

## Define print behavior for catDB
setMethod(f="show", signature="catDB", 
	definition=function(object) {    
	cat("An instance of '", class(object), "' containing:", "\n\t", 
             "catmap: ", length(object@catmap), " mapping data frames (", paste(names(object@catmap), collapse=", "), ") ", "\n\t",
	     "catlist: ", length(object@catlist), " category lists (", paste(names(object@catlist), collapse=", "), ")", "\n\t",
	     "idconv: ", length(object@idconv), " id conversion data frames (", paste(names(object@idconv), collapse=", "), ")", 
	     "\n", sep="")
})

## Extend names() method
setMethod(f="names", signature="catDB",
    definition=function(x) {
    	return(slotNames(x))
})
	
############################
## Construct catDB object ##
############################
## (A.1) Generate sample data frames with assigned gene-to-GO mappings,
## one for MF, one for BP and one for CC mappings
## custom mappings can be used here, but need to have the same format as GO_XX_DFs in the following examples
## (A.1.1) Obtain mappings from geneontology.org
.readGOorg <- function(myfile, colno, org) {
	go_org <- read.delim(myfile, na.strings = "", header=FALSE, comment.char = "!", sep="\t")
        go_org <- go_org[ , colno]
        names(go_org) <- c("GOID", "GeneID", "GOCAT")
	if(org == "Arabidopsis") {
		go_org[,"GeneID"] <- gsub(".*(AT.G\\d\\d\\d\\d\\d).*", "\\1", as.character(go_org[,2]), perl=TRUE)
        	go_org <- go_org[grep("^AT.G\\d\\d\\d\\d\\d", as.character(go_org$GeneID), perl=TRUE),]
        	go_org <- go_org[!duplicated(paste(go_org[,"GOID"], gsub("\\.\\d{1,}", "", as.character(go_org[,"GeneID"]), perl=TRUE), sep="_")),]
        }
        go_org <- na.omit(go_org)
	## Make GO cat labels consistent, e.g. from BioMart
	gocat <- as.character(go_org$GOCAT)
	gocat[gocat=="molecular_function"] <- "F"
	gocat[gocat=="biological_process"] <- "P"
	gocat[gocat=="cellular_component"] <- "C"
	go_org[,"GOCAT"] <- gocat
	## Split data frame by category
        GO_MF_DF <- go_org[go_org[,3]=="F",]
        GO_BP_DF <- go_org[go_org[,3]=="P",]
        GO_CC_DF <- go_org[go_org[,3]=="C",]
        
	## Generates data frame (go_df) containing the commonly used components for all GO nodes: GOID, GO Term and Ontology Type. This step is only required if "go_df" hasn't been imported with the above load() function.
	go_df <- data.frame(GOID=names(Term(GOTERM)), Term=Term(GOTERM), Ont=Ontology(GOTERM))
        go_df <- na.omit(go_df) 
	
	## Return results
	dflist <- list(D_BP=GO_BP_DF, D_CC=GO_CC_DF, D_MF=GO_MF_DF, GO_DF=go_df)
	return(dflist)
}
## Usage:
# catdf <- .readGOorg(myfile = "data/GO/GOannotationsBiomart_mod.txt", org="", colno = c(1,2,3))

## (A.1.2) Obtain mappings from BioC
.sampleDFgene2GO <- function(lib) {
        require(lib, character.only=TRUE)
	mylibbase <- gsub(".db", "", lib) 
        GOMF <- eapply(get(paste(mylibbase, "GO", sep="")), getOntology, "MF") # generates list with GeneID components containing MFGOs
        GO_MF_DF <- data.frame(GOID=unlist(GOMF), GeneID=rep(names(GOMF), as.vector(sapply(GOMF, length))), Count=rep(as.vector(sapply(GOMF, length)), as.vector(sapply(GOMF, length))))
        GOBP <- eapply(get(paste(mylibbase, "GO", sep="")), getOntology, "BP") # generates list with GeneID components containing BPGOs
        GO_BP_DF <- data.frame(GOID=unlist(GOBP), GeneID=rep(names(GOBP), as.vector(sapply(GOBP, length))), Count=rep(as.vector(sapply(GOBP, length)), as.vector(sapply(GOBP, length))))
        GOCC <- eapply(get(paste(mylibbase, "GO", sep="")), getOntology, "CC") # generates list with GeneID components containing CCGOs
        GO_CC_DF <- data.frame(GOID=unlist(GOCC), GeneID=rep(names(GOCC), as.vector(sapply(GOCC, length))), Count=rep(as.vector(sapply(GOCC, length)), as.vector(sapply(GOCC, length))))
	## Generates data frame (go_df) containing the commonly used components for all GO nodes: GOID, GO Term and Ontology Type. This step is only required if "go_df" hasn't been imported with the above load() function.
	go_df <- data.frame(GOID=names(Term(GOTERM)), Term=Term(GOTERM), Ont=Ontology(GOTERM))
        go_df <- na.omit(go_df) 
	## Return results
	dflist <- list(D_BP=GO_BP_DF, D_CC=GO_CC_DF, D_MF=GO_MF_DF, GO_DF=go_df)
	return(dflist)
}
## Usage:
# catdf <- .sampleDFgene2GO(lib="ath1121501.db")

## (A.2) Generate list containing gene-to-GO-OFFSPRING associations including assiged nodes
## This is slow (3 minutes), but needs to be done only once! 
.gene2GOlist <- function(catdf, rootUK=FALSE) { # If the argument 'rootUK' is set to TRUE then the root nodes are treated as terminal nodes to account for the new unknown terms 
        ## Import required data frames
	GO_MF_DF <- catdf$D_MF
        GO_BP_DF <- catdf$D_BP
        GO_CC_DF <- catdf$D_CC
       
	## Populate each GO node with associated gene ids 
	for(i in c("MF","BP","CC")) {
                if(i=="MF") {
                        go_offspr_list <- as.list(GOMFOFFSPRING) }
                if(i=="BP") {
                        go_offspr_list <- as.list(GOBPOFFSPRING) }
                if(i=="CC") {
                        go_offspr_list <- as.list(GOCCOFFSPRING) }
                go_offspr_list <- lapply(go_offspr_list, unlist); go_offspr_list <- lapply(go_offspr_list, as.vector) # clean-up step for the list
                go_offspr_list_temp <- lapply(names(go_offspr_list), function(x) c(x, go_offspr_list[[x]]) ) # include list component (GOID) names in corresponding (GOID) vectors
                names(go_offspr_list_temp) <- names(go_offspr_list) # names list components after go_offspr_list
                go_offspr_list <- go_offspr_list_temp
                go_offspr_list <- lapply(go_offspr_list, function(x) x[!is.na(x)]) # remove NAs in vectors
                
		## Treat root nodes as terminal nodes to account for the new unknown terms. This step removes the offspring information from the root nodes.
		if(rootUK==TRUE) {
			if(i=="MF") { go_offspr_list[["GO:0003674"]] <- c("GO:0003674") } 
			if(i=="BP") { go_offspr_list[["GO:0008150"]] <- c("GO:0008150") }
			if(i=="CC") { go_offspr_list[["GO:0005575"]] <- c("GO:0005575") }
		}

		## Retrieve gene (affy) IDs for GOID vectors
		if(i=="MF") {
                        MF_node_gene_list <- lapply(go_offspr_list, function(x) unique(as.vector(GO_MF_DF[GO_MF_DF$GOID %in% x, 2]))) 
                }
		if(i=="BP") {
                        BP_node_gene_list <- lapply(go_offspr_list, function(x) unique(as.vector(GO_BP_DF[GO_BP_DF$GOID %in% x, 2]))) 
                }
		if(i=="CC") {
                        CC_node_gene_list <- lapply(go_offspr_list, function(x) unique(as.vector(GO_CC_DF[GO_CC_DF$GOID %in% x, 2])))  
                }
		cat("\n", paste("'", i, "_node_gene_list'", sep=""), "with gene-to-GO-OFFSPRING associations created.", "\n") 
        }
	node_gene_list <-list(L_BP=BP_node_gene_list, L_CC=CC_node_gene_list, L_MF=MF_node_gene_list)
	return(node_gene_list)
}
## Usage:
# catlist <- .gene2GOlist(catdf=catdf, rootUK=FALSE)

## (A.3) Generate AffyID-to-GeneID mappings when working with chip feature IDs
## This function creates a AffyID-to-GeneID mapping data frame using by default the TAIR mappings for the Arabidopsis ATH1 chip. 
## Once the decoding data frame 'affy2locusDF' is created, the function returns for a query set of AffyIDs the corresponding GeneIDs.
## To use the function for the mappings of other chips, one needs to create the corresponding decoding data frame 'affy2locusDF'.
.AffyID2GeneID <- function(map = "ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2008-5-29.txt", download=FALSE, catdb=NULL, affyIDs, probe2gene=1) {
        if(download==TRUE) {
                cat("\n", "Downloading AffyID-to-GeneID mappings", "\n")
                affy2locus <- read.delim(map, na.strings = "", fill=TRUE, header=TRUE, sep="\t")[,-c(2:4,7:9)]
                names(affy2locus) <- c("AffyID", "AGI", "Desc")
                row.names(affy2locus) <- as.vector(affy2locus[,1])
                my_list <- apply(affy2locus[,-c(3)], 1, list); my_list <- lapply(my_list, unlist)
                my_list <- lapply(my_list, function(x) as.vector(x[-1]))
                my_list <- lapply(my_list, strsplit, ";"); my_list <- lapply(my_list, unlist)
                affy2locusDF <- data.frame(unlist(my_list))
                affy2locusDF <- data.frame(rep(names(unlist(lapply(my_list, length))), as.vector(unlist(lapply(my_list, length)))), affy2locusDF)
                names(affy2locusDF) <- c("AffyID", "GeneID")
                return(affy2locusDF)
        }
        if(!missing(affyIDs)) {
		affy2locusDF <- idconv(catdb)[[1]]
		if(probe2gene==1) { # For probe sets that match several loci, only the first locus ID will be used
			affy2locusDF <- affy2locusDF[!duplicated(affy2locusDF$AffyID),]
		}	
        	GeneIDs <- unique(as.vector(affy2locusDF[affy2locusDF[,1] %in% affyIDs, 2]))
        	return(GeneIDs)
        }
}
## Usage:
# affy2locusDF <- systemPipeR:::.AffyID2GeneID(map = "ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt", download=TRUE)
# catdb_conv <- makeCATdb(myfile="data/GO/GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=list(affy=affy2locusDF))
# systemPipeR:::.AffyID2GeneID(catdb=catdb_conv, affyIDs=c("244901_at", "244902_at"))

## Constructor function to generate catDB object
makeCATdb <- function(myfile, lib=NULL, org="", colno = c(1,2,3), idconv=NULL, rootUK=FALSE) {
	if(!is.null(lib) & !is.null(myfile)) stop("Arguments myfile and lib are exclusive. One of them needs to be assigned NULL.")
	if(!is.null(lib)) {
		catdf <- .sampleDFgene2GO(lib=lib)
	}
	if(!is.null(myfile)) {
		catdf <- .readGOorg(myfile=myfile, org=org, colno=colno)
	}
	catlist <- .gene2GOlist(catdf=catdf, rootUK=rootUK)
	catdblist <- c(catmap=list(catdf), catlist=list(catlist), idconv=list(idconv))
	catdb <- as(catdblist, "catDB")
	return(catdb)
}
## Usage:
# catdb <- makeCATdb(myfile="data/GO/GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
# catdb <- makeCATdb(myfile=NULL, lib="ath1121501.db", org="", colno=c(1,2,3), idconv=NULL)
# save(catdb, file="data/GO/catdb.RData") 
# load("data/GO/catdb.RData")

######################################################################################################
### GOHyperGAll: Global Hypergeometric Test Using Custom Gene-to-GO Mappings Plus GO Slim Analysis ###
######################################################################################################
## Utility: To test a sample population of genes for over-representation of GO terms, the 
## function 'GOHyperGAll' computes for all GO nodes a hypergeometric distribution test and 
## returns the corresponding raw and Bonferroni corrected p-values. A subsequent filter function 
## performs a GO Slim analysis using default or custom GO Slim categories. 
## The associated publication is available in Plant Physiol (2008) 147, 41-57.
## Note: GOHyperGAll provides similar utilities as the GOHyperG function in the GOstats package 
## from BioConductor. The main difference is that GOHyperGAll simplifies the usage of custom 
## chip-to-gene and gene-to-GO mappings.
##
## How it works:
## (A) Generate the required data objects (slow, but needs to be done only once)
## (B) Define GOhyperG_All function
## (C) Subsetting and plotting of results by assigned nodes or goSlim categories

############################
## (B) GOhyperG_All function
############################
## (B.1) Define GOhyperG_All function
GOHyperGAll <- function(catdb, gocat="MF", sample, Nannot=2) {
	## go_df containing the commonly used components for all GO nodes: GOID, GO Term and Ontology Type. This step is only required if "go_df" hasn't been imported with the above load() function.
        go_df <- catmap(catdb)$GO_DF
	
	## (m): Obtain for every node in GO tree their number of associated genes or chip features
        if(gocat=="MF") {node_list <- catlist(catdb)$L_MF}
        if(gocat=="BP") {node_list <- catlist(catdb)$L_BP}
        if(gocat=="CC") {node_list <- catlist(catdb)$L_CC}
        node_stats_df <- data.frame(NodeSize=sapply(node_list, length))
        node_stats_df <- data.frame(GOID=row.names(node_stats_df), node_stats_df)
        row.names(node_stats_df) <- 1:length(node_stats_df[,1])
        m <- as.vector(node_stats_df$NodeSize)        
                
	## (x): Obtain for every node in GO tree the number of matching genes in sample set
        node_sample_stats <- sapply(node_list, function(x) { sum(unlist(x) %in% sample) } )
        node_sample_stats <- as.vector(node_sample_stats)
        x <- node_sample_stats        

	## (n): Obtain the number of unique genes at GO nodes with direct annotations
        if(gocat=="MF") { GO_DF <- catmap(catdb)$D_MF }
        if(gocat=="BP") { GO_DF <- catmap(catdb)$D_BP }
        if(gocat=="CC") { GO_DF <- catmap(catdb)$D_CC }
        n <- length(unique(GO_DF[, 2]))

	## (k): Obtain number of unique genes in test sample that have GO mappings
        k <- length(unique(GO_DF[GO_DF[,2] %in% sample, 2]))

	## Obtain gene/chip keys matching at GO nodes
        match_key <- sapply(node_list, function(x) { x[unlist(x) %in% sample] } )
        match_key <- sapply(match_key, function(x) { paste(x, collapse=" ") } )
	match_key <- as.vector(match_key)
	key <- match_key; key[key==""] <- "NA"

	## Apply phyper function
        phyp_v <- phyper(x-1, m, n-m , k, lower.tail = FALSE)

	## P-value correction according to Bioinformatics, 20, 3710-3715
	Ncorrect <- table(GO_DF[GO_DF$GeneID %in% sample, 1]) # Obtain the GO nodes with direct annotations from sample set 
	Ncorrect <- sum(Ncorrect >= Nannot) # Count only those that have 2 or more annotations from sample set
	if(Ncorrect<=1) { 
		adj_phyp_v <- phyp_v # no adjustment necessary if Ncorrect <= 1
	} else {
		adj_phyp_v <- phyp_v * Ncorrect # Calculates simple Bonferroni correction. 
		adj_phyp_v[adj_phyp_v >= 1] <- 1
		# adj_phyp_v <- sapply(phyp_v, p.adjust, method=Padj, n = Ncorrect) # Runs p.adjust(). This is disabled because most adjustment methods require that the length of the p-value vector is >= n.  
	}
	
	## Generate output data format
	result_df <- data.frame(node_stats_df, SampleMatch=x, Phyper=phyp_v, Padj=adj_phyp_v, SampleKeys=key)
        result_df <- merge(result_df, go_df, x.by="GOID", y.by="GOID", all.x=TRUE)
        result_df <- result_df[order(result_df$Phyper), ]
	result_df <- result_df[,c(1:5,7:8,6)]
        return(result_df)
}
# test_sample <- unique(as.character(catmap(catdb)$D_MF[1:100,"GeneID"]))
# GOHyperGAll(catdb=catdb, gocat="MF", sample=test_sample, Nannot=2)[1:20,]

####################################################################################
## (C) Subsetting of results from GOHyperGAll by assigned nodes or goSlim categories
####################################################################################
## (C.1) Define subsetting function
GOHyperGAll_Subset <- function(catdb, GOHyperGAll_result, sample=test_sample, type="goSlim", myslimv) { # type: "goSlim" or "assigned"; optional argument "myslimv" to privde custom goSlim vector
        if(type=="goSlim") {
                if(missing(myslimv)) {
                        slimv <- c("GO:0003674", "GO:0008150", "GO:0005575", "GO:0030246","GO:0008289","GO:0003676","GO:0000166","GO:0019825","GO:0005515","GO:0003824","GO:0030234","GO:0003774","GO:0004871","GO:0005198","GO:0030528","GO:0045182","GO:0005215","GO:0006519","GO:0007154","GO:0016043","GO:0006412","GO:0006464","GO:0006810","GO:0007275","GO:0007049","GO:0005975","GO:0006629","GO:0006139","GO:0019748","GO:0015979","GO:0005618","GO:0005829","GO:0005783","GO:0005768","GO:0005794","GO:0005739","GO:0005777","GO:0009536","GO:0005840","GO:0005773","GO:0005764","GO:0005856","GO:0005634","GO:0005886","GO:0005576") # contains new unknown terms: "GO:0003674", "GO:0008150", "GO:0005575"
                        # slimv <- c("GO:0005554", "GO:0000004", "GO:0008372", "GO:0030246","GO:0008289","GO:0003676","GO:0000166","GO:0019825","GO:0005515","GO:0003824","GO:0030234","GO:0003774","GO:0004871","GO:0005198","GO:0030528","GO:0045182","GO:0005215","GO:0006519","GO:0007154","GO:0016043","GO:0006412","GO:0006464","GO:0006810","GO:0007275","GO:0007049","GO:0005975","GO:0006629","GO:0006139","GO:0019748","GO:0015979","GO:0005618","GO:0005829","GO:0005783","GO:0005768","GO:0005794","GO:0005739","GO:0005777","GO:0009536","GO:0005840","GO:0005773","GO:0005764","GO:0005856","GO:0005634","GO:0005886","GO:0005576") # contains old unknown terms: "GO:0005554", "GO:0000004", "GO:0008372" 
                        } else { 
                        slimv <- myslimv } 
                GOHyperGAll_subset <- GOHyperGAll_result[GOHyperGAll_result[,1] %in% slimv, ]
        }
        if(type=="assigned") {
        	## Import required data frames
		GO_MF_DF <- catmap(catdb)$D_MF
		GO_BP_DF <- catmap(catdb)$D_BP
		GO_CC_DF <- catmap(catdb)$D_CC
                termGO <- c(as.vector(GO_MF_DF[GO_MF_DF$GeneID %in% sample, 1]), 
                            as.vector(GO_BP_DF[GO_BP_DF$GeneID %in% sample, 1]), 
                            as.vector(GO_CC_DF[GO_CC_DF$GeneID %in% sample, 1]))
                subset_v <- unique(termGO)
                GOHyperGAll_subset <- GOHyperGAll_result[GOHyperGAll_result[,1] %in% subset_v, ]
        }
        GOHyperGAll_subset
}
## Usage:
# GOHyperGAll_result <- GOHyperGAll(catdb=catdb, gocat="MF", sample=test_sample, Nannot=2)
# GOHyperGAll_Subset(catdb, GOHyperGAll_result, sample=test_sample, type="goSlim") 

#########################################################
## (D) Reduce GO Term Redundancy in 'GOHyperGAll_results' 
#########################################################
## (D.1) The function 'GOHyperGAll_Simplify' subsets the data frame 'GOHyperGAll_result' by a user 
## specified adjusted p-value cutoff and removes from it all GO nodes with overlapping children sets 
## (OFFSPRING). Only the best scoring nodes remain in the data frame. 
## The argument 'correct' is experimental. It aims to favor the selection of distal (information rich) 
## GO terms that have at the same time a large number of sample matches. The following calculation is used 
## for this adjustment: phyper x Number_of_children / SampleMatch
## Define GOHyperGAll_Simplify()
GOHyperGAll_Simplify <- function(GOHyperGAll_result, gocat="MF", cutoff=0.001, correct=TRUE) { # gocat: "MF", "BP" or "CC"; cutoff: p-value cutoff; correct: TRUE or FALSE 
	if(gocat!=as.vector(GOHyperGAll_result$Ont[!is.na(GOHyperGAll_result$Ont)])[1]) { stop("The GO categories in GOHyperGAll_Simplify() and GOHyperGAll_result need to match") }
	testDF <- GOHyperGAll_result[GOHyperGAll_result$Padj<=cutoff,]
	testDF <- data.frame(testDF, test=rep(0, times=length(testDF[,1])))
	testDF <- testDF[!is.na(testDF$Ont),]
	GOIDv <- NULL
	GO_OL_Matchv <- NULL
	while(sum(testDF$test==0)>0) {
		clusterv <- NULL
		test <- as.vector(testDF[,1])
		for(j in 1:length(test)) {
			if(gocat=="MF") { mymatch <- sum(unique(na.omit(c(test[j], as.list(GOMFOFFSPRING)[[test[j]]])) %in% na.omit(c(test[1], as.list(GOMFOFFSPRING)[[test[1]]]))))
				if(mymatch==1) { mymatch <- length(as.list(GOMFOFFSPRING)[[test[j]]])  }
			}
			if(gocat=="BP") { mymatch <- sum(unique(na.omit(c(test[j], as.list(GOBPOFFSPRING)[[test[j]]])) %in% na.omit(c(test[1], as.list(GOBPOFFSPRING)[[test[1]]]))))
				if(mymatch==1) { mymatch <- length(as.list(GOBPOFFSPRING)[[test[j]]])  }
			}
			if(gocat=="CC") { mymatch <- sum(unique(na.omit(c(test[j], as.list(GOCCOFFSPRING)[[test[j]]])) %in% na.omit(c(test[1], as.list(GOCCOFFSPRING)[[test[1]]]))))
				if(mymatch==1) { mymatch <- length(as.list(GOCCOFFSPRING)[[test[j]]])  }
			}
			clusterv <- c(clusterv, mymatch)
		}
		clusterv[clusterv==0] <- NA
		testDF <- data.frame(testDF[,-9], test=clusterv)
		if(correct==TRUE) { 
			testDF <- data.frame(testDF, decide=testDF$Padj * (testDF$test/testDF$SampleMatch)) 
			} else {
			testDF <- data.frame(testDF, decide=testDF$Padj) }
		GOIDv <- c(GOIDv, as.vector(testDF[order(testDF[,10]),][1,1]))
		GO_OL_Matchv <- c(GO_OL_Matchv, length(unique(unlist(strsplit(as.vector(testDF[!is.na(testDF$test),8]), " ")))))
		testDF <- testDF[is.na(testDF$test),]
		testDF <- testDF[order(testDF[,5]),-c(9,10)]
		testDF <- data.frame(testDF, test=rep(0, times=length(testDF[,1])))
		cat(GOIDv, "\n")
	}
	simplifyDF <- data.frame(GOID=GOIDv, GO_OL_Match=GO_OL_Matchv)
	simplifyDF
}

## Apply GOHyperGAll_Simplify
## simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result, gocat="MF", cutoff=0.001, correct=TRUE)
## data.frame(GOHyperGAll_result[GOHyperGAll_result[,1] %in% simplifyDF[,1], -8], GO_OL_Match=simplifyDF[,2])

########################################
## (D.2) Batch Analysis of Gene Clusters
########################################
## The function 'GOCluster_Report' performs the three GO analyses in batch mode: 'GOHyperGAll', 
## 'GOHyperGAll_Subset' or 'GOHyperGAll_Simplify'. It processes many groups of genes (e.g. 
## gene expression clusters) and organizes the results in a single data frame.

GOCluster_Report <- function(catdb, setlist, id_type="affy", method="all", CLSZ=10, cutoff=0.001, gocats=c("MF", "BP", "CC"), myslimv="default", correct=TRUE, recordSpecGO=NULL, ...) { # CLSZ: minimum cluster size; method: "all", "slim" or "simplify"; gocat: "MF", "BP" or "CC"; cutoff: adjusted p-value cutoff; recordSpecGO: argument to include one specific GOID in each of the 3 ontologies, e.g: recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575")
    	CL_DF <- data.frame(geneID=unlist(setlist), CLID=rep(names(setlist), sapply(setlist, length)), ClusterSize=rep(sapply(setlist, length), sapply(setlist, length)))
	cluster_loop <- unique(as.vector(CL_DF[CL_DF[,3]>=CLSZ,2]))
        if(length(cluster_loop[grep("CL", cluster_loop)])>0) {
                cluster_loop <- paste("CL", sort(as.numeric(gsub("CL","", as.character(cluster_loop)))), sep="") 
        }
	if(method=="all") {
        	containerDF <- data.frame(CLID=NULL, CLSZ=NULL, GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
		for(i in cluster_loop) {
			cat("\n", "Processing cluster no", i, "with method: \"all\" (GOHyperGAll) \n")
			if(id_type=="affy") {
				affy_sample <- CL_DF[CL_DF[,2]==i, 1]
				test_sample <- .AffyID2GeneID(catdb=catdb, affyIDs=affy_sample, probe2gene=1)
			}
			if(id_type=="gene") {
				test_sample <- as.vector(CL_DF[CL_DF[,2]==i, 1])
			}	
			containerDF2 <- data.frame(GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
			count <- 0
			for(j in gocats) {
				count <- count+1
				GOHyperGAll_result <- GOHyperGAll(catdb=catdb, gocat=j, sample=test_sample, ...)
				tempDF <- GOHyperGAll_result[GOHyperGAll_result$Padj <= cutoff, ]
				if(length(tempDF[,1])==0) { # If filter returns empty data frame, then include at least the first two best scoring GO entries 
                                        tempDF <- GOHyperGAll_result[1:2,]
                                } 
				containerDF2 <- rbind(containerDF2, tempDF)
				if(length(recordSpecGO)>0) {
					containerDF2 <- rbind(containerDF2, GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],])
				}
				no_annot <- test_sample[!test_sample %in% eval(parse(text=paste("catmap(catdb)$D_", j, sep="")))[,2]]
				no_annot <- no_annot[no_annot!="no_match"]
				containerDF2 <- rbind(containerDF2, data.frame(GOID=paste("no_annot_", j, sep=""), NodeSize=NA, SampleMatch=length(no_annot), Phyper=NA, Padj=NA, Term=NA, Ont=j, SampleKeys=paste(no_annot, collapse=", ")))
			}
			tempDF2 <- data.frame(CLID=rep(i, times=length(containerDF2[,1])), CLSZ=rep(unique(as.vector(CL_DF[CL_DF[,2]==i,3])), times=length(containerDF2[,1])), containerDF2)
			containerDF <- rbind(containerDF, tempDF2)
		}
		return(containerDF)
	}	
	if(method=="slim") {
        	containerDF <- data.frame(CLID=NULL, CLSZ=NULL, GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
		for(i in cluster_loop) {
			cat("\n", "Processing cluster no", i, "with method: \"slim\" (GOHyperGAll_Subset) \n")
			if(id_type=="affy") {
				affy_sample <- CL_DF[CL_DF[,2]==i, 1]
				test_sample <- .AffyID2GeneID(catdb=catdb, affyIDs=affy_sample, probe2gene=1)
			}
			if(id_type=="gene") {
				test_sample <- as.vector(CL_DF[CL_DF[,2]==i, 1])
			}	
			containerDF2 <- data.frame(GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
			count <- 0
			for(j in gocats) {
				count <- count+1
				GOHyperGAll_result <- GOHyperGAll(catdb, gocat=j, sample=test_sample, ...)
                        	if(any(myslimv == "default")) {
                               		slimv <- c("GO:0003674", "GO:0008150", "GO:0005575", "GO:0030246","GO:0008289","GO:0003676","GO:0000166","GO:0019825","GO:0005515","GO:0003824","GO:0030234","GO:0003774","GO:0004871","GO:0005198","GO:0030528","GO:0045182","GO:0005215","GO:0006519","GO:0007154","GO:0016043","GO:0006412","GO:0006464","GO:0006810","GO:0007275","GO:0007049","GO:0005975","GO:0006629","GO:0006139","GO:0019748","GO:0015979","GO:0005618","GO:0005829","GO:0005783","GO:0005768","GO:0005794","GO:0005739","GO:0005777","GO:0009536","GO:0005840","GO:0005773","GO:0005764","GO:0005856","GO:0005634","GO:0005886","GO:0005576") # contains new unknown terms: "GO:0003674", "GO:0008150", "GO:0005575"
                                } else { 
                                	slimv <- myslimv 
				} 
				tempDF <- GOHyperGAll_Subset(catdb, GOHyperGAll_result, sample=test_sample, type="goSlim", myslimv=slimv)
				containerDF2 <- rbind(containerDF2, tempDF)
				if(length(recordSpecGO)>0) {
					containerDF2 <- rbind(containerDF2, GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],])
				}
				no_annot <- test_sample[!test_sample %in% eval(parse(text=paste("catmap(catdb)$D_", j, sep="")))[,2]]
				no_annot <- no_annot[no_annot!="no_match"]
				containerDF2 <- rbind(containerDF2, data.frame(GOID=paste("no_annot_", j, sep=""), NodeSize=NA, SampleMatch=length(no_annot), Phyper=NA, Padj=NA, Term=NA, Ont=j, SampleKeys=paste(no_annot, collapse=", ")))
			}
			tempDF2 <- data.frame(CLID=rep(i, times=length(containerDF2[,1])), CLSZ=rep(unique(as.vector(CL_DF[CL_DF[,2]==i,3])), times=length(containerDF2[,1])), containerDF2)
			containerDF <- rbind(containerDF, tempDF2)
		}
		return(containerDF)
	}	
	if(method=="simplify") {
        	containerDF <- data.frame(CLID=NULL, CLSZ=NULL, GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL, GO_OL_Match=NULL)
		for(i in cluster_loop) {
			cat("\n", "Processing cluster no", i, "with method: \"simplify\" (GOHyperGAll_Simplify) \n")
			if(id_type=="affy") {
				affy_sample <- CL_DF[CL_DF[,2]==i, 1]
				test_sample <- .AffyID2GeneID(catdb=catdb, affyIDs=affy_sample, probe2gene=1)
			}
			if(id_type=="gene") {
				test_sample <- as.vector(CL_DF[CL_DF[,2]==i, 1])
			}	
			containerDF2 <- data.frame(GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL, GO_OL_Match=NULL)
			count <- 0
			for(j in gocats) {
				count <- count+1
				GOHyperGAll_result <- GOHyperGAll(catdb, gocat=j, sample=test_sample, ...)
				simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result, gocat=j, cutoff=cutoff, correct=correct)
				if(length(simplifyDF)==0) { # If simplifyDF() returns empty data frame, then include at least the first two best scoring GO entries 
					simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result[1:2,], gocat=j, cutoff=1, correct=TRUE) 
				}
				tempDF <- data.frame(GOHyperGAll_result[GOHyperGAll_result[,1] %in% simplifyDF[,1], ], GO_OL_Match=simplifyDF[,2])
				containerDF2 <- rbind(containerDF2, tempDF)
				if(length(recordSpecGO)>0) {
					containerDF2 <- rbind(containerDF2, data.frame(GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],], GO_OL_Match=GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],3]))
				}
				no_annot <- test_sample[!test_sample %in% eval(parse(text=paste("catmap(catdb)$D_", j, sep="")))[,2]]
				no_annot <- no_annot[no_annot!="no_match"]
				containerDF2 <- rbind(containerDF2, data.frame(GOID=paste("no_annot_", j, sep=""), NodeSize=NA, SampleMatch=length(no_annot), Phyper=NA, Padj=NA, Term=NA, Ont=j, SampleKeys=paste(no_annot, collapse=", "), GO_OL_Match=length(no_annot)))
			}
			tempDF2 <- data.frame(CLID=rep(i, times=length(containerDF2[,1])), CLSZ=rep(unique(as.vector(CL_DF[CL_DF[,2]==i,3])), times=length(containerDF2[,1])), containerDF2)
			containerDF <- rbind(containerDF, tempDF2)
		}
		containerDF <- containerDF[, c(1:9,11,10)]
		return(containerDF)
	}	
}

## Apply GOCluster_Report
# testlist <- list(Set1=test_sample)
# GOBatchResult <- GOCluster_Report(catdb=catdb, setlist=testlist, method="all", id_type="gene", CLSZ=10, cutoff=0.001, gocats=c("MF", "BP", "CC"), recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575"))

######################
## (D.3) GO Barplot ##
######################
goBarplot <- function(GOBatchResult, gocat) {
	x <- GOBatchResult
	x <- x[,c("SampleMatch", "Term", "CLID", "Ont")]
	x <- x[x$Ont==gocat,1:3]
	colnames(x)[grep("CLID", colnames(x))] <- "Sample"
	term <- as.character(x$Term); term[is.na(term)] <- "no GO assignment"; x[,"Term"] <- term
	ontnames <- c(CC="Cellular Component", BP="Biological Process", MF="Molecular Function")
	x[,2] <- factor(x[,2], levels=unique(x[,2]), ordered=TRUE) # Defines plotting order of bars!!!
	p <- ggplot(x, aes(Term, SampleMatch, fill=Sample)) + 
		geom_bar(position="dodge", stat="identity") +
		coord_flip() + 
		theme(axis.text.y=element_text(angle=0, hjust=1)) + 
		xlab("GO Term") + 
		ylab("Gene Count") + 
		ylim(0, max(x$SampleMatch)) + 
		ggtitle(ontnames[gocat])
	print(p)
}
## Usage:
# goBarplot(GOBatchResult, gocat="MF")

