##############################################
## Set Intersect and Venn Diagram Functions ##
##############################################
## Author: Thomas Girke
## Last update: Jan 5, 2015

##############################################
## Class and Method Definitions for VENNset ##
##############################################
## Define VENNset class
setClass("VENNset", representation(setlist="list", intersectmatrix="matrix", complexitylevels="integer", vennlist="list"))

## Methods to return VENNset components
setGeneric(name="setlist", def=function(x) standardGeneric("setlist"))
setMethod(f="setlist", signature="VENNset", definition=function(x) {return(x@setlist)})
setGeneric(name="intersectmatrix", def=function(x) standardGeneric("intersectmatrix"))
setMethod(f="intersectmatrix", signature="VENNset", definition=function(x) {return(x@intersectmatrix)})
setGeneric(name="complexitylevels", def=function(x) standardGeneric("complexitylevels"))
setMethod(f="complexitylevels", signature="VENNset", definition=function(x) {return(x@complexitylevels)})
setGeneric(name="vennlist", def=function(x) standardGeneric("vennlist"))
setMethod(f="vennlist", signature="VENNset", definition=function(x) {return(x@vennlist)})
setGeneric(name="as.list", def=function(x) standardGeneric("as.list"))
setMethod(f="as.list", signature="VENNset", definition=function(x) {
	mylist <- list(setlist=x@setlist, 
                       intersectmatrix=x@intersectmatrix, 
                       complexitylevels=x@complexitylevels, 
                       vennlist=x@vennlist)
	return(mylist)
})

## Constructor methods
## List to VENNset with: as(mylist, "VENNset")
setAs(from="list", to="VENNset",  
        def=function(from) {
		new("VENNset", setlist=from$setlist,
		             intersectmatrix=from$intersectmatrix,
		             complexitylevels=from$complexitylevels,
			     vennlist=from$vennlist)
})

## Define print behavior for VENNset
setMethod(f="show", signature="VENNset", 
	definition=function(object) {    
	cat("An instance of '", class(object), "' with ", length(object@setlist), " label sets ", "\n", sep="")
})

## Extend names() method
setMethod(f="names", signature="VENNset",
    definition=function(x) {
    	return(slotNames(x))
})

## Extend length() method
setMethod(f="length", signature="VENNset",
    definition=function(x) {
        return(length(x@setlist))
})

###################################################
## Class and Method Definitions for INTERSECTset ##
###################################################
## Define INTERSECTset class
setClass("INTERSECTset", representation(setlist="list", intersectmatrix="matrix", complexitylevels="integer", intersectlist="list"))

## Methods to return INTERSECTset components
setMethod(f="setlist", signature="INTERSECTset", definition=function(x) {return(x@setlist)})
setMethod(f="intersectmatrix", signature="INTERSECTset", definition=function(x) {return(x@intersectmatrix)})
setMethod(f="complexitylevels", signature="INTERSECTset", definition=function(x) {return(x@complexitylevels)})
setGeneric(name="intersectlist", def=function(x) standardGeneric("intersectlist"))
setMethod(f="intersectlist", signature="INTERSECTset", definition=function(x) {return(x@intersectlist)})
setMethod(f="as.list", signature="INTERSECTset", definition=function(x) {
	mylist <- list(setlist=x@setlist, 
                       intersectmatrix=x@intersectmatrix, 
                       complexitylevels=x@complexitylevels, 
                       intersectlist=x@intersectlist)
	return(mylist)
})

## Constructor methods
## List to INTERSECTset with: as(mylist, "INTERSECTset")
setAs(from="list", to="INTERSECTset",  
        def=function(from) {
		new("INTERSECTset", setlist=from$setlist,
		             intersectmatrix=from$intersectmatrix,
		             complexitylevels=from$complexitylevels,
			     intersectlist=from$intersectlist)
})

## Define print behavior for INTERSECTset
setMethod(f="show", signature="INTERSECTset", 
	definition=function(object) {    
	cat("An instance of '", class(object), "' with ", length(object@setlist), " label sets and 'complexity = ", paste0(unique(complexitylevels(vennset)), collapse=", "), "'", "\n", sep="")
})

## Extend names() method
setMethod(f="names", signature="INTERSECTset",
    definition=function(x) {
    	return(slotNames(x))
})

## Extend length() method
setMethod(f="length", signature="INTERSECTset",
    definition=function(x) {
        return(length(x@setlist))
})

################################
## Generic Intersect Function ##
################################
## Computation of (1) Venn Intersects and (2) Pairwise Intersects
overLapper <- function(setlist, complexity="default", sep="_", cleanup=FALSE, keepdups=FALSE, type) {
	## Default complexity levels
	if(complexity[1]=="default") complexity <- 1:length(setlist)
	
	## Check validity of inputs
	if(class(setlist)!="list" | length(names(setlist))==0) stop("Unexpected input.
             The input 'setlist' needs to be of class 'list' where each list component stores a
             label set as 'vector' and the name of each label set is provided under the
             name slot of each list component.")
	if(!all(sapply(setlist, is.vector)) | length(setlist) < 2) stop("Unexpected input. 
	     The input 'setlist' needs to be a list with at least 2 components each 
             containig a vector of set labels.")
	if(length(type)!=1 & all(c("vennsets", "intersects") %in% type)) stop("Argument 'type' needs to be assigned 'vennsets' or 'intersects'.")
	if(type=="vennsets" & !identical(complexity, 1:length(setlist))) stop("When assigning 'vennsets' to 'type', then 'complexity' needs to be assigned 'default' or '1:length(setlist)'.")
	
	## Clean up of sample sets to minimize formatting issues 
	if(cleanup==TRUE) {
		## Set all characters to upper case 
		setlist <- sapply(setlist, function(x) gsub("([A-Z])", "\\U\\1", x, perl=TRUE, ignore.case=TRUE))
		## Remove leading and trailing spaces
		setlist <- sapply(setlist, function(x) gsub("^ {1,}| {1,}$", "", x, perl=TRUE, ignore.case=TRUE))
	}
	
	## Append object counter to retain duplicates 
	if(keepdups==TRUE) {
		dupCount <- function(setlist=setlist) {
			count <- table(setlist)
			paste(rep(names(count), count), unlist(sapply(count, function(x) seq(1, x))), sep=".")
		}
		mynames <- names(setlist)
		setlist <- lapply(setlist, function(x) dupCount(x)) # lapply necessary for numeric data!
		names(setlist) <- mynames
	}	

	## Create intersect matrix (removes duplicates!)
	setunion <- sort(unique(unlist(setlist)))
	setmatrix <- sapply(names(setlist), function(x) setunion %in% unique(setlist[[x]])) 
	rownames(setmatrix) <- setunion
	storage.mode(setmatrix) <- "numeric"

	## Create all possible sample combinations within requested complexity levels
	labels <- names(setlist)
	allcombl <- lapply(complexity, function(x) combn(labels, m=x, simplify=FALSE))
	allcombl <- unlist(allcombl, recursive=FALSE)
	complevels <- sapply(allcombl, length)
	
	## Return intersect list for generated sample combinations 
	if(type=="intersects") {
		OLlist <- sapply(seq(along=allcombl), function(x) setunion[rowSums(setmatrix[, rep(allcombl[[x]], 2)]) == 2 * length(allcombl[[x]])])
		names(OLlist) <- sapply(allcombl, paste, collapse=sep)
		OLlist <- list(setlist=setlist, intersectmatrix=setmatrix, complexitylevels=complevels, intersectlist=OLlist)
		intersectset <- as(OLlist, "INTERSECTset")
		return(intersectset)
	}	

	## Return Venn intersect list for generated sample combinations 
	if(type=="vennsets") {
		vennSets <- function(setmatrix=setmatrix, allcombl=allcombl, index=1) {
			mycol1 <- which(colnames(setmatrix) %in% allcombl[[index]])
			mycol2 <- which(!colnames(setmatrix) %in% allcombl[[index]])
			cond1 <- rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
			cond2 <- rowSums(setmatrix[, rep(mycol2, 2)]) == 0
			return(setunion[cond1 & cond2])
		}
		vennOLlist <- sapply(seq(along=allcombl), function(x) vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x))
		names(vennOLlist) <- sapply(allcombl, paste, collapse=sep)
		OLlist <- list(setlist=setlist, intersectmatrix=setmatrix, complexitylevels=complevels, vennlist=vennOLlist)
		vennset <- as(OLlist, "VENNset")
		return(vennset)
	}
}

####################################
## Venn Diagram Plotting Function ##
####################################
vennPlot <- function(x, mymain="Venn Diagram", mysub="default", setlabels="default", yoffset=seq(0,10, by=0.34), 
                     ccol=rep(1,31), colmode=1, lcol=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), 
                     lines=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), mylwd=3, diacol=1, type="ellipse", 
                     ccex=1.0, lcex=1.0, sepsplit="_", ...) {
	
	## Check validity of inputs and 
	if(!any(c(class(x)=="VENNset", class(x)=="list", is.numeric(x)))) {
		stop("x needs to be one of: VENNset, list of VENNsets, named numeric vector, or list of named numeric vectors.")
	}
	if(class(x)=="list") {
		if(length(unique(sapply(x, length))) != 1) stop("List components need to have identical length.") 
	}
	
	## Compute count set(s) 
	# If x is VENNset
	if(class(x)=="VENNset") { 
	 	counts <- list(sapply(vennlist(x), length))
		myclass <- "VENNset"
	# If x is list of VENNsets
	} else if(class(x)=="list" & all(sapply(x, class)=="VENNset")) { 
		counts <- lapply(x, function(y) sapply(vennlist(y), length))	
		myclass <- "VENNset"
	## If x is count set (named numeric vector) 
	} else if(is.numeric(x) & is.list(x)==FALSE) {
		counts <- list(x)
		myclass <- "numeric"
	## If x is list of count sets
	} else if(class(x)=="list" & all(sapply(x, is.numeric))) {
		counts <- x
		myclass <- "numeric"
	} else {
		stop("x needs to be one of: VENNset, list of VENNsets, named numeric vector, or list of named numeric vectors.")
	}
	
	## Check for supported number of Venn counts: 3, 7, 15 and 31
	if(!length(counts[[1]]) %in%  c(3,7,15,31)) stop("Only 2-5 way venn comparisons are supported.")

        ## Function to return for a set label the index of matches in the name field of a counts object
        grepLabel <- function(label, x=names(counts[[1]])) {
                x <- strsplit(x, sepsplit)
                as.numeric(which(sapply(x, function(y) any(y==label))))
        }
	
	## 2-way Venn diagram
	if(length(counts[[1]])==3) {
		## Define subtitle
		if(mysub=="default") {
			if(myclass=="numeric") {
                        	n <- names(counts[[1]])[1:2]
                        	if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        	sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        	if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
				} else {
                                	sample_counts <- sapply(n, function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
                        	}
				mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), sep="")
			} else if(myclass=="VENNset") {
				if(class(x)=="list") x <- x[[1]]
				sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
                        	mysub <- paste(paste("Unique objects: All =", length(unique(unlist(setlist(x))))), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), sep="")
			} else {
				mysub <- mysub
			}
		}
		## Plot venn shapes
		symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
		
		## Add counts
		for(i in seq(along=counts)) {
			olDF <- data.frame(x=c(3.1, 7.0, 5.0), 
                                           y=c(6.0, 6.0, 6.0), 
                                           counts=counts[[i]])
                        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                        if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) } # For coloring several numbers per intersect differently. ccol can needs to be list to color each field differently..
		}
                
		## Add sample labels
		if(length(setlabels)==1 & setlabels[1]=="default") { 
			setlabels <- names(counts[[1]][1:2])
		} else {
			setlabels <- setlabels
		}
		text(c(2.0, 8.0), c(8.8, 8.8), labels=setlabels, col=lcol, cex=lcex, ...)	
	}
 
	## 3-way Venn diagram
	if(length(counts[[1]])==7) { 
		## Define subtitle
		if(mysub=="default") {
			if(myclass=="numeric") {
                        	n <- names(counts[[1]])[1:3]
                        	if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        	sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        	if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
				} else {
				        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
				}
                        	mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), sep="")
			} else if(myclass=="VENNset") {
				if(class(x)=="list") x <- x[[1]]
				sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
                        	mysub <- paste(paste("Unique objects: All =", length(unique(unlist(setlist(x))))), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), sep="")
			} else {
				mysub <- mysub
			}
		}
		## Plot venn shapes
		symbols(x=c(4, 6, 5), y=c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)
		
		## Add counts
		for(i in seq(along=counts)) {
			olDF <- data.frame(x=c(3.0, 7.0, 5.0, 5.0, 3.8, 6.3, 5.0), 
                                           y=c(6.5, 6.5, 3.0, 7.0, 4.6, 4.6, 5.3), 
                                           counts=counts[[i]])
	        	 if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                         if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }

		}

                ## Add sample labels
		if(length(setlabels)==1 & setlabels[1]=="default") { 
			setlabels <- names(counts[[1]][1:3])
		} else {
			setlabels <- setlabels
		}
		text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels=setlabels, col=lcol, cex=lcex, ...)	
	}
	
	## 4-way Venn diagram with ellipses
	if(length(counts[[1]])==15 & type=="ellipse") {
		## Define subtitle
		if(mysub=="default") {
			if(myclass=="numeric") {
                        	n <- names(counts[[1]])[1:4]
                        	if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        	sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        	if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
				} else {
			        	sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
				}
                        	mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
			} else if(myclass=="VENNset") {
				if(class(x)=="list") x <- x[[1]]
				sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
                        	mysub <- paste(paste("Unique objects: All =", length(unique(unlist(setlist(x))))), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
			} else { 
				mysub <- mysub
			}
		}
		## Plot ellipse
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments  
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		## Plot ellipse as 4-way venn diagram
		ellipseVenn <- function(...) {
			split.screen(c(1,1))
			plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
			## Add counts
			for(i in seq(along=counts)) {
				olDF <- data.frame(x=c(1.5, 3.5, 6.5, 8.5, 2.9, 3.1, 5.0, 5.0, 6.9, 7.1, 3.6, 5.8, 4.2, 6.4, 5.0), 
                                                   y=c(4.8, 7.2, 7.2, 4.8, 5.9, 2.2, 0.7, 6.0, 2.2, 5.9, 4.0, 1.4, 1.4, 4.0, 2.8), 
                                                   counts=counts[[i]])
	        	        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                                if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }
			}
			## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") { 
				setlabels <- names(counts[[1]][1:4])
			} else {
				setlabels <- setlabels
			}
			text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels=setlabels, col=lcol, cex=lcex, ...)
			close.screen(all=TRUE) 
		}
		ellipseVenn(...)
	} 

	## 4-way Venn diagram with circles (pseudo-venn diagram that misses two overlap sectors) 
	if(length(counts[[1]])==15 & type=="circle") {
		## Define subtitle
		if(mysub=="default") {
			if(myclass=="numeric") {
                        	n <- names(counts[[1]])[1:4]
                        	if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        	sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        	if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
				} else {
			        	sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
				}
                        	mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
			} else if(myclass=="VENNset") {
				if(class(x)=="list") x <- x[[1]]
				sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
                        	mysub <- paste(paste("Unique objects: All =", length(unique(unlist(setlist(x))))), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
			} else { 
				mysub <- mysub
			}
		}
		## Plot venn shapes
		symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)
		
		## Add counts
		for(i in seq(along=counts)) {
		        olDF <- data.frame(x=c(3.0, 6.5, 3.0, 6.5, 4.8, 3.0, 4.8, 4.8, 6.5, 4.8, 3.9, 5.7, 3.9, 5.7, 4.8), 
                                           y=c(7.2, 7.2, 3.2, 3.2, 7.2, 5.2, 0.4, 0.4, 5.2, 3.2, 6.3, 6.3, 4.2, 4.2, 5.2), 
                                           counts=counts[[i]])
			 if(colmode==1) { text(olDF$x[-c(7,8)], olDF$y[-c(7,8)] + yoffset[i], olDF$counts[-c(7,8)], col=ccol, cex=ccex, ...) } # rows 14-15 of olDF are printed in next step
			 if(colmode==2) { text(olDF$x[-c(7,8)], olDF$y[-c(7,8)] + yoffset[i], olDF$counts[-c(7,8)], col=ccol[[i]], cex=ccex[i], ...) }
			 text(c(4.8), c(0.8) + yoffset[i], paste("Only in ", names(counts[[1]][1]), " & ", names(counts[[1]][4]), ": ", olDF$counts[7], "; Only in ", names(counts[[1]][2]), " & ", names(counts[[1]][3]), ": ", olDF$counts[8], sep=""), col=diacol, cex=ccex, ...)
		}

                ## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") { 
				setlabels <- names(counts[[1]][1:4])
			} else {
				setlabels <- setlabels
			}
		text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels=setlabels, col=lcol, cex=lcex, ...)
	} 
	
	## 5-way Venn diagram
	if(length(counts[[1]])==31) {
		## Define subtitle
		if(mysub=="default") {
			if(myclass=="numeric") {
                        	n <- names(counts[[1]])[1:5]
                        	if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        	sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        	if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  } 
				} else {
			        	sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
				}
                        	mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), paste("; S5 =", sample_counts[5]), sep="")
			} else if(myclass=="VENNset") {
				if(class(x)=="list") x <- x[[1]]
				sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
                        	mysub <- paste(paste("Unique objects: All =", length(unique(unlist(setlist(x))))), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), paste("; S5 =", sample_counts[5]), sep="")
			} else { 
				mysub <- mysub
			}
	 	}
		## Plot ellipse
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments  
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		## Plot ellipse as 5-way venn diagram
		ellipseVenn <- function(...) {
			split.screen(c(1,1))
			screen(1, new=FALSE)
			plotellipse(center=c(4.83,6.2), radius=c(1.43,4.11), rotate=0, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.25,5.4), radius=c(1.7,3.6), rotate=66, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.1,3.5), radius=c(1.55,3.9), rotate=150, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.48,3.15), radius=c(1.55,3.92), rotate=210, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(3.7,4.8), radius=c(1.7,3.6), rotate=293.5, segments=360, xlab="", ylab="", col=lines[5], axes=FALSE, lwd=mylwd, ...)

			## Add counts
			for(i in seq(along=counts)) {
				olDF <- data.frame(x=c(4.85, 8.0, 7.1, 3.5, 2.0, 5.90, 4.4, 4.60, 3.60, 7.1, 6.5, 3.2, 5.4, 6.65, 3.40, 5.00, 6.02, 3.60, 5.20, 4.03, 4.20, 6.45, 6.8, 3.39, 6.03, 5.74, 4.15, 3.95, 5.2, 6.40, 5.1), 
                                                   y=c(8.30, 6.2, 1.9, 1.6, 5.4, 6.85, 6.6, 2.45, 6.40, 4.3, 6.0, 4.6, 2.1, 3.40, 3.25, 6.43, 6.38, 5.10, 2.49, 6.25, 3.08, 5.30, 4.0, 3.80, 3.20, 5.95, 5.75, 3.75, 3.0, 4.50, 4.6),
					counts=counts[[i]]) 
	        	        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                                if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }
			}
			## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") { 
				setlabels <- names(counts[[1]][1:5])
			} else {
				setlabels <- setlabels
			}
			text(c(5.7, 7.9, 8.5, 4.2, 0.8), c(9.9, 7.9, 1.9, 0.0, 7.3), adj=c(0, 0.5), labels=setlabels, col=lcol, cex=lcex, ...)
			close.screen(all=TRUE) 
		}
		ellipseVenn(...)
	} 
}

##########################################
## Bar Plot Function for Intersect Sets ##
##########################################
## Plots the counts of VENNset/INTERSECTset generated by overLapper function
olBarplot <- function(x, mincount=0, complexity="default", myxlabel="default", myylabel="Counts", mytitle="default", ...) {
	## Check validity of inputs 
	if(!any(c("VENNset", "INTERSECTset") %in% class(x))) stop("x needs to be of class VENNset or INTERSECTset")

	## Generate counts 
	if(class(x)=="VENNset") counts <- sapply(vennlist(x), length)
	if(class(x)=="INTERSECTset") counts <- sapply(intersectlist(x), length)
	
	## Complexity filter
	if(complexity[1]!="default") {
		complfilter <- complexitylevels(x) %in% complexity
	} else {
		complfilter <- complexitylevels(x) %in% complexitylevels(x)
	}

	## Min count filter
	mincountfilter <- counts >= mincount
	
	## Apply filters
	myfilter <- complfilter & mincountfilter
	counts <- counts[myfilter]
	
	## Color bars by default by complexity levels 
	mycol <- complexitylevels(x)
	mycol <- mycol[myfilter] 

	## Define x-axis label
	if(myxlabel=="default") {
		myxlabel <- paste("Intersect Sets (min count cutoff = ", mincount, ")", sep="")
	} else {
		myxlabel <- myxlabel
	}
	
	## Define main title
	if(mytitle=="default") {
		mytitle <- paste("Intersect Plot of", class(x), "Object")
	} else {
		mytitle < mytitle	
	}
	
	## Generate bar plot with ggplot2
	df_plot <- data.frame(Intersect_Sets=names(counts), Counts=counts, Level=as.character(mycol))
	df_plot[,1] <- factor(df_plot[,1], levels=unique(df_plot[,1]), ordered=TRUE) # Defines plotting order of bars!!!	
	ggplot(df_plot, aes(Intersect_Sets, Counts, fill = Level)) + 
		geom_bar(position="stack", stat="identity", ...) + 
		coord_flip() + 
		theme(legend.position="none") +
		theme(axis.text.y=element_text(angle=0, hjust=1)) + 
		labs(x = myxlabel, y = myylabel) +
		ggtitle(mytitle) 
}


