######################
## Filter VCF files ##
######################
filterVars <- function (files, filter, varcaller="gatk", organism, out_dir="results") {
  stopifnot(is.character(files))
  stopifnot(is.character(filter) && length(filter) == 1)
  stopifnot(is.character(organism) && length(organism) == 1)
  stopifnot(is.character(out_dir) && length(out_dir) == 1)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(out_dir)) stop("Cannot create output directory", out_dir)
  ## global functions or variables
  totalDepth<- refDepth <- altDepth <- totalDepth <- refDepth <- altDepth <- NULL
  if (class(files) %in% c("SYSfiles", "SYSfiles2")) {
    retrun(warning('filterVars: New version of SPR no longer accept "SYSfiles", "SYSfiles2" objects as inputs.\n',
                   'Use `getColumn` to get a vector of paths instead.'))
  }
  if(!all(check_files <- file.exists(files))) stop("Some files are missing:\n", paste0(files[!check_files], collapse = ",\n"))
  if(!all(check_ext <- stringr::str_detect(files, "\\.vcf$"))) stop("filterVars: All files need to end with .vcf\n", paste0(files[!check_ext], collapse = ",\n"))
  outfiles <- file.path(out_dir, basename(gsub("\\.vcf", "_filter.vcf", files)))
  for (i in seq(along = files)) {
    vcf <- readVcf(files[i], organism)
    vr <- as(vcf, "VRanges")
    if (varcaller == "gatk") {
      vrfilt <- vr[eval(parse(text = filter)), ]
    }
    if (varcaller == "bcftools") {
      vrsambcf <- vr
      vr <- unlist(values(vr)$DP4)
      vr <- matrix(vr, ncol = 4, byrow = TRUE)
      totalDepth(vrsambcf) <- as.integer(values(vrsambcf)$DP)
      refDepth(vrsambcf) <- rowSums(vr[, 1:2])
      altDepth(vrsambcf) <- rowSums(vr[, 3:4])
      vrfilt <- vrsambcf[eval(parse(text = filter)), ]
    }
    vcffilt <- asVCF(vrfilt)
    writeVcf(vcffilt, outfiles[i], index = TRUE)
    print(paste("Generated file", i, gsub(".*/", "", paste0(outfiles[i], ".bgz"))))
  }
  out_paths <- paste0(outfiles, ".bgz")
  names(out_paths) <- names(files)
  out_paths
}
## Usage for GATK:
# filter <- "totalDepth(vr) >= 20 & (altDepth(vr) / totalDepth(vr) >= 0.8) & rowSums(softFilterMatrix(vr))==6"
# filterVars(args, filter, varcaller="gatk", organism="Pinfest")
## Usage for BCFtools:
# filter <- "rowSums(vr) >= 20 & (rowSums(vr[,3:4])/rowSums(vr[,1:4]) >= 0.8)"
# filterVars(args, filter, varcaller="bcftools", organism="Pinfest")

#################
## VAR Reports ##
#################
## Report for locatVariants() where data for each variant is collapsed to a single line.
.allAnnot <- function(x, vcf) {
  rd <- rowRanges(vcf)
  ## Make variant calls in rd unique by collapsing duplicated ones
  VARID <- VARID <- unique(names(rd))
  REF <- tapply(as.character(values(rd)$REF), factor(names(rd)), function(i) paste(unique(i), collapse=" "))
  ALT <- tapply(as.character(unlist(values(rd)$ALT)), factor(names(rd)), function(i) paste(unique(i), collapse=" "))
  QUAL <- tapply(values(rd)$QUAL, factor(names(rd)), function(i) paste(unique(i), collapse=" "))

  ## fix names field in x if incomplete
  if(any(names(x)=="")) {
    index <- unique(names(rd)); names(index) <- gsub("_.*", "", index)
    names(x) <- index[paste(as.character(seqnames(x)), ":", start(x), sep="")]
  }

  ## Make annotated variant calls in x unique by collapsing duplicated ones
  LOCATION <- tapply(as.character(values(x)$LOCATION), as.factor(names(x)), function(i) paste(i, collapse=" "))
  GENEID <- tapply(values(x)$GENEID, factor(names(x)), function(i) paste(unique(i), collapse=" "))

  ## Assemble results in data.frame
  df <- data.frame(VARID=VARID,
                   REF=REF[VARID],
                   ALT=ALT[VARID],
                   QUAL=QUAL[VARID],
                   LOCATION=LOCATION[VARID],
                   GENEID=GENEID[VARID])
  df[,"LOCATION"] <- gsub("NA", "", df$LOCATION)
  df[,"GENEID"] <- gsub("NA", "", df$GENEID)
  return(df)
}
## Usage:
# allvar <- locateVariants(rd, txdb, AllVariants())
# varreport <- .allAnnot(allvar, vcf)

## Report for predictCoding() where data for each variant if collapsed to one line.
.codingReport <- function(x, txdb) {
  txids <- values(transcripts(txdb))$tx_name; names(txids) <- values(transcripts(txdb))$tx_id
  #myf <- as.factor(names(values(x)$CDSLOC))
  myf <- as.factor(names(x))
  if(length(myf)>0) {
    df <- data.frame(VARID=tapply(as.character(myf), myf, unique),
                     Strand=tapply(as.character(strand(x)), myf, unique),
                     Consequence=tapply(as.character(values(x)$CONSEQUENCE), myf, function(i) paste(unique(i), collapse=" ")),
                     Codon=tapply(paste(start(values(x)$CDSLOC), "_", as.character(values(x)$REFCODON), "/", as.character(values(x)$VARCODON), sep=""), myf, paste, collapse=" "),
                     AA=tapply(paste(sapply(values(x)$PROTEINLOC, paste, collapse="_"), "_", as.character(values(x)$REFAA), "/", as.character(values(x)$VARAA), sep=""), myf, paste, collapse=" "),
                     TXIDs=tapply(txids[values(x)$TXID], myf, paste, collapse=" "),
                     GENEIDcode=tapply(values(x)$GENEID, myf, function(i) paste(unique(i), collapse=" ")))
  } else {
    df <- data.frame(VARID=NA, Strand=NA, Consequence=NA, Codon=NA, AA=NA, TXIDs=NA, GENEIDcode=NA)[FALSE,,drop=FALSE]
  }
  return(df)
}
## Usage:
# codereport <- predictCoding(vcf, txdb, seqSource=fa)
# codereport <- .codingReport(coderport, txdb)

####################
## Variant Report ##
####################
variantReport <- function (files, txdb, fa, organism, out_dir="results") {
  if (class(files) %in% c("SYSfiles", "SYSfiles2")) {
    retrun(warning('filterVars: New version of SPR no longer accept "SYSfiles", "SYSfiles2" objects as inputs.\n',
                   'Use `getColumn` to get a vector of paths instead.'))
  }

  stopifnot(is.character(files))
  stopifnot(is.character(organism) && length(organism) == 1)
  stopifnot(is.character(out_dir) && length(out_dir) == 1)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(out_dir)) stop("Cannot create output directory", out_dir)

  if(!all(check_files <- file.exists(files))) stop("Some files are missing:\n", paste0(files[!check_files], collapse = ",\n"))
  if(!all(check_ext <- stringr::str_detect(files, "(\\.vcf|\\.bgz)$"))) stop("filterVars: All files need to end with .vcf or .bgz\n", paste0(files[!check_ext], collapse = ",\n"))
  outfiles <- file.path(out_dir, basename(gsub("(\\.vcf$|\\.vcf.bgz$|\\.bgz$)", "_anno.tsv", files)))

  for (i in seq(along = files)) {
    vcf <- readVcf(files[i], organism)
    allvar <- locateVariants(vcf, txdb, AllVariants())
    varreport <- .allAnnot(allvar, vcf)
    coding <- predictCoding(vcf, txdb, seqSource = fa)
    codereport <- .codingReport(coding, txdb)
    vr <- as(vcf, "VRanges")
    varid <- paste(as.character(seqnames(vr)), ":", start(vr),
                   "_", ref(vr), "/", alt(vr), sep = "")
    vrdf <- data.frame(row.names = varid, as.data.frame(vr))
    vrdf <- vrdf[, c("totalDepth", "refDepth", "altDepth")]
    fullreport <- cbind(varreport, codereport[rownames(varreport), -1])
    fullreport <- cbind(VARID = as.character(fullreport[, 1]),
                        vrdf[as.character(rownames(fullreport)), ],
                        fullreport[, -1])
    fullreport <- data.frame(lapply(fullreport, as.character),
                             stringsAsFactors = FALSE)
    write.table(fullreport, file = outfiles[i], row.names = FALSE,
                quote = FALSE, sep = "\t", na = "")
    print(paste("Generated file", i, gsub(".*/", "", outfiles[i])))
  }
  names(outfiles) <- names(files)
  outfiles
}
## Usage:
# variantReport(args=args, txdb=txdb, fa=fa, organism="Pinfest")

#############################
## Combine Variant Reports ##
#############################
combineVarReports <- function(files, filtercol, ncol=15) {
  if (class(files) %in% c("SYSfiles", "SYSfiles2")) {
    retrun(warning('filterVars: New version of SPR no longer accept "SYSfiles", "SYSfiles2" objects as inputs.\n',
                   'Use `getColumn` to get a vector of paths instead.'))
  }

  stopifnot(is.character(files))
  if(!all(check_files <- file.exists(files))) stop("Some files are missing:\n", paste0(files[!check_files], collapse = ",\n"))
  if(!all(check_ext <- stringr::str_detect(files, "(\\.tsv)$"))) stop("filterVars: All files need to end with .tsv\n", paste0(files[!check_ext], collapse = ",\n"))

  samples <- names(files)
  for(i in seq(along=samples)) {
    if(i==1) {
      varDF <- read.delim(files[i], colClasses=rep("character", ncol))
      varDF <- cbind(Sample=samples[i], varDF)
      if(filtercol[1]!="All") varDF <- varDF[varDF[,names(filtercol)]==filtercol,]
    } else {
      tmpDF <- read.delim(files[i])
      tmpDF <- read.delim(files[i], colClasses=rep("character", ncol))
      tmpDF <- cbind(Sample=samples[i], tmpDF)
      if(filtercol[1]!="All") tmpDF <- tmpDF[tmpDF[,names(filtercol)]==filtercol,]
      varDF <- rbind(as.data.frame(as.matrix(varDF)), as.data.frame(as.matrix(tmpDF)))
    }
    varDF <- varDF[order(varDF$VARID),]
  }
  return(varDF)
}
## Usage:
# args <- systemArgs(sysma="annotate_vars.param", mytargets="targets_gatk_filtered.txt")
# combineDF <- combineVarReports(args, filtercol=c(Consequence="nonsynonymous"))

###########################################
## Create summary statistics of variants ##
###########################################
varSummary <- function(files) {
  if (class(files) %in% c("SYSfiles", "SYSfiles2")) {
    retrun(warning('filterVars: New version of SPR no longer accept "SYSfiles", "SYSfiles2" objects as inputs.\n',
                   'Use `getColumn` to get a vector of paths instead.'))
  }

  stopifnot(is.character(files))
  if(!all(check_files <- file.exists(files))) stop("Some files are missing:\n", paste0(files[!check_files], collapse = ",\n"))
  if(!all(check_ext <- stringr::str_detect(files, "(\\.tsv)$"))) stop("filterVars: All files need to end with .tsv\n", paste0(files[!check_ext], collapse = ",\n"))

  for(i in seq(along=files)) {
    annotDF <- read.delim(files[i])
    count <- c(all=length(annotDF[,1]),
               table(unlist(strsplit(as.character(annotDF$LOCATION), " "))),
               table(unlist(strsplit(as.character(annotDF$Consequence), " "))))
    if(i==1) {
      countDF <- data.frame(count)
    } else {
      countDF <- cbind(countDF, count[rownames(countDF)])
    }
  }
  countDF[is.na(countDF)] <- 0
  colnames(countDF) <- names(files)
  return(countDF)
}
## Usage:
# args <- systemArgs(sysma="annotate_vars.param", mytargets="targets_gatk_filtered.txt")
# varSummaryDF <- varSummary(args)
