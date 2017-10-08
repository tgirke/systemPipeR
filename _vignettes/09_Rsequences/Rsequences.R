## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(limma) 
    library(ggplot2) }) 

## ----package_requrirements, eval=FALSE-----------------------------------
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("Biostrings", "GenomicRanges", "GenomicRanges", "rtracklayer", "systemPipeR", "seqLogo", "ShortRead"))

## ----string_matching_base1, eval=TRUE------------------------------------
myseq <- c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC")

## ----string_matching_base2, eval=TRUE------------------------------------
myseq[grep("ATG", myseq)] 

## ----string_matching_base3, eval=TRUE------------------------------------
pos1 <- regexpr("AT", myseq) 
as.numeric(pos1); attributes(pos1)$match.length # Returns position information of matches

## ----string_matching_base4, eval=TRUE------------------------------------
pos2 <- gregexpr("AT", myseq) 
as.numeric(pos2[[1]]); attributes(pos2[[1]])$match.length # Returns positions of matches in first sequence

## ----string_matching_base5, eval=TRUE------------------------------------
gsub("^ATG", "atg", myseq) 

## ----positional_parsing_base, eval=TRUE----------------------------------
nchar(myseq) # Computes length of strings
substring(myseq[1], c(1,3), c(2,5)) # Positional parsing of several fragments from one string
substring(myseq, c(1,4,7), c(2,6,10)) # Positional parsing of many strings

## ----random_sequences_base, eval=TRUE------------------------------------
rand <- sapply(1:100, function(x) paste(sample(c("A","T","G","C"), sample(10:20), replace=T), collapse=""))
rand[1:3]

## ----count_sequences_base, eval=TRUE-------------------------------------
table(c(rand[1:4], rand[1]))

## ----parse_from_ref, eval=TRUE, message=FALSE, warnings=FALSE------------
library(Biostrings)
ref <- DNAString(paste(sample(c("A","T","G","C"), 100000, replace=T), collapse=""))
randstart <- sample(1:(length(ref)-15), 1000)
randreads <- Views(ref, randstart, width=15)
rand_set <- DNAStringSet(randreads)
unlist(rand_set)

## ----download_sequences, eval=TRUE, message=FALSE, warnings=FALSE--------
dir.create("data", showWarnings = FALSE)
# system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn")
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn", "data/AE004437.ffn")

## ----import_sequences1, eval=TRUE----------------------------------------
myseq <- readDNAStringSet("data/AE004437.ffn")
myseq[1:3]

## ----import_sequences2, eval=TRUE----------------------------------------
sub <- myseq[grep("99.*", names(myseq))]
length(sub)

## ----export_sequences, eval=TRUE-----------------------------------------
writeXStringSet(sub, file="./data/AE004437sub.ffn", width=80)

## ----xstring_container, eval=TRUE, message=FALSE, warnings=FALSE---------
library(Biostrings)
d <- DNAString("GCATAT-TAC")
d
d[1:4]

## ----rnastring_container, eval=TRUE, message=FALSE, warnings=FALSE-------
r <- RNAString("GCAUAU-UAC") 
r <- RNAString(d) # Converts d to RNAString object
r

## ----aastring_container, eval=TRUE, message=FALSE, warnings=FALSE--------
p <- AAString("HCWYHH")
p

## ----bstring_container, eval=TRUE, message=FALSE, warnings=FALSE---------
b <- BString("I store any set of characters. Other XString objects store only the IUPAC characters.")
b

## ----xtringset_container, eval=TRUE--------------------------------------
dset <- DNAStringSet(c("GCATATTAC", "AATCGATCC", "GCATATTAC")) 
names(dset) <- c("seq1", "seq2", "seq3") # Assigns names
dset[1:2]

## ----xtringset_container_utilities, eval=TRUE----------------------------
width(dset) # Returns the length of each sequences
d <- dset[[1]] # The [[ subsetting operator returns a single entry as XString object
dset2 <- c(dset, dset) # Appends/concatenates two XStringSet objects
dsetchar <- as.character(dset) # Converts XStringSet to named vector 
dsetone <- unlist(dset) # Collapses many sequences to a single one stored in a DNAString container

## ----xtringset_subsetting, eval=TRUE-------------------------------------
DNAStringSet(dset, start=c(1,2,3), end=c(4,8,5)) 

## ----xmultiple_alignment, eval=TRUE--------------------------------------
origMAlign <- readDNAMultipleAlignment(filepath = system.file("extdata",
              "msx2_mRNA.aln", package = "Biostrings"), format = "clustal")
origMAlign

## ----rev_comp, eval=TRUE-------------------------------------------------
randset <- DNAStringSet(rand)
complement(randset[1:2])
reverse(randset[1:2])
reverseComplement(randset[1:2])

## ----translate, eval=TRUE, message=FALSE, warnings=FALSE-----------------
translate(randset[1:2])

## ----pattern_matching1, eval=TRUE----------------------------------------
myseq1 <- readDNAStringSet("./data/AE004437.ffn") 
mypos <- matchPattern("ATGGTG", myseq1[[1]], max.mismatch=1) 

## ----pattern_matching2, eval=TRUE----------------------------------------
countPattern("ATGGCT", myseq1[[1]], max.mismatch=1) 

## ----pattern_matching3, eval=TRUE----------------------------------------
vcountPattern("ATGGCT", myseq1, max.mismatch=1)[1:20]

## ----pattern_matching4, eval=TRUE----------------------------------------
tmp <- c(DNAStringSet("ATGGTG"), DNAStringSet(mypos)) 

## ----pattern_matching5, eval=TRUE----------------------------------------
consensusMatrix(tmp)[1:4,] 

## ----pattern_matching6, eval=TRUE----------------------------------------
myvpos <- vmatchPattern("ATGGCT", myseq1, max.mismatch=1) 
myvpos # The results are stored as MIndex object.
Views(myseq1[[1]], start(myvpos[[1]]), end(myvpos[[1]])) # Retrieves the result for single entry

## ----all_matches, eval=FALSE---------------------------------------------
## sapply(seq(along=myseq1), function(x)
##        as.character(Views(myseq1[[x]], start(myvpos[[x]]), end(myvpos[[x]]))))[1:4]

## ----regex_pattern_matching, eval=TRUE-----------------------------------
myseq <- DNAStringSet(c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC"))
myseq[grep("^ATG", myseq, perl=TRUE)] # String searching with regular expression support
pos1 <- regexpr("AT", myseq) # Searches 'myseq' for first match of pattern "AT"
as.numeric(pos1); attributes(pos1)$match.length # Returns position information of matches
pos2 <- gregexpr("AT", myseq) # Searches 'myseq' for all matches of pattern "AT"
as.numeric(pos2[[1]]); attributes(pos2[[1]])$match.length # Match positions in first sequence
DNAStringSet(gsub("^ATG", "NNN", myseq)) # String substitution with regular expression support

## ----pwm_logo, eval=TRUE-------------------------------------------------
library(seqLogo) 
pwm <- PWM(DNAStringSet(c("GCT", "GGT", "GCA"))) 
pwm
seqLogo(t(t(pwm) * 1/colSums(pwm)))

## ----pwm_logo2, eval=TRUE------------------------------------------------
library(ggplot2); library(ggseqlogo)
pwm <- PWM(DNAStringSet(c("GCT", "GGT", "GCA"))) 
ggseqlogo(pwm)

## ----pwm_search, eval=TRUE-----------------------------------------------
chr <- DNAString("AAAGCTAAAGGTAAAGCAAAA") 
matchPWM(pwm, chr, min.score=0.9) 

## ----raw_to_char, eval=TRUE----------------------------------------------
phred <- 1:9
phreda <- paste(sapply(as.raw((phred)+33), rawToChar), collapse="")
phreda
as.integer(charToRaw(phreda))-33 

## ----raw_to_char2, eval=TRUE---------------------------------------------
dset <- DNAStringSet(sapply(1:100, function(x) paste(sample(c("A","T","G","C"), 20, replace=T), collapse=""))) # Creates random sample sequence.
myqlist <- lapply(1:100, function(x) sample(1:40, 20, replace=T)) # Creates random Phred score list.
myqual <- sapply(myqlist, function(x) toString(PhredQuality(x))) # Converts integer scores into ASCII characters.
myqual <- PhredQuality(myqual) # Converts to a PhredQuality object.
dsetq1 <- QualityScaledDNAStringSet(dset, myqual) # Combines DNAStringSet and quality data in QualityScaledDNAStringSet object.
dsetq1[1:2]

## ----read_fastq1, eval=TRUE, message=FALSE, warnings=FALSE---------------
library(ShortRead)
download.file("http://faculty.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_6_10_2012/Rsequences/data.zip", "data.zip")
unzip("data.zip")

## ----read_fastq2, eval=TRUE, message=FALSE, warnings=FALSE---------------
fastq <- list.files("data", "*.fastq$"); fastq <- paste("data/", fastq, sep="")
names(fastq) <- paste("flowcell6_lane", 1:length(fastq), sep="_") 
(fq <- readFastq(fastq[1])) # Imports first FASTQ file
countLines(dirPath="./data", pattern=".fastq$")/4 # Counts numbers of reads in FASTQ files
id(fq)[1] # Returns ID field
sread(fq)[1] # Returns sequence
quality(fq)[1] # Returns Phred scores 
as(quality(fq), "matrix")[1:4,1:12] # Coerces Phred scores to numeric matrix
ShortReadQ(sread=sread(fq), quality=quality(fq), id=id(fq)) # Constructs a ShortReadQ from components

## ----see_fastq, eval=TRUE, message=FALSE, warnings=FALSE, fig.height=12, fig.width=14----
library(systemPipeR)
fqlist <- seeFastq(fastq=fastq, batchsize=800, klength=8) # For real data set batchsize to at least 10^5 
seeFastqPlot(fqlist)

## ----shortread_fastq_qc, eval=FALSE--------------------------------------
## sp <- SolexaPath(system.file('extdata', package='ShortRead'))
## fl <- file.path(analysisPath(sp), "s_1_sequence.txt")
## fls <- c(fl, fl)
## coll <- QACollate(QAFastqSource(fls), QAReadQuality(), QAAdapterContamination(),
## 	    QANucleotideUse(), QAQualityUse(), QASequenceUse(), QAFrequentSequence(n=10),
## 		QANucleotideByCycle(), QAQualityByCycle())
## x <- qa2(coll, verbose=TRUE)
## res <- report(x)
## if(interactive())
## browseURL(res)

## ----adaptor_trimming, eval=TRUE-----------------------------------------
fqtrim <- trimLRPatterns(Rpattern="GCCCGGGTAA", subject=fq)
sread(fq)[1:2] # Before trimming
sread(fqtrim)[1:2] # After trimming

## ----read_counting, eval=TRUE--------------------------------------------
tables(fq)$distribution # Counts read occurences
sum(srduplicated(fq)) # Identifies duplicated reads
fq[!srduplicated(fq)]

## ----trim_tails, eval=TRUE-----------------------------------------------
cutoff <- 30
cutoff <- rawToChar(as.raw(cutoff+33))
sread(trimTails(fq, k=2, a=cutoff, successive=FALSE))[1:2]

## ----remove_low_quality_reads, eval=TRUE---------------------------------
cutoff <- 30
qcount <- rowSums(as(quality(fq), "matrix") <= 20) 
fq[qcount == 0] # Number of reads where all Phred scores >= 20

## ----remove_N_reads, eval=TRUE-------------------------------------------
filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) # Removes reads with >=20 of one nucleotide
filter <- compose(filter1, filter2)
fq[filter(fq)]

## ----stream_fastq1, eval=TRUE--------------------------------------------
fq <- yield(FastqStreamer(fastq[1], 50)) # Imports first 50 reads 
fq <- yield(FastqSampler(fastq[1], 50)) # Random samples 50 reads 

## ----stream_fastq2, eval=TRUE, message=FALSE, warnings=FALSE-------------
f <- FastqStreamer(fastq[1], 50) 
while(length(fq <- yield(f))) {
	fqsub <- fq[grepl("^TT", sread(fq))] 
	writeFastq(fqsub, paste(fastq[1], "sub", sep="_"), mode="a", compress=FALSE)
}
close(f)

## ----genomicranges1, eval=TRUE, message=FALSE, warnings=FALSE------------
library(GenomicRanges); library(rtracklayer)
gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), ranges = IRanges(1:10, end = 7:16, names = head(letters, 10)), strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)), score = 1:10, GC = seq(1, 0, length = 10)) # Example of creating a GRanges object with its constructor function.

## ----genomicranges2, eval=TRUE, message=FALSE, warnings=FALSE------------
gff <- import.gff("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/gff3.gff") # Imports a simplified GFF3 genome annotation file.
seqlengths(gff) <- end(ranges(gff[which(values(gff)[,"type"]=="chromosome"),])) 
names(gff) <- 1:length(gff) # Assigns names to corresponding slot
gff[1:4,]

## ----genomicranges3, eval=TRUE, message=FALSE, warnings=FALSE------------
as.data.frame(gff)[1:4, 1:7]

## ----genomicranges4, eval=TRUE, message=FALSE, warnings=FALSE------------
gff_rd <- as(gff, "RangedData") 
gff_gr <- as(gff_rd, "GRanges") 

## ----genomicranges_access2, eval=TRUE------------------------------------
gff[1:4]
gff[1:4, c("type", "ID")] 
gff[2] <- gff[3] 

## ----genomicranges_access3, eval=TRUE------------------------------------
c(gff[1:2], gff[401:402]) 

## ----genomicranges_access4, eval=TRUE------------------------------------
seqnames(gff)
ranges(gff)
strand(gff)
seqlengths(gff) 
start(gff[1:4])
end(gff[1:4])
width(gff[1:4]) 

## ----genomicranges_access5, eval=TRUE------------------------------------
values(gff) # or elementMetadata(gff)
values(gff)[, "type"][1:20] 
gff[values(gff)[ ,"type"] == "gene"] 

## ----genomicranges_utilities1, eval=TRUE---------------------------------
gff <- gff[values(gff)$type != "chromosome"] 

## ----genomicranges_utilities2, eval=TRUE---------------------------------
strand(gff) <- "*" 

## ----genomicranges_utilities3, eval=TRUE---------------------------------
reduce(gff) 

## ----genomicranges_utilities4, eval=TRUE---------------------------------
gaps(gff) 

## ----genomicranges_utilities4b, eval=TRUE--------------------------------
setdiff(as(seqinfo(gff), "GRanges"), gff) 

## ----genomicranges_utilities5, eval=TRUE---------------------------------
disjoin(gff)

## ----genomicranges_utilities6, eval=TRUE---------------------------------
coverage(gff)

## ----genomicranges_utilities7, eval=TRUE---------------------------------
findOverlaps(gff, gff[1:4]) 

## ----genomicranges_utilities8, eval=TRUE---------------------------------
countOverlaps(gff, gff[1:4])[1:40]

## ----genomicranges_utilities9, eval=TRUE---------------------------------
subsetByOverlaps(gff, gff[1:4]) 

## ----genomicrangeslist_objects, eval=TRUE--------------------------------
sp <- split(gff, seq(along=gff)) # Stores every range in separate component of a GRangesList object
split(gff, seqnames(gff)) # Stores ranges of each chromosome in separate component.
unlist(sp) # Returns data as GRanges object
sp[1:4, "type"] # Subsetting of GRangesList objects is similar to GRanges objects.
lapply(sp[1:4], length) # Looping over GRangesList objects similar to lists

## ----txdb_objects1, eval=TRUE, message=FALSE, warnings=FALSE-------------
library(GenomicFeatures)
download.file("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/gff3.gff", "data/gff3.gff")
txdb <- makeTxDbFromGFF(file="data/gff3.gff", format="gff", dataSource="TAIR", organism="Arabidopsis thaliana")
saveDb(txdb, file="./data/TAIR10.sqlite")
txdb <- loadDb("./data/TAIR10.sqlite")
transcripts(txdb)
transcriptsBy(txdb, by = "gene")
exonsBy(txdb, by = "gene")

## ----txdb_objects1_biomart, eval=FALSE, message=FALSE, warnings=FALSE----
## library(GenomicFeatures); library("biomaRt")
## txdb <- makeTxDbFromBiomart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org")

## ----biomart_basics, eval=FALSE, message=FALSE, warnings=FALSE-----------
## listMarts() # Lists BioMart databases
## listMarts(host="plants.ensembl.org")
## mymart <- useMart("plants_mart", host="plants.ensembl.org") # Select one, here plants_mart_25
## listDatasets(mymart) # List datasets available in the selected BioMart database
## mymart <- useMart("plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")
## listAttributes(mymart) # List available features
## getBM(attributes=c("ensembl_gene_id", "description"), mart=mymart)[1:4,]

## ----getseq_gff, eval=TRUE, message=FALSE, warnings=FALSE----------------
gff <- gff[values(gff)$type != "chromosome"] # Remove chromosome ranges
rand <- DNAStringSet(sapply(unique(as.character(seqnames(gff))), function(x) paste(sample(c("A","T","G","C"), 200000, replace=T), collapse="")))
writeXStringSet(DNAStringSet(rand), "./data/test")
getSeq(FaFile("./data/test"), gff)

## ----extractTranscritpSeqs, eval=FALSE, message=FALSE, warnings=FALSE----
## library(GenomicFeatures); library(Biostrings); library(Rsamtools)
## extractTranscriptSeqs(FaFile("mygenome.fasta"), exonsBy(txdb, "tx", use.names=TRUE))

## ----hw6a_demultiplex, eval=FALSE----------------------------------------
## demultiplex <- function(x, barcode, nreads) {
## 	f <- FastqStreamer(x, nreads)
## 	while(length(fq <- yield(f))) {
## 		for(i in barcode) {
## 			pattern <- paste("^", i, sep="")
## 			fqsub <- fq[grepl(pattern, sread(fq))]
## 			if(length(fqsub) > 0) {
## 				writeFastq(fqsub, paste(x, i, sep="_"), mode="a", compress=FALSE)
## 			}
## 		}
## 	}
## 	close(f)
## }
## demultiplex(x=fastq[1], barcode=c("TT", "AA", "GG"), nreads=50)

## ----hw6b_sequence_parsing, eval=FALSE-----------------------------------
## download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.gff", "data/AE004437.gff")
## download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.fna", "data/AE004437.fna")
## chr <- readDNAStringSet("data/AE004437.fna")
## gff <- import("data/AE004437.gff")
## gffgene <- gff[values(gff)[,"type"]=="gene"]
## gene <- DNAStringSet(Views(chr[[1]], IRanges(start(gffgene), end(gffgene))))
## names(gene) <- values(gffgene)[,"locus_tag"]
## pos <- values(gffgene[strand(gffgene) == "+"])[,"locus_tag"]
## p1 <- translate(gene[names(gene) %in% pos])
## names(p1) <- names(gene[names(gene) %in% pos])
## neg <- values(gffgene[strand(gffgene) == "-"])[,"locus_tag"]
## p2 <- translate(reverseComplement(gene[names(gene) %in% neg]))
## names(p2) <- names(gene[names(gene) %in% neg])
## writeXStringSet(c(p1, p2), "./data/mypep.fasta")

## ----sessionInfo, eval=TRUE----------------------------------------------
sessionInfo()

