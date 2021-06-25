## It will check the VENNset class and methods
library(systemPipeR)
skip_on_bioc()
skip_on_ci()
setwd(file.path(tempdir(), "rnaseq"))

test_that("check_features", {
    ## Create paths
    file <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
    txdb <- makeTxDbFromGFF(file=file, format="gff3", organism="Arabidopsis")
    
    targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
    dir_path <- system.file("extdata/cwl", package="systemPipeR")
    args <- loadWorkflow(targets=targetspath, wf_file="hisat2/hisat2-mapping-se.cwl", 
                         input_file="hisat2/hisat2-mapping-se.yml", dir_path=dir_path)
    args <- renderWF(args, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
    expect_s4_class(args, "SYSargs2")
    args <- runCommandline(args, make_bam = TRUE, dir = TRUE)
    outpaths <- subsetWF(args , slot="output", subset=1, index=1)
    file.exists(outpaths)
    
    ## genFeatures
    file <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
    txdb <- makeTxDbFromGFF(file=file, format="gff3", organism="Arabidopsis")
    feat <- genFeatures(txdb, featuretype="all", reduce_ranges=TRUE, upstream=1000, downstream=0, verbose=TRUE)
    
    ## featuretypeCounts
    featureCounts <- featuretypeCounts(bfl=BamFileList(outpaths, yieldSize=50000), grl=feat, singleEnd=TRUE, readlength=c(74:76,99:102), type="data.frame")
    write.table(featureCounts, "results/featureCounts.xls", quote=FALSE, row.names=FALSE, sep="\t")
    featureCounts2 <- featuretypeCounts(bfl=BamFileList(outpaths), grl=feat, singleEnd=TRUE, readlength=NULL, type="data.frame")
    write.table(featureCounts2, "results/featureCounts2.xls", quote=FALSE, row.names=FALSE, sep="\t")
    
    ## plotfeaturetypeCounts
    featureCounts <- read.delim("results/featureCounts.xls", check.names=FALSE)
    myplots <- plotfeaturetypeCounts(x=featureCounts, graphicsfile="results/featureCounts.pdf", graphicsformat="pdf", scales="fixed", anyreadlength=FALSE, drop_N_total_aligned=TRUE, scale_count_val=10^6, scale_length_val=10^3)
    featureCounts2 <- read.delim("results/featureCounts2.xls", check.names=FALSE)
    myplots <- plotfeaturetypeCounts(x=featureCounts2, graphicsfile="results/featureCounts2.pdf", graphicsformat="pdf", scales="fixed", anyreadlength=TRUE, drop_N_total_aligned=TRUE, scale_count_val=10^6, scale_length_val=10^3)
    ## featureCoverage'
    grl <- cdsBy(txdb, "tx", use.names=TRUE)
    fcov <- featureCoverage(bfl=BamFileList(outpaths[1]), grl=grl[1:4], resizereads=NULL, readlengthrange=NULL, Nbins=NULL, method=mean, fixedmatrix=TRUE, resizefeatures=TRUE, upstream=20, downstream=20, outfile="results/zzz.xls", overwrite=TRUE)
    ## plotfeatureCoverage
    grl <- cdsBy(txdb, "tx", use.names=TRUE)
    fcov <- featureCoverage(bfl=BamFileList(outpaths[1:2]), grl=grl[1:4], resizereads=NULL, readlengthrange=NULL, Nbins=20, method=mean, fixedmatrix=TRUE, resizefeatures=TRUE, upstream=20, downstream=20, outfile="results/zzz.xls", overwrite=TRUE)
    plotfeatureCoverage(covMA=fcov, method=mean, scales="fixed", extendylim=2, scale_count_val=10^6)
    ## predORF
    file <- system.file("extdata", "someORF.fa", package="Biostrings")
    dna <- readDNAStringSet(file)
    orf <- predORF(dna[1:4], n=1, type="df", mode="orf", strand="antisense", startcodon="ATG", stopcodon=c("TAA", "TAG", "TGA"))
    ## scaleRanges
    gff <- system.file("extdata/annotation", "tair10.gff", package="systemPipeRdata")
    txdb <- makeTxDbFromGFF(file=gff, format="gff3", organism="Arabidopsis")
    futr <- fiveUTRsByTranscript(txdb, use.names=TRUE)
    genome <- system.file("extdata/annotation", "tair10.fasta", package="systemPipeRdata")
    dna <- extractTranscriptSeqs(FaFile(genome), futr)
    uorf <- predORF(dna, n="all", mode="orf", longest_disjoint=TRUE, strand="sense")
    grl_scaled <- scaleRanges(subject=futr, query=uorf, type="uORF", verbose=TRUE)
    rtracklayer::export.gff3(unlist(grl_scaled), "uorf.gff")
    })
