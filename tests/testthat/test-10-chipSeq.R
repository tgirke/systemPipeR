## It will check the ChipSeq.R functions
library(systemPipeR)
library(systemPipeRdata)
skip_on_bioc()
skip_on_ci()

genWorkenvir("rnaseq", mydirname = file.path(tempdir(), "chipseq"))
setwd(file.path(tempdir(), "chipseq"))

test_that("check_chipSeq_fnc", {
    dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-idx", package = "systemPipeR")
    idx <- loadWorkflow(targets = NULL, wf_file = "bowtie2-index.cwl", 
                        input_file = "bowtie2-index.yml", dir_path = dir_path)
    idx <- renderWF(idx)

    ## Run in single machine
    runCommandline(idx, make_bam = FALSE, dir=FALSE)
    
    targets <- system.file("extdata", "targetsPE_chip.txt", package = "systemPipeR")
    dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-pe", package = "systemPipeR")
    args <- loadWF(targets = targets, wf_file = "bowtie2-mapping-pe.cwl", 
                   input_file = "bowtie2-mapping-pe.yml", dir_path = dir_path)
    args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
                                         FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
    args <- args[1:4]
    args <- runCommandline(args, make_bam = TRUE, dir=FALSE, force = F)
    
    writeTargetsout(x = args, file = "targets_bam.txt", step = 1, 
                    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE, 
                    remove = TRUE)
    outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
    file.exists(outpaths)
    
    dir_path <- system.file("extdata/cwl/mergeBamByFactor", package = "systemPipeR")
    args <- loadWF(targets = "targets_bam.txt", wf_file = "merge-bam.cwl", 
                   input_file = "merge-bam.yml", dir_path = dir_path)
    args <- renderWF(args, inputvars = c(FileName = "_BAM_PATH_", 
                                         SampleName = "_SampleName_"))
    
    args_merge <- mergeBamByFactor(args = args, overwrite = TRUE)
    writeTargetsout(x = args_merge, file = "targets_mergeBamByFactor.txt", 
                    step = 1, new_col = "FileName", new_col_output_index = 1, 
                    overwrite = TRUE, remove = TRUE)
    
    dir_path <- system.file("extdata/cwl/MACS2/MACS2-noinput/", package = "systemPipeR")
    args <- loadWF(targets = "targets_mergeBamByFactor.txt", wf_file = "macs2.cwl", 
                   input_file = "macs2.yml", dir_path = dir_path)
    args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
                                         SampleName = "_SampleName_"))
    
    args <- runCommandline(args, make_bam = FALSE, dir=FALSE, force = TRUE)
    outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
    file.exists(outpaths)
    writeTargetsout(x = args, file = "targets_macs.txt", step = 1, 
                    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
    
    writeTargetsRef(infile = "targets_mergeBamByFactor.txt", outfile = "targets_bam_ref.txt", 
                    silent = FALSE, overwrite = TRUE)
    dir_path <- system.file("extdata/cwl/MACS2/MACS2-input/", package = "systemPipeR")
    args_input <- loadWF(targets = "targets_bam_ref.txt", wf_file = "macs2-input.cwl", 
                         input_file = "macs2.yml", dir_path = dir_path)
    args_input <- renderWF(args_input, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
                                                     FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
    cmdlist(args_input)[1]
    ## Run
    args_input <- runCommandline(args_input,make_bam = FALSE, dir=FALSE, force = T)
    outpaths_input <- subsetWF(args_input, slot = "output", subset = 1, 
                               index = 1)
    file.exists(outpaths_input)
    writeTargetsout(x = args_input, file = "targets_macs_input.txt", 
                    step = 1, new_col = "FileName", new_col_output_index = 1, 
                    overwrite = TRUE)
    
    library(ChIPpeakAnno)
    library(GenomicFeatures)
    dir_path <- system.file("extdata/cwl/annotate_peaks", package = "systemPipeR")
    args <- loadWF(targets = "targets_macs.txt", wf_file = "annotate_peaks.cwl", 
                   input_file = "annotate_peaks.yml", dir_path = dir_path)
    args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
                                         SampleName = "_SampleName_"))
    
    expect_warning(txdb <- makeTxDbFromGFF(file = "data/tair10.gff", format = "gff", 
                            dataSource = "TAIR", organism = "Arabidopsis thaliana"))
    ge <- genes(txdb, columns = c("tx_name", "gene_id", "tx_type"))
    for (i in seq(along = args)) {
        peaksGR <- as(read.delim(infile1(args)[i], comment = "#"), 
                      "GRanges")
        annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData = genes(txdb))
        df <- data.frame(as.data.frame(annotatedPeak), as.data.frame(values(ge[values(annotatedPeak)$feature, 
        ])))
        outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
        write.table(df, outpaths[i], quote = FALSE, row.names = FALSE, 
                    sep = "\t")
    }
    writeTargetsout(x = args, file = "targets_peakanno.txt", step = 1, 
                    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
    
    library(GenomicRanges)
    dir_path <- system.file("extdata/cwl/count_rangesets", package = "systemPipeR")
    args <- loadWF(targets = "targets_macs.txt", wf_file = "count_rangesets.cwl", 
                   input_file = "count_rangesets.yml", dir_path = dir_path)
    args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
                                         SampleName = "_SampleName_"))
    
    ## Bam Files
    targets <- system.file("extdata", "targetsPE_chip.txt", package = "systemPipeR")
    dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-pe", package = "systemPipeR")
    args_bam <- loadWF(targets = targets, wf_file = "bowtie2-mapping-pe.cwl", 
                       input_file = "bowtie2-mapping-pe.yml", dir_path = dir_path)
    args_bam <- renderWF(args_bam, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
                                                 SampleName = "_SampleName_"))
    args_bam <- output_update(args_bam, dir = FALSE, replace = TRUE, 
                              extension = c(".sam", ".bam"))
    outpaths <- subsetWF(args_bam, slot = "output", subset = 1, index = 1)
    
    bfl <- BamFileList(outpaths[1:4], yieldSize = 50000, index = character())
    countDFnames <- countRangeset(bfl, args, mode = "Union", ignore.strand = TRUE)
    writeTargetsout(x = args, file = "targets_countDF.txt", step = 1, 
                    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
    
    dir_path <- system.file("extdata/cwl/rundiff", package = "systemPipeR")
    args_diff <- loadWF(targets = "targets_countDF.txt", wf_file = "rundiff.cwl", 
                        input_file = "rundiff.yml", dir_path = dir_path)
    args_diff <- renderWF(args_diff, inputvars = c(FileName = "_FASTQ_PATH1_", 
                                                   SampleName = "_SampleName_"))
    
    cmp <- readComp(file = args_bam, format = "matrix")
    cmp[[1]] <- cmp[[1]][1,]
    args_bam <- args_bam[1:4]
    dbrlist <- runDiff(args = args_diff, diffFct = run_edgeR, targets = targets.as.df(targets(args_bam)), 
                       cmp = cmp[[1]], independent = TRUE, dbrfilter = c(Fold = 2, 
                                                                         FDR = 1))
    writeTargetsout(x = args_diff, file = "targets_rundiff.txt", 
                    step = 1, new_col = "FileName", new_col_output_index = 1, 
                    overwrite = TRUE)
    
    library(Biostrings)
    library(seqLogo)
    library(BCRANK)
    dir_path <- system.file("extdata/cwl/annotate_peaks", package = "systemPipeR")
    args <- loadWF(targets = "targets_macs.txt", wf_file = "annotate_peaks.cwl", 
                   input_file = "annotate_peaks.yml", dir_path = dir_path)
    args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
                                         SampleName = "_SampleName_"))
    
    rangefiles <- infile1(args)
    for (i in seq(along = rangefiles)) {
        df <- read.delim(rangefiles[i], comment = "#")
        peaks <- as(df, "GRanges")
        names(peaks) <- paste0(as.character(seqnames(peaks)), "_", 
                               start(peaks), "-", end(peaks))
        peaks <- peaks[order(values(peaks)$X.log10.pvalue., decreasing = TRUE)]
        pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
        names(pseq) <- names(peaks)
        writeXStringSet(pseq, paste0(rangefiles[i], ".fasta"))
    }
    
})

