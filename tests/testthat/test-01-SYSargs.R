## It will check the SYSargs class and methods
skip_on_bioc()
skip_on_ci()

library(systemPipeR)
library(systemPipeRdata)

genWorkenvir("chipseq", mydirname = file.path(tempdir(), "chipseq"))
setwd(file.path(tempdir(), "chipseq"))

test_that("check_SYSargs", {
    ## build instance 
    ## systemArgs() // Methods
    param <- system.file("extdata", "rsubread.param", package="systemPipeR")
    targets <- system.file("extdata", "targets.txt", package="systemPipeR")
    args <- systemArgs(sysma=param, mytargets=targets)
    expect_type(names(args), "character")
    expect_type(modules(args), "character")
    expect_type(cores(args), "double")
    expect_type(outpaths(args), "character")
    expect_type(sysargs(args), "character")
    expect_type(reference(args), "character")
    expect_type(infile1(args), "character")
    expect_type(infile2(args), "character")
    expect_type(outfile1(args), "character")
    expect_type(software(args), "character")
    expect_type(results(args), "character")
    expect_type(targetsin(args), "list")
    expect_type(targetsout(args), "list")
    expect_type(targetsheader(args), "character")
    expect_s4_class(args, "SYSargs")
    ## reference // infile1 // outfile1
    library(Rsubread)
    buildindex(basename = reference(args), reference = reference(args))
    args <- args[1]
    align(index = reference(args), readfile1 = infile1(args), input_format = "FASTQ",
          output_file = outfile1(args), output_format = "SAM", nthreads = 8, indels = 1, TH1 = 2)
    for(i in seq(along=outfile1(args))) asBam(file=outfile1(args)[i],
                                              destination=gsub(".sam", "", outfile1(args)[i]), overwrite=TRUE, indexDestination=TRUE)
    expect_true(file.exists(outpaths(args)))
    # requires Hisat2 installed... 
    dir_path <- system.file("extdata/cwl", package = "systemPipeR")
    idx <- loadWorkflow(targets = NULL, wf_file = "hisat2/hisat2-index.cwl", input_file = "hisat2/hisat2-index.yml",
                        dir_path = dir_path)
    idx <- renderWF(idx)
    cmdlist(idx)
    ## runCommandline() 
    runCommandline(idx, make_bam = FALSE, dir=FALSE, force = TRUE)
    ## Construct SYSargs object from param and targets files 
    param <- system.file("extdata", "hisat2.param", package="systemPipeR")
    targets <- system.file("extdata", "targets.txt", package="systemPipeR")
    args <- systemArgs(sysma=param, mytargets=targets)
    runCommandline(args[2], make_bam = TRUE, dir=FALSE)
    expect_s4_class(args, "SYSargs")
    expect_true(file.exists(outpaths(args[2])))
    ## alignStats()
    read_statsDF <- alignStats(args[2]) 
    expect_s3_class(read_statsDF, "data.frame")
    write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
})
