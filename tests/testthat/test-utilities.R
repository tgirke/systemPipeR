library(systemPipeR)
skip_on_bioc()
skip_on_ci()

test_that("check_evalCode", {
    rmdPath <- system.file("extdata", "SPRtest.Rmd", package="systemPipeR")
    output <- evalCode(rmdPath, eval = TRUE, output="test.Rmd")
    expect_equal(output, "test.Rmd")
})

## It will check the utilities.R functions
library(systemPipeRdata)

genWorkenvir("rnaseq", mydirname = file.path(tempdir(), "rnaseq_test"))
setwd(file.path(tempdir(), "rnaseq_test"))
test_that("check_SYSargs2_test", {
    skip_on_bioc()
    ## Preprocessing of single-end reads
    dir_path <- system.file("extdata/cwl/preprocessReads/trim-se", package="systemPipeR")
    targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
    trim <- loadWorkflow(targets=targetspath, wf_file="trim-se.cwl", input_file="trim-se.yml", dir_path=dir_path)
    trim <- renderWF(trim, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))
    preprocessReads(args=trim[1], Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)
    expect_true(file.exists(trim$output[[1]][[1]]))
    ## Preprocessing of paired-end reads
    dir_path <- system.file("extdata/cwl/preprocessReads/trim-pe", package="systemPipeR")
    targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
    trim <- loadWorkflow(targets=targetspath, wf_file="trim-pe.cwl", input_file="trim-pe.yml", dir_path=dir_path)
    trim <- renderWF(trim, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
    preprocessReads(args=trim[1], Fct="trimLRPatterns(Rpattern='GCCCGGGTAA', subject=fq)", batchsize=100000, overwrite=TRUE, compress=TRUE)
    expect_true(file.exists(trim$output[[1]][[1]][1]))
    expect_true(file.exists(trim$output[[1]][[1]][2]))
})