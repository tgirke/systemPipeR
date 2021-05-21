#devtools::test()

library(systemPipeR)
library(systemPipeRdata)

genWorkenvir("rnaseq", mydirname = file.path(tempdir(), "rnaseq"))
setwd(file.path(tempdir(), "rnaseq"))

# test_that("check_SYSargs2", {
#     ## build instance
#     ## loadWorkflow() // renderWF()
#     targets <- system.file("extdata", "targets.txt", package = "systemPipeR")
#     dir_path <- system.file("extdata/cwl/rsubread/rsubread-se", package = "systemPipeR")
#     args <- loadWorkflow(targets = targets, wf_file = "rsubread-mapping-se.cwl",
#                              input_file = "rsubread-mapping-se.yml", dir_path = dir_path)
#     args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
#     args <- args[1:2]
#     expect_s4_class(args, "SYSargs2")
#     ## build index
#     dir_path <- system.file("extdata/cwl/rsubread/rsubread-idx", package = "systemPipeR")
#     idx <- loadWorkflow(targets = NULL, wf_file = "rsubread-index.cwl", input_file = "rsubread-index.yml",
#                         dir_path = dir_path)
#     idx <- renderWF(idx)
#     runCommandline(args = idx, make_bam = FALSE)
#     ## Run alignment
#     ## runCommandline() //check.output()
#     args <- runCommandline(args=args)
#     out <- check.output(args)
#     expect_type(out, "logical")
#     # symLink2bam()
#     symLink2bam(sysargs=args, command="ln -s", htmldir=c(tempdir(), "/rnaseq/somedir/"),
#                 ext=c(".bam", ".bai"), urlbase="http://myserver.edu/~username/",
#                 urlfile="IGVurl.txt")
#     expect_true(file.exists("somedir/M1A.bam"))
# })

test_that("check_SYSargs2_test", {
    ## build instance 
    ## loadWorkflow() // renderWF()
    dir_path <- system.file("extdata/cwl/test/", package = "systemPipeR")
    dir.create("param/docopt.R/test")
    file.copy(system.file("extdata/docopt.R/test/test.doc.R", package = "systemPipeR"), to = "param/docopt.R/test", recursive = TRUE)
    args <- loadWorkflow(targets = NULL, wf_file = "test.cwl",
                         input_file = "test.yml", dir_path = dir_path)
    args <- renderWF(args)
    expect_s4_class(args, "SYSargs2")
    ## runCommandline() //check.output()
    args <- runCommandline(args=args)
    out <- check.output(args)
    expect_type(out, "logical")
    ## symLink2bam()
    symLink2bam(sysargs=args, command="ln -s", htmldir=c(tempdir(), "/rnaseq/somedir/"), 
                ext=c(".bam", ".bai"), urlbase="http://myserver.edu/~username/", 
                urlfile="IGVurl.txt")
    expect_true(file.exists("somedir/.bam"))
})
