## It will check the SYSargs2 class and methods
library(systemPipeR)
library(systemPipeRdata)

genWorkenvir("rnaseq", mydirname = file.path(tempdir(), "rnaseq"))
setwd(file.path(tempdir(), "rnaseq"))

test_that("check_SYSargs2_test", {
    skip_on_bioc()
    ## build instance 
    ## loadWorkflow() // renderWF()
    dir_path <- system.file("extdata/cwl/test/", package = "systemPipeR")
    dir.create("param/docopt.R/test", recursive = TRUE)
    file.copy(system.file("extdata/docopt.R/test/test.doc.R", package = "systemPipeR"), to = "param/docopt.R/test", recursive = TRUE)
    args <- loadWorkflow(targets = NULL, wf_file = "test.cwl",
                         input_file = "test.yml", dir_path = dir_path)
    args <- renderWF(args)
    ## Methods
    expect_type(names(args), "character")
    expect_type(targets(args), "list")
    expect_type(targetsheader(args), "list")
    expect_type(modules(args), "character")
    expect_type(wf(args), "list")
    expect_type(clt(args), "list")
    expect_type(yamlinput(args), "list")
    expect_type(cmdlist(args), "list")
    expect_type(input(args), "list")
    expect_type(output(args), "list")
    expect_type(files(args), "list")
    expect_type(inputvars(args), "list")
    expect_type(cmdToCwl(args), "list")
    expect_s4_class(args, "SYSargs2")
    expect_equal(baseCommand(args), "Rscript")
    expect_message(SampleName(args))
    ## runCommandline() //check.output()
    args <- runCommandline(args=args, force=TRUE, make_bam = FALSE)
    #args2 <- runCommandline(args=args, dir = TRUE, force=TRUE)
    out <- check.output(args)
    expect_setequal(out$Existing_Files, "1")
    ## symLink2bam()
    symLink2bam(sysargs=args, command="ln -s", htmldir=c(tempdir(), "/rnaseq/somedir/"), 
                ext=c(".bam", ".bai"), urlbase="http://myserver.edu/~username/", 
                urlfile="IGVurl.txt")
    expect_true(file.exists("somedir/.bam"))
    ## readComp()
    targetspath <- system.file("extdata", "targets.txt", package="systemPipeR")
    cmp <- readComp(targetspath, format="vector", delim="-")
    expect_type(cmp, "list")
    expect_error(readComp(args, format="vector", delim="-"))
})


# requires Hisat2 installed... 
test_that("check_SYSargs2_hisat2", {
    skip_on_bioc()
    skip_on_ci()
    ## build instance
    ## loadWorkflow() // renderWF()
    targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
    dir_path <- system.file("extdata/cwl", package="systemPipeR")
    WF <- loadWF(targets=targetspath, wf_file="hisat2/hisat2-mapping-pe.cwl",
                 input_file="hisat2/hisat2-mapping-pe.yml",
                 dir_path=dir_path)
    WF <- renderWF(WF, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
    WF <- WF[1:2]
    expect_s4_class(WF, "SYSargs2")
    ## Replacement Methods
    yamlinput(WF, "thread") <- 5L
    expect_equal(yamlinput(WF)$thread, 5)
    ## build index
    idx <- loadWorkflow(targets = NULL, wf_file = "hisat2/hisat2-index.cwl", input_file = "hisat2/hisat2-index.yml",
                        dir_path = dir_path)
    idx <- renderWF(idx)
    runCommandline(args = idx, make_bam = FALSE, dir=FALSE)
    ## Run alignment
    ## runCommandline() //check.output()
    WF <- runCommandline(args=WF, make_bam = TRUE, dir=FALSE)
    out <- check.output(WF)
    expect_equal(out$Existing_Files, c(2,2))
    ## alignStats()
    read_statsDF <- alignStats(WF, subset = "FileName") 
    expect_equal(read_statsDF$FileName, c("M1A", "M1B"))
    write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
})
    