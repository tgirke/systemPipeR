## It will check the SYSargs2 class and methods
library(systemPipeR)
library(systemPipeRdata)

genWorkenvir("rnaseq", mydirname = file.path(tempdir(), "rnaseq2"))
setwd(file.path(tempdir(), "rnaseq2"))

test_that("check_SYSargsList_test", {
    ## build instance
    dir_path <- system.file("extdata/cwl/example", package = "systemPipeR")
    sal <- SPRproject(overwrite = TRUE)
    targetspath <- system.file("extdata/cwl/example/targets_example.txt", package = "systemPipeR")
    ## appendStep <-  //SYSargsList() // cmdlist() // length(sal)
    appendStep(sal) <- SYSargsList(
        targets = targetspath, step_name = "example",
        wf_file = "example.cwl", input_file = "example.yml", dir_path = dir_path,
        inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_")
    )
    ## Methods
    expect_type(names(sal), "character")
    expect_type(length(sal), "integer")
    expect_type(stepsWF(sal), "list")
    expect_type(statusWF(sal), "list")
    expect_type(dependency(sal), "list")
    expect_type(projectInfo(sal), "list")
    expect_type(targetsWF(sal), "list")
    expect_type(SEobj(sal), "list")
    expect_type(outfiles(sal), "list")
    expect_type(cmdlist(sal), "list")
    expect_type(runInfo(sal), "list")
    expect_type(sysargslist(sal), "list")
    expect_error(SampleName(sal))
    baseCommand(sal)
    ## Class
    expect_s4_class(sal, "SYSargsList")
    expect_length(sal, 1)
    expect_length(cmdlist(sal)[[1]], 3)
    ## runWF()
    sal <- runWF(sal)
    expect_setequal(sal$statusWF[[1]]$status.summary, "Success")
    check <- check.output(sal)
    expect_setequal(check$example$Existing_Files, "1")
    ## replacement methods
    renameStep(sal, 1) <- "newStep"
    expect_error(appendStep(sal) <- sal)
    appendStep(sal) <- LineWise(
        {
            a <- log(-1)
        },
        step_name = "R_code",
        dependency = "newStep"
    )
    sal <- runWF(sal)
    ## `+` method
    sal <- sal[1] + sal[2]
    expect_length(sal, 2)
})

# requires trimmomatic/Hisat2 installed...
test_that("check_sal_hisat2", {
    skip_on_bioc()
    skip_on_ci()
    ## build instance
    sal2 <- SPRproject(overwrite = TRUE)
    targetspath <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
    appendStep(sal2) <- SYSargsList(
        targets = targetspath,
        step_name = "Quality",
        wf_file = "trimmomatic/workflow_trimmomatic-pe.cwl",
        input_file = "trimmomatic/trimmomatic-pe.yml",
        dir_path = system.file("extdata/cwl", package = "systemPipeR"),
        inputvars = c(FileName1 = "_FASTQ_PATH1_", FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_")
    )
    appendStep(sal2) <- SYSargsList(
        step_name = "Index", dir = FALSE,
        targets = NULL,
        wf_file = "hisat2/workflow_hisat2-index.cwl",
        input_file = "hisat2/hisat2-index.yml",
        dir_path = system.file("extdata/cwl", package = "systemPipeR"),
        dependency = "Quality"
    )
    appendStep(sal2) <- SYSargsList(
        targets = "Quality",
        step_name = "Mapping", dir = TRUE,
        wf_file = "workflow-hisat2/workflow_hisat2-pe.cwl",
        input_file = "workflow-hisat2/workflow_hisat2-pe.yml",
        dir_path = system.file("extdata/cwl", package = "systemPipeR"),
        inputvars = c(
            trimmomatic_1_paired = "_FASTQ_PATH1_",
            trimmomatic_2_paired = "_FASTQ_PATH2_", SampleName = "_SampleName_"
        ),
        rm_targets_col = c("FileName1", "FileName2"),
        dependency = c("Quality", "Index")
    )
    sal2 <- subset(sal2, subset_steps = c(1, 3), input_targets = 1:4, keep_steps = TRUE)
    expect_s4_class(sal2, "SYSargsList")
    ## Replacement Methods
    yamlinput(sal2, step = 1, paramName = "thread") <- 5L
    expect_equal(yamlinput(sal2, 1)$thread, 5)
    ## Run alignment
    ## runWF() // check.output()
    sal2 <- runWF(sal2)
    out <- check.output(sal2)
    expect_equal(sum(out$Quality$Existing_Files), 16)
})
