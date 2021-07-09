## It will check the SYSargs2 class and methods
library(systemPipeR)
setwd(file.path(tempdir(), "rnaseq"))

test_that("check_SYSargsList_test", {
    ## build instance 
    dir_path <- system.file("extdata/cwl/example", package="systemPipeR")
    sal <- SPRproject(overwrite = TRUE)
    targetspath <- system.file("extdata/cwl/example/targets_example.txt", package="systemPipeR")
    ## appendStep <-  //SYSargsList() // cmdlist() // length(sal)
    appendStep(sal) <- SYSargsList(targets=targetspath, step_name="example", 
                                   wf_file="example.cwl", input_file="example.yml", dir_path = dir_path, 
                                   inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
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
    expect_type(sysargslist(sal), "list")
    ## Class
    expect_s4_class(sal, "SYSargsList")
    expect_length(sal, 1)
    expect_length(cmdlist(sal)[[1]], 3)
    ##runWF()
   # sal <- runWF(sal)
   # expect_setequal(statusWF(sal)[[1]][[1]], "DONE")
   # check <- check.output(sal)
   # expect_setequal(check$Step_1$Existing_Files, "1")
   # 
   ## replacement methods
    renameStep(sal, 1) <- "newStep"
   expect_error(appendStep(sal) <- sal)
   appendStep(sal) <- LineWise({
                                a <- log(-1)
                                }, step_name = "R_code", 
                               dependency= "newStep")

   ## `+` method
   sal <- sal[1] + sal[2] 
   expect_length(sal, 2)
})

