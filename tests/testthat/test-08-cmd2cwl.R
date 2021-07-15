## It will check the SYSargs class and methods
skip_on_bioc()
skip_on_ci()

## It will check the cmd2cwl class and methods
library(systemPipeR)

test_that("check_LineWise_test", {
    ## build simple example
    command <- "
    hisat2 \
    -S <F, out: ./results/M1A.sam> \
    -x <F: ./data/tair10.fasta> \
     -k <int: 1> \
    -min-intronlen <int: 30> \
    -max-intronlen <int: 3000> \
    -threads <int: 4> \
    -U <F: ./data/SRR446027_1.fastq.gz>
"
    cmd <- createParamFiles(command, writeParamFiles = FALSE)
    expect_s4_class(cmd, "SYSargs2")
    writeParamFiles(cmd, overwrite = TRUE)
    expect_true(file.exists("param/cwl/hisat2/hisat2.yml"))
    expect_true(file.exists("param/cwl/hisat2/hisat2.cwl"))
    ## printParam
    print_cmd <- printParam(cmd, position = "baseCommand")
    expect_s3_class(print_cmd, "cwlParse")
    ## subsetParam
    cmd2 <- subsetParam(cmd, position = "inputs", index = 1, trim = TRUE)
    expect_equal(cmdlist(cmd2)[[1]][[1]], "hisat2 -S ./results/M1A.sam")
    ## renameParam
    cmd3 <- renameParam(cmd, "input", index = 1, rename = "sam")
    expect_equal(names(cmd3$cmdToCwl$inputs)[1], "sam")
    ## replaceParam
    new_inputs <- new_inputs <- list(
        "new_input1" = list(type = "File", preF = "-b", yml = "myfile"),
        "new_input2" = "-L <int: 4>"
    )
    cmd4 <- replaceParam(cmd, "inputs", index = 1:2, replace = new_inputs)
    expect_equal(names(cmd4$cmdToCwl$inputs)[1], "new_input1")
    ## appendParam
    newIn <- new_inputs <- list(
        "new_input1" = list(type = "File", preF = "-b1", yml = "myfile1"),
        "new_input2" = list(type = "File", preF = "-b2", yml = "myfile2"),
        "new_input3" = "-b3 <F: myfile3>"
    )
    cmd5 <- appendParam(cmd, "inputs", index = 1:2, append = new_inputs)
    expect_equal(names(cmd5$cmdToCwl$inputs)[10], "new_input3")
})
