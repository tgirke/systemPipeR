skip_on_bioc()
skip_on_ci()

library(systemPipeR)
test_that("check_evalCode", {
    rmdPath <- system.file("extdata/", "systemPipeTEST.Rmd", package="systemPipeR")
    evalCode(rmdPath, eval = TRUE, output="test.Rmd")
})