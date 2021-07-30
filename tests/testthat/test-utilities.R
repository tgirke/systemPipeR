library(systemPipeR)
skip_on_bioc()
skip_on_ci()

test_that("check_evalCode", {
    rmdPath <- system.file("extdata", "SPRtest.Rmd", package="systemPipeR")
    output <- evalCode(rmdPath, eval = TRUE, output="test.Rmd")
    expect_equal(output, "test.Rmd")
})