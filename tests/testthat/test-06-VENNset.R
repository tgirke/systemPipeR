## It will check the VENNset class and methods
library(systemPipeR)
skip_on_bioc()

test_that("check_VENNset", {
    ## Sample data
    setlist <- list(A=sample(letters, 18), B=sample(letters, 16),
                    C=sample(letters, 20), D=sample(letters, 22),
                    E=sample(letters, 18), F=sample(letters, 22))
    
    ## Create VENNset
    vennset <- overLapper(setlist[1:5], type="vennsets")
    ## Class
    expect_s4_class(vennset, "VENNset")
    expect_length(vennset, 5)
    ## Accessor methods for VENNset/INTERSECTset objects
    expect_type(names(vennset), "character")
    expect_type(setlist(vennset), "list")
    expect_type(intersectmatrix(vennset), "double")
    expect_type(complexitylevels(vennset), "integer")
    expect_type(vennlist(vennset), "list")
    ## Coerce VENNset/INTERSECTset object to list
    expect_type(as.list(vennset), "list")
    ## Bar plot of standard intersect counts
    interset <- overLapper(setlist, type="intersects")
    ## Class
    expect_s4_class(interset, "INTERSECTset")
    expect_length(interset, 6)
    plot <- olBarplot(interset, mincount=1)
    expect_s3_class(plot, "ggplot")
})
