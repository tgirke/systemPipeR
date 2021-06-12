## It will check the catDB class and methods
library(systemPipeR)
skip_on_bioc()
skip_on_ci()

test_that("check_catDB", {
    ## Obtain annotations from BioMart
    library("biomaRt")
    listMarts()  # To choose BioMart database
    listMarts(host = "plants.ensembl.org")
    m <- useMart("plants_mart", host = "plants.ensembl.org")
    listDatasets(m)
    m <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "plants.ensembl.org")
    listAttributes(m)  # Choose data types you want to download
    go <- getBM(attributes = c("go_id", "tair_locus", "namespace_1003"), mart = m)
    go <- go[go[, 3] != "", ]
    go[, 3] <- as.character(go[, 3])
    go[go[, 3] == "molecular_function", 3] <- "F"
    go[go[, 3] == "biological_process", 3] <- "P"
    go[go[, 3] == "cellular_component", 3] <- "C"
    go[1:4, ]
    dir.create("./data2/GO", recursive = TRUE)
    write.table(go, "data2/GO/GOannotationsBiomart_mod.txt", quote = FALSE, row.names = FALSE,
                col.names = FALSE, sep = "\t")
    ## Create catDB instance (takes a while but needs to be done only once)
    catdb <- makeCATdb(myfile = "data2/GO/GOannotationsBiomart_mod.txt", lib = NULL, org = "",
                       colno = c(1, 2, 3), idconv = NULL)
    catdb
    ## Class
    expect_s4_class(catdb, "catDB")
    expect_length(catdb, 1)
    ## Methods
     expect_type(names(catdb), "character")
     expect_type(catmap(catdb), "list")
     expect_type(catlist(catdb), "list")
     expect_type(idconv(catdb), "NULL")
})
