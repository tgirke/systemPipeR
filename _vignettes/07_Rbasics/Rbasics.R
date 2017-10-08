## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")),
    warning=FALSE, message=FALSE)

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(limma) 
    library(ggplot2) }) 

## ----install_cran, eval=FALSE--------------------------------------------
## install.packages(c("pkg1", "pkg2"))
## install.packages("pkg.zip", repos=NULL)

## ----install_bioc, eval=FALSE--------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## library(BiocInstaller)
## BiocVersion()
## biocLite()
## biocLite(c("pkg1", "pkg2"))

## ----closing_r, eval=FALSE-----------------------------------------------
## q()

## ----r_assignment, eval=FALSE--------------------------------------------
## object <- ...

## ----r_ls, eval=FALSE----------------------------------------------------
## ls()

## ----r_dirshow, eval=FALSE-----------------------------------------------
## dir()

## ----r_dirpath, eval=FALSE-----------------------------------------------
## getwd()

## ----r_setwd, eval=FALSE-------------------------------------------------
## setwd("/home/user")

## ----r_syntax, eval=FALSE------------------------------------------------
## object <- function_name(arguments)
## object <- object[arguments]

## ----r_find_help, eval=FALSE---------------------------------------------
## ?function_name

## ----r_package_load, eval=FALSE------------------------------------------
## library("my_library")

## ----r_package_functions, eval=FALSE-------------------------------------
## library(help="my_library")

## ----r_load_vignette, eval=FALSE-----------------------------------------
## vignette("my_library")

## ----r_execute_script, eval=FALSE----------------------------------------
## source("my_script.R")

## ----sh_execute_script, eval=FALSE, engine="sh"--------------------------
## $ Rscript my_script.R
## $ R CMD BATCH my_script.R
## $ R --slave < my_script.R

## ----r_numeric_data, eval=TRUE-------------------------------------------

x <- c(1, 2, 3)
x
is.numeric(x)
as.character(x)

## ----r_character_data, eval=TRUE-----------------------------------------
x <- c("1", "2", "3")
x
is.character(x)
as.numeric(x)

## ----r_complex_data, eval=TRUE-------------------------------------------
c(1, "b", 3)

## ----r_logical_data, eval=TRUE-------------------------------------------
x <- 1:10 < 5
x  
!x
which(x) # Returns index for the 'TRUE' values in logical vector

## ----r_vector_object, eval=TRUE------------------------------------------
myVec <- 1:10; names(myVec) <- letters[1:10]  
myVec[1:5]
myVec[c(2,4,6,8)]
myVec[c("b", "d", "f")]

## ----r_factor_object, eval=TRUE------------------------------------------
factor(c("dog", "cat", "mouse", "dog", "dog", "cat"))

## ----r_matrix_object, eval=TRUE------------------------------------------
myMA <- matrix(1:30, 3, 10, byrow = TRUE) 
class(myMA)
myMA[1:2,]
myMA[1, , drop=FALSE]

## ----r_dataframe_object, eval=TRUE---------------------------------------
myDF <- data.frame(Col1=1:10, Col2=10:1) 
myDF[1:2, ]

## ----r_list_object, eval=TRUE--------------------------------------------
myL <- list(name="Fred", wife="Mary", no.children=3, child.ages=c(4,7,9)) 
myL
myL[[4]][1:2] 

## ----r_function_object, eval=FALSE---------------------------------------
## myfct <- function(arg1, arg2, ...) {
## 	function_body
## }

## ----r_subset_by_index, eval=TRUE----------------------------------------
myVec <- 1:26; names(myVec) <- LETTERS 
myVec[1:4]

## ----r_subset_by_logical, eval=TRUE--------------------------------------
myLog <- myVec > 10
myVec[myLog] 

## ----r_subset_by_names, eval=TRUE----------------------------------------
myVec[c("B", "K", "M")]

## ----r_subset_by_dollar, eval=TRUE---------------------------------------
iris$Species[1:8]

## ----r_combine_vectors, eval=TRUE----------------------------------------
c(1, 2, 3)
x <- 1:3; y <- 101:103
c(x, y)
iris$Species[1:8]

## ----r_cbind_rbind, eval=TRUE--------------------------------------------
ma <- cbind(x, y)
ma
rbind(ma, ma)

## ----r_length_dim, eval=TRUE---------------------------------------------
length(iris$Species)
dim(iris)

## ----col_row_names, eval=TRUE--------------------------------------------
rownames(iris)[1:8]
colnames(iris)

## ----name_slots, eval=TRUE-----------------------------------------------
names(myVec)
names(myL)

## ----sort_objects, eval=TRUE---------------------------------------------
sort(10:1)

## ----order_objects, eval=TRUE--------------------------------------------
sortindex <- order(iris[,1], decreasing = FALSE)
sortindex[1:12]
iris[sortindex,][1:2,]
sortindex <- order(-iris[,1]) # Same as decreasing=TRUE

## ----order_columns, eval=TRUE--------------------------------------------
iris[order(iris$Sepal.Length, iris$Sepal.Width),][1:2,]

## ----comparison_operators, eval=TRUE-------------------------------------
1==1

## ----logical_operators, eval=TRUE----------------------------------------
x <- 1:10; y <- 10:1
x > y & x > 5

## ----logical_calculations, eval=TRUE-------------------------------------
x + y
sum(x)
mean(x)
apply(iris[1:6,1:3], 1, mean) 

## ----read_delim, eval=FALSE----------------------------------------------
## myDF <- read.delim("myData.xls", sep="\t")

## ----read_excel, eval=FALSE----------------------------------------------
## library(gdata)
## myDF <- read.xls"myData.xls")

## ----read_gs, eval=FALSE-------------------------------------------------
## library("googlesheets"); library("dplyr"); library(knitr)
## gs_auth() # Creates authorizaton token (.httr-oauth) in current directory if not present
## sheetid <-"1U-32UcwZP1k3saKeaH1mbvEAOfZRdNHNkWK2GI1rpPM"
## gap <- gs_key(sheetid)
## mysheet <- gs_read(gap, skip=4)
## myDF <- as.data.frame(mysheet)
## myDF

## ----write_table, eval=FALSE---------------------------------------------
## write.table(myDF, file="myfile.xls", sep="\t", quote=FALSE, col.names=NA)

## ----readlines, eval=FALSE-----------------------------------------------
## myDF <- readLines("myData.txt")

## ----writelines, eval=FALSE----------------------------------------------
## writeLines(month.name, "myData.txt")

## ----paste_windows, eval=FALSE-------------------------------------------
## read.delim("clipboard")

## ----paste_osx, eval=FALSE-----------------------------------------------
## read.delim(pipe("pbpaste"))

## ----copy_windows, eval=FALSE--------------------------------------------
## write.table(iris, "clipboard", sep="\t", col.names=NA, quote=F)

## ----copy_osx, eval=FALSE------------------------------------------------
## zz <- pipe('pbcopy', 'w')
## write.table(iris, zz, sep="\t", col.names=NA, quote=F)
## close(zz)

## ----unique, eval=TRUE---------------------------------------------------
length(iris$Sepal.Length)
length(unique(iris$Sepal.Length))

## ----table, eval=TRUE----------------------------------------------------
table(iris$Species)

## ----aggregate, eval=TRUE------------------------------------------------
aggregate(iris[,1:4], by=list(iris$Species), FUN=mean, na.rm=TRUE)

## ----intersect, eval=TRUE------------------------------------------------
month.name %in% c("May", "July")

## ----merge, eval=TRUE----------------------------------------------------
frame1 <- iris[sample(1:length(iris[,1]), 30), ]
frame1[1:2,]
dim(frame1)
my_result <- merge(frame1, iris, by.x = 0, by.y = 0, all = TRUE)
dim(my_result)

## ----tidyverse_install, eval=FALSE---------------------------------------
## install.packages("tidyverse")

## ----data_frame_tbl1, eval=TRUE------------------------------------------
library(tidyverse)
as_data_frame(iris) # coerce data.frame to data frame tbl

## ----data_frame_tbl2, eval=FALSE-----------------------------------------
## as_tibble(iris) # newer function provided by tibble package
## tbl_df(iris) # this alternative exists for historical reasons

## ----tabular_sample, eval=TRUE-------------------------------------------
write_tsv(iris, "iris.txt") # Creates sample file

## ----tabular_import1, eval=TRUE------------------------------------------
iris_df <- read_tsv("iris.txt") # Import with read_tbv from readr package
iris_df

## ----tabular_import2, eval=TRUE------------------------------------------
library(data.table)
iris_df <- as_data_frame(fread("iris.txt")) # Import with fread and conversion to tibble
iris_df

## ----tabular_import_ignore, eval=FALSE-----------------------------------
## fread("grep -v '^#' iris.txt")

## ----tabular_export_readr, eval=FALSE------------------------------------
## write_tsv(iris_df, "iris.txt")

## ----dplyr_bind, eval=TRUE-----------------------------------------------
bind_cols(iris_df, iris_df)
bind_rows(iris_df, iris_df)

## ----plyr_get_cols, eval=TRUE--------------------------------------------
iris_df[[5]][1:12]
iris_df$Species[1:12]

## ----plyr_filter, eval=TRUE----------------------------------------------
filter(iris_df, Sepal.Length > 7.5, Species=="virginica")

## ----plyr_filter_base, eval=TRUE-----------------------------------------
iris_df[iris_df[, "Sepal.Length"] > 7.5 & iris_df[, "Species"]=="virginica", ]

## ----plyr_filter_boolean, eval=TRUE--------------------------------------
filter(iris_df, Sepal.Length > 7.5 | Sepal.Length < 5.5, Species=="virginica")

## ----plyr_subset, eval=TRUE----------------------------------------------
slice(iris_df, 1:2)

## ----plyr_subset_base, eval=TRUE-----------------------------------------
iris_df[1:2,]

## ----plyr_sample_set2, eval=TRUE-----------------------------------------
df1 <- bind_cols(data_frame(ids1=paste0("g", 1:10)), as_data_frame(matrix(1:40, 10, 4, dimnames=list(1:10, paste0("CA", 1:4)))))
df1

## ----plyr_subset_names, eval=TRUE----------------------------------------
slice(df1, match(c("g10", "g4", "g4"), df1$ids1))

## ----plyr_subset_names_base, eval=TRUE-----------------------------------
df1_old <- as.data.frame(df1)
rownames(df1_old) <- df1_old[,1]
df1_old[c("g10", "g4", "g4"),]

## ----plyr_order1, eval=TRUE----------------------------------------------
arrange(iris_df, Species, Sepal.Length, Sepal.Width)

## ----plyr_order2, eval=TRUE----------------------------------------------
arrange(iris_df, desc(Species), Sepal.Length, Sepal.Width)

## ----plyr_order_base, eval=TRUE------------------------------------------
iris_df[order(iris_df$Species, iris_df$Sepal.Length, iris_df$Sepal.Width), ]
iris_df[order(iris_df$Species, decreasing=TRUE), ] 

## ----plyr_col_select1, eval=TRUE-----------------------------------------
select(iris_df, Species, Petal.Length, Sepal.Length)

## ----plyr_col_select2, eval=TRUE-----------------------------------------
select(iris_df, Sepal.Length : Petal.Width)

## ----plyr_col_drop, eval=TRUE--------------------------------------------
select(iris_df, -(Sepal.Length : Petal.Width))

## ----plyr_col_rename, eval=TRUE------------------------------------------
rename(iris_df, new_col_name = Species)

## ----baser_col_rename, eval=FALSE----------------------------------------
## colnames(iris_df)[colnames(iris_df)=="Species"] <- "new_col_names"

## ----plyr_unique, eval=TRUE----------------------------------------------
distinct(iris_df, Species, .keep_all=TRUE)

## ----baser_unique, eval=TRUE---------------------------------------------
iris_df[!duplicated(iris_df$Species),]

## ----plyr_mutate, eval=TRUE----------------------------------------------
mutate(iris_df, Ratio = Sepal.Length / Sepal.Width, Sum = Sepal.Length + Sepal.Width)

## ----plyr_transmute, eval=TRUE-------------------------------------------
transmute(iris_df, Ratio = Sepal.Length / Sepal.Width, Sum = Sepal.Length + Sepal.Width)

## ----plyr_bind_cols, eval=TRUE-------------------------------------------
bind_cols(iris_df, iris_df)

## ----plyr_summarize1, eval=TRUE------------------------------------------
summarize(iris_df, mean(Petal.Length))

## ----plyr_summarize2, eval=TRUE------------------------------------------
summarize_all(iris_df[,1:4], mean)

## ----plyr_summarize, eval=TRUE-------------------------------------------
summarize(group_by(iris_df, Species), mean(Petal.Length))

## ----plyr_summarize3, eval=TRUE------------------------------------------
summarize_all(group_by(iris_df, Species), mean) 

## ----plyr_join_sample, eval=TRUE-----------------------------------------
df1 <- bind_cols(data_frame(ids1=paste0("g", 1:10)), as_data_frame(matrix(1:40, 10, 4, dimnames=list(1:10, paste0("CA", 1:4)))))
df1
df2 <- bind_cols(data_frame(ids2=paste0("g", c(2,5,11,12))), as_data_frame(matrix(1:16, 4, 4, dimnames=list(1:4, paste0("CB", 1:4)))))
df2

## ----plyr_inner_join, eval=TRUE------------------------------------------
inner_join(df1, df2, by=c("ids1"="ids2"))

## ----plyr_left_join, eval=TRUE-------------------------------------------
left_join(df1, df2, by=c("ids1"="ids2"))

## ----plyr_right_join, eval=TRUE------------------------------------------
right_join(df1, df2, by=c("ids1"="ids2"))

## ----plyr_full_join, eval=TRUE-------------------------------------------
full_join(df1, df2, by=c("ids1"="ids2"))

## ----plyr_anti_join, eval=TRUE-------------------------------------------
anti_join(df1, df2, by=c("ids1"="ids2"))

## ----plyr_chaining1, eval=TRUE-------------------------------------------
iris_df %>% # Declare data frame to use 
    select(Sepal.Length:Species) %>% # Select columns
    filter(Species=="setosa") %>% # Filter rows by some value
    arrange(Sepal.Length) %>% # Sort by some column
    mutate(Subtract=Petal.Length - Petal.Width) # Calculate and append
    # write_tsv("iris.txt") # Export to file, omitted here to show result 

## ----plyr_chaining2, eval=TRUE-------------------------------------------
iris_df %>% # Declare data frame to use 
    group_by(Species) %>% # Group by species
    summarize(Mean_Sepal.Length=mean(Sepal.Length), 
              Max_Sepal.Length=max(Sepal.Length),
              Min_Sepal.Length=min(Sepal.Length),
              SD_Sepal.Length=sd(Sepal.Length),
              Total=n()) 

## ----plyr_chaining3, eval=TRUE-------------------------------------------
iris_df %>% 
    group_by(Species) %>% 
    summarize_all(mean) %>% 
    reshape2::melt(id.vars=c("Species"), variable.name = "Samples", value.name="Values") %>%
    ggplot(aes(Samples, Values, fill = Species)) + 
        geom_bar(position="dodge", stat="identity")

## ----load_sqlite, eval=TRUE----------------------------------------------
library(RSQLite)
mydb <- dbConnect(SQLite(), "test.db") # Creates database file test.db
mydf1 <- data.frame(ids=paste0("id", seq_along(iris[,1])), iris)
mydf2 <- mydf1[sample(seq_along(mydf1[,1]), 10),]
dbWriteTable(mydb, "mydf1", mydf1)
dbWriteTable(mydb, "mydf2", mydf2)

## ----list_tables, eval=TRUE----------------------------------------------
dbListTables(mydb)

## ----import_sqlite_tables, eval=TRUE-------------------------------------
dbGetQuery(mydb, 'SELECT * FROM mydf2')

## ----query_sqlite_tables, eval=TRUE--------------------------------------
dbGetQuery(mydb, 'SELECT * FROM mydf1 WHERE "Sepal.Length" < 4.6')

## ----join_sqlite_tables, eval=TRUE---------------------------------------
dbGetQuery(mydb, 'SELECT * FROM mydf1, mydf2 WHERE mydf1.ids = mydf2.ids')

## ----sample_data, eval=TRUE----------------------------------------------
set.seed(1410)
y <- matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))

## ----basic_scatter_plot, eval=TRUE---------------------------------------
plot(y[,1], y[,2]) 

## ----pairs_scatter_plot, eval=TRUE---------------------------------------
pairs(y) 

## ----labels_scatter_plot, eval=TRUE--------------------------------------
plot(y[,1], y[,2], pch=20, col="red", main="Symbols and Labels")
text(y[,1]+0.03, y[,2], rownames(y))

## ----row_scatter_plot, eval=TRUE-----------------------------------------
plot(y[,1], y[,2], type="n", main="Plot of Labels")
text(y[,1], y[,2], rownames(y)) 

## ----plot_usage, eval=FALSE----------------------------------------------
## grid(5, 5, lwd = 2)
## op <- par(mar=c(8,8,8,8), bg="lightblue")
## plot(y[,1], y[,2], type="p", col="red", cex.lab=1.2, cex.axis=1.2,
##      cex.main=1.2, cex.sub=1, lwd=4, pch=20, xlab="x label",
##      ylab="y label", main="My Main", sub="My Sub")
## par(op)

## ----plot_regression, eval=TRUE------------------------------------------
plot(y[,1], y[,2])
myline <- lm(y[,2]~y[,1]); abline(myline, lwd=2) 
summary(myline) 

## ----plot_regression_log, eval=TRUE--------------------------------------
plot(y[,1], y[,2], log="xy") 

## ----plot_regression_math, eval=TRUE-------------------------------------
plot(y[,1], y[,2]); text(y[1,1], y[1,2], expression(sum(frac(1,sqrt(x^2*pi)))), cex=1.3) 

## ----plot_line_single, eval=TRUE-----------------------------------------
plot(y[,1], type="l", lwd=2, col="blue") 

## ----plot_line_many, eval=TRUE-------------------------------------------
split.screen(c(1,1)) 
plot(y[,1], ylim=c(0,1), xlab="Measurement", ylab="Intensity", type="l", lwd=2, col=1)
for(i in 2:length(y[1,])) { 
	screen(1, new=FALSE)
	plot(y[,i], ylim=c(0,1), type="l", lwd=2, col=i, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n") 
}
close.screen(all=TRUE) 

## ----plot_bar_simple, eval=TRUE------------------------------------------
barplot(y[1:4,], ylim=c(0, max(y[1:4,])+0.3), beside=TRUE, legend=letters[1:4]) 
text(labels=round(as.vector(as.matrix(y[1:4,])),2), x=seq(1.5, 13, by=1) + sort(rep(c(0,1,2), 4)), y=as.vector(as.matrix(y[1:4,]))+0.04) 

## ----plot_bar_error, eval=TRUE-------------------------------------------
bar <- barplot(m <- rowMeans(y) * 10, ylim=c(0, 10))
stdev <- sd(t(y))
arrows(bar, m, bar, m + stdev, length=0.15, angle = 90)

## ----plot_hist, eval=TRUE------------------------------------------------
hist(y, freq=TRUE, breaks=10)

## ----plot_dens, eval=TRUE------------------------------------------------
plot(density(y), col="red")

## ----plot_pie, eval=TRUE-------------------------------------------------
pie(y[,1], col=rainbow(length(y[,1]), start=0.1, end=0.8), clockwise=TRUE)
legend("topright", legend=row.names(y), cex=1.3, bty="n", pch=15, pt.cex=1.8, 
col=rainbow(length(y[,1]), start=0.1, end=0.8), ncol=1) 

## ----color_palette, eval=TRUE--------------------------------------------
palette()
palette(rainbow(5, start=0.1, end=0.2))
palette()
palette("default")

## ----color_grey, eval=TRUE-----------------------------------------------
gray(seq(0.1, 1, by= 0.2))

## ----color_gradient, eval=TRUE-------------------------------------------
library(gplots)
colorpanel(5, "darkblue", "yellow", "white")

## ----save_graphics, eval=FALSE-------------------------------------------
## pdf("test.pdf")
## plot(1:10, 1:10)
## dev.off()

## ----save_graphics_svg, eval=FALSE---------------------------------------
## library("RSvgDevice")
## devSVG("test.svg")
## plot(1:10, 1:10)
## dev.off()

## ----import_data1, eval=FALSE--------------------------------------------
## my_mw <- read.delim(file="MolecularWeight_tair7.xls", header=T, sep="\t")
## my_mw[1:2,]

## ----import_data2, eval=FALSE--------------------------------------------
## my_target <- read.delim(file="TargetP_analysis_tair7.xls", header=T, sep="\t")
## my_target[1:2,]

## ----import_data1b, eval=TRUE--------------------------------------------
my_mw <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/MolecularWeight_tair7.xls", header=T, sep="\t") 
my_mw[1:2,]

## ----import_data2b, eval=TRUE--------------------------------------------
my_target <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/TargetP_analysis_tair7.xls", header=T, sep="\t") 
my_target[1:2,]

## ----col_names_uni, eval=TRUE--------------------------------------------
colnames(my_target)[1] <- "ID"
colnames(my_mw)[1] <- "ID" 

## ----merge_tables, eval=TRUE---------------------------------------------
my_mw_target <- merge(my_mw, my_target, by.x="ID", by.y="ID", all.x=T)

## ----merge_tables_shorten, eval=TRUE-------------------------------------
my_mw_target2a <- merge(my_mw, my_target[1:40,], by.x="ID", by.y="ID", all.x=T)  # To remove non-matching rows, use the argument setting 'all=F'.
my_mw_target2 <- na.omit(my_mw_target2a) # Removes rows containing "NAs" (non-matching rows).

## ----homework3D_solutions, eval=FALSE, echo=FALSE------------------------
## my_mw_target_tmp <- merge(my_mw, my_target[1:40,], by.x="ID", by.y="ID", all=FALSE)
## all(my_mw_target2 == my_mw_target_tmp)
## my_mw_target2a <- as.matrix(my_mw_target2a)
## my_mw_target2a[is.na(my_mw_target2a)] <- 0
## my_mw_target2a <- as.data.frame(my_mw_target2a)

## ----filter_tables1, eval=TRUE-------------------------------------------
query <- my_mw_target[my_mw_target[, 2] > 100000 & my_mw_target[, 4] == "C", ] 
query[1:4, ]
dim(query)

## ----string_sub, eval=TRUE-----------------------------------------------
my_mw_target3 <- data.frame(loci=gsub("\\..*", "", as.character(my_mw_target[,1]), perl = TRUE), my_mw_target)
my_mw_target3[1:3,1:8]

## ----calcul_1, eval=TRUE-------------------------------------------------
mycounts <- table(my_mw_target3[,1])[my_mw_target3[,1]]
my_mw_target4 <- cbind(my_mw_target3, Freq=mycounts[as.character(my_mw_target3[,1])]) 

## ----calcul_2, eval=TRUE-------------------------------------------------
data.frame(my_mw_target4, avg_AA_WT=(my_mw_target4[,3] / my_mw_target4[,4]))[1:2,5:11] 

## ----calcul_3, eval=TRUE-------------------------------------------------
mymean <- apply(my_mw_target4[,6:9], 1, mean)
mystdev <- apply(my_mw_target4[,6:9], 1, sd, na.rm=TRUE)
data.frame(my_mw_target4, mean=mymean, stdev=mystdev)[1:2,5:12] 

## ----plot_example, eval=TRUE---------------------------------------------
plot(my_mw_target4[1:500,3:4], col="red")

## ----export_example, eval=TRUE-------------------------------------------
write.table(my_mw_target4, file="my_file.xls", quote=F, sep="\t", col.names = NA) 

## ----source_example, eval=FALSE------------------------------------------
## source("exerciseRbasics.R")

## ----install_rmarkdown, eval=FALSE---------------------------------------
## install.packages("rmarkdown")

## ----render_rmarkdown, eval=FALSE, message=FALSE-------------------------
## rmarkdown::render("sample.Rmd", clean=TRUE, output_format="html_document")

## ----render_commandline, eval=FALSE, message=FALSE, engine="sh"----------
## $ echo "rmarkdown::render('sample.Rmd', clean=TRUE)" | R --slave
## $ Rscript -e "rmarkdown::render('sample.Rmd', clean=TRUE)"

## ----render_makefile, eval=FALSE, message=FALSE, engine="sh"-------------
## $ make -B

## ----kable---------------------------------------------------------------
library(knitr)
kable(iris[1:12,])

## ----some_jitter_plot, eval=TRUE-----------------------------------------
library(ggplot2)
dsmall <- diamonds[sample(nrow(diamonds), 1000), ]
ggplot(dsmall, aes(color, price/carat)) + geom_jitter(alpha = I(1 / 2), aes(color=color))

## ----some_custom_inserted_plot, eval=TRUE, warning=FALSE, message=FALSE----
png("myplot.png")
ggplot(dsmall, aes(color, price/carat)) + geom_jitter(alpha = I(1 / 2), aes(color=color))
dev.off()

## ----rmarkdown_symbolic_link, eval=FALSE---------------------------------
## cd ~/.html
## ln -s ~/bigdata/today/rmarkdown/sample.html sample.html

## ----fluidpage, eval=FALSE-----------------------------------------------
## ui <- fluidPage()

## ----server, eval=FALSE--------------------------------------------------
## server <- function(input, output) {}

## ----shinyapp, eval=FALSE------------------------------------------------
## shinyApp(ui = ui, server = server)

## ----runshinyapp1, eval=FALSE--------------------------------------------
## library(shiny)
## runApp("myappdir") # To show code in app, add argument: display.mode="showcase"

## ----deployshinyapp1, eval=FALSE-----------------------------------------
## setwd("myappdir")
## library(rsconnect)
## deployApp()

## ----embedshiny, eval=FALSE----------------------------------------------
## <iframe src="https://tgirke.shinyapps.io/diamonds/" style="border: none; width: 880px; height: 900px"></iframe>

## ----learnshiny, eval=FALSE----------------------------------------------
## mydir <- system.file("examples", package="shiny")
## dir.create('my_shiny_test_dir')
## file.copy(mydir, "my_shiny_test_dir", recursive=TRUE)
## setwd("my_shiny_test_dir/examples")
## runApp("01_hello") # Runs first example app in directory
## dir() # Lists available Shiny examples (directories).

## ----sessionInfo---------------------------------------------------------
sessionInfo()

