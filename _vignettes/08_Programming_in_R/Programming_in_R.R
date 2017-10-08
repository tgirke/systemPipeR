## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
suppressPackageStartupMessages({
    library(limma) 
    library(ggplot2) }) 

## ----if_statement, eval=FALSE--------------------------------------------
## if(TRUE) {
## 	statements_1
## } else {
## 	statements_2
## }

## ----if_statement_example, eval=TRUE-------------------------------------
if(1==0) { 
	print(1) 
} else { 
	print(2) 
}

## ----ifelse_statement, eval=FALSE----------------------------------------
## ifelse(test, true_value, false_value)

## ----ifelse_statement_example, eval=TRUE---------------------------------
x <- 1:10 
ifelse(x<5, x, 0)

## ----for_loop, eval=FALSE------------------------------------------------
## for(variable in sequence) {
## 	statements
## }

## ----for_loop_example, eval=TRUE-----------------------------------------
mydf <- iris
myve <- NULL
for(i in seq(along=mydf[,1])) {
	myve <- c(myve, mean(as.numeric(mydf[i,1:3])))
}
myve[1:8]

## ----for_loop_inject_example, eval=TRUE----------------------------------
myve <- numeric(length(mydf[,1]))
for(i in seq(along=myve)) {
	myve[i] <- mean(as.numeric(mydf[i,1:3]))
}
myve[1:8]

## ----for_loop_stop_example, eval=FALSE-----------------------------------
## x <- 1:10
## z <- NULL
## for(i in seq(along=x)) {
## 	if(x[i] < 5) {
## 		z <- c(z, x[i]-1)
## 	} else {
## 		stop("values need to be < 5")
## 	}
## }

## ----while_loop, eval=FALSE----------------------------------------------
## while(condition) {
## 	statements
## }

## ----while_loop_example, eval=TRUE---------------------------------------
z <- 0
while(z<5) { 
	z <- z + 2
	print(z)  
}

## ----apply_loop, eval=FALSE----------------------------------------------
## apply(X, MARGIN, FUN, ARGs)

## ----apply_loop_example, eval=TRUE---------------------------------------
apply(iris[1:8,1:3], 1, mean)

## ----tapply_loop, eval=FALSE---------------------------------------------
## tapply(vector, factor, FUN)

## ----tapply_loop_example, eval=TRUE--------------------------------------
iris[1:2,]
tapply(iris$Sepal.Length, iris$Species, mean)

## ----lapply_loop_example, eval=TRUE--------------------------------------
x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
lapply(x, mean)
sapply(x, mean)

## ----lapply_loop_fct_example, eval=FALSE---------------------------------
## lapply(names(x), function(x) mean(x))
## sapply(names(x), function(x) mean(x))

## ----function_def_syntax, eval=FALSE-------------------------------------
## myfct <- function(arg1, arg2, ...) {
## 	function_body
## }

## ----function_call_syntax, eval=FALSE------------------------------------
## myfct(arg1=..., arg2=...)

## ----define_function_example, eval=TRUE----------------------------------
myfct <- function(x1, x2=5) { 
	z1 <- x1 / x1
	z2 <- x2 * x2
        myvec <- c(z1, z2) 
        return(myvec)
} 

## ----usage_function_example1, eval=TRUE----------------------------------
myfct(x1=2, x2=5) 

## ----usage_function_example2, eval=TRUE----------------------------------
myfct(2, 5) 

## ----usage_function_example3, eval=TRUE----------------------------------
myfct(x1=2) 

## ----usage_function_example4, eval=TRUE----------------------------------
myfct 

## ----grep_fct, eval=TRUE-------------------------------------------------
month.name[grep("A", month.name)] 

## ----gsub_fct, eval=TRUE-------------------------------------------------
gsub('(i.*a)', 'xxx_\\1', "virginica", perl = TRUE) 

## ----ls_fct, eval=TRUE---------------------------------------------------
mylist <- ls()
mylist[1] 

## ----eval_expr, eval=TRUE------------------------------------------------
get(mylist[1])

## ----eval_expr2, eval=TRUE-----------------------------------------------
eval(parse(text=mylist[1])) 

## ----back_ref, eval=TRUE-------------------------------------------------
x <- gsub("(a)","\\1_", month.name[1], perl=T) 
x

## ----split_string, eval=TRUE---------------------------------------------
strsplit(x,"_")

## ----reverse_string, eval=TRUE-------------------------------------------
paste(rev(unlist(strsplit(x, NULL))), collapse="") 

## ----sys_time, eval=TRUE-------------------------------------------------
system.time(ls()) 

## ----sys_date, eval=TRUE-------------------------------------------------
date() 

## ----sys_sleep, eval=TRUE------------------------------------------------
Sys.sleep(1) 

## ----read_lines, eval=TRUE-----------------------------------------------
cat(month.name, file="zzz.txt", sep="\n")
x <- readLines("zzz.txt")
x[1:6] 
x <- x[c(grep("^J", as.character(x), perl = TRUE))]
t(as.data.frame(strsplit(x, "u")))

## ----system_blast, eval=TRUE---------------------------------------------
system("blastall -p blastp -i seq.fasta -d uniprot -o seq.blastp")

## ----r_script1, eval=FALSE-----------------------------------------------
## source("my_script.R")

## ----r_cmd_script1, eval=FALSE, engine="sh"------------------------------
## Rscript my_script.R # or just ./myscript.R after making it executable
## R CMD BATCH my_script.R # Alternative way 1
## R --slave < my_script.R # Alternative way 2

## ----r_cmd_script2, eval=FALSE, engine="sh"------------------------------
## myarg <- commandArgs()
## print(iris[1:myarg[6], ])

## ----r_cmd_script3, eval=FALSE, engine="sh"------------------------------
## Rscript test.R 10

## ----package_skeleton1, eval=FALSE---------------------------------------
## package.skeleton(name="mypackage", code_files=c("script1.R", "script2.R"))

## ----r_build_package, eval=FALSE-----------------------------------------
## R CMD build mypackage
## R CMD check mypackage_1.0.tar.gz

## ----install_package, eval=FALSE-----------------------------------------
## install.packages("mypackage_1.0.tar.gz", repos=NULL)

## ----exercise1_for, eval=TRUE--------------------------------------------
myMA <- matrix(rnorm(500), 100, 5, dimnames=list(1:100, paste("C", 1:5, sep="")))
myve_for <- NULL
for(i in seq(along=myMA[,1])) {
	myve_for <- c(myve_for, mean(as.numeric(myMA[i, ])))
}
myResult <- cbind(myMA, mean_for=myve_for)
myResult[1:4, ]

## ----exercise1_while, eval=TRUE------------------------------------------
z <- 1
myve_while <- NULL
while(z <= length(myMA[,1])) {
	myve_while <- c(myve_while, mean(as.numeric(myMA[z, ])))
	z <- z + 1
}
myResult <- cbind(myMA, mean_for=myve_for, mean_while=myve_while)
myResult[1:4, -c(1,2)]

## ----exercise1_confirm, eval=TRUE----------------------------------------
all(myResult[,6] == myResult[,7])

## ----exercise1_apply, eval=TRUE------------------------------------------
myve_apply <- apply(myMA, 1, mean)
myResult <- cbind(myMA, mean_for=myve_for, mean_while=myve_while, mean_apply=myve_apply)
myResult[1:4, -c(1,2)]

## ----exercise1_noloops, eval=TRUE----------------------------------------
mymean <- rowMeans(myMA)
myResult <- cbind(myMA, mean_for=myve_for, mean_while=myve_while, mean_apply=myve_apply, mean_int=mymean)
myResult[1:4, -c(1,2,3)]

## ----exercise2_fct, eval=TRUE--------------------------------------------
myMA <- matrix(rnorm(100000), 10000, 10, dimnames=list(1:10000, paste("C", 1:10, sep="")))
myMA[1:2,]
myList <- tapply(colnames(myMA), c(1,1,1,2,2,2,3,3,4,4), list) 
names(myList) <- sapply(myList, paste, collapse="_")
myMAmean <- sapply(myList, function(x) apply(myMA[,x], 1, mean))
myMAmean[1:4,] 

## ----exercise2_fct_solution, eval=FALSE, echo=FALSE, keep.source=TRUE----
## myMAcomp <- function(myMA=myMA, group=c(1,1,1,2,2,2,3,3,4,4), myfct=mean) {
## 	myList <- tapply(colnames(myMA), group, list)
## 	names(myList) <- sapply(myList, paste, collapse="_")
## 	myMAmean <- sapply(myList, function(x) apply(myMA[, x, drop=FALSE], 1, myfct))
## 	return(myMAmean)
## }
## myMAcomp(myMA=myMA, group=c(1,1,1,2,2,2,3,3,4,4), myfct=mean)[1:2,]

## ----nested_loops1, eval=TRUE--------------------------------------------
setlist <- lapply(11:30, function(x) sample(letters, x, replace=TRUE))
names(setlist) <- paste("S", seq(along=setlist), sep="") 
setlist[1:6]

## ----nested_loops2, eval=TRUE--------------------------------------------
setlist <- sapply(setlist, unique)
olMA <- sapply(names(setlist), function(x) sapply(names(setlist), 
               function(y) sum(setlist[[x]] %in% setlist[[y]])))
olMA[1:12,] 

## ----nested_loops3, eval=TRUE--------------------------------------------
image(olMA)

## ----package_skeleton2, eval=FALSE---------------------------------------
## package.skeleton(name="mypackage", code_files=c("script1.R"))

## ----build_package_tar, eval=FALSE---------------------------------------
## system("R CMD build mypackage")

## ----install_package_tar, eval=FALSE-------------------------------------
## install.packages("mypackage_1.0.tar.gz", repos=NULL, type="source")
## library(mypackage)
## ?myMAcomp # Opens help for function defined by mypackage

## ----hw_revcomp, eval=TRUE-----------------------------------------------
x <- c("ATGCATTGGACGTTAG")  
x
x <- substring(x, 1:nchar(x), 1:nchar(x)) 
x
x <- rev(x) 
x
x <- paste(x, collapse="")
x
chartr("ATGC", "TACG", x) 

## ----hw_translate, eval=TRUE---------------------------------------------
AAdf <- read.table(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/AA.txt", header=TRUE, sep="\t") 
AAdf[1:4,]
AAv <- as.character(AAdf[,2]) 
names(AAv) <- AAdf[,1] 
AAv
y <- gsub("(...)", "\\1_", x) 
y <- unlist(strsplit(y, "_")) 
y <- y[grep("^...$", y)] 
AAv[y] 

## ----sessionInfo---------------------------------------------------------
sessionInfo()

