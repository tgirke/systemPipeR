#############################################################
## R code for creating dynamic programming matrix for HW4B ##
#############################################################

## Define input sequences
x <- "FIPFSAGPRNCIGQK"
y <- "PFGFGKRSCMGRRLA"

## Vectorize input sequences
x <- substring(x, 1:nchar(x), 1:nchar(x))
y <- substring(y, 1:nchar(y), 1:nchar(y))

## Define gap penality
gp <- 8 # Gap penalty

## Create dynamic programming matrix based on x, y and gp
ma <- matrix(NA, length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))
ma[1,] <- seq(0, -(length(ma[1,])-1) * gp, -gp)
ma[,1] <- seq(0, -(length(ma[,1])-1) * gp, -gp)
ma

## If desired, write ma to tabular file and upload it to Google Sheets 
write.table(ma, file="ma.xls", quote=FALSE, na = "", col.names = NA, sep="\t")

