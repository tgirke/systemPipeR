#########################################################
## Makefile to execute systemPipeR and generate Report ##
#########################################################
## Usage:
## (1) Change value assigned to 'MAIN" to base name of your *.Rmd (or *.Rnw) source file 
## (2) To run workflow, execute from command-line 
##		$ make # or make all
## (3) To clean up directory, execute from command-line 
## 		$ make clean

## Define suffixes used for file interconversions
# .SUFFIXES: .tex .pdf .Rnw .R
.SUFFIXES: .Rmd .html .md .R

MAIN = sample

#######################################
## Build PDF report with knitr/Latex ##
#######################################
# all: $(MAIN).pdf 
all: $(MAIN).html 

# .Rnw.pdf:
# 	R CMD Sweave --engine=knitr::knitr --pdf $<
.Rmd.html:
	Rscript -e "rmarkdown::render('$<')"

########################################
## Build PDF report with Sweave/Latex ##
########################################
# all: $(MAIN).tex $(MAIN).R $(MAIN).pdf 
# 
# .Rnw.R:
# 	R CMD Stangle $<
# 
# .Rnw.tex:
# 	R CMD Sweave $<
# 
# .tex.pdf:
# 	pdflatex $<
# 	bibtex $*
# 	pdflatex $<
# 	pdflatex $<
 
########################
## Clean-up directory ##
########################
clean:
	rm -fv $(MAIN).aux $(MAIN).log $(MAIN).bbl $(MAIN).blg $(MAIN).out $(MAIN).toc $(MAIN).tex
	rm -fv $(MAIN)-*.pdf $(MAIN)-*.png $(MAIN)-*.jpg
	rm -fv *~ 
