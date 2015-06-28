#############################################################
## Makefile to execute systemPipeR and generate PDF Report ##
#############################################################
## Usage:
## (1) Change string assigned to 'MAIN" to base name of your *.Rnw source file 
## (2) To run workflow, execute from command-line 
##		$ make # or make all
## (3) To clean up directory, execute from command-line 
## 		$ make clean

## Define suffixes used for file interconversions
.SUFFIXES: .tex .pdf .Rnw .R

MAIN = systemPipeRNAseq # change assiged string to base name of your *.Rnw source file

all: $(MAIN).tex $(MAIN).R $(MAIN).pdf 

.Rnw.R:
	R CMD Stangle $<

.Rnw.tex:
	R CMD Sweave $<

.tex.pdf:
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

clean:
	rm -fv $(MAIN).aux $(MAIN).log $(MAIN).bbl $(MAIN).blg $(MAIN).toc $(MAIN).out
	rm -fv $(MAIN)-*.pdf $(MAIN)-*.png $(MAIN)-*.jpg
	rm -fv *~ 
