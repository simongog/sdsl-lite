include ../../../Make.helper

CONFIG_FILES= ../test_case.config 

all: lcp.pdf

lcp.pdf: lcp.tex
	@echo "Use pdflatex to generate lcp.pdf"
	@pdflatex lcp.tex >> LaTeX.Log 2>&1

lcp.tex: ../results/all.txt ../../basic_functions.R lcp.R $(CONFIG_FILES)
	@echo "Use R to generate lcp.tex"
	@R --vanilla < lcp.R > R.log 2>&1

clean:
	rm -f lcp.pdf lcp.aux lcp.tex fig* \
	lcp.log R.log LaTeX.log
