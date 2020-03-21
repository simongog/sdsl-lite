include ../../../Make.helper

CONFIG_FILES= ../test_case.config

all: gm.pdf

gm.pdf: gm.tex
	@echo "Use pdflatex to generate gm.pdf"
	@pdflatex gm.tex >> LaTeX.Log 2>&1

gm.tex: ../results/all.txt ../../basic_functions.R gm.R $(CONFIG_FILES)
	@echo "Use R to generate gm.tex"
	@R --vanilla < gm.R > R.log 2>&1

clean:
	rm -f gm.pdf gm.aux gm.tex fig-gm* \
	      gm.log R.log LaTeX.log
