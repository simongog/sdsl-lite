include ../../../Make.helper

CONFIG_FILES= ../test_case.config

all: rrr.pdf

rrr.pdf: rrr.tex
	@echo "Use pdflatex to generate rrr.pdf"
	@pdflatex rrr.tex >> LaTeX.Log 2>&1

rrr.tex: ../results/all.txt ../../basic_functions.R rrr.R $(CONFIG_FILES)
	@echo "Use R to generate rrr.tex"
	@R --vanilla < rrr.R > R.log 2>&1

clean:
	rm -f rrr.pdf rrr.aux rrr.tex fig-rrr* \
	      rrr.log R.log LaTeX.log
