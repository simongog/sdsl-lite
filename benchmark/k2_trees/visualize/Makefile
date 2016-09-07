include ../../../Make.helper

CONFIG_FILES= ../test_case.config 

all: wt.pdf

wt.pdf: wt.tex
	@echo "Use pdflatex to generate wt.pdf"
	@pdflatex wt.tex >> LaTeX.Log 2>&1

wt.tex: ../results/all.txt ../../basic_functions.R wt.R $(CONFIG_FILES)
	@echo "Use R to generate wt.tex"
	@R --vanilla < wt.R > R.log 2>&1

clean:
	rm -f wt.pdf wt.aux wt.tex fig* \
	wt.log R.log LaTeX.log
