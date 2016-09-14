include ../../../Make.helper

CONFIG_FILES= ../test_case.config

all: k2.pdf

k2.pdf: k2.tex
	@echo "Use pdflatex to generate k2.pdf"
	@pdflatex k2.tex >> LaTeX.Log 2>&1

k2.tex: ../results/all.txt ../../basic_functions.R k2.R $(CONFIG_FILES)
	@echo "Use R to generate k2.tex"
	@R --vanilla < k2.R > R.log 2>&1

clean:
	rm -f k2.pdf k2.aux k2.tex fig* \
	k2.log R.log LaTeX.log
