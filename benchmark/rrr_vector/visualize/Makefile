include ../../Make.helper

CONFIG_FILES=../index.config ../test_case.config

all: rrr.pdf

rrr.pdf: rrr.tex
	@echo "Use pdflatex to generate count.pdf"
	@pdflatex rrr.tex >> LaTeX.Log 2>&1

count.tex: ../results/all.txt ../../basic_functions.R rrr.R $(CONFIG_FILES)
	@echo "Use R to generate count.tex"
	@R --vanilla < rrr.R > R.log 2>&1

clean:
	rm -f rrr.pdf rrr.aux rrr.tex tbl-rrr* \
	      rrr.log R.log LaTeX.log
