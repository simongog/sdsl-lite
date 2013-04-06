include ../../../Make.helper
CONFIG_FILES=index-filter.config ../index.config ../test_case.config ../sample.config

all: extract.pdf

extract.pdf: extract.tex tbl-extract.tex fig-extract.tex
	@echo "Use pdflatex to generate extract.pdf"
	@pdflatex extract.tex >> LaTeX.Log 2>&1

tbl-extract.tex fig-extract.tex: ../../basic_functions.R extract.R $(CONFIG_FILES) ../results/all.txt
	@echo "Use R to generate fig-extract.tex and tbl-extract.tex"
	@R --vanilla < extract.R > R.log 2>&1

clean:
	rm -f extract.pdf extract.aux extract.log R.log LaTeX.log \
          tbl-extract.tex fig-extract.tex
