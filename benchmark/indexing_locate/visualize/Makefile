include ../../../Make.helper
CONFIG_FILES=index-filter.config ../index.config ../test_case.config ../sample.config

all: locate.pdf

locate.pdf: locate.tex tbl-locate.tex fig-locate.tex
	@echo "Use pdflatex to generate locate.pdf"
	@pdflatex locate.tex >> LaTeX.Log 2>&1

tbl-locate.tex fig-locate.tex: ../../basic_functions.R locate.R $(CONFIG_FILES) ../results/all.txt
	@echo "Use R to generate tbl-locate.tex and fig-locate.tex"
	@R --vanilla < locate.R > R.log 2>&1

clean:
	rm -f locate.pdf locate.aux fig-locate.tex tbl-locate.tex \
	   locate.log R.log LaTeX.log
