include ../../Make.helper
CONFIG_FILES=index-filter.config ../index.config ../test_case.config ../sample.config

all: locate.pdf

locate.pdf: locate.tex 
	@echo "Use pdflatex to generate locate.pdf"
	@pdflatex locate.tex >> LaTeX.Log 2>&1

locate.tex: ../../basic_functions.R locate.R $(CONFIG_FILES) ../results/all.txt
	@echo "Use R to generate locate.tex"
	@R --vanilla < locate.R > R.log 2>&1

clean:
	rm -f locate.pdf locate.aux locate.tex \
	   locate.log R.log LaTeX.log
