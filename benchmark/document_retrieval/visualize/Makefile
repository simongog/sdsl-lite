include ../../../Make.helper
CONFIG_FILES=../index.config ../test_case.config ../pattern_length.config

all: doc_re_time.pdf

doc_re_time.pdf: doc_re_time.tex fig-runtime.tex tbl-indexes.tex
	@echo "Use pdflatex to generate doc_re_time.pdf"
#	@pdflatex doc_re_time.tex >> LaTeX.Log 2>&1
	@pdflatex doc_re_time.tex 

#tbl-locate.tex fig-locate.tex: ../../basic_functions.R locate.R $(CONFIG_FILES) ../results/all.txt
#	@echo "Use R to generate tbl-locate.tex and fig-locate.tex"
fig-runtime.tex tbl-indexes.tex: ../../basic_functions.R runtime.R $(CONFIG_FILES) ../results/all.txt
	@R --vanilla < runtime.R > R.log 2>&1

clean:
	rm -f doc_re_time.pdf doc_re_time.aux fig-runtime.tex tbl-indexes.tex \
	   doc_re_time.log R.log LaTeX.log
