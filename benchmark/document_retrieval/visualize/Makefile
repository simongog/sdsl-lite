include ../../../Make.helper
CONFIG_FILES=../index.config ../test_case.config ../pattern_length.config \
			 ../index_int.config ../test_case_int.config ../pattern_length_int.config

all: doc_re_time.pdf

doc_re_time.pdf: doc_re_time.tex \
                 fig-runtime.tex fig-runtime_int.tex \
				 tbl-sizes.tex tbl-sizes_int.tex \
		         tbl-indexes.tex tbl-indexes_int.tex \
				 tbl-collections.tex tbl-collections_int.tex 
	@echo "Use pdflatex to generate doc_re_time.pdf"
	@pdflatex doc_re_time.tex >> LaTeX.log 2>&1

#tbl-locate.tex fig-locate.tex: ../../basic_functions.R locate.R $(CONFIG_FILES) ../results/all.txt
#	@echo "Use R to generate tbl-locate.tex and fig-locate.tex"
fig-runtime.tex fig-runtime_int.tex tbl-indexes.tex tbl-indexes_int.tex tbl-sizes.tex tbl-sizes_int.tex tbl-collections.tex tbl-collections_int.tex: ../../basic_functions.R runtime.R \
		                         $(CONFIG_FILES) ../results/all.txt \
								 ../info/sizes.txt ../info/sizes_int.txt
	@R --vanilla < runtime.R > R.log 2>&1

clean:
	rm -f doc_re_time.pdf doc_re_time.aux \
		  fig-runtime.tex tbl-indexes.tex \
		  fig-runtime_int.tex tbl-indexes_int.tex \
		  tbl-sizes.tex tbl-sizes_int.tex \
		  tbl-collections.tex tbl-collections_int.tex \
	      doc_re_time.log R.log LaTeX.log
