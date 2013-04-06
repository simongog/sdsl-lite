include ../../../Make.helper
COMPILE_IDS:=$(call config_ids,../compile_options.config)
TABLES = $(foreach COMPILE_ID,$(COMPILE_IDS),tbl-count-$(COMPILE_ID).tex)

CONFIG_FILES=index-filter.config ../index.config ../test_case.config

all: count.pdf

count.pdf: count.tex
	@echo "Use pdflatex to generate count.pdf"
	@pdflatex count.tex >> LaTeX.Log 2>&1

count.tex: ../results/all.txt ../../basic_functions.R count.R $(CONFIG_FILES)
	@echo "Use R to generate count.tex"
	@R --vanilla < count.R > R.log 2>&1

clean:
	rm -f $(TABLES) count.pdf count.aux count.tex \
	   count.log R.log LaTeX.log
