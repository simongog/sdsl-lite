include ../../Make.helper
COMPILE_IDS:=$(call config_ids,../compile_options.config)

TABLES = $(foreach COMPILE_ID,$(COMPILE_IDS),tbl-count-$(COMPILE_ID).tex)

all: count.pdf

count.pdf: $(TABLES)	
	@echo "Use pdflatex to generate count.pdf"
	@pdflatex count.tex >> LaTeX.Log 2>&1

tbl-count-%: ../results/all.txt ../index.config ../../basic_functions.R count.R count.tex
		

count.tex:
	R --vanilla < count.R > R.log 2>&1

clean:
	rm -f $(TABLES) count.pdf count.aux count.tex \
	   count.log R.log LaTeX.log
