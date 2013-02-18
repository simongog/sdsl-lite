TABLES = tbl-count-NOOPT.tex tbl-count-NOSSE.tex tbl-count-SSE.tex \
		   tbl-count-HP.tex 


all: count.pdf

count.pdf: $(TABLES)	
	pdflatex count.tex

tbl-count-%: ../results/count.txt  basic_functions.R count.R count.tex
	R --vanilla < count.R > R.log 2>&1

clean:
	rm -f $(TABLES) count.pdf count.aux count.log R.log
