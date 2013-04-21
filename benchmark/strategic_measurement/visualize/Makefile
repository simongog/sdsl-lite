count.pdf: count.tex
	pdflatex count.tex

count.tex: ../results/count.txt count.R
	R --vanilla < count.R

locate.pdf: locate.tex
	pdflatex locate.tex

locate.tex: ../results/locate.txt locate.R
	R --vanilla < locate.R

cand.pdf: cand.tex
	pdflatex cand.tex

cand.tex: ../results/cand.txt cand.R
	R --vanilla < cand.R

occ.pdf: occ.tex tbl-occ.tex
	@pdflatex occ.tex

tbl-occ.tex: occ.R ../results/occ.txt
	@R --vanilla < occ.R 
