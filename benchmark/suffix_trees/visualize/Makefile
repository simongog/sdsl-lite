include ../../../Make.helper

CONFIG_FILES=../test_case.config

all: suffix_trees.pdf

suffix_trees.pdf: suffix_trees.tex
	@echo "Use pdflatex to generate suffix_trees.pdf"
	@pdflatex suffix_trees.tex >> LaTeX.log 2>&1

suffix_trees.tex: ../results/all.txt ../../basic_functions.R suffix_trees.R $(CONFIG_FILES)
	@echo "Use R to generate tex files."
	@R --vanilla < suffix_trees.R > R.log 2>&1

graph_data: ../results/all.txt ../../basic_functions.R graph_data.R $(CONFIG_FILES)
	@echo "Use R to generate raw graph data."
	@R --vanilla < graph_data.R > R.log 2>&1

clean:
	rm -f suffix_trees.pdf suffix_trees.aux suffix_trees.log fig-space.tex fig-time.tex R.log LaTeX.log
