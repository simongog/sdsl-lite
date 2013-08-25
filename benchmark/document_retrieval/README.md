# Benchmarking operation top-k search on simple document search implementations

## Methodology

The benchmark is close to the one used in the ESA 2010 article
of Culpepper, Navarro, Puglisi and Turpin.

Explored dimensions:

  * text type
  * instance size (just adjust the `test_case.config` file for this)
  * index implementations
    - [Sadakane's method](../../tutorial/document_listing/doc_list_index_sada.hpp)
    - [Wavelet tree greedy traversal](../../tutorial/document_listing/doc_list_index_greedy.hpp)
    - [Wavelet tree quantile probing](../../tutorial/document_listing/doc_list_index_qprobing.hpp)
 
## Directory structure
  * [bin](./bin): Contains the executables of the project.
    - `doc_list_index_*` generates and queries and indexes
  * [indexes](./indexes): Contains the 
  * [results](./results): Contains the results of the experiments.
  * [src](./src): Contains the source code of the benchmark.
  * [visualize](./visualize): Contains a `R`-script which generates
                              a report.

## Prerequisites
  * For the visualization you need the following software:
    - [R][RPJ] with package `xtable`. You can install the
      package by calling `install.packages("xtable")` in R.
    - [pdflatex][LT] to generate the pdf reports.

## Usage
    TODO
	

[RPJ]: http://www.r-project.org/ "R"
[LT]: http://www.tug.org/applications/pdftex/ "pdflatex"
