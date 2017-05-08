# Benchmarking k2 trees

## Methodology

Explored dimensions:

  * k2 tree implementations
  * test cases
  * methods (`adj`, `neighbors`, `reverse_neighbors`)

## Data

  * The data input for the benchmarking are arc files. An arc file is a text
      file where each line represents a directed edge between two nodes, the
      first column is the origin node and the second the target node.

### Test cases

  * EXAMPLE test case uses an small file with a little more than 100 nodes and
      roughly the same number of edges from web data commons.
  * HOSTGRAPH is a test case where the data comes from the Web Corpus released
      by the Common Crawl Foundation in April 2014. The file aggregates the
      page graph by subdomain/host. It has 123.660.351 edges.

## Directory structure

  * [bin](./bin): Contains the executables of the project.
    * `build_*` generates the binary file with the graph from the arc files.
    * `gen_*` executes the experiments.
  * [results](./results): Contains the results of the experiments.
  * [src](./src):  Contains the source code of the benchmark.
  * [visualize](./visualize): Contains a `R`-script which generates
    a report in LaTeX format.

## Prerequisites

  * For the visualization you need the following software:
    - [R][RPJ] with package `tikzDevice`. You can install the
      package by calling
      `install.packages("filehash", repos="http://cran.r-project.org")`
      and
      `install.packages("tikzDevice", repos="http://R-Forge.R-project.org")`
      in `R`.
    - [pdflatex][LT] to generate the pdf reports.

## Usage

 * `make timing` compiles the programs, downloads or generates
    the test instances, builds the k2 trees,
    runs the performance tests and generated a report located at
   `visualize/k2.pdf`. The raw numbers of the timings
   can be found in the `results/all.txt`. The default benchmark
   took 75 minutes on my machine (MacBookPro Retina 2.6Ghz Intel
   Core i5 16GB 1600 Mhz DDR3, SSD). Have a look at the
   [complete report][RES].
 * All created binaries and test results can be deleted
   by calling `make cleanall`.

## Customization of the benchmark

The project contains several configuration files:

  * [k2tree.config][K2CONFIG]: Specify different k2 tree implementations.
  * [test_case.config][TCCONF]: Specify test instances by ID, path, LaTeX-name
                                for the report, and download URL.
  * [compile_options.config][CCONF]: Specify compile options by option string.

Note that the benchmark will execute every combination of k2 trees and test cases.

[RPJ]: http://www.r-project.org/ "R"
[LT]: http://www.tug.org/applications/pdftex/ "pdflatex"
[K2CONFIG]: ./k2tree.config "k2tree.config"
[TCCONF]: ./test_case.config "test_case.config"
[CCONF]: ./compile_options.config "compile_options.config"
[RES]: https://users.dcc.uchile.cl/~fmontoto/static/k2.pdf "k2.pdf"
