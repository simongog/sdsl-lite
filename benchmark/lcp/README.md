# Benchmarking LCP algorithms 

## Methodology

Explored dimensions:
  
  * lcp algorithms
  * test cases

## Directory structure

  * [bin](./bin): Contains the executables of the project.
  * [results](./results): Contains the results of the experiments.
  * [src](./src):  Contains the source code of the benchmark.
  * [visualize](./visualize): Contains a `R`-script which generates
							  a report in LaTeX format.

## Prerequisites

  * For the visualization you need the following software:
    - [R][RPJ] with package `xtable`. You can install the
      package by calling `install.packages("xtable")` in R.
    - [pdflatex][LT] to generate the pdf reports.
		
## Usage

 * `make timing` compiles the programs, downloads
    the test instances, builds the LCP arrays and generates a report located at
   `visualize/lcp.pdf`. The raw numbers of the timings 
   can be found in the `results/all.txt`. The execution of the
   default benchmark took 66 minutes on my machine (MacBookPro Retina
   2.6GHz Intel Core i7, 16GB 1600 Mhz DDR3, SSD). 
   Have a look at the [generated report][RES].
 * All created binaries and test results can be deleted
   by calling `make cleanall`.

## Customization of the benchmark

  The project contains several configuration files:
 
  * [wt.config][LCPCONFIG]: Specify different LCP algorithms.
  * [test_case.config][TCCONF]: Specify test instances by ID, path, LaTeX-name 
    for the report, and download URL.
  * [compile_options.config][CCONF]: Specify compile options by option string.

  Note that the benchmark will execute every combination of lcp algorithms and test cases.

[RPJ]: http://www.r-project.org/ "R"
[LT]: http://www.tug.org/applications/pdftex/ "pdflatex"
[LCPCONFIG]: ./lcp.config "lcp.config"
[TCCONF]: ./test_case.config "test_case.config"
[CCONF]: ./compile_options.config "compile_options.config"
[RES]: https://github.com/simongog/simongog.github.com/raw/master/assets/images/lcp.pdf "lcp.pdf"
