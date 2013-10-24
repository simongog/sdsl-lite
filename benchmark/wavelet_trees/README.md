# Benchmarking wavelet trees 

## Methodology

Explored dimensions:
  
  * block size K (K in [5..255])
  * wavelet tree implementations
  * test cases
  * methods (`access`, `rank`, `select`) 

## Directory structure

  * [bin](./bin): Contains the executables of the project.
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

 * `make all` compiles the programs, downloads or generates
    the test instances, builds the wavelet trees, 
    runs the performance tests and generated a report located at
   `visualize/wt.pdf`. The raw numbers of the timings 
   can be found in the `results/all.txt`.

## Customization of the benchmark
  The project contains several configuration files:
 
  * [wt.config][WTCONFIG]: Specify different wavelet tree implementations.
  * [test_case.config][TCCONF]: Specify test instances by
       ID, path, LaTeX-name for the report, and download URL.
  * [compile_options.config][CCONF]: Specify compile
    options by option string.

  Note that the benchmark will execute every combination of your
  choices. 

[RPJ]: http://www.r-project.org/ "R"
[LT]: http://www.tug.org/applications/pdftex/ "pdflatex"
[WTCONFIG]: ./wt.config "wt.config"
[TCCONF]: ./test_case.config "test_case.config"
[CCONF]: ./compile_options.config "compile_options.config"
