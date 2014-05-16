# Benchmarking wavelet trees 

## Methodology

Explored dimensions:
  
  * self - delimiting code implementations
  * test cases
  * methods (`encoding`, `decoding`) 

## Directory structure

  * [bin](./bin): Contains the executables of the project.
  * [results](./results): Contains the results of the experiments.
  * [src](./src):  Contains the source code of the benchmark.
  * [visualize](./visualize): Contains LaTex files and a makefile for generating a report

## Prerequisites

  * To run the test on larger test cases (>= 200 MB), you should have at least 2 GB 
    of free memory (some vectors have very poor compression).
  * For the visualization you need the following software:
    - [pdflatex][LT] to generate the pdf reports.
    - [pgfplots][PGFP] installed in [LT] to generate plots in pdf reports.
		
## Usage

 * `make timing` compiles the programs, downloads or generates
    the test instances, builds the compression vectors, 
    runs the performance tests and generated a report located at
   `visualize/self_delimiting_codes.pdf`. The raw numbers of the encoding / decoding
   rates and compression can be found in the file `results/result.csv`.
   The used test cases can be found in file `results/tc.csv`.
   The tested vectors can be found in file `results/vat.csv`.
   The default benchmark took 14 minutes on my machine (Asus P50IJ
   Pentium(R) Dual-Core CPU T4500 @ 2.30GHz 2GB).
 * All created binaries and test results can be deleted
   by calling `make cleanall`.

## Customization of the benchmark

The project contains several configuration files:
 
  * [vectors.config][VCONFIG]: Specify different compression vectors and their used coders.
  * [test_case.config][TCCONFIG]: Specify test instances by ID, path, LaTeX-name 
                                for the report, and download URL.
  * [compile_options.config][CCONFIG]: Specify compile options by option string.

Note that the benchmark will execute every combination of vectors and test cases.

[LT]: http://www.tug.org/applications/pdftex/ "pdflatex"
[PGFP]: http://www.ctan.org/pkg/pgfplots "pgfplots"
[VCONFIG]: ./vectors.config "vectors.config"
[TCCONFIG]: ./test_case.config "test_case.config"
[CCONFIG]: ./compile_options.config "compile_options.config"
