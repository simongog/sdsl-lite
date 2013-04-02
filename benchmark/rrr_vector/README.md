# Benchmarking bitvector [rrr_vector](../../include/sdsl/rrr_vector.hpp)

## Methodology

Explored dimensions:
  
  * block size K (K in [5..255])
  * instance type (artificial, WT based)
  * instance size (MB and GB range)
  * methods (`access`, `rank`, `select`) 
  * compile options

## Directory structure

  * [bin](./bin): Contains the executables of the project.
    * `rrr_time_and_space_*` generates `rrr_vector`s, answers
               queries and outputs space information.
    * `generate_rnd_bitvector` generates evenly distributed bitvectors.
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
  * The processing of the 1GB bitvectors requires about
    1.5GB of RAM.
		
## Usage

 * `make timing` compiles the programs, downloads or generates
    the test instances, builds the compressed bitvectors, 
    runs the performance tests and generated a report located at
   `visualize/rrr.pdf`. The raw numbers of the timings 
   can be found in the `results/all.txt`.
   On my machine (MacBookPro Retina 2.6GHz Intel Core i7,
   16GB 1600 Mhz DDR3, SSD) the benchmark took about 58 minutes
   (excluding the time to download the test instances).
   Have a look at the [generated report][RES].
 * All created indexes and test results can be deleted
   by calling `make cleanall`.

## Customization of the benchmark
  The project contains several configuration files:
 
  * [block_size.config][KCONFIG]: Specify different block sizes.
  * [test_case.config][TCCONF]: Specify test instances by
       ID, path, LaTeX-name for the report, and download URL.
  * [compile_options.config][CCONF]: Specify compile
    options by ID and option string.

  Note that the benchmark will execute every combination of your
  choices. 

[RPJ]: http://www.r-project.org/ "R"
[LT]: http://www.tug.org/applications/pdftex/ "pdflatex"
[RES]: https://github.com/simongog/simongog.github.com/raw/master/assets/images/rrr.pdf "rrr.pdf"
[KCONFIG]: ./block_size.config "block_size.config"
[TCCONF]: ./test_case.config "test_case.config"
[CCONF]: ./compile_options.config "compile_options.config"
[VCONF]: ./visualize/index-filter.config "index-filter.config"
