# Benchmarking operation `count` on FM-indexes



## Methodology

Explored Dimensions:
  
  * text type
  * text size (just adjust the test_case.config file for this)
  * compile options
  * index implementations

Pattern selection:

We use the benchmark code including the random pattern selection
and test cases form the [Pizza&Chili][pz] website.

## Directory structure

  * bin: Contains the executables of the project.
    * `build_idx_*` generates indexes
    * `query_idx_*` executes the count experiments 
    * `info_*` outputs the space breakdown of an index.
    * `genpattern` pattern generation form [Pizza&Chili][pz] website.
  * src:  Contains the source code of the benchmark.
  * visualize: `R` and `pdflatex` scripts to generate
               reports of the collected data.
  * results: Contains the results of the experiments.

	Files included in this archive form the Pizza&Chili website:
	  * src/genpatterns.c
      * src/run_quries_sdsl.cpp is a adapted version of the
	    Pizza&Chili file run_queries.c .

## Prerequisites
  * For the visualization you need the following software:
    - [R][RPJ] with package `xtable`. You can install the
      package by calling `install.packages("xtable")` in R.
    - [pdflatex][LT] to generate the pdf reports.

		
## Usage

 * Populate the data directory with some test files. E.g. you
   can download the 200MB [Pizza&Chili][pz] test cases by 
   executing `../data/get_corpus.sh`. 
 * `make timing`  compiles the programs, builds the indexes and
   runs the performance tests and generated a report located as
   `visualization/count.pdf`. The raw numbers of the timing
   can be found in the `results/all.txt`. 
   Indexes and temporary files are stored in the
   directory `data` and `tmp`. For the 1GB of
   [Pizza&Chili][pz] data the project will produce about
   12GB of additional data. On a machine equipped with
   2GB of main memory and a recent processor the benchmark
   will take about 45 minutes.
 * All created indexes and test results can be deleted
   by calling `make cleanall`.

## Customization of the benchmark
  The project contains several configuration files:
 
  * `index.config`: Here you can specify the data structures. Select a
                    ID, sdsl-class and a LaTeX name for the report.
  * `test_case.config`: Specify ID, path to the test case and LaTeX
                        name for the report.
  * `compile_options.config`: Specify ID and compile options string.

  Note that the benchmark will execute every combination of your
  choices. 

  Finally, the visualization can also be configured:

  * `visualization/count.config` Here you can specify which
  indexes should be listed in the report. 

