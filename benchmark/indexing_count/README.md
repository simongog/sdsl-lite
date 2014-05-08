# Benchmarking operation `count` on FM-indexes

## Methodology

Explored dimensions:
  
  * text type 
  * instance size (just adjust the test_case.config file for this)
  * compile options
  * index implementations

Pattern selection:

We use the benchmark code including the random pattern selection
and test cases from the [Pizza&Chili][pz] website.

## Directory structure

  * [bin](./bin): Contains the executables of the project.
    * `build_idx_*` generates indexes
    * `query_idx_*` executes the count experiments 
    * `info_*` outputs the space breakdown of an index.
    * `genpattern` pattern generation from [Pizza&Chili][pz] website.
  * [indexes](./indexes): Contains the generated indexes.
  * [results](./visualize): Contains the results of the experiments.
  * [src](./src):  Contains the source code of the benchmark.
  * [visualize](./visualize): Contains a `R`-script which generates
               a report.

	Files included in this archive from the Pizza&Chili website:
	  * [src/genpatterns.c](./src/genpatterns.c)
      * [src/run_quries_sdsl.cpp](./src/run_quries_sdsl.cpp) 
	    is a customized version of the Pizza&Chili file run_queries.c .

## Prerequisites
  * For the visualization you need the following software:
    - [R][RPJ] with packages `xtable`,`plyr`. You can install the
      package by calling `install.packages("xtable", "plyr")` in R.
    - [pdflatex][LT] to generate the pdf reports.
  * The construction of the 200MB indexes requires about 1GB
    of RAM.
		
## Usage

 * `make timing`  compiles the programs, downloads test the
   200MB [Pizza&Chili][pz] test cases, builds the indexes,
   runs the performance tests, and generated a report located at
   `visualize/count.pdf`. The raw numbers of the timings
   can be found in `results/all.txt`. 
   Indexes and temporary files are stored in the
   directory `indexes` and `tmp`. For the 5 x 200 MB of
   [Pizza&Chili][pz] data the project will produce about
   7.2 GB of additional data. On my machine (MacBookPro Retina
   2.6GHz Intel Core i7, 16GB 1600 Mhz DDR3, SSD) the
   benchmark, invoced by `make timing`, took about 11 minutes
   (excluding the time to download the test instances).
   Have a look at the [generated report][RES].
 * All created indexes and test results can be deleted
   by calling `make cleanall`.

## Customization of the benchmark
  The project contains several configuration files:
 
  * [index.config](./index.config): Specify data structures'
			ID, sdsl-class, and LaTeX-name for the report.
  * [test_case.config](./test_case.config): Specify test cases' 
			ID, path, LaTeX-name for the report, and download URL.
  * [compile_options.config](./compile_options.config): Specify 
			compile options' ID and string.

  Note that the benchmark will execute every combination of your choices. 

  Finally, the visualization can also be configured:

  * [visualize/index-filter.config](./visualize/index-filter.config): 
	  Specify which indexes should be listed in the report. 

[sdsl]: https://github.com/simongog/sdsl "sdsl"
[pz]: http://pizzachili.di.unipi.it "Pizza&Chili"
[RPJ]: http://www.r-project.org/ "R"
[LT]: http://www.tug.org/applications/pdftex/ "pdflatex"
[RES]: https://github.com/simongog/simongog.github.com/raw/master/assets/images/count.pdf "count.pdf"
